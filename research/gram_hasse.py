#!/usr/bin/env python3
"""Exact Hasse--Minkowski obstruction filter for order-23 Gram snapshots.

If an invertible rational matrix R satisfies G = R R^T, then G is rationally
congruent to the identity.  This script checks that necessary condition using
exact LDL^T diagonalization and exact Hilbert symbols at every relevant place.

No-obstruction is deliberately not reported as a decomposition: rational
congruence neither constructs R nor requires its entries to be in {-1, +1}.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tempfile
from fractions import Fraction
from pathlib import Path
from typing import Any


ORDER = 23
BINARY_NORMALIZATION = "G=24I-J+4A"
MULTILEVEL_NORMALIZATION = "diag=23;offdiag={-5,-1,3}"
FOUR_LEVEL_NORMALIZATION = "diag=23;offdiag={-5,-1,3,7}"
SUPPORTED_SNAPSHOTS = {
    "gram-tabu": {BINARY_NORMALIZATION},
    "gram-gm-switch": {BINARY_NORMALIZATION},
    "gram-multilevel-tabu": {
        MULTILEVEL_NORMALIZATION,
        FOUR_LEVEL_NORMALIZATION,
    },
}
MAX_FACTOR_BITS = 64
POLLARD_ATTEMPTS = 32
POLLARD_STEPS_PER_ATTEMPT = 250_000
SMALL_PRIMES = (
    2,
    3,
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
    47,
    53,
    59,
    61,
    67,
    71,
    73,
    79,
    83,
    89,
    97,
)


class InputError(ValueError):
    """The snapshot or selected hit is malformed or internally inconsistent."""


class FactorizationError(RuntimeError):
    """Bounded deterministic factorization did not finish."""


def parse_decimal(value: Any, field: str) -> int:
    if not isinstance(value, str) or not value or not value.isascii():
        raise InputError(f"{field} must be a non-empty ASCII decimal string")
    if value[0] == "-":
        digits = value[1:]
    else:
        digits = value
    if not digits or not digits.isdigit():
        raise InputError(f"{field} must be an ASCII decimal string")
    return int(value)


def parse_edge_set(
    hit: dict[str, Any], edges_field: str, count_field: str
) -> set[tuple[int, int]]:
    edges = hit.get(edges_field)
    if not isinstance(edges, list):
        raise InputError(f"{edges_field} must be a list")
    edge_count = hit.get(count_field)
    if type(edge_count) is not int or edge_count != len(edges):
        raise InputError(f"{count_field} does not match {edges_field}")

    edge_set: set[tuple[int, int]] = set()
    for position, edge in enumerate(edges):
        if (
            not isinstance(edge, list)
            or len(edge) != 2
            or isinstance(edge[0], bool)
            or isinstance(edge[1], bool)
            or not isinstance(edge[0], int)
            or not isinstance(edge[1], int)
        ):
            raise InputError(
                f"{edges_field}[{position}] must contain two integer vertices"
            )
        left, right = edge
        if not (1 <= left < right <= ORDER):
            raise InputError(
                f"{edges_field}[{position}] must satisfy "
                f"1 <= left < right <= {ORDER}"
            )
        pair = (left - 1, right - 1)
        if pair in edge_set:
            raise InputError(f"{edges_field}[{position}] is duplicated")
        edge_set.add(pair)
    return edge_set


def build_gram(
    hit: dict[str, Any], normalization: str
) -> list[list[int]]:
    if normalization == BINARY_NORMALIZATION:
        level_edges = {3: parse_edge_set(hit, "edges", "edge_count")}
    elif normalization in (
        MULTILEVEL_NORMALIZATION,
        FOUR_LEVEL_NORMALIZATION,
    ):
        minus5_edges = parse_edge_set(
            hit, "minus5_edges", "minus5_count"
        )
        plus3_edges = parse_edge_set(hit, "plus3_edges", "plus3_count")
        plus7_edges = (
            parse_edge_set(hit, "plus7_edges", "plus7_count")
            if normalization == FOUR_LEVEL_NORMALIZATION
            else set()
        )
        overlap = (
            (minus5_edges & plus3_edges)
            | (minus5_edges & plus7_edges)
            | (plus3_edges & plus7_edges)
        )
        if overlap:
            raise InputError(
                "non-default level edge arrays must be pairwise disjoint"
            )
        level_edges = {
            -5: minus5_edges,
            3: plus3_edges,
            7: plus7_edges,
        }
    else:
        raise InputError(f"unsupported normalization {normalization!r}")

    gram: list[list[int]] = []
    for row in range(ORDER):
        values: list[int] = []
        for column in range(ORDER):
            if row == column:
                values.append(23)
            else:
                edge = (min(row, column), max(row, column))
                values.append(
                    next(
                        (
                            level
                            for level, edges in level_edges.items()
                            if edge in edges
                        ),
                        -1,
                    )
                )
        gram.append(values)
    return gram


def bareiss_determinant(matrix: list[list[int]]) -> int:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise InputError("determinant input must be a non-empty square matrix")
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1

    for pivot_column in range(size - 1):
        pivot_row = next(
            (
                row
                for row in range(pivot_column, size)
                if work[row][pivot_column] != 0
            ),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != pivot_column:
            work[pivot_column], work[pivot_row] = (
                work[pivot_row],
                work[pivot_column],
            )
            sign = -sign

        pivot = work[pivot_column][pivot_column]
        for row in range(pivot_column + 1, size):
            for column in range(pivot_column + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_column] * work[pivot_column][column]
                )
                quotient, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("non-exact division in Bareiss elimination")
                work[row][column] = quotient
            work[row][pivot_column] = 0
        previous = pivot

    return sign * work[-1][-1]


def positive_ldl_diagonal(matrix: list[list[int]]) -> list[Fraction] | None:
    """Return exact positive LDL^T pivots, or None if the form is not PD."""

    size = len(matrix)
    lower = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    diagonal: list[Fraction] = []
    for row in range(size):
        lower[row][row] = Fraction(1)
        pivot = Fraction(matrix[row][row]) - sum(
            lower[row][column] * lower[row][column] * diagonal[column]
            for column in range(row)
        )
        if pivot <= 0:
            return None
        diagonal.append(pivot)
        for next_row in range(row + 1, size):
            residual = Fraction(matrix[next_row][row]) - sum(
                lower[next_row][column]
                * lower[row][column]
                * diagonal[column]
                for column in range(row)
            )
            lower[next_row][row] = residual / pivot
    return diagonal


def is_prime_64(value: int) -> bool:
    """Deterministic Miller--Rabin primality test for unsigned 64-bit inputs."""

    if value < 2:
        return False
    for prime in SMALL_PRIMES:
        if value % prime == 0:
            return value == prime

    exponent = value - 1
    shifts = 0
    while exponent % 2 == 0:
        shifts += 1
        exponent //= 2

    # Jim Sinclair's deterministic witness set for n < 2^64.
    for base in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if base % value == 0:
            continue
        witness = pow(base, exponent, value)
        if witness in (1, value - 1):
            continue
        for _ in range(shifts - 1):
            witness = witness * witness % value
            if witness == value - 1:
                break
        else:
            return False
    return True


def deterministic_parameter(value: int, attempt: int, label: str) -> int:
    digest = hashlib.sha256(f"{value}:{attempt}:{label}".encode("ascii")).digest()
    return int.from_bytes(digest[:8], "big")


def pollard_brent(value: int) -> int:
    """Return a nontrivial factor, with deterministic and bounded attempts."""

    for prime in SMALL_PRIMES:
        if value % prime == 0:
            return prime

    for attempt in range(POLLARD_ATTEMPTS):
        current = 1 + deterministic_parameter(value, attempt, "y") % (value - 1)
        constant = 1 + deterministic_parameter(value, attempt, "c") % (value - 1)
        batch = 32 + deterministic_parameter(value, attempt, "m") % 225
        factor = 1
        radius = 1
        steps = 0
        saved = current
        anchor = current

        while factor == 1 and steps < POLLARD_STEPS_PER_ATTEMPT:
            anchor = current
            for _ in range(radius):
                current = (current * current + constant) % value
                steps += 1
                if steps >= POLLARD_STEPS_PER_ATTEMPT:
                    break

            offset = 0
            product = 1
            while (
                offset < radius
                and factor == 1
                and steps < POLLARD_STEPS_PER_ATTEMPT
            ):
                saved = current
                count = min(batch, radius - offset)
                for _ in range(count):
                    current = (current * current + constant) % value
                    product = product * abs(anchor - current) % value
                    steps += 1
                    if steps >= POLLARD_STEPS_PER_ATTEMPT:
                        break
                factor = math.gcd(product, value)
                offset += count
            radius *= 2

        if factor == value:
            factor = 1
            while factor == 1 and steps < POLLARD_STEPS_PER_ATTEMPT:
                saved = (saved * saved + constant) % value
                factor = math.gcd(abs(anchor - saved), value)
                steps += 1
        if 1 < factor < value:
            return factor

    raise FactorizationError(
        "deterministic Pollard-Brent budget exhausted "
        f"({POLLARD_ATTEMPTS} attempts, "
        f"{POLLARD_STEPS_PER_ATTEMPT} steps each)"
    )


_FACTOR_CACHE: dict[int, dict[int, int]] = {}


def factor_integer(value: int) -> dict[int, int]:
    """Factor a positive unsigned 64-bit integer exactly."""

    if value <= 0:
        raise FactorizationError("factorization input must be positive")
    if value.bit_length() > MAX_FACTOR_BITS:
        raise FactorizationError(
            f"factorization input exceeds the supported {MAX_FACTOR_BITS}-bit bound"
        )
    if value in _FACTOR_CACHE:
        return _FACTOR_CACHE[value].copy()

    factors: list[int] = []
    pending = [value]
    while pending:
        part = pending.pop()
        if part == 1:
            continue
        if is_prime_64(part):
            factors.append(part)
            continue
        divisor = pollard_brent(part)
        if divisor in (1, part) or part % divisor:
            raise FactorizationError("Pollard-Brent returned an invalid divisor")
        pending.extend((divisor, part // divisor))

    result: dict[int, int] = {}
    for prime in sorted(factors):
        result[prime] = result.get(prime, 0) + 1
    product = math.prod(prime**exponent for prime, exponent in result.items())
    if product != value or any(not is_prime_64(prime) for prime in result):
        raise FactorizationError("factorization failed its exact product check")
    _FACTOR_CACHE[value] = result.copy()
    return result


def valuation_and_unit(value: Fraction, prime: int, modulus: int) -> tuple[int, int]:
    numerator = value.numerator
    denominator = value.denominator
    valuation = 0
    while numerator % prime == 0:
        numerator //= prime
        valuation += 1
    while denominator % prime == 0:
        denominator //= prime
        valuation -= 1
    unit = numerator % modulus * pow(denominator % modulus, -1, modulus) % modulus
    return valuation, unit


def hilbert_symbol(left: Fraction, right: Fraction, prime: int) -> int:
    """Return the exact Hilbert symbol (left, right)_prime in {-1, +1}."""

    if left == 0 or right == 0:
        raise ArithmeticError("Hilbert symbols require nonzero arguments")
    if prime == 2:
        left_value, left_unit = valuation_and_unit(left, 2, 8)
        right_value, right_unit = valuation_and_unit(right, 2, 8)
        exponent = (
            ((left_unit - 1) // 2) * ((right_unit - 1) // 2)
            + left_value * ((right_unit * right_unit - 1) // 8)
            + right_value * ((left_unit * left_unit - 1) // 8)
        )
        return -1 if exponent % 2 else 1

    left_value, left_unit = valuation_and_unit(left, prime, prime)
    right_value, right_unit = valuation_and_unit(right, prime, prime)
    exponent = left_value * right_value * ((prime - 1) // 2)
    if right_value % 2 and pow(left_unit, (prime - 1) // 2, prime) == prime - 1:
        exponent += 1
    if left_value % 2 and pow(right_unit, (prime - 1) // 2, prime) == prime - 1:
        exponent += 1
    return -1 if exponent % 2 else 1


def hasse_invariant(diagonal: list[Fraction], prime: int) -> int:
    invariant = 1
    for right in range(len(diagonal)):
        for left in range(right):
            invariant *= hilbert_symbol(diagonal[left], diagonal[right], prime)
    return invariant


def analyze_hit(
    hit: dict[str, Any], hit_index: int, normalization: str
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "hit_index": hit_index,
        "status": "error",
        "bad_primes": [],
        "claim": "necessary_condition_only",
    }
    try:
        gram = build_gram(hit, normalization)
        determinant = bareiss_determinant(gram)
        claimed_determinant = parse_decimal(hit.get("determinant"), "determinant")
        if determinant != claimed_determinant:
            raise InputError(
                "claimed determinant does not match the reconstructed Gram matrix"
            )
        if determinant <= 0:
            result.update(
                {
                    "determinant": str(determinant),
                    "status": "rejected",
                    "rejection_reasons": ["nonpositive_determinant"],
                }
            )
            return result

        square_root = math.isqrt(determinant)
        if square_root * square_root != determinant:
            result.update(
                {
                    "determinant": str(determinant),
                    "status": "rejected",
                    "rejection_reasons": ["determinant_not_a_square"],
                }
            )
            return result
        claimed_root = parse_decimal(hit.get("square_root"), "square_root")
        if square_root != claimed_root:
            raise InputError(
                "claimed square_root does not match the exact Gram determinant"
            )

        diagonal = positive_ldl_diagonal(gram)
        positive_definite = diagonal is not None
        if hit.get("positive_definite") is not positive_definite:
            raise InputError(
                "claimed positive_definite does not match exact LDL^T elimination"
            )
        if not positive_definite:
            result.update(
                {
                    "determinant": str(determinant),
                    "square_root": str(square_root),
                    "status": "rejected",
                    "rejection_reasons": ["not_positive_definite"],
                }
            )
            return result
        assert diagonal is not None
        if math.prod(diagonal, start=Fraction(1)) != determinant:
            raise ArithmeticError("LDL^T pivots disagree with the exact determinant")

        divisible = square_root % (1 << 22) == 0
        if hit.get("divisible_by_2_22") is not divisible:
            raise InputError(
                "claimed divisible_by_2_22 does not match the exact square root"
            )
        qualified = divisible and positive_definite
        if hit.get("qualified") is not qualified:
            raise InputError("claimed qualified flag does not match exact checks")

        factors = factor_integer(square_root)
        tested_primes = sorted({2, *factors})
        invariants = {
            prime: hasse_invariant(diagonal, prime) for prime in tested_primes
        }
        bad_primes = [
            prime for prime in tested_primes if invariants[prime] == -1
        ]

        # For an integral form, odd primes outside det(G) have invariant +1.
        # Hilbert reciprocity therefore gives a useful independent consistency
        # check on the finite list above (the positive-definite real invariant
        # is +1).
        product_formula_consistent = math.prod(invariants.values()) == 1
        if not product_formula_consistent:
            raise ArithmeticError(
                "Hilbert reciprocity check failed; refusing to classify the hit"
            )

        result.update(
            {
                "determinant": str(determinant),
                "square_root": str(square_root),
                "square_root_factorization": {
                    str(prime): exponent for prime, exponent in factors.items()
                },
                "positive_definite": True,
                "tested_primes": tested_primes,
                "hasse_invariants": {
                    str(prime): invariants[prime] for prime in tested_primes
                },
                "bad_primes": bad_primes,
                "hilbert_reciprocity_checked": True,
            }
        )
        if bad_primes:
            result.update(
                {
                    "status": "rejected",
                    "rejection_reasons": ["finite_hasse_invariant"],
                }
            )
        else:
            result.update(
                {
                    "status": "no_obstruction",
                    "rejection_reasons": [],
                    "warning": (
                        "No rational obstruction is not a {-1,+1} "
                        "decomposition and does not construct a factor."
                    ),
                }
            )
    except (ArithmeticError, FactorizationError, InputError) as error:
        result["error"] = str(error)
    return result


def select_indices(arguments: argparse.Namespace, hit_count: int) -> list[int]:
    if arguments.all_hits:
        return list(range(hit_count))
    indices: list[int] = []
    seen: set[int] = set()
    for index in arguments.hit_index:
        if not 0 <= index < hit_count:
            raise InputError(
                f"hit index {index} is outside the available range 0..{hit_count - 1}"
            )
        if index not in seen:
            indices.append(index)
            seen.add(index)
    return indices


def atomic_write(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as output:
            output.write(payload)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Reject order-23 Gram candidates with an exact rational "
            "Hasse--Minkowski obstruction."
        )
    )
    parser.add_argument("--snapshot", required=True, type=Path)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument(
        "--hit-index",
        action="append",
        type=int,
        default=[],
        help="zero-based hit index; repeat to select multiple hits",
    )
    selection.add_argument(
        "--all-hits",
        action="store_true",
        help="analyze every hit in the snapshot",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="atomically write JSON here instead of standard output",
    )
    return parser


def run(arguments: argparse.Namespace) -> tuple[dict[str, Any], int]:
    try:
        snapshot_bytes = arguments.snapshot.read_bytes()
        snapshot = json.loads(snapshot_bytes)
        if not isinstance(snapshot, dict):
            raise InputError("snapshot root must be an object")
        if snapshot.get("challenge_id") != "maxdet-23-v1":
            raise InputError("snapshot challenge_id must be maxdet-23-v1")
        engine = snapshot.get("engine")
        if not isinstance(engine, str) or engine not in SUPPORTED_SNAPSHOTS:
            raise InputError(
                "snapshot engine must be one of "
                f"{sorted(SUPPORTED_SNAPSHOTS)!r}"
            )
        if type(snapshot.get("schema_version")) is not int:
            raise InputError("snapshot schema_version must be an integer")
        if snapshot.get("schema_version") != 1:
            raise InputError("snapshot schema_version must be 1")
        normalization = snapshot.get("normalization")
        if normalization not in SUPPORTED_SNAPSHOTS[engine]:
            raise InputError(
                f"unsupported snapshot normalization "
                f"{normalization!r} for engine {engine!r}"
            )
        hits = snapshot.get("hits")
        if not isinstance(hits, list):
            raise InputError("snapshot hits must be a list")
        statistics = snapshot.get("statistics")
        parameters = snapshot.get("parameters")
        if not isinstance(statistics, dict) or not isinstance(parameters, dict):
            raise InputError("snapshot statistics and parameters must be objects")

        def source_unsigned(container: dict[str, Any], field: str) -> int:
            value = container.get(field)
            if type(value) is not int or value < 0:
                raise InputError(f"snapshot {field} must be an unsigned integer")
            return value

        complete = snapshot.get("complete")
        termination = snapshot.get("termination")
        if type(complete) is not bool:
            raise InputError("snapshot complete must be boolean")
        if not isinstance(termination, str):
            raise InputError("snapshot termination must be a string")
        indices = select_indices(arguments, len(hits))
        results = []
        for index in indices:
            hit = hits[index]
            if not isinstance(hit, dict):
                results.append(
                    {
                        "hit_index": index,
                        "status": "error",
                        "bad_primes": [],
                        "claim": "necessary_condition_only",
                        "error": "hit must be an object",
                    }
                )
            else:
                results.append(analyze_hit(hit, index, normalization))

        counts = {
            status: sum(result["status"] == status for result in results)
            for status in ("rejected", "no_obstruction", "error")
        }
        report = {
            "schema_version": 1,
            "engine": "gram-hasse",
            "method": "exact-rational-hasse-minkowski-obstruction",
            "snapshot": {
                "path": str(arguments.snapshot),
                "sha256": hashlib.sha256(snapshot_bytes).hexdigest(),
                "challenge_id": snapshot["challenge_id"],
                "engine": engine,
                "normalization": normalization,
                "complete": complete,
                "termination": termination,
                "stored_hit_count": len(hits),
                "max_stored_hits": source_unsigned(
                    parameters, "max_stored_hits"
                ),
                "exact_square_observations": source_unsigned(
                    statistics, "exact_squares"
                ),
                "qualified_survivor_observations": source_unsigned(
                    statistics, "qualified_survivors"
                ),
                "unrecorded_unique_square_observations": source_unsigned(
                    statistics, "unrecorded_square_observations"
                ),
            },
            "selection": {"hit_indices": indices},
            "selection_scope": (
                "all selected indices refer only to the snapshot's stored hits; "
                "search observations omitted by deduplication or capacity are "
                "not classified"
            ),
            "claim_boundary": (
                "rejected proves only that no rational R can satisfy G=RR^T; "
                "no_obstruction does not construct a factor and does not imply "
                "a {-1,+1} decomposition"
            ),
            "results": results,
            "summary": {
                "selected": len(results),
                **counts,
            },
        }
        return report, 1 if counts["error"] else 0
    except (OSError, UnicodeError, json.JSONDecodeError, InputError) as error:
        return (
            {
                "schema_version": 1,
                "engine": "gram-hasse",
                "status": "error",
                "error": str(error),
                "claim_boundary": "no mathematical classification was made",
            },
            2,
        )


def main() -> int:
    arguments = make_parser().parse_args()
    if (
        arguments.output is not None
        and arguments.output.resolve() == arguments.snapshot.resolve()
    ):
        print(
            json.dumps(
                {
                    "schema_version": 1,
                    "engine": "gram-hasse",
                    "status": "error",
                    "error": "--output must not overwrite --snapshot",
                },
                sort_keys=True,
            )
        )
        return 2
    report, return_code = run(arguments)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(payload, end="")
    else:
        try:
            atomic_write(arguments.output, payload)
        except OSError as error:
            print(
                json.dumps(
                    {
                        "schema_version": 1,
                        "engine": "gram-hasse",
                        "status": "error",
                        "error": str(error),
                    },
                    sort_keys=True,
                )
            )
            return 2
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
