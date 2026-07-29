#!/usr/bin/env python3
"""Exact local screens of the order-23 Gram defect graph.

This is an experimental research tool, not part of the trusted verifier.  It
reproduces the reference one-edge toggle and one-for-one edge-swap screens,
supports an explicit two-for-two reference edge-swap screen, and exhausts
additions of up to four cross-block edges to the Ehlich K3 + 5 K4 candidate
Gram graph.

The large screens are exact.  They use the matrix determinant lemma through an
integer adjugate correction and memoize only update patterns having identical
entries in the exact inverse correction matrix.  Every cached pattern is
independently checked by full integer Bareiss elimination.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
import time
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


ORDER = 23
FLOOR = 2_779_447_296_000_000
REQUIRED_DIVISOR = 1 << 22
REFERENCE_SQUARED = 7_725_327_271_241_711_616_000_000_000_000
EHLICH_SQUARED = 8_894_085_385_420_800_000_000_000_000_000
IDEAL_BLOCK_SIZES = (3, 4, 4, 4, 4, 4)

Edge = Tuple[int, int]
Update = Tuple[int, int, int]
Matrix = List[List[int]]


@dataclass
class Counters:
    expected_cases: int = 0
    shard_expected_cases: int = 0
    examined: int = 0
    exact_pattern_evaluations: int = 0
    exact_positive: int = 0
    exact_perfect_square: int = 0
    square_above_floor: int = 0
    square_divisible_by_2_22: int = 0
    positive_definite_after_numeric_gates: int = 0
    survivors: int = 0


@dataclass
class ModeResult:
    mode: str
    base: str
    base_defect_edge_count: int
    counters: Counters
    square_hits: List[dict] = field(default_factory=list)
    exact_pattern_cache_size: int = 0
    elapsed_seconds: float = 0.0
    complete: bool = True
    termination: str = "completed"
    last_case_index: Optional[int] = None


class Heartbeat:
    def __init__(self, mode: str, interval: float, expected: int) -> None:
        self.mode = mode
        self.interval = interval
        self.expected = expected
        self.started = time.monotonic()
        self.next = self.started + interval

    def maybe(self, counters: Counters) -> None:
        if self.interval <= 0 or counters.examined & 4095:
            return
        now = time.monotonic()
        if now < self.next:
            return
        event = {
            "elapsed_seconds": round(now - self.started, 3),
            "event": "heartbeat",
            "exact_pattern_evaluations": counters.exact_pattern_evaluations,
            "examined": counters.examined,
            "expected": self.expected,
            "mode": self.mode,
            "square_hits": counters.exact_perfect_square,
        }
        print(json.dumps(event, sort_keys=True, separators=(",", ":")),
              file=sys.stderr, flush=True)
        self.next = now + self.interval

    def elapsed(self) -> float:
        return time.monotonic() - self.started


def exact_determinant(matrix: Sequence[Sequence[int]]) -> int:
    """Fraction-free Bareiss determinant with exact integer divisions."""

    order = len(matrix)
    if order == 0:
        return 1
    if any(len(row) != order for row in matrix):
        raise ValueError("determinant requires a square matrix")
    work = [list(row) for row in matrix]
    previous_pivot = 1
    sign = 1
    for column in range(order - 1):
        pivot_row = column
        while pivot_row < order and work[pivot_row][column] == 0:
            pivot_row += 1
        if pivot_row == order:
            return 0
        if pivot_row != column:
            work[pivot_row], work[column] = (
                work[column],
                work[pivot_row],
            )
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, order):
            left = work[row][column]
            for inner in range(column + 1, order):
                numerator = (
                    work[row][inner] * pivot
                    - left * work[column][inner]
                )
                if column:
                    quotient, remainder = divmod(numerator, previous_pivot)
                    if remainder:
                        raise ArithmeticError(
                            "fraction-free determinant division failed"
                        )
                    work[row][inner] = quotient
                else:
                    work[row][inner] = numerator
            work[row][column] = 0
        previous_pivot = pivot
    return sign * work[-1][-1]


def fraction_determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    order = len(matrix)
    if order == 0:
        return Fraction(1)
    work = [list(row) for row in matrix]
    determinant = Fraction(1)
    sign = 1
    for column in range(order):
        pivot_row = column
        while pivot_row < order and not work[pivot_row][column]:
            pivot_row += 1
        if pivot_row == order:
            return Fraction(0)
        if pivot_row != column:
            work[pivot_row], work[column] = (
                work[column],
                work[pivot_row],
            )
            sign = -sign
        pivot = work[column][column]
        determinant *= pivot
        for row in range(column + 1, order):
            if not work[row][column]:
                continue
            factor = work[row][column] / pivot
            for inner in range(column + 1, order):
                work[row][inner] -= factor * work[column][inner]
            work[row][column] = Fraction(0)
    return sign * determinant


def exact_inverse(matrix: Sequence[Sequence[int]]) -> List[List[Fraction]]:
    order = len(matrix)
    augmented: List[List[Fraction]] = []
    for row_index, row in enumerate(matrix):
        augmented.append(
            [Fraction(value) for value in row]
            + [Fraction(int(row_index == column)) for column in range(order)]
        )
    for column in range(order):
        pivot_row = column
        while pivot_row < order and not augmented[pivot_row][column]:
            pivot_row += 1
        if pivot_row == order:
            raise ValueError("cannot invert singular matrix")
        if pivot_row != column:
            augmented[pivot_row], augmented[column] = (
                augmented[column],
                augmented[pivot_row],
            )
        pivot = augmented[column][column]
        augmented[column] = [value / pivot for value in augmented[column]]
        for row in range(order):
            if row == column:
                continue
            factor = augmented[row][column]
            if not factor:
                continue
            augmented[row] = [
                left - factor * right
                for left, right in zip(augmented[row], augmented[column])
            ]
    return [row[order:] for row in augmented]


def exact_positive_definite(matrix: Sequence[Sequence[int]]) -> bool:
    return all(
        exact_determinant([row[:size] for row in matrix[:size]]) > 0
        for size in range(1, len(matrix) + 1)
    )


def read_reference(path: Path) -> Matrix:
    rows = [
        [int(token) for token in line.split()]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(rows) != ORDER or any(len(row) != ORDER for row in rows):
        raise ValueError("reference must be a 23 by 23 matrix")
    if any(value not in (-1, 1) for row in rows for value in row):
        raise ValueError("reference contains an entry outside {-1,+1}")
    return rows


def gram(matrix: Sequence[Sequence[int]]) -> Matrix:
    return [
        [
            sum(left * right for left, right in zip(row, other))
            for other in matrix
        ]
        for row in matrix
    ]


def ideal_gram() -> Matrix:
    block_of: List[int] = []
    for block, size in enumerate(IDEAL_BLOCK_SIZES):
        block_of.extend([block] * size)
    return [
        [
            ORDER
            if row == column
            else 3
            if block_of[row] == block_of[column]
            else -1
            for column in range(ORDER)
        ]
        for row in range(ORDER)
    ]


def all_edges() -> List[Edge]:
    return [
        (first, second)
        for first in range(ORDER)
        for second in range(first + 1, ORDER)
    ]


ALL_EDGES = all_edges()
EDGE_NUMBERS = {edge: index + 1 for index, edge in enumerate(ALL_EDGES)}


def edges_with_value(matrix: Matrix, value: int) -> List[Edge]:
    return [edge for edge in ALL_EDGES if matrix[edge[0]][edge[1]] == value]


def apply_updates(base: Matrix, updates: Sequence[Update]) -> Matrix:
    candidate = [row[:] for row in base]
    for first, second, delta in updates:
        candidate[first][second] += delta
        candidate[second][first] += delta
    return candidate


def inverse_identity_check(
    matrix: Matrix, inverse: Sequence[Sequence[Fraction]]
) -> None:
    for row in range(ORDER):
        for column in range(ORDER):
            value = sum(
                Fraction(matrix[row][inner]) * inverse[inner][column]
                for inner in range(ORDER)
            )
            if value != int(row == column):
                raise ArithmeticError("exact inverse identity check failed")


def exact_adjugate(
    determinant: int,
    inverse: Sequence[Sequence[Fraction]],
) -> Matrix:
    adjugate: Matrix = []
    for row in inverse:
        integer_row = []
        for value in row:
            scaled = determinant * value
            if scaled.denominator != 1:
                raise ArithmeticError("inverse does not yield an integer adjugate")
            integer_row.append(scaled.numerator)
        adjugate.append(integer_row)
    return adjugate


def integer_update_correction(
    base_determinant: int,
    adjugate: Matrix,
    updates: Sequence[Update],
) -> Matrix:
    u_indices: List[int] = []
    v_indices: List[int] = []
    coefficients: List[int] = []
    for first, second, delta in updates:
        u_indices.extend((first, second))
        v_indices.extend((second, first))
        coefficients.extend((delta, delta))
    rank = len(u_indices)
    return [
        [
            base_determinant * int(row == column)
            + coefficients[row] * adjugate[v_indices[row]][u_indices[column]]
            for column in range(rank)
        ]
        for row in range(rank)
    ]


def determinant_from_updates(
    base_determinant: int,
    adjugate: Matrix,
    updates: Sequence[Update],
) -> int:
    rank = 2 * len(updates)
    numerator = exact_determinant(
        integer_update_correction(base_determinant, adjugate, updates)
    )
    denominator = base_determinant ** (rank - 1)
    value, remainder = divmod(numerator, denominator)
    if remainder:
        raise ArithmeticError("determinant lemma produced a non-integer")
    return value


def inverse_code_table(
    inverse: Sequence[Sequence[Fraction]],
) -> Tuple[List[List[int]], List[Fraction]]:
    values = sorted({value for row in inverse for value in row})
    codes = {value: index for index, value in enumerate(values)}
    return [[codes[value] for value in row] for row in inverse], values


def update_pattern_key(
    inverse_codes: Sequence[Sequence[int]],
    updates: Sequence[Update],
) -> Tuple[int, ...]:
    """Exact correction-matrix key without expensive Fraction hashing."""

    u_indices: List[int] = []
    v_indices: List[int] = []
    coefficients: List[int] = []
    for first, second, delta in updates:
        u_indices.extend((first, second))
        v_indices.extend((second, first))
        coefficients.extend((delta, delta))
    rank = len(u_indices)
    return (
        rank,
        *coefficients,
        *(
            inverse_codes[v_indices[row]][u_indices[column]]
            for row in range(rank)
            for column in range(rank)
        ),
    )


def shard_size(total: int, shard_index: int, shard_count: int) -> int:
    if shard_index >= total:
        return 0
    return (total - 1 - shard_index) // shard_count + 1


def edge_json(edge: Edge) -> dict:
    return {
        "edge_number": EDGE_NUMBERS[edge],
        "first": edge[0] + 1,
        "second": edge[1] + 1,
    }


def evaluate_determinant(
    *,
    determinant: int,
    candidate_factory,
    case_index: int,
    defect_edge_count: int,
    added: Sequence[Edge],
    removed: Sequence[Edge],
    counters: Counters,
    square_hits: List[dict],
) -> None:
    if determinant <= 0:
        return
    counters.exact_positive += 1
    root = math.isqrt(determinant)
    if root * root != determinant:
        return

    # A square is rare enough to justify an independent full Bareiss check.
    candidate = candidate_factory()
    if exact_determinant(candidate) != determinant:
        raise ArithmeticError(
            "determinant-lemma result disagrees with full Bareiss"
        )
    counters.exact_perfect_square += 1
    above_floor = root > FLOOR
    divisible = root % REQUIRED_DIVISOR == 0
    if above_floor:
        counters.square_above_floor += 1
    if divisible:
        counters.square_divisible_by_2_22 += 1
    positive_definite = False
    if above_floor and divisible:
        positive_definite = exact_positive_definite(candidate)
        if positive_definite:
            counters.positive_definite_after_numeric_gates += 1
            counters.survivors += 1
    square_hits.append(
        {
            "above_floor": above_floor,
            "added_edges": [edge_json(edge) for edge in added],
            "case_index": case_index,
            "defect_edge_count": defect_edge_count,
            "determinant": str(determinant),
            "divisible_by_2_22": divisible,
            "positive_definite": positive_definite,
            "removed_edges": [edge_json(edge) for edge in removed],
            "square_root": str(root),
        }
    )


def run_reference_toggle(
    reference: Matrix, args: argparse.Namespace
) -> ModeResult:
    counters = Counters(
        expected_cases=len(ALL_EDGES),
        shard_expected_cases=shard_size(
            len(ALL_EDGES), args.shard_index, args.shard_count
        ),
    )
    result = ModeResult(
        mode="reference-toggle",
        base="published-reference",
        base_defect_edge_count=len(edges_with_value(reference, 3)),
        counters=counters,
    )
    heartbeat = Heartbeat(
        result.mode, args.heartbeat_seconds, counters.shard_expected_cases
    )
    for case_index, edge in enumerate(ALL_EDGES):
        if case_index % args.shard_count != args.shard_index:
            continue
        removing = reference[edge[0]][edge[1]] == 3
        updates = [(edge[0], edge[1], -4 if removing else 4)]
        candidate = apply_updates(reference, updates)
        determinant = exact_determinant(candidate)
        counters.examined += 1
        counters.exact_pattern_evaluations += 1
        evaluate_determinant(
            determinant=determinant,
            candidate_factory=lambda candidate=candidate: candidate,
            case_index=case_index,
            defect_edge_count=result.base_defect_edge_count
            + (-1 if removing else 1),
            added=[] if removing else [edge],
            removed=[edge] if removing else [],
            counters=counters,
            square_hits=result.square_hits,
        )
        heartbeat.maybe(counters)
    result.exact_pattern_cache_size = counters.exact_pattern_evaluations
    result.elapsed_seconds = heartbeat.elapsed()
    return result


def run_reference_swap(
    reference: Matrix, args: argparse.Namespace
) -> ModeResult:
    present = edges_with_value(reference, 3)
    absent = edges_with_value(reference, -1)
    total = len(present) * len(absent)
    counters = Counters(
        expected_cases=total,
        shard_expected_cases=shard_size(
            total, args.shard_index, args.shard_count
        ),
    )
    result = ModeResult(
        mode="reference-swap",
        base="published-reference",
        base_defect_edge_count=len(present),
        counters=counters,
    )
    heartbeat = Heartbeat(
        result.mode, args.heartbeat_seconds, counters.shard_expected_cases
    )
    case_index = 0
    for removed in present:
        for added in absent:
            current = case_index
            case_index += 1
            if current % args.shard_count != args.shard_index:
                continue
            updates = [
                (removed[0], removed[1], -4),
                (added[0], added[1], 4),
            ]
            candidate = apply_updates(reference, updates)
            determinant = exact_determinant(candidate)
            counters.examined += 1
            counters.exact_pattern_evaluations += 1
            evaluate_determinant(
                determinant=determinant,
                candidate_factory=lambda candidate=candidate: candidate,
                case_index=current,
                defect_edge_count=result.base_defect_edge_count,
                added=[added],
                removed=[removed],
                counters=counters,
                square_hits=result.square_hits,
            )
            heartbeat.maybe(counters)
    if case_index != total:
        raise AssertionError("reference swap completion count mismatch")
    result.exact_pattern_cache_size = counters.exact_pattern_evaluations
    result.elapsed_seconds = heartbeat.elapsed()
    return result


def reference_double_swap_cases(
    present: Sequence[Edge],
    absent: Sequence[Edge],
    shard_index: int = 0,
    shard_count: int = 1,
) -> Iterator[Tuple[int, Tuple[Edge, Edge], Tuple[Edge, Edge]]]:
    """Yield one shard of two removals and two additions.

    Global indices use the lexicographic rank of the removed pair as the major
    coordinate and the lexicographic rank of the added pair as the minor
    coordinate.  ``islice`` advances over other shards in C rather than
    yielding and rejecting their combinations in Python.
    """

    added_pair_count = math.comb(len(absent), 2)
    for removed_index, removed in enumerate(
        itertools.combinations(present, 2)
    ):
        block_offset = removed_index * added_pair_count
        local_start = (shard_index - block_offset) % shard_count
        owned_additions = itertools.islice(
            itertools.combinations(absent, 2),
            local_start,
            None,
            shard_count,
        )
        for owned_index, added in enumerate(owned_additions):
            added_index = local_start + owned_index * shard_count
            if added_index >= added_pair_count:
                raise AssertionError("double-swap added-pair index overflow")
            case_index = block_offset + added_index
            yield case_index, removed, added


def run_reference_double_swap(
    reference: Matrix, args: argparse.Namespace
) -> ModeResult:
    present = edges_with_value(reference, 3)
    absent = edges_with_value(reference, -1)
    total = math.comb(len(present), 2) * math.comb(len(absent), 2)
    counters = Counters(
        expected_cases=total,
        shard_expected_cases=shard_size(
            total, args.shard_index, args.shard_count
        ),
    )
    result = ModeResult(
        mode="reference-double-swap",
        base="published-reference",
        base_defect_edge_count=len(present),
        counters=counters,
    )
    heartbeat = Heartbeat(
        result.mode, args.heartbeat_seconds, counters.shard_expected_cases
    )

    base_determinant = exact_determinant(reference)
    inverse = exact_inverse(reference)
    inverse_identity_check(reference, inverse)
    adjugate = exact_adjugate(base_determinant, inverse)
    inverse_codes, _ = inverse_code_table(inverse)
    determinant_cache: Dict[Tuple[int, ...], int] = {}
    emitted = 0

    for case_index, removed, added in reference_double_swap_cases(
        present,
        absent,
        args.shard_index,
        args.shard_count,
    ):
        if (
            args.seconds
            and (counters.examined & 1023) == 0
            and heartbeat.elapsed() >= args.seconds
        ):
            result.complete = False
            result.termination = "time-limit"
            break
        emitted += 1
        updates = [
            (removed[0][0], removed[0][1], -4),
            (removed[1][0], removed[1][1], -4),
            (added[0][0], added[0][1], 4),
            (added[1][0], added[1][1], 4),
        ]
        pattern = update_pattern_key(inverse_codes, updates)
        determinant = determinant_cache.get(pattern)
        if determinant is None:
            determinant = determinant_from_updates(
                base_determinant, adjugate, updates
            )
            determinant_cache[pattern] = determinant
            counters.exact_pattern_evaluations += 1

            # Validate every structural cache insertion independently.  Cache
            # reuse is therefore an exact symmetry reduction, never a
            # heuristic determinant filter.
            if (
                exact_determinant(apply_updates(reference, updates))
                != determinant
            ):
                raise ArithmeticError(
                    "double-swap determinant lemma self-check failed"
                )

        counters.examined += 1
        result.last_case_index = case_index
        evaluate_determinant(
            determinant=determinant,
            candidate_factory=lambda updates=tuple(updates): apply_updates(
                reference, updates
            ),
            case_index=case_index,
            defect_edge_count=result.base_defect_edge_count,
            added=added,
            removed=removed,
            counters=counters,
            square_hits=result.square_hits,
        )
        heartbeat.maybe(counters)

    if result.complete and emitted != counters.shard_expected_cases:
        raise AssertionError("double-swap shard completion count mismatch")
    result.exact_pattern_cache_size = len(determinant_cache)
    result.elapsed_seconds = heartbeat.elapsed()
    return result


def ideal_addition_cases(
    absent: Sequence[Edge],
    max_added: int,
    shard_index: int = 0,
    shard_count: int = 1,
) -> Iterator[Tuple[int, Tuple[Edge, ...]]]:
    """Yield only this shard's combinations with their unsharded case index.

    ``itertools.islice`` performs the stepping in C.  A shard therefore does
    not enter Python merely to construct and reject every combination owned by
    another shard.
    """

    cardinality_offset = 0
    for count in range(1, max_added + 1):
        cardinality_cases = math.comb(len(absent), count)
        local_start = (shard_index - cardinality_offset) % shard_count
        owned = itertools.islice(
            itertools.combinations(absent, count),
            local_start,
            None,
            shard_count,
        )
        for owned_index, selected in enumerate(owned):
            case_index = (
                cardinality_offset
                + local_start
                + owned_index * shard_count
            )
            if case_index >= cardinality_offset + cardinality_cases:
                raise AssertionError("sharded combination index overflow")
            yield case_index, selected
        cardinality_offset += cardinality_cases


def run_ideal_add(ideal: Matrix, args: argparse.Namespace) -> ModeResult:
    absent = edges_with_value(ideal, -1)
    total = sum(math.comb(len(absent), count)
                for count in range(1, args.max_added + 1))
    counters = Counters(
        expected_cases=total,
        shard_expected_cases=shard_size(
            total, args.shard_index, args.shard_count
        ),
    )
    result = ModeResult(
        mode="ideal-add",
        base="ehlich-k3-plus-five-k4",
        base_defect_edge_count=len(edges_with_value(ideal, 3)),
        counters=counters,
    )
    heartbeat = Heartbeat(
        result.mode, args.heartbeat_seconds, counters.shard_expected_cases
    )

    base_determinant = exact_determinant(ideal)
    inverse = exact_inverse(ideal)
    inverse_identity_check(ideal, inverse)
    adjugate = exact_adjugate(base_determinant, inverse)
    inverse_codes, _ = inverse_code_table(inverse)
    determinant_cache: Dict[Tuple[int, ...], int] = {}
    emitted = 0

    for case_index, selected in ideal_addition_cases(
        absent,
        args.max_added,
        args.shard_index,
        args.shard_count,
    ):
        if (
            args.seconds
            and (counters.examined & 1023) == 0
            and heartbeat.elapsed() >= args.seconds
        ):
            result.complete = False
            result.termination = "time-limit"
            break
        emitted += 1
        updates = [(first, second, 4) for first, second in selected]
        pattern = update_pattern_key(inverse_codes, updates)
        determinant = determinant_cache.get(pattern)
        if determinant is None:
            determinant = determinant_from_updates(
                base_determinant, adjugate, updates
            )
            determinant_cache[pattern] = determinant
            counters.exact_pattern_evaluations += 1

            # There are only a few hundred symmetry patterns.  Check every
            # cache insertion against the independent full Bareiss path so a
            # bad cache key can never turn the exhaustive screen into an
            # unsafe filter.
            if exact_determinant(apply_updates(ideal, updates)) != determinant:
                raise ArithmeticError("determinant lemma self-check failed")

        counters.examined += 1
        result.last_case_index = case_index
        evaluate_determinant(
            determinant=determinant,
            candidate_factory=lambda updates=tuple(updates): apply_updates(
                ideal, updates
            ),
            case_index=case_index,
            defect_edge_count=result.base_defect_edge_count + len(selected),
            added=selected,
            removed=[],
            counters=counters,
            square_hits=result.square_hits,
        )
        heartbeat.maybe(counters)

    if result.complete and emitted != counters.shard_expected_cases:
        raise AssertionError("ideal shard completion count mismatch")
    result.exact_pattern_cache_size = len(determinant_cache)
    result.elapsed_seconds = heartbeat.elapsed()
    return result


def validate_bases(reference: Matrix, ideal: Matrix) -> None:
    if exact_determinant(reference) != REFERENCE_SQUARED:
        raise ValueError("reference Gram determinant mismatch")
    if exact_determinant(ideal) != EHLICH_SQUARED:
        raise ValueError("ideal Gram determinant mismatch")
    if not exact_positive_definite(reference):
        raise ValueError("reference Gram is not positive definite")
    if not exact_positive_definite(ideal):
        raise ValueError("ideal Gram is not positive definite")
    reference_edges = set(edges_with_value(reference, 3))
    ideal_edges = set(edges_with_value(ideal, 3))
    if len(reference_edges) != 45 or len(ideal_edges) != 33:
        raise ValueError("unexpected defect edge count")
    if not ideal_edges <= reference_edges:
        raise ValueError(
            "reference is not the ideal defect graph plus bridges"
        )


def result_json(
    args: argparse.Namespace, results: Sequence[ModeResult]
) -> dict:
    return {
        "challenge_id": "maxdet-23-v1",
        "claim_boundary": (
            "Gram screening supplies necessary conditions only; no matrix "
            "construction or optimality claim."
        ),
        "floor": str(FLOOR),
        "ideal_gram_determinant": str(EHLICH_SQUARED),
        "max_added": args.max_added,
        "mode": args.mode,
        "reference_gram_determinant": str(REFERENCE_SQUARED),
        "reference_matrix": args.reference.as_posix(),
        "results": [
            {
                "base": result.base,
                "base_defect_edge_count": result.base_defect_edge_count,
                "complete": result.complete,
                "counters": vars(result.counters),
                "elapsed_seconds": round(result.elapsed_seconds, 6),
                "exact_pattern_cache_size": result.exact_pattern_cache_size,
                "last_case_index": result.last_case_index,
                "mode": result.mode,
                "square_hits": result.square_hits,
                "termination": result.termination,
            }
            for result in results
        ],
        "schema_version": 1,
        "seconds": args.seconds,
        "shard_count": args.shard_count,
        "shard_index": args.shard_index,
    }


def atomic_write_json(path: Path, document: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(document, stream, sort_keys=True, separators=(",", ":"))
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=(
            "all",
            "reference-toggle",
            "reference-swap",
            "reference-double-swap",
            "ideal-add",
        ),
        default="all",
        help=(
            "'all' runs the established lightweight screens; the much larger "
            "reference-double-swap mode must be selected explicitly."
        ),
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("references/orrick-et-al-2003/matrix.txt"),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--max-added", type=int, choices=(1, 2, 3, 4), default=3
    )
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--heartbeat-seconds", type=float, default=10.0)
    parser.add_argument(
        "--seconds",
        type=float,
        default=0.0,
        help=(
            "Stop a bounded exhaustive shard after this search time; "
            "0 means complete."
        ),
    )
    args = parser.parse_args()
    if args.shard_count <= 0:
        parser.error("--shard-count must be positive")
    if not 0 <= args.shard_index < args.shard_count:
        parser.error("--shard-index must be in [0, shard-count)")
    if not math.isfinite(args.heartbeat_seconds) or args.heartbeat_seconds < 0:
        parser.error("--heartbeat-seconds must be finite and non-negative")
    if not math.isfinite(args.seconds) or args.seconds < 0:
        parser.error("--seconds must be finite and non-negative")
    return args


def main() -> int:
    args = parse_arguments()
    reference_matrix = read_reference(args.reference)
    reference = gram(reference_matrix)
    ideal = ideal_gram()
    validate_bases(reference, ideal)

    results: List[ModeResult] = []
    if args.mode in ("all", "reference-toggle"):
        results.append(run_reference_toggle(reference, args))
    if args.mode in ("all", "reference-swap"):
        results.append(run_reference_swap(reference, args))
    if args.mode == "reference-double-swap":
        results.append(run_reference_double_swap(reference, args))
    if args.mode in ("all", "ideal-add"):
        results.append(run_ideal_add(ideal, args))

    atomic_write_json(args.output, result_json(args, results))
    survivors = 0
    for result in results:
        counters = result.counters
        survivors += counters.survivors
        print(
            f"{result.mode} examined={counters.examined} "
            f"patterns={counters.exact_pattern_evaluations} "
            f"squares={counters.exact_perfect_square} "
            f"survivors={counters.survivors} "
            f"elapsed={result.elapsed_seconds:.3f}s"
        )
    print(f"finished survivors={survivors} output={args.output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ArithmeticError, OSError, ValueError) as error:
        print(f"gram_edge_search: {error}", file=sys.stderr)
        raise SystemExit(2)
