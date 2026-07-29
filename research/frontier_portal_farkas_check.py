#!/usr/bin/env python3
"""Exactly check sparse linear certificates for frontier shell slices.

No LP/MILP solver or floating-point arithmetic is used.  The checker reads
the complete 1,382-column shell bound to the published frontier Gram and
verifies two eight-term integer linear functionals.  Both functionals
annihilate every shell column outside the size-six orbit.  The target Gram
therefore forces two exact sums on the selected size-six columns.

Given the independently derived exact orbit count of three selected columns,
the missing representative triples (0,1,2), (0,1,3), and (0,2,5) each miss a
forced sum by 16.  This is a direct integer contradiction, not a solver
infeasibility claim.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
ORDER = 23
FRONTIER = 2_779_447_296_000_000
SMALL_ORBIT = (8161, 26521, 8356323, 8357403, 8363525, 8380805)
DEFAULT_SHELL = (
    ROOT / "runs/direct-search/gram-shell-reference-columns-29760.json"
)
DEFAULT_ORBIT_PROOF = (
    ROOT
    / (
        "runs/direct-search/frontier-factor-class-expansion-20260728/"
        "calibration/orbit-proof-v4.json"
    )
)
DEFAULT_OUTPUT = (
    ROOT
    / (
        "runs/direct-search/frontier-portal-harvest-20260729/"
        "exact-farkas-report.json"
    )
)

# Each tuple is (row, column, coefficient), with zero-based row coordinates.
# The factored forms are:
#   L_cross_9_10  = -(s_9+s_10)(s_15+s_16+s_17+s_18)
#   L_cross_3_4   =  (s_3+s_4)(s_15+s_16+s_17+s_18)
FUNCTIONALS = {
    "negative_cross_rows_9_10_to_15_18": tuple(
        (row, column, -1)
        for row in (9, 10)
        for column in (15, 16, 17, 18)
    ),
    "positive_cross_rows_3_4_to_15_18": tuple(
        (row, column, 1)
        for row in (3, 4)
        for column in (15, 16, 17, 18)
    ),
}
REPRESENTATIVE_CERTIFICATES = (
    {
        "triple": (0, 1, 2),
        "functional": "negative_cross_rows_9_10_to_15_18",
    },
    {
        "triple": (0, 1, 3),
        "functional": "positive_cross_rows_3_4_to_15_18",
    },
    {
        "triple": (0, 2, 5),
        "functional": "negative_cross_rows_9_10_to_15_18",
    },
)

Matrix = list[list[int]]
Term = tuple[int, int, int]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def display_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def atomic_write(path: Path, contents: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(contents)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def read_matrix(path: Path) -> Matrix:
    matrix = [
        [int(token) for token in line.split()]
        for line in path.read_text(encoding="ascii").splitlines()
        if line.strip()
    ]
    if len(matrix) != ORDER or any(
        len(row) != ORDER or any(value not in (-1, 1) for value in row)
        for row in matrix
    ):
        raise ValueError(f"{path}: expected exactly {ORDER}x{ORDER} signs")
    return matrix


def gram(matrix: Matrix) -> Matrix:
    return [
        [
            sum(
                matrix[row][index] * matrix[column][index]
                for index in range(ORDER)
            )
            for column in range(ORDER)
        ]
        for row in range(ORDER)
    ]


def determinant(matrix: Matrix) -> int:
    """Exact fraction-free Bareiss determinant."""
    work = [row[:] for row in matrix]
    sign = 1
    denominator = 1
    for pivot_index in range(len(work) - 1):
        pivot_row = next(
            (
                row
                for row in range(pivot_index, len(work))
                if work[row][pivot_index]
            ),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = (
                work[pivot_row],
                work[pivot_index],
            )
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, len(work)):
            for column in range(pivot_index + 1, len(work)):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                if numerator % denominator:
                    raise ArithmeticError("Bareiss division was not exact")
                work[row][column] = numerator // denominator
            work[row][pivot_index] = 0
        denominator = pivot
    return sign * work[-1][-1]


def sign(mask: int, row: int) -> int:
    return 1 if mask >> row & 1 else -1


def functional_on_mask(mask: int, terms: Iterable[Term]) -> int:
    return sum(
        coefficient * sign(mask, row) * sign(mask, column)
        for row, column, coefficient in terms
    )


def functional_on_gram(target: Matrix, terms: Iterable[Term]) -> int:
    return sum(
        coefficient * target[row][column]
        for row, column, coefficient in terms
    )


def shell_inputs(path: Path) -> tuple[list[int], Path, dict[str, Any]]:
    report = json.loads(path.read_text(encoding="utf-8"))
    results = report.get("results")
    if (
        not isinstance(results, list)
        or len(results) != 1
        or not isinstance(results[0], dict)
    ):
        raise ValueError("shell report must contain exactly one result")
    result = results[0]
    masks = result.get("shell_sign_masks")
    if (
        not isinstance(masks, list)
        or len(masks) != 1382
        or len(set(masks)) != 1382
        or not set(SMALL_ORBIT).issubset(masks)
        or any(
            type(mask) is not int
            or mask < 0
            or mask >= 1 << ORDER
            or mask & 1 == 0
            for mask in masks
        )
    ):
        raise ValueError("shell report is not the complete 1,382-column shell")
    source_value = result.get("source")
    if not isinstance(source_value, str):
        raise ValueError("shell report lacks a matrix source")
    source = Path(source_value)
    if not source.is_absolute():
        source = ROOT / source
    source = source.resolve()
    if sha256_file(source) != result.get("source_sha256"):
        raise ValueError("shell source SHA-256 mismatch")
    if result.get("determinant") != str(FRONTIER * FRONTIER):
        raise ValueError("shell Gram determinant is not the frontier square")
    return [int(mask) for mask in masks], source, report


def terms_payload(terms: Iterable[Term]) -> list[dict[str, int]]:
    return [
        {
            "row_zero_based": row,
            "column_zero_based": column,
            "row_one_based": row + 1,
            "column_one_based": column + 1,
            "coefficient": coefficient,
        }
        for row, column, coefficient in terms
    ]


def check_certificate(
    triple: tuple[int, int, int],
    name: str,
    terms: tuple[Term, ...],
    masks: list[int],
    target: Matrix,
) -> dict[str, Any]:
    small = set(SMALL_ORBIT)
    free_masks = [mask for mask in masks if mask not in small]
    free_values = [
        functional_on_mask(mask, terms) for mask in free_masks
    ]
    if any(value != 0 for value in free_values):
        raise ArithmeticError(
            f"{name} does not exactly annihilate the free shell"
        )
    selected_masks = [SMALL_ORBIT[index] for index in triple]
    target_value = functional_on_gram(target, terms)
    fixed_value = sum(
        functional_on_mask(mask, terms) for mask in selected_masks
    )
    residual_value = target_value - fixed_value
    if residual_value == 0:
        raise ArithmeticError(f"{triple}: certificate does not contradict")
    return {
        "triple": list(triple),
        "selected_masks": selected_masks,
        "functional": name,
        "functional_terms": terms_payload(terms),
        "target_gram_value": target_value,
        "fixed_triple_value": fixed_value,
        "target_minus_fixed": residual_value,
        "free_shell_column_count": len(free_masks),
        "free_shell_value_histogram": {
            str(value): count
            for value, count in sorted(Counter(free_values).items())
        },
        "all_free_shell_columns_annihilated_exactly": True,
        "contradiction": (
            f"the free columns contribute exactly 0, but would have to "
            f"contribute {residual_value}"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shell-report", type=Path, default=DEFAULT_SHELL)
    parser.add_argument("--orbit-proof", type=Path, default=DEFAULT_ORBIT_PROOF)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--stdout",
        action="store_true",
        help="print the exact report instead of writing --output",
    )
    arguments = parser.parse_args()

    shell_path = arguments.shell_report.expanduser().resolve()
    orbit_path = arguments.orbit_proof.expanduser().resolve()
    masks, source, _ = shell_inputs(shell_path)
    target_matrix = read_matrix(source)
    target = gram(target_matrix)
    source_determinant = abs(determinant(target_matrix))
    target_determinant = determinant(target)
    if (
        source_determinant != FRONTIER
        or target_determinant != FRONTIER * FRONTIER
    ):
        raise ArithmeticError("source matrix or Gram determinant mismatch")

    orbit_proof = json.loads(orbit_path.read_text(encoding="utf-8"))
    orbit_certificate = orbit_proof.get("orbit_count_exact_certificate")
    if (
        orbit_proof.get("shell_report_sha256") != sha256_file(shell_path)
        or orbit_proof.get("source_matrix_sha256") != sha256_file(source)
        or orbit_proof.get("shell_size") != 1382
        or not isinstance(orbit_certificate, dict)
        or orbit_certificate.get("rank") != 4
        or orbit_certificate.get("shell_orbit_sizes") != [6, 432, 432, 512]
        or orbit_certificate.get("unique_counts") != [3, 6, 6, 8]
    ):
        raise ValueError("orbit proof is not the exact 3/6/6/8 certificate")

    functional_checks: dict[str, dict[str, Any]] = {}
    for name, terms in FUNCTIONALS.items():
        values = [
            functional_on_mask(mask, terms)
            for mask in masks
            if mask not in set(SMALL_ORBIT)
        ]
        if values != [0] * 1376:
            raise ArithmeticError(f"{name}: non-small shell is not a kernel")
        functional_checks[name] = {
            "factored_form_zero_based": (
                "-(s[9]+s[10])*(s[15]+s[16]+s[17]+s[18])"
                if name.startswith("negative")
                else "(s[3]+s[4])*(s[15]+s[16]+s[17]+s[18])"
            ),
            "terms": terms_payload(terms),
            "target_gram_value": functional_on_gram(target, terms),
            "non_small_shell_column_count": len(values),
            "non_small_shell_value_histogram": {
                "0": len(values)
            },
            "exact_kernel_check": True,
        }

    representative_checks = []
    for certificate in REPRESENTATIVE_CERTIFICATES:
        name = str(certificate["functional"])
        representative_checks.append(
            check_certificate(
                certificate["triple"],
                name,
                FUNCTIONALS[name],
                masks,
                target,
            )
        )

    target_values = {
        name: functional_on_gram(target, terms)
        for name, terms in FUNCTIONALS.items()
    }
    triple_table = []
    valid_triples = []
    invalid_triples = []
    for triple in combinations(range(len(SMALL_ORBIT)), 3):
        values = {
            name: sum(
                functional_on_mask(SMALL_ORBIT[index], terms)
                for index in triple
            )
            for name, terms in FUNCTIONALS.items()
        }
        satisfies = values == target_values
        record = {
            "triple": list(triple),
            "functional_values": values,
            "satisfies_both_forced_equalities": satisfies,
        }
        triple_table.append(record)
        (valid_triples if satisfies else invalid_triples).append(
            list(triple)
        )
    expected_valid = [
        [0, 1, 4],
        [0, 2, 4],
        [0, 4, 5],
        [1, 2, 3],
        [1, 3, 5],
        [2, 3, 5],
    ]
    if valid_triples != expected_valid or len(invalid_triples) != 14:
        raise ArithmeticError("unexpected exact small-triple partition")

    script = Path(__file__).resolve()
    report = {
        "schema_version": 1,
        "claim": (
            "Given the complete exported frontier-Gram shell and the exact "
            "forced size-six-orbit count of three, the representative slices "
            "012, 013, and 025 are exactly infeasible. In fact the two sparse "
            "integer equalities reject 14 of the 20 possible triples."
        ),
        "claim_boundary": (
            "This checker verifies the contradiction entirely in Python "
            "integers. Completeness of the exported shell and the independent "
            "3/6/6/8 orbit-count derivation are bound by hashes rather than "
            "re-derived here."
        ),
        "method": (
            "two sparse eight-term integer linear combinations of off-diagonal "
            "Gram equations; every non-small shell column has coefficient zero"
        ),
        "arithmetic": "Python arbitrary-precision integers; no floating point",
        "frontier": str(FRONTIER),
        "source_matrix": {
            "path": display_path(source),
            "sha256": sha256_file(source),
            "exact_absolute_determinant": str(source_determinant),
        },
        "target_gram": {
            "exact_determinant": str(target_determinant),
            "frontier_squared": str(FRONTIER * FRONTIER),
        },
        "shell": {
            "path": display_path(shell_path),
            "sha256": sha256_file(shell_path),
            "size": len(masks),
            "small_orbit_masks": list(SMALL_ORBIT),
            "non_small_column_count": len(masks) - len(SMALL_ORBIT),
        },
        "orbit_count_certificate": {
            "path": display_path(orbit_path),
            "sha256": sha256_file(orbit_path),
            "rank": orbit_certificate["rank"],
            "shell_orbit_sizes": orbit_certificate["shell_orbit_sizes"],
            "unique_counts": orbit_certificate["unique_counts"],
            "forced_small_orbit_count": orbit_certificate["unique_counts"][0],
        },
        "functionals": functional_checks,
        "representative_certificates": representative_checks,
        "all_size_three_small_orbit_subsets": triple_table,
        "valid_triples": valid_triples,
        "invalid_triples": invalid_triples,
        "valid_triple_count": len(valid_triples),
        "invalid_triple_count": len(invalid_triples),
        "all_checks_passed": True,
        "checker": {
            "path": display_path(script),
            "sha256": sha256_file(script),
        },
    }
    serialized = (
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    if arguments.stdout:
        print(serialized.decode("utf-8"), end="")
    else:
        output = arguments.output.expanduser().resolve()
        if output.exists():
            parser.error(f"refusing to overwrite {output}")
        atomic_write(output, serialized)
        print(
            f"exact certificates passed: representatives=3 "
            f"invalid_triples={len(invalid_triples)} "
            f"report={display_path(output)}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
