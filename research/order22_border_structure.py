#!/usr/bin/env python3
"""Exact structural audit for order-23 matrices built from order-22 borders.

The last row and column are treated as the border.  The report binds the raw
matrix, checks the exact Schur/adjugate determinant identity, fingerprints the
order-22 core and order-23 Gram matrices, and records the complete one-entry
core-determinant sensitivity profile.  No floating-point arithmetic is used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence


Matrix = list[list[int]]


def parse_labeled_path(text: str) -> tuple[str, Path]:
    if "=" not in text:
        raise argparse.ArgumentTypeError("matrix must be LABEL=PATH")
    label, raw_path = text.split("=", 1)
    if not label or not raw_path:
        raise argparse.ArgumentTypeError("matrix must be LABEL=PATH")
    return label, Path(raw_path)


def read_sign_matrix(path: Path, order: int) -> tuple[bytes, Matrix]:
    payload = path.read_bytes()
    try:
        rows = [
            [int(token) for token in line.split()]
            for line in payload.decode("ascii").splitlines()
            if line.strip()
        ]
    except (UnicodeDecodeError, ValueError) as error:
        raise ValueError(f"{path}: not an ASCII integer matrix") from error
    if len(rows) != order or any(len(row) != order for row in rows):
        raise ValueError(f"{path}: expected exactly {order} by {order}")
    if any(value not in (-1, 1) for row in rows for value in row):
        raise ValueError(f"{path}: entries must be -1 or 1")
    return payload, rows


def bareiss_determinant(matrix: Sequence[Sequence[int]]) -> int:
    order = len(matrix)
    if order == 0:
        return 1
    work = [list(map(int, row)) for row in matrix]
    previous_pivot = 1
    sign = 1
    for column in range(order - 1):
        pivot_row = next(
            (row for row in range(column, order) if work[row][column]),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[pivot_row], work[column] = work[column], work[pivot_row]
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
                        raise ArithmeticError("non-exact Bareiss division")
                    work[row][inner] = quotient
                else:
                    work[row][inner] = numerator
            work[row][column] = 0
        previous_pivot = pivot
    return sign * work[-1][-1]


def exact_adjugate(matrix: Matrix, determinant: int) -> Matrix:
    order = len(matrix)
    augmented: list[list[Fraction]] = [
        [Fraction(value) for value in matrix[row]]
        + [Fraction(int(row == column)) for column in range(order)]
        for row in range(order)
    ]
    for column in range(order):
        pivot_row = next(
            (
                row
                for row in range(column, order)
                if augmented[row][column]
            ),
            None,
        )
        if pivot_row is None:
            raise ArithmeticError("singular core has no adjugate inverse")
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
                value - factor * pivot_value
                for value, pivot_value in zip(
                    augmented[row], augmented[column]
                )
            ]
    adjugate: Matrix = []
    for row in range(order):
        integer_row: list[int] = []
        for value in augmented[row][order:]:
            scaled = value * determinant
            if scaled.denominator != 1:
                raise ArithmeticError("inverse did not yield integer adjugate")
            integer_row.append(scaled.numerator)
        adjugate.append(integer_row)
    expected = [
        [determinant if row == column else 0 for column in range(order)]
        for row in range(order)
    ]
    if multiply(matrix, adjugate) != expected:
        raise ArithmeticError("core times adjugate identity failed")
    return adjugate


def multiply(left: Matrix, right: Matrix) -> Matrix:
    return [
        [
            sum(left[row][inner] * right[inner][column]
                for inner in range(len(right)))
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def matrix_vector(matrix: Matrix, vector: Sequence[int]) -> list[int]:
    return [
        sum(value * vector[column] for column, value in enumerate(row))
        for row in matrix
    ]


def gram(matrix: Matrix) -> Matrix:
    return [
        [
            sum(
                matrix[row][column] * matrix[other][column]
                for column in range(len(matrix[0]))
            )
            for other in range(len(matrix))
        ]
        for row in range(len(matrix))
    ]


def transpose(matrix: Matrix) -> Matrix:
    return [list(column) for column in zip(*matrix)]


def canonical_json_sha256(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def histogram(values: Iterable[int]) -> dict[str, int]:
    counts = Counter(values)
    return {str(value): counts[value] for value in sorted(counts)}


def off_diagonal(matrix: Matrix) -> list[int]:
    return [
        matrix[row][column]
        for row in range(len(matrix))
        for column in range(row + 1, len(matrix))
    ]


def gram_summary(matrix: Matrix) -> dict[str, object]:
    exact = gram(matrix)
    off = off_diagonal(exact)
    levels = sorted(set(abs(value) for value in off))
    degree_profiles = {
        str(level): sorted(
            sum(
                abs(exact[row][column]) == level
                for column in range(len(exact))
                if column != row
            )
            for row in range(len(exact))
        )
        for level in levels
    }
    return {
        "sha256": canonical_json_sha256(exact),
        "off_diagonal_histogram": histogram(off),
        "off_diagonal_absolute_histogram": histogram(map(abs, off)),
        "absolute_level_degree_profiles": degree_profiles,
    }


def matrix_summary(matrix: Matrix) -> dict[str, object]:
    return {
        "row_gram": gram_summary(matrix),
        "column_gram": gram_summary(transpose(matrix)),
    }


def nonzero_gcd(values: Iterable[int]) -> int:
    result = 0
    for value in values:
        result = math.gcd(result, abs(value))
    return result


def adjugate_summary(adjugate: Matrix) -> dict[str, object]:
    flat = [value for row in adjugate for value in row]
    quantum = nonzero_gcd(flat)
    return {
        "sha256": canonical_json_sha256(adjugate),
        "entry_gcd": quantum,
        "zero_entries": sum(value == 0 for value in flat),
        "absolute_entry_histogram": histogram(map(abs, flat)),
        "gcd_scaled_absolute_entry_histogram": histogram(
            abs(value) // quantum for value in flat
        ),
        "row_l1_histogram": histogram(
            sum(abs(value) for value in row) for row in adjugate
        ),
        "column_l1_histogram": histogram(
            sum(abs(adjugate[row][column]) for row in range(len(adjugate)))
            for column in range(len(adjugate))
        ),
    }


def core_flip_profile(
    core: Matrix,
    determinant: int,
    adjugate: Matrix,
) -> dict[str, object]:
    records = []
    for row in range(len(core)):
        for column in range(len(core)):
            flipped_determinant = (
                determinant
                - 2 * core[row][column] * adjugate[column][row]
            )
            records.append(
                {
                    "row_1_based": row + 1,
                    "column_1_based": column + 1,
                    "flipped_determinant": flipped_determinant,
                    "absolute_flipped_determinant": abs(flipped_determinant),
                    "cofactor": adjugate[column][row],
                }
            )
    ranked = sorted(
        records,
        key=lambda record: (
            -int(record["absolute_flipped_determinant"]),
            int(record["row_1_based"]),
            int(record["column_1_based"]),
        ),
    )
    neutral = [
        record
        for record in records
        if record["absolute_flipped_determinant"] == abs(determinant)
    ]
    return {
        "absolute_determinant_histogram": histogram(
            int(record["absolute_flipped_determinant"])
            for record in records
        ),
        "maximal_determinant_flip_count": len(neutral),
        "maximal_determinant_flips": neutral,
        "top_32_by_absolute_core_determinant": ranked[:32],
    }


def sibling_border_log(path: Path) -> dict[str, object] | None:
    candidates = [path.with_name("run.jsonl")]
    if path.name == "best.matrix.txt":
        candidates.extend(
            [path.with_name("screen.jsonl"), path.with_name("screen-rerun.jsonl")]
        )
    for candidate in candidates:
        if not candidate.is_file():
            continue
        finished = None
        for line in candidate.read_text(encoding="utf-8").splitlines():
            record = json.loads(line)
            if record.get("event") == "finished":
                finished = record
        if finished is not None:
            return {
                "path": str(candidate),
                "raw_sha256": hashlib.sha256(candidate.read_bytes()).hexdigest(),
                "assignments_completed": finished.get("assignments_completed"),
                "best_border_columns_up_to_global_sign": finished.get(
                    "best_border_columns_up_to_global_sign"
                ),
                "complete": finished.get("complete"),
            }
    return None


def audit_matrix(label: str, path: Path) -> dict[str, object]:
    resolved = path.resolve()
    payload, matrix = read_sign_matrix(resolved, 23)
    determinant = bareiss_determinant(matrix)
    core = [row[:22] for row in matrix[:22]]
    core_determinant = bareiss_determinant(core)
    if not core_determinant:
        raise ArithmeticError(f"{path}: singular last-border core")
    adjugate = exact_adjugate(core, core_determinant)
    border_column = [matrix[row][22] for row in range(22)]
    border_row = matrix[22][:22]
    corner = matrix[22][22]
    adjugate_column = matrix_vector(adjugate, border_column)
    schur_determinant = (
        corner * core_determinant
        - sum(
            border_row[index] * adjugate_column[index]
            for index in range(22)
        )
    )
    if schur_determinant != determinant:
        raise ArithmeticError(f"{path}: Schur determinant identity failed")
    fixed_column_optimum = abs(core_determinant) + sum(
        abs(value) for value in adjugate_column
    )
    if abs(determinant) != fixed_column_optimum:
        raise ArithmeticError(f"{path}: stored border row is not optimal")
    border_quantum = nonzero_gcd(adjugate_column)
    return {
        "label": label,
        "path": str(resolved),
        "raw_sha256": hashlib.sha256(payload).hexdigest(),
        "determinant": determinant,
        "absolute_determinant": abs(determinant),
        "matrix_structure": matrix_summary(matrix),
        "core": {
            "determinant": core_determinant,
            "absolute_determinant": abs(core_determinant),
            "raw_sha256": hashlib.sha256(
                "".join(
                    " ".join(map(str, row)) + "\n"
                    for row in core
                ).encode("ascii")
            ).hexdigest(),
            "matrix_structure": matrix_summary(core),
            "adjugate": adjugate_summary(adjugate),
            "one_entry_flip_profile": core_flip_profile(
                core, core_determinant, adjugate
            ),
        },
        "border": {
            "corner": corner,
            "column": border_column,
            "row": border_row,
            "adjugate_times_column": adjugate_column,
            "adjugate_times_column_gcd": border_quantum,
            "gcd_scaled_absolute_histogram": histogram(
                abs(value) // border_quantum for value in adjugate_column
            ),
            "l1_norm": sum(abs(value) for value in adjugate_column),
            "fixed_column_exact_optimum": fixed_column_optimum,
            "stored_border_is_exactly_optimal": True,
            "sibling_exhaustive_log": sibling_border_log(resolved),
        },
    }


def comparison_groups(
    records: Sequence[dict[str, object]],
) -> dict[str, object]:
    selectors = {
        "order23_row_gram": ("matrix_structure", "row_gram", "sha256"),
        "order23_column_gram": (
            "matrix_structure",
            "column_gram",
            "sha256",
        ),
        "order22_row_gram": (
            "core",
            "matrix_structure",
            "row_gram",
            "sha256",
        ),
        "order22_column_gram": (
            "core",
            "matrix_structure",
            "column_gram",
            "sha256",
        ),
        "order22_adjugate": ("core", "adjugate", "sha256"),
    }
    result: dict[str, object] = {}
    for name, selector in selectors.items():
        groups: dict[str, list[str]] = {}
        for record in records:
            value: object = record
            for key in selector:
                value = value[key]  # type: ignore[index]
            groups.setdefault(str(value), []).append(str(record["label"]))
        result[name] = [
            {"sha256": fingerprint, "labels": sorted(labels)}
            for fingerprint, labels in sorted(groups.items())
        ]
    return result


def atomic_write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix",
        action="append",
        required=True,
        type=parse_labeled_path,
        metavar="LABEL=PATH",
    )
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    labels = [label for label, _ in arguments.matrix]
    if len(labels) != len(set(labels)):
        parser.error("matrix labels must be unique")
    records = [
        audit_matrix(label, path) for label, path in arguments.matrix
    ]
    report = {
        "engine": "order22-border-structure-v1",
        "arithmetic": "exact-integer-and-rational-only",
        "last_row_and_column_are_border": True,
        "records": records,
        "comparison_groups": comparison_groups(records),
    }
    if arguments.output:
        atomic_write_json(arguments.output, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
