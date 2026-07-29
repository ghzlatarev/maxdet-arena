#!/usr/bin/env python3
"""Enumerate closed-quadruple switches of a square sign matrix.

A closed quadruple is a set of four rows whose componentwise product is
constant.  Its columns split into four fields, after identifying a four-sign
pattern with its negation.  Negating the four selected entries in any union of
fields gives a switched sign matrix.  Column quadruples are handled by
transposition.

The script computes every determinant with exact Bareiss elimination and can
write the nontrivial switched matrices for use as independent search seeds.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import tempfile


Matrix = list[list[int]]


def read_matrix(path: Path) -> Matrix:
    rows = [
        [int(token) for token in line.split()]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("input must be a nonempty square matrix")
    if any(value not in (-1, 1) for row in rows for value in row):
        raise ValueError("input entries must all be +1 or -1")
    return rows


def determinant(matrix: Matrix) -> int:
    """Return the exact determinant using fraction-free Bareiss elimination."""
    work = [row[:] for row in matrix]
    size = len(work)
    sign = 1
    denominator = 1
    for pivot_index in range(size - 1):
        pivot_row = next(
            (
                row
                for row in range(pivot_index, size)
                if work[row][pivot_index] != 0
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
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
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


def transpose(matrix: Matrix) -> Matrix:
    return [list(column) for column in zip(*matrix)]


def canonical_pattern(pattern: tuple[int, ...]) -> tuple[int, ...]:
    negated = tuple(-value for value in pattern)
    return min(pattern, negated)


def closed_quadruples(matrix: Matrix) -> list[tuple[int, int, int, int]]:
    result: list[tuple[int, int, int, int]] = []
    for quadruple in itertools.combinations(range(len(matrix)), 4):
        products = {
            matrix[quadruple[0]][column]
            * matrix[quadruple[1]][column]
            * matrix[quadruple[2]][column]
            * matrix[quadruple[3]][column]
            for column in range(len(matrix))
        }
        if len(products) == 1:
            result.append(quadruple)
    return result


def field_columns(
    matrix: Matrix, quadruple: tuple[int, int, int, int]
) -> list[list[int]]:
    by_pattern: dict[tuple[int, ...], list[int]] = {}
    for column in range(len(matrix)):
        pattern = tuple(matrix[row][column] for row in quadruple)
        by_pattern.setdefault(canonical_pattern(pattern), []).append(column)
    if len(by_pattern) != 4:
        raise ValueError(
            f"closed quadruple {quadruple} produced {len(by_pattern)} fields, not 4"
        )
    return [by_pattern[key] for key in sorted(by_pattern)]


def switch_fields(
    matrix: Matrix,
    quadruple: tuple[int, int, int, int],
    fields: list[list[int]],
    field_mask: int,
) -> Matrix:
    switched = [row[:] for row in matrix]
    for field_index, columns in enumerate(fields):
        if not field_mask & (1 << field_index):
            continue
        for row in quadruple:
            for column in columns:
                switched[row][column] *= -1
    return switched


def matrix_text(matrix: Matrix) -> str:
    return "\n".join(" ".join(str(value) for value in row) for row in matrix) + "\n"


def matrix_sha256(matrix: Matrix) -> str:
    return hashlib.sha256(matrix_text(matrix).encode("utf-8")).hexdigest()


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent, text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            output.write(content)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def enumerate_orientation(
    matrix: Matrix,
    orientation: str,
    output_dir: Path | None,
) -> list[dict[str, object]]:
    working = matrix if orientation == "rows" else transpose(matrix)
    records: list[dict[str, object]] = []
    for quadruple in closed_quadruples(working):
        fields = field_columns(working, quadruple)
        for field_mask in range(1, 1 << len(fields)):
            switched_working = switch_fields(working, quadruple, fields, field_mask)
            switched = (
                switched_working
                if orientation == "rows"
                else transpose(switched_working)
            )
            record: dict[str, object] = {
                "orientation": orientation,
                "quadruple_zero_based": list(quadruple),
                "field_sizes": [len(field) for field in fields],
                "field_mask": field_mask,
                "absolute_determinant": abs(determinant(switched)),
                "matrix_sha256": matrix_sha256(switched),
            }
            if output_dir is not None:
                stem = "-".join(str(index) for index in quadruple)
                output_path = (
                    output_dir
                    / f"{orientation}-{stem}-mask-{field_mask:02x}.matrix.txt"
                )
                atomic_write(output_path, matrix_text(switched))
                record["output"] = str(output_path)
            records.append(record)
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="write every nontrivial switched matrix to this directory",
    )
    args = parser.parse_args()

    matrix = read_matrix(args.matrix)
    original_determinant = abs(determinant(matrix))
    records = [
        *enumerate_orientation(matrix, "rows", args.output_dir),
        *enumerate_orientation(matrix, "columns", args.output_dir),
    ]
    summary = {
        "schema_version": 1,
        "input": str(args.matrix),
        "order": len(matrix),
        "input_sha256": matrix_sha256(matrix),
        "input_absolute_determinant": original_determinant,
        "closed_row_quadruples": len(closed_quadruples(matrix)),
        "closed_column_quadruples": len(closed_quadruples(transpose(matrix))),
        "switched_variants": len(records),
        "same_score_variants": sum(
            record["absolute_determinant"] == original_determinant
            for record in records
        ),
        "best_absolute_determinant": max(
            [original_determinant]
            + [int(record["absolute_determinant"]) for record in records]
        ),
        "records": records,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
