#!/usr/bin/env python3
"""Recover and exhaust a six-generator neutral-switch compound cube.

The input is a five-node prefix

    base --A0--> a0
      |
      B0
      v
      b0 --A1--> a1 --B1--> b1

where each transition flips twelve entries on one of two alternating
9-by-9 support families.  Within a family, A0/A1 (or B0/B1) must be
disjoint cyclic assignments of the same degree pattern
(2,2,2,1,1,1,1,1,1) on both sides.  The third assignment is the exact
complement in the two allowed arm rectangles.  The program derives A2
and B2, exhausts all 2^6 compound states with exact Bareiss
determinants, and atomically archives every state tied with the base.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile


Matrix = tuple[tuple[int, ...], ...]
Coordinate = tuple[int, int]


def read_matrix(path: Path) -> tuple[Matrix, bytes]:
    raw = path.read_bytes()
    rows = tuple(
        tuple(int(token) for token in line.split())
        for line in raw.decode("utf-8").splitlines()
        if line.strip()
    )
    if len(rows) != 23 or any(len(row) != 23 for row in rows):
        raise ValueError(f"{path} must contain exactly a 23x23 matrix")
    if any(value not in (-1, 1) for row in rows for value in row):
        raise ValueError(f"{path} entries must all be +1 or -1")
    return rows, raw


def matrix_text(matrix: Matrix) -> str:
    return (
        "\n".join(
            " ".join(str(value) for value in row)
            for row in matrix
        )
        + "\n"
    )


def atomic_write(path: Path, content: str) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent, text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            output.write(content)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary_name, path)
        directory = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def determinant(matrix: Matrix) -> int:
    work = [list(row) for row in matrix]
    order = len(work)
    sign = 1
    denominator = 1
    for pivot_index in range(order - 1):
        pivot_row = next(
            (
                row
                for row in range(pivot_index, order)
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
        for row in range(pivot_index + 1, order):
            for column in range(pivot_index + 1, order):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index]
                    * work[pivot_index][column]
                )
                if numerator % denominator:
                    raise ArithmeticError("Bareiss division was not exact")
                work[row][column] = numerator // denominator
            work[row][pivot_index] = 0
        denominator = pivot
    return sign * work[-1][-1]


def difference(first: Matrix, second: Matrix) -> frozenset[Coordinate]:
    return frozenset(
        (row, column)
        for row in range(23)
        for column in range(23)
        if first[row][column] != second[row][column]
    )


def apply_mask(base: Matrix, mask: frozenset[Coordinate]) -> Matrix:
    return tuple(
        tuple(
            -base[row][column]
            if (row, column) in mask
            else base[row][column]
            for column in range(23)
        )
        for row in range(23)
    )


def xor_masks(
    masks: list[frozenset[Coordinate]],
) -> frozenset[Coordinate]:
    result: set[Coordinate] = set()
    for mask in masks:
        result.symmetric_difference_update(mask)
    return frozenset(result)


def degree_profile(
    mask: frozenset[Coordinate],
) -> tuple[dict[int, int], dict[int, int]]:
    row_degrees: dict[int, int] = {}
    column_degrees: dict[int, int] = {}
    for row, column in mask:
        row_degrees[row] = row_degrees.get(row, 0) + 1
        column_degrees[column] = column_degrees.get(column, 0) + 1
    return row_degrees, column_degrees


def derive_third_assignment(
    name: str,
    first: frozenset[Coordinate],
    second: frozenset[Coordinate],
) -> tuple[frozenset[Coordinate], dict[str, object]]:
    if len(first) != 12 or len(second) != 12:
        raise ValueError(f"{name} transitions must each have 12 flips")
    if first & second:
        raise ValueError(f"{name} first two assignments must be disjoint")

    first_rows, first_columns = degree_profile(first)
    second_rows, second_columns = degree_profile(second)
    expected_degrees = [1] * 6 + [2] * 3
    if (
        sorted(first_rows.values()) != expected_degrees
        or sorted(first_columns.values()) != expected_degrees
        or first_rows != second_rows
        or first_columns != second_columns
    ):
        raise ValueError(
            f"{name} assignments must share the "
            "(2x3,1x6) row/column degree profile"
        )

    double_rows = sorted(
        row for row, degree in first_rows.items() if degree == 2
    )
    single_rows = sorted(
        row for row, degree in first_rows.items() if degree == 1
    )
    double_columns = sorted(
        column
        for column, degree in first_columns.items()
        if degree == 2
    )
    single_columns = sorted(
        column
        for column, degree in first_columns.items()
        if degree == 1
    )
    allowed = frozenset(
        (row, column)
        for row in double_rows
        for column in single_columns
    ) | frozenset(
        (row, column)
        for row in single_rows
        for column in double_columns
    )
    if not first <= allowed or not second <= allowed:
        raise ValueError(
            f"{name} assignments contain edges outside the arm rectangles"
        )
    third = allowed - first - second
    third_rows, third_columns = degree_profile(third)
    if (
        len(third) != 12
        or third_rows != first_rows
        or third_columns != first_columns
    ):
        raise ValueError(
            f"{name} arm complement is not a third valid assignment"
        )
    metadata: dict[str, object] = {
        "double_rows_one_based": [value + 1 for value in double_rows],
        "single_rows_one_based": [value + 1 for value in single_rows],
        "double_columns_one_based": [
            value + 1 for value in double_columns
        ],
        "single_columns_one_based": [
            value + 1 for value in single_columns
        ],
    }
    return third, metadata


def mask_record(mask: frozenset[Coordinate]) -> list[list[int]]:
    return [
        [row + 1, column + 1]
        for row, column in sorted(mask)
    ]


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--a0-endpoint", type=Path, required=True)
    parser.add_argument("--b0-endpoint", type=Path, required=True)
    parser.add_argument("--a1-endpoint", type=Path, required=True)
    parser.add_argument("--b1-endpoint", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    inputs: dict[str, tuple[Matrix, bytes, Path]] = {}
    for name, path in (
        ("base", args.base),
        ("a0_endpoint", args.a0_endpoint),
        ("b0_endpoint", args.b0_endpoint),
        ("a1_endpoint", args.a1_endpoint),
        ("b1_endpoint", args.b1_endpoint),
    ):
        matrix, raw = read_matrix(path)
        inputs[name] = (matrix, raw, path)

    base = inputs["base"][0]
    a0 = difference(base, inputs["a0_endpoint"][0])
    b0 = difference(base, inputs["b0_endpoint"][0])
    a1 = difference(
        inputs["b0_endpoint"][0], inputs["a1_endpoint"][0]
    )
    b1 = difference(
        inputs["a1_endpoint"][0], inputs["b1_endpoint"][0]
    )
    a2, a_metadata = derive_third_assignment("A", a0, a1)
    b2, b_metadata = derive_third_assignment("B", b0, b1)
    generators = [a0, a1, a2, b0, b1, b2]
    generator_names = ["A0", "A1", "A2", "B0", "B1", "B2"]

    for left in range(len(generators)):
        for right in range(left):
            if generators[left] & generators[right]:
                raise ValueError(
                    f"{generator_names[left]} and "
                    f"{generator_names[right]} are not disjoint"
                )

    expected_prefix = {
        "a0_endpoint": 0b000001,
        "b0_endpoint": 0b001000,
        "a1_endpoint": 0b001010,
        "b1_endpoint": 0b011010,
    }
    for name, state in expected_prefix.items():
        mask = xor_masks(
            [
                generators[index]
                for index in range(6)
                if state & (1 << index)
            ]
        )
        if apply_mask(base, mask) != inputs[name][0]:
            raise ValueError(f"{name} does not match the inferred prefix")

    if args.output_dir.exists():
        raise FileExistsError(
            f"fresh output directory already exists: {args.output_dir}"
        )
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir()

    base_score = abs(determinant(base))
    states: list[dict[str, object]] = []
    matrices: dict[int, Matrix] = {}
    scores: dict[int, int] = {}
    for state in range(1 << 6):
        selected = [
            generators[index]
            for index in range(6)
            if state & (1 << index)
        ]
        matrix = apply_mask(base, xor_masks(selected))
        score = abs(determinant(matrix))
        matrices[state] = matrix
        scores[state] = score
        states.append(
            {
                "state_hex": f"{state:02x}",
                "generators": [
                    generator_names[index]
                    for index in range(6)
                    if state & (1 << index)
                ],
                "hamming_from_base": 12 * len(selected),
                "absolute_determinant": str(score),
            }
        )

    ties = sorted(
        state for state, score in scores.items() if score == base_score
    )
    promotions = sorted(
        state for state, score in scores.items() if score > base_score
    )
    tie_edges: list[dict[str, str]] = []
    tie_set = set(ties)
    for state in ties:
        for index, name in enumerate(generator_names):
            neighbor = state ^ (1 << index)
            if neighbor in tie_set and state < neighbor:
                tie_edges.append(
                    {
                        "first_state_hex": f"{state:02x}",
                        "generator": name,
                        "second_state_hex": f"{neighbor:02x}",
                    }
                )

    artifacts: list[dict[str, str]] = []
    for state in ties:
        text = matrix_text(matrices[state])
        filename = f"tie-{state:02x}.matrix.txt"
        atomic_write(args.output_dir / filename, text)
        artifacts.append(
            {
                "state_hex": f"{state:02x}",
                "path": filename,
                "raw_sha256": sha256(text.encode("utf-8")),
            }
        )

    report = {
        "schema_version": 1,
        "event": "finished",
        "complete": True,
        "method": "six-generator-neutral-compound-cube-v1",
        "base_absolute_determinant": str(base_score),
        "best_absolute_determinant": str(max(scores.values())),
        "compound_states": len(states),
        "tie_states": len(ties),
        "promotion_states": len(promotions),
        "inputs": {
            name: {
                "path": str(path),
                "raw_sha256": sha256(raw),
            }
            for name, (_, raw, path) in inputs.items()
        },
        "families": {
            "A": a_metadata,
            "B": b_metadata,
        },
        "generators": {
            name: mask_record(mask)
            for name, mask in zip(generator_names, generators)
        },
        "states": states,
        "tie_edges": tie_edges,
        "tie_artifacts": artifacts,
        "promotion_state_hexes": [
            f"{state:02x}" for state in promotions
        ],
    }
    atomic_write(
        args.output_dir / "report.json",
        json.dumps(report, indent=2, sort_keys=True) + "\n",
    )
    print(
        json.dumps(
            {
                "best_absolute_determinant": str(max(scores.values())),
                "compound_states": len(states),
                "promotion_states": len(promotions),
                "tie_edges": len(tie_edges),
                "tie_states": len(ties),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
