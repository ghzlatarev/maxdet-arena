#!/usr/bin/env python3
"""Reproduce the second published order-22 D-optimal design.

This is Example 4 of Lin and Phoa, "Constructing near-Hadamard designs
with (almost) D-optimality by General Supplementary Difference Sets",
Statistica Sinica 26 (2016), 413-427, DOI 10.5705/ss.202014.0115.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import tempfile
from pathlib import Path


ORDER = 22
BLOCK_ORDER = 10
EXPECTED_DETERMINANT = -409_600_000_000_000
EXPECTED_SHA256 = (
    "68d673e47304b076acfc4361aee752efa68e203651132eee2860591c15dd34f8"
)

A_FIRST = (-1, -1, 1, -1, 1, -1, 1, -1, 1, 1)
B_FIRST = (-1, -1, -1, -1, 1, 1, -1, 1, 1, 1)


def circulant(first: tuple[int, ...]) -> list[list[int]]:
    """Return the circulant whose successive rows are right shifts."""

    if len(first) != BLOCK_ORDER or any(value not in (-1, 1) for value in first):
        raise ValueError("invalid first row")
    return [
        [first[(column - row) % BLOCK_ORDER] for column in range(BLOCK_ORDER)]
        for row in range(BLOCK_ORDER)
    ]


def transpose(matrix: list[list[int]]) -> list[list[int]]:
    return [list(column) for column in zip(*matrix)]


def build_matrix() -> list[list[int]]:
    ones = [1] * BLOCK_ORDER
    minus_ones = [-1] * BLOCK_ORDER
    a = circulant(A_FIRST)
    b = circulant(B_FIRST)
    at = transpose(a)
    bt = transpose(b)

    matrix = [
        [1, 1, *ones, *ones],
        [1, -1, *ones, *minus_ones],
    ]
    matrix.extend([1, 1, *a[row], *b[row]] for row in range(BLOCK_ORDER))
    matrix.extend(
        [1, -1, *bt[row], *[-value for value in at[row]]]
        for row in range(BLOCK_ORDER)
    )
    if len(matrix) != ORDER or any(
        len(row) != ORDER or any(value not in (-1, 1) for value in row)
        for row in matrix
    ):
        raise ArithmeticError("constructed matrix is not a 22x22 sign matrix")
    return matrix


def bareiss_determinant(source: list[list[int]]) -> int:
    work = [row[:] for row in source]
    previous = 1
    sign = 1
    for column in range(len(work) - 1):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[pivot], work[column] = work[column], work[pivot]
            sign = -sign
        pivot_value = work[column][column]
        for row in range(column + 1, len(work)):
            for inner in range(column + 1, len(work)):
                numerator = (
                    work[row][inner] * pivot_value
                    - work[row][column] * work[column][inner]
                )
                if column and numerator % previous:
                    raise ArithmeticError("non-exact Bareiss division")
                work[row][inner] = numerator if not column else numerator // previous
            work[row][column] = 0
        previous = pivot_value
    return sign * work[-1][-1]


def payload(matrix: list[list[int]]) -> bytes:
    return "".join(" ".join(map(str, row)) + "\n" for row in matrix).encode("ascii")


def atomic_write(path: Path, contents: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not path.is_file():
        raise ValueError("output path is not a regular file")
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(contents)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()

    matrix = build_matrix()
    determinant = bareiss_determinant(matrix)
    contents = payload(matrix)
    digest = hashlib.sha256(contents).hexdigest()
    if determinant != EXPECTED_DETERMINANT:
        raise ArithmeticError(
            f"determinant changed: {determinant} != {EXPECTED_DETERMINANT}"
        )
    if digest != EXPECTED_SHA256:
        raise ArithmeticError(f"payload SHA-256 changed: {digest}")
    atomic_write(arguments.output.expanduser().absolute(), contents)
    print(
        f"wrote order-{ORDER} GSDS matrix; "
        f"det={determinant} sha256={digest}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
