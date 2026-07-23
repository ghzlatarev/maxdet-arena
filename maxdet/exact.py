"""Dependency-free exact integer linear algebra.

The trusted path intentionally avoids NumPy, BLAS, floating point, and symbolic
math packages. Candidate ranking is based only on Python arbitrary-precision
integers.
"""

from __future__ import annotations

from collections.abc import Sequence

from .errors import VerificationError

Matrix = Sequence[Sequence[int]]


def _copy_square(matrix: Matrix) -> list[list[int]]:
    rows = [list(row) for row in matrix]
    size = len(rows)
    if size == 0:
        raise ValueError("matrix must not be empty")
    if any(len(row) != size for row in rows):
        raise ValueError("matrix must be square")
    return rows


def bareiss_determinant(matrix: Matrix) -> int:
    """Return det(matrix) exactly with fraction-free Bareiss elimination.

    Row pivoting is deterministic: the first non-zero row is selected.
    """

    work = _copy_square(matrix)
    size = len(work)
    if size == 1:
        return work[0][0]

    sign = 1
    previous_pivot = 1

    for column in range(size - 1):
        pivot_row = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign

        pivot = work[column][column]
        for row in range(column + 1, size):
            for inner_column in range(column + 1, size):
                numerator = (
                    work[row][inner_column] * pivot
                    - work[row][column] * work[column][inner_column]
                )
                if column:
                    quotient, remainder = divmod(numerator, previous_pivot)
                    if remainder:
                        raise VerificationError(
                            "Bareiss exact-division invariant failed"
                        )
                    numerator = quotient
                work[row][inner_column] = numerator
            work[row][column] = 0

        previous_pivot = pivot

    return sign * work[-1][-1]


def modular_determinant(matrix: Matrix, prime: int) -> int:
    """Return det(matrix) modulo a prime using modular elimination."""

    if prime <= 2:
        raise ValueError("prime must be greater than two")
    work = _copy_square(matrix)
    size = len(work)
    for row in work:
        for column, value in enumerate(row):
            row[column] = value % prime

    determinant = 1
    for column in range(size):
        pivot_row = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            determinant = -determinant

        pivot = work[column][column]
        determinant = (determinant * pivot) % prime
        inverse = pow(pivot, prime - 2, prime)
        for row in range(column + 1, size):
            if work[row][column] == 0:
                continue
            factor = work[row][column] * inverse % prime
            for inner_column in range(column, size):
                work[row][inner_column] = (
                    work[row][inner_column]
                    - factor * work[column][inner_column]
                ) % prime

    return determinant % prime


def gram_matrix(matrix: Matrix) -> list[list[int]]:
    """Return A A^T using integer arithmetic."""

    work = _copy_square(matrix)
    size = len(work)
    return [
        [
            sum(work[row][k] * work[column][k] for k in range(size))
            for column in range(size)
        ]
        for row in range(size)
    ]


def normalize_signs(matrix: Matrix) -> list[list[int]]:
    """Normalize sign-equivalent matrices to positive first row and column.

    This deliberately does not claim canonicalization under row or column
    permutations. The resulting hash is a sign-normalized identity, not a full
    Hadamard-equivalence class identifier.
    """

    work = _copy_square(matrix)
    size = len(work)

    for column in range(size):
        if work[0][column] == -1:
            for row in range(size):
                work[row][column] *= -1

    for row in range(size):
        if work[row][0] == -1:
            for column in range(size):
                work[row][column] *= -1

    return work


def matrix_text(matrix: Matrix) -> str:
    """Serialize a matrix in the only accepted arena text format."""

    work = _copy_square(matrix)
    return "".join(" ".join(str(value) for value in row) + "\n" for row in work)
