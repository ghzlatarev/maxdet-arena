#!/usr/bin/env python3
"""Align a sign-matrix elite to a target under Hadamard monomial actions.

Alternating Hungarian assignments optimize signed row and column permutations
to minimize Hamming distance.  This is a heuristic alignment, not an
equivalence test: a nonzero final distance says nothing about inequivalence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np

try:
    from scipy.optimize import linear_sum_assignment
except ImportError as error:  # pragma: no cover - dependency diagnostic
    raise SystemExit("align_elite.py requires scipy") from error


def read_matrix(path: Path) -> tuple[np.ndarray, bytes]:
    raw = path.read_bytes()
    rows = [
        [int(token) for token in line.split()]
        for line in raw.decode("utf-8").splitlines()
        if line.strip()
    ]
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError(f"{path} must contain a nonempty square matrix")
    if any(value not in (-1, 1) for row in rows for value in row):
        raise ValueError(f"{path} entries must all be +1 or -1")
    return np.asarray(rows, dtype=np.int8), raw


def determinant(matrix: np.ndarray) -> int:
    work = [[int(value) for value in row] for row in matrix.tolist()]
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


def matrix_text(matrix: np.ndarray) -> str:
    return (
        "\n".join(
            " ".join(str(int(value)) for value in row)
            for row in matrix.tolist()
        )
        + "\n"
    )


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


def signed_assignment(correlations: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    targets, sources = linear_sum_assignment(-np.abs(correlations))
    permutation = sources[np.argsort(targets)]
    diagonal = correlations[np.arange(correlations.shape[0]), permutation]
    signs = np.where(diagonal < 0, -1, 1).astype(np.int8)
    return permutation, signs


def materialize(
    source: np.ndarray,
    row_permutation: np.ndarray,
    row_signs: np.ndarray,
    column_permutation: np.ndarray,
    column_signs: np.ndarray,
) -> np.ndarray:
    return (
        row_signs[:, None]
        * column_signs[None, :]
        * source[row_permutation][:, column_permutation]
    )


def align(
    target: np.ndarray,
    source: np.ndarray,
    seed: int,
    restarts: int,
    maximum_iterations: int,
) -> tuple[np.ndarray, dict[str, object]]:
    randomizer = np.random.default_rng(seed)
    best_distance = target.size + 1
    best: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None = None
    improvements = 0
    target_wide = target.astype(np.int16)

    for restart in range(restarts):
        if restart == 0:
            column_permutation = np.arange(target.shape[1])
            column_signs = np.ones(target.shape[1], dtype=np.int8)
        else:
            column_permutation = randomizer.permutation(target.shape[1])
            column_signs = randomizer.choice(
                np.asarray([-1, 1], dtype=np.int8), target.shape[1]
            )

        previous_distance = -1
        for _ in range(maximum_iterations):
            columns_aligned = (
                source[:, column_permutation] * column_signs[None, :]
            )
            row_correlations = target_wide @ columns_aligned.astype(np.int16).T
            row_permutation, row_signs = signed_assignment(row_correlations)

            rows_aligned = source[row_permutation] * row_signs[:, None]
            column_correlations = target_wide.T @ rows_aligned.astype(np.int16)
            column_permutation, column_signs = signed_assignment(
                column_correlations
            )

            candidate = materialize(
                source,
                row_permutation,
                row_signs,
                column_permutation,
                column_signs,
            )
            distance = int(np.count_nonzero(candidate != target))
            if distance < best_distance:
                best_distance = distance
                best = (
                    row_permutation.copy(),
                    row_signs.copy(),
                    column_permutation.copy(),
                    column_signs.copy(),
                )
                improvements += 1
            if distance == previous_distance:
                break
            previous_distance = distance

    if best is None:
        raise RuntimeError("alignment produced no candidate")
    row_permutation, row_signs, column_permutation, column_signs = best
    result = materialize(
        source,
        row_permutation,
        row_signs,
        column_permutation,
        column_signs,
    )
    metadata: dict[str, object] = {
        "hamming_distance": best_distance,
        "improvements": improvements,
        "row_permutation_zero_based": row_permutation.tolist(),
        "row_signs": row_signs.tolist(),
        "column_permutation_zero_based": column_permutation.tolist(),
        "column_signs": column_signs.tolist(),
    }
    return result, metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument("--restarts", type=int, default=10_000)
    parser.add_argument("--max-iterations", type=int, default=30)
    args = parser.parse_args()

    if args.restarts <= 0 or args.max_iterations <= 0:
        parser.error("--restarts and --max-iterations must be positive")
    resolved = {
        path.resolve()
        for path in (args.target, args.source, args.output, args.metadata)
    }
    if len(resolved) != 4:
        parser.error("target, source, output, and metadata must be distinct")

    target, target_bytes = read_matrix(args.target)
    source, source_bytes = read_matrix(args.source)
    if target.shape != source.shape:
        raise ValueError("target and source orders differ")
    source_score = abs(determinant(source))
    raw_distance = int(np.count_nonzero(source != target))
    aligned, alignment = align(
        target,
        source,
        args.seed,
        args.restarts,
        args.max_iterations,
    )
    aligned_score = abs(determinant(aligned))
    if aligned_score != source_score:
        raise ArithmeticError("alignment did not preserve exact determinant")

    aligned_text = matrix_text(aligned)
    record = {
        "schema_version": 1,
        "event": "finished",
        "complete": True,
        "target": str(args.target),
        "target_raw_sha256": hashlib.sha256(target_bytes).hexdigest(),
        "source": str(args.source),
        "source_raw_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "output": str(args.output),
        "output_raw_sha256": hashlib.sha256(
            aligned_text.encode("utf-8")
        ).hexdigest(),
        "order": int(target.shape[0]),
        "seed": args.seed,
        "restarts": args.restarts,
        "maximum_iterations": args.max_iterations,
        "source_absolute_determinant": str(source_score),
        "aligned_absolute_determinant": str(aligned_score),
        "raw_hamming_distance": raw_distance,
        **alignment,
    }
    atomic_write(args.output, aligned_text)
    atomic_write(
        args.metadata,
        json.dumps(record, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
