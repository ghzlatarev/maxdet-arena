#!/usr/bin/env python3
"""Classify order-23 sign matrices under Hadamard equivalence.

For every pivot (r, c), the matrix is dephased so pivot row r and pivot
column c are positive.  Deleting those lines leaves a 22 by 22 binary
incidence matrix.  Its color-preserving bipartite-graph certificate is
invariant under row and column permutations.  The minimum certificate over
all 23^2 pivots is therefore a complete invariant for signed row/column
permutation (H-) equivalence.

This is a research audit, not part of the trusted verifier.  It requires the
pinned package ``pynauty==2.8.8.1``.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import sys
from typing import Iterable

try:
    import pynauty
except ImportError:
    print(
        "error: install pynauty==2.8.8.1 in an isolated environment",
        file=sys.stderr,
    )
    raise SystemExit(2)


Matrix = list[list[int]]
ROOT = Path(__file__).resolve().parents[1]
FRONTIER = 2_779_447_296_000_000

LOCAL_CORPUS_PATHS = (
    "references/orrick-et-al-2003/matrix.txt",
    "runs/direct-search/reference-data/orrick-pre-april2003.matrix.txt",
    "runs/direct-search/best-below/frontier-class14-late-28751.matrix.txt",
    (
        "runs/direct-search/h24-deletion-remaining-classes/elites/"
        "class-9-r10-c24.matrix.txt"
    ),
    "runs/direct-search/best-below/frontier-class51-28752.matrix.txt",
    (
        "runs/direct-search/multiflip/class9-depth12-tie-replay-29831/"
        "tie.matrix.txt"
    ),
    (
        "runs/direct-search/multiflip/depth12-tie-harvest-29851/"
        "tie.matrix.txt"
    ),
    (
        "runs/direct-search/multiflip/tie2-radius12-harvest-29862/"
        "tie.matrix.txt"
    ),
    "runs/direct-search/gram-shell-milp-reference-29761.matrix.txt",
    "runs/direct-search/gram-shell-milp-reference-seed29762.matrix.txt",
    "runs/direct-search/gram-shell-milp-reference-seed29768.matrix.txt",
    "runs/direct-search/gram-shell-milp-reference-seed29773.matrix.txt",
    (
        "runs/direct-search/hamming-sphere-long/tie1-r24-29901/"
        "tie.matrix.txt"
    ),
    (
        "runs/direct-search/hamming-sphere-long/class9-r32-29902/"
        "tie.matrix.txt"
    ),
    (
        "runs/direct-search/multiflip/sphere32-radius12-harvest-29940/"
        "tie.matrix.txt"
    ),
)
LOCAL_CORPUS_GLOBS = ("runs/direct-search/gram-factors/*/*.matrix.txt",)


def read_matrix(path: Path) -> tuple[bytes, Matrix]:
    raw = path.read_bytes()
    try:
        rows = [
            [int(token) for token in line.split()]
            for line in raw.decode("utf-8").splitlines()
            if line.strip()
        ]
    except (UnicodeDecodeError, ValueError) as error:
        raise ValueError(f"{path}: invalid integer matrix") from error
    if len(rows) != 23 or any(len(row) != 23 for row in rows):
        raise ValueError(f"{path}: expected exactly 23 rows and 23 columns")
    if any(value not in (-1, 1) for row in rows for value in row):
        raise ValueError(f"{path}: entries must all be +1 or -1")
    return raw, rows


def determinant(matrix: Matrix) -> int:
    """Return the exact determinant by fraction-free Bareiss elimination."""
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


def transpose(matrix: Matrix) -> Matrix:
    return [list(column) for column in zip(*matrix)]


def graph_certificate(
    vertex_count: int,
    adjacency: dict[int, list[int]],
    coloring: list[set[int]] | None = None,
) -> bytes:
    graph = pynauty.Graph(
        number_of_vertices=vertex_count,
        directed=False,
        adjacency_dict=adjacency,
        vertex_coloring=coloring or [],
    )
    return pynauty.certificate(graph)


def h_certificate(matrix: Matrix) -> bytes:
    """Return the exact H-equivalence certificate for one sign matrix."""
    size = len(matrix)
    side = size - 1
    best: bytes | None = None
    for pivot_row in range(size):
        kept_rows = [row for row in range(size) if row != pivot_row]
        for pivot_column in range(size):
            kept_columns = [
                column for column in range(size) if column != pivot_column
            ]
            adjacency = {vertex: [] for vertex in range(2 * side)}
            pivot = matrix[pivot_row][pivot_column]
            for row_vertex, row in enumerate(kept_rows):
                for column_vertex, column in enumerate(kept_columns):
                    dephased = (
                        matrix[row][column]
                        * matrix[pivot_row][column]
                        * matrix[row][pivot_column]
                        * pivot
                    )
                    if dephased == -1:
                        other = side + column_vertex
                        adjacency[row_vertex].append(other)
                        adjacency[other].append(row_vertex)
            certificate = graph_certificate(
                2 * side,
                adjacency,
                [set(range(side)), set(range(side, 2 * side))],
            )
            if best is None or certificate < best:
                best = certificate
    assert best is not None
    return best


def normalized_gram_graph(
    matrix: Matrix,
) -> tuple[bytes, int, int, tuple[int, ...]]:
    """Return certificate, automorphism order, edge count, and degrees.

    This binary defect-graph encoding applies when row switching puts every
    off-diagonal Gram entry in {-1, 3}, as it does for this frontier corpus.
    """
    size = len(matrix)
    gram = [
        [
            sum(matrix[row][k] * matrix[column][k] for k in range(size))
            for column in range(size)
        ]
        for row in range(size)
    ]
    signs = [1]
    for index in range(1, size):
        value = gram[0][index]
        signs.append(1 if value % 4 == 3 else -1)
    normalized = [
        [signs[row] * signs[column] * gram[row][column] for column in range(size)]
        for row in range(size)
    ]
    values = {
        normalized[row][column]
        for row in range(size)
        for column in range(row)
    }
    if values != {-1, 3}:
        raise ValueError(
            "normalized Gram is outside the {-1, 3} binary defect model: "
            f"{sorted(values)}"
        )
    adjacency = {vertex: [] for vertex in range(size)}
    for row in range(size):
        for column in range(row):
            if normalized[row][column] == 3:
                adjacency[row].append(column)
                adjacency[column].append(row)
    graph = pynauty.Graph(
        number_of_vertices=size,
        directed=False,
        adjacency_dict=adjacency,
    )
    certificate = pynauty.certificate(graph)
    automorphisms = pynauty.autgrp(graph)
    group_order = int(round(automorphisms[1] * (10 ** automorphisms[2])))
    degrees = tuple(sorted(len(neighbors) for neighbors in adjacency.values()))
    edge_count = sum(degrees) // 2
    return certificate, group_order, edge_count, degrees


def display_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def expand_inputs(
    explicit: Iterable[Path],
    patterns: Iterable[str],
    local_corpus: bool,
) -> list[Path]:
    paths = list(explicit)
    globs = list(patterns)
    if local_corpus:
        paths.extend(ROOT / relative for relative in LOCAL_CORPUS_PATHS)
        globs.extend(LOCAL_CORPUS_GLOBS)
    for pattern in globs:
        base = ROOT if local_corpus and not Path(pattern).is_absolute() else Path(".")
        matches = sorted(base.glob(pattern))
        if not matches:
            raise ValueError(f"glob matched no files: {pattern}")
        paths.extend(matches)
    unique: dict[Path, None] = {}
    for path in paths:
        candidate = path if path.is_absolute() else Path.cwd() / path
        unique[candidate.resolve()] = None
    if not unique:
        raise ValueError("provide matrix paths, --glob, or --local-corpus")
    missing = [path for path in unique if not path.is_file()]
    if missing:
        raise ValueError("missing input: " + ", ".join(map(str, missing)))
    return list(unique)


def sha256_hex(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def audit(paths: list[Path], expected_determinant: int) -> dict[str, object]:
    unique_matrices: dict[str, dict[str, object]] = {}
    duplicates: list[dict[str, str]] = []
    for path in paths:
        raw, matrix = read_matrix(path)
        raw_hash = sha256_hex(raw)
        if raw_hash in unique_matrices:
            duplicates.append(
                {
                    "matrix_sha256": raw_hash,
                    "first": str(unique_matrices[raw_hash]["path"]),
                    "duplicate": display_path(path),
                }
            )
            continue
        absolute_determinant = abs(determinant(matrix))
        if absolute_determinant != expected_determinant:
            raise ValueError(
                f"{path}: |det|={absolute_determinant}, "
                f"expected {expected_determinant}"
            )
        direct = h_certificate(matrix)
        transposed = h_certificate(transpose(matrix))
        gram, automorphism_order, edges, degrees = normalized_gram_graph(matrix)
        unique_matrices[raw_hash] = {
            "path": display_path(path),
            "matrix_sha256": raw_hash,
            "h_certificate": direct,
            "ht_certificate": min(direct, transposed),
            "self_dual": direct == transposed,
            "gram_certificate": gram,
            "gram_automorphism_order": automorphism_order,
            "gram_edges": edges,
            "gram_degrees": degrees,
        }

    h_classes: dict[bytes, list[dict[str, object]]] = defaultdict(list)
    ht_classes: dict[bytes, list[dict[str, object]]] = defaultdict(list)
    gram_classes: dict[bytes, list[dict[str, object]]] = defaultdict(list)
    for entry in unique_matrices.values():
        h_classes[entry["h_certificate"]].append(entry)
        ht_classes[entry["ht_certificate"]].append(entry)
        gram_classes[entry["gram_certificate"]].append(entry)

    def classes_payload(
        classes: dict[bytes, list[dict[str, object]]],
    ) -> list[dict[str, object]]:
        payload = []
        for certificate, members in sorted(
            classes.items(), key=lambda item: sha256_hex(item[0])
        ):
            payload.append(
                {
                    "certificate_sha256": sha256_hex(certificate),
                    "count": len(members),
                    "members": [
                        {
                            "path": member["path"],
                            "matrix_sha256": member["matrix_sha256"],
                        }
                        for member in sorted(members, key=lambda item: item["path"])
                    ],
                }
            )
        return payload

    gram_payload = classes_payload(gram_classes)
    for group, (_, members) in zip(
        gram_payload,
        sorted(gram_classes.items(), key=lambda item: sha256_hex(item[0])),
    ):
        representative = members[0]
        group.update(
            {
                "automorphism_group_order": representative[
                    "gram_automorphism_order"
                ],
                "edge_count": representative["gram_edges"],
                "degree_multiset": list(representative["gram_degrees"]),
            }
        )

    return {
        "schema_version": 1,
        "method": (
            "minimum color-preserving nauty certificate over all 23^2 "
            "dephased pivots"
        ),
        "pynauty_version": getattr(pynauty, "__version__", "unknown"),
        "expected_absolute_determinant": expected_determinant,
        "input_file_count": len(paths),
        "duplicate_file_count": len(duplicates),
        "unique_matrix_count": len(unique_matrices),
        "h_class_count": len(h_classes),
        "ht_class_count": len(ht_classes),
        "gram_class_count": len(gram_classes),
        "h_classes": classes_payload(h_classes),
        "ht_classes": classes_payload(ht_classes),
        "gram_classes": gram_payload,
        "duplicates": duplicates,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", nargs="*", type=Path)
    parser.add_argument(
        "--glob",
        action="append",
        default=[],
        help="add paths matching a glob (repeatable)",
    )
    parser.add_argument(
        "--local-corpus",
        action="store_true",
        help="audit the frozen 2026-07-28 local frontier-corpus path set",
    )
    parser.add_argument(
        "--expected-absolute-determinant",
        type=int,
        default=FRONTIER,
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="write the JSON audit to a fresh path instead of stdout",
    )
    arguments = parser.parse_args()
    try:
        paths = expand_inputs(
            arguments.matrix,
            arguments.glob,
            arguments.local_corpus,
        )
        result = audit(paths, arguments.expected_absolute_determinant)
    except (OSError, ValueError, ArithmeticError) as error:
        parser.error(str(error))
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(serialized, end="")
    else:
        if arguments.output.exists():
            parser.error(f"refusing to overwrite --output: {arguments.output}")
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(serialized, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
