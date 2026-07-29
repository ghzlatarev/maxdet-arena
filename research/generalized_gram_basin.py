#!/usr/bin/env python3
"""Exactly classify local order-23 matrices by signed row/column Gram basins.

For a sign matrix H, row sign changes act on G = H H^T by signed
permutation congruence.  The parity of Hamming distances gives a canonical
switching: after choosing any anchor row, switch every off-diagonal entry of
G to 3 modulo 4.  The remaining equivalence is ordinary permutation of the
23 rows.

The normalized complete edge-labeled graph is encoded as a colored incidence
graph and certified with pinned ``pynauty==2.8.8.1``.  Row and column
certificates are sorted to form an unordered H/transpose (HT) Gram-basin key.
This is an exact classifier of the supplied Gram matrices.  It is not a
factor-equivalence classifier and makes only local-corpus novelty claims.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import glob as glob_module
import hashlib
import json
from pathlib import Path
import struct
import sys
from typing import Sequence

try:
    import pynauty
except ImportError:
    print(
        "error: install pynauty==2.8.8.1 in an isolated environment",
        file=sys.stderr,
    )
    raise SystemExit(2)


ORDER = 23
PINNED_PYNAUTY = "2.8.8.1"
ROOT = Path(__file__).resolve().parents[1]
GRAM_LABEL_PALETTE = tuple(range(-21, 24, 4))
CERTIFICATE_DOMAIN = b"maxdet-23-general-signed-gram-certificate-v1\0"
BASIN_DOMAIN = b"maxdet-23-unordered-row-column-gram-basin-v1\0"
DESCRIPTOR_DOMAIN = b"maxdet-23-signed-gram-descriptor-v1\0"

Matrix = list[list[int]]
Gram = list[list[int]]


def sha256_bytes(contents: bytes) -> str:
    return hashlib.sha256(contents).hexdigest()


def display_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def read_matrix(path: Path) -> tuple[bytes, Matrix]:
    try:
        raw = path.read_bytes()
        text = raw.decode("ascii")
        matrix = [
            [int(token) for token in line.split()]
            for line in text.splitlines()
            if line.strip()
        ]
    except (OSError, UnicodeDecodeError, ValueError) as error:
        raise ValueError(f"{path}: cannot parse sign matrix: {error}") from error
    if len(matrix) != ORDER or any(len(row) != ORDER for row in matrix):
        raise ValueError(f"{path}: expected exactly {ORDER}x{ORDER} entries")
    if any(value not in (-1, 1) for row in matrix for value in row):
        raise ValueError(f"{path}: every entry must be -1 or +1")
    return raw, matrix


def determinant(matrix: Sequence[Sequence[int]]) -> int:
    """Return the exact integer determinant by Bareiss elimination."""

    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("Bareiss input must be a nonempty square matrix")
    work = [list(row) for row in matrix]
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


def transpose(matrix: Sequence[Sequence[int]]) -> Matrix:
    return [list(column) for column in zip(*matrix)]


def gram(matrix: Sequence[Sequence[int]]) -> Gram:
    size = len(matrix)
    if size != ORDER or any(len(row) != ORDER for row in matrix):
        raise ValueError(f"Gram input must be exactly {ORDER}x{ORDER}")
    return [
        [
            sum(matrix[row][index] * matrix[column][index] for index in range(ORDER))
            for column in range(ORDER)
        ]
        for row in range(ORDER)
    ]


def normalize_signed_gram(matrix: Gram, anchor: int = 0) -> Gram:
    """Switch an order-23 sign Gram so every off-diagonal is 3 modulo 4."""

    if len(matrix) != ORDER or any(len(row) != ORDER for row in matrix):
        raise ValueError(f"signed Gram must be exactly {ORDER}x{ORDER}")
    if not 0 <= anchor < ORDER:
        raise ValueError("Gram normalization anchor is out of range")
    for row in range(ORDER):
        if matrix[row][row] != ORDER:
            raise ValueError("signed Gram diagonal must equal 23")
        for column in range(row):
            if matrix[row][column] != matrix[column][row]:
                raise ValueError("signed Gram must be symmetric")

    signs = [1 if matrix[anchor][index] % 4 == 3 else -1 for index in range(ORDER)]
    signs[anchor] = 1
    normalized = [
        [
            signs[row] * signs[column] * matrix[row][column]
            for column in range(ORDER)
        ]
        for row in range(ORDER)
    ]
    for row in range(ORDER):
        for column in range(row):
            value = normalized[row][column]
            if value % 4 != 3:
                raise ValueError(
                    "input is not an order-23 sign Gram: mod-4 switching "
                    f"failed at ({row}, {column}) with value {value}"
                )
            if value not in GRAM_LABEL_PALETTE:
                raise ValueError(
                    f"normalized off-diagonal Gram value is out of range: {value}"
                )
    return normalized


def gram_descriptor(normalized: Gram) -> dict[str, object]:
    """Return a useful but deliberately noncanonical structural descriptor."""

    label_to_index = {
        label: index for index, label in enumerate(GRAM_LABEL_PALETTE)
    }
    global_histogram = [0] * len(GRAM_LABEL_PALETTE)
    profiles: list[list[int]] = []
    for row in range(ORDER):
        profile = [0] * len(GRAM_LABEL_PALETTE)
        for column in range(ORDER):
            if row == column:
                continue
            index = label_to_index[normalized[row][column]]
            profile[index] += 1
            if column < row:
                global_histogram[index] += 1
        profiles.append(profile)
    profiles.sort()
    descriptor = {
        "palette": list(GRAM_LABEL_PALETTE),
        "global_label_histogram": global_histogram,
        "vertex_incident_label_profiles": profiles,
    }
    canonical = json.dumps(
        descriptor, separators=(",", ":"), sort_keys=True
    ).encode("ascii")
    descriptor["descriptor_sha256"] = sha256_bytes(
        DESCRIPTOR_DOMAIN + canonical
    )
    descriptor["claim"] = (
        "Search descriptor only; equal descriptors do not prove Gram "
        "equivalence. The pynauty certificate is the exact classifier."
    )
    return descriptor


def incidence_graph(
    normalized: Gram,
) -> tuple[pynauty.Graph, tuple[int, ...]]:
    """Encode the normalized complete labeled graph as a colored incidence graph."""

    edge_count = ORDER * (ORDER - 1) // 2
    adjacency = {vertex: [] for vertex in range(ORDER + edge_count)}
    color_vertices = {label: set() for label in GRAM_LABEL_PALETTE}
    edge_vertex = ORDER
    for row in range(ORDER):
        for column in range(row):
            adjacency[row].append(edge_vertex)
            adjacency[column].append(edge_vertex)
            adjacency[edge_vertex].extend((row, column))
            color_vertices[normalized[row][column]].add(edge_vertex)
            edge_vertex += 1
    histogram = tuple(
        len(color_vertices[label]) for label in GRAM_LABEL_PALETTE
    )
    # Empty pynauty cells are semantically invisible.  Keep the nonempty cells
    # in fixed palette order, and bind the complete palette histogram into the
    # outer certificate payload below.
    coloring = [set(range(ORDER))]
    coloring.extend(
        color_vertices[label]
        for label in GRAM_LABEL_PALETTE
        if color_vertices[label]
    )
    graph = pynauty.Graph(
        number_of_vertices=ORDER + edge_count,
        directed=False,
        adjacency_dict=adjacency,
        vertex_coloring=coloring,
    )
    return graph, histogram


def certificate_payload(
    normalized: Gram,
) -> tuple[bytes, bytes, tuple[int, ...]]:
    """Return complete domain-separated bytes and the raw pynauty certificate."""

    graph, histogram = incidence_graph(normalized)
    raw_certificate = pynauty.certificate(graph)
    header = bytearray(CERTIFICATE_DOMAIN)
    header.extend(struct.pack(">H", ORDER))
    for label, count in zip(GRAM_LABEL_PALETTE, histogram):
        header.extend(struct.pack(">hH", label, count))
    header.extend(struct.pack(">I", len(raw_certificate)))
    payload = bytes(header) + raw_certificate
    return payload, raw_certificate, histogram


def certificate_record(matrix: Gram, anchor: int = 0) -> dict[str, object]:
    normalized = normalize_signed_gram(matrix, anchor=anchor)
    payload, raw_certificate, histogram = certificate_payload(normalized)
    return {
        "certificate_bytes": payload,
        "certificate_sha256": sha256_bytes(payload),
        "pynauty_certificate_sha256": sha256_bytes(raw_certificate),
        "label_histogram": list(histogram),
        "descriptor": gram_descriptor(normalized),
    }


def basin_payload(row_certificate: bytes, column_certificate: bytes) -> bytes:
    ordered = sorted((row_certificate, column_certificate))
    payload = bytearray(BASIN_DOMAIN)
    for certificate in ordered:
        payload.extend(struct.pack(">I", len(certificate)))
        payload.extend(certificate)
    return bytes(payload)


def classify_matrix(label: str, path: Path) -> dict[str, object]:
    raw, matrix = read_matrix(path)
    exact_determinant = determinant(matrix)
    row = certificate_record(gram(matrix))
    column = certificate_record(gram(transpose(matrix)))
    basin = basin_payload(
        row.pop("certificate_bytes"),
        column.pop("certificate_bytes"),
    )
    unordered = sorted(
        (str(row["certificate_sha256"]), str(column["certificate_sha256"]))
    )
    return {
        "label": label,
        "path": display_path(path),
        "matrix_sha256": sha256_bytes(raw),
        "determinant": str(exact_determinant),
        "absolute_determinant": str(abs(exact_determinant)),
        "row_gram": row,
        "column_gram": column,
        "unordered_row_column_gram_certificate_sha256s": unordered,
        "ht_gram_basin_key_sha256": sha256_bytes(basin),
    }


def group_records(
    records: Sequence[dict[str, object]],
    kind: str,
) -> list[dict[str, object]]:
    groups: dict[str, list[dict[str, object]]] = defaultdict(list)

    def key(record: dict[str, object]) -> str:
        if kind in ("row_gram", "column_gram"):
            metadata = record[kind]
            if not isinstance(metadata, dict):
                raise ArithmeticError(f"malformed {kind} record")
            return str(metadata["certificate_sha256"])
        if kind == "ht_gram_basin":
            return str(record["ht_gram_basin_key_sha256"])
        raise ValueError(f"unknown grouping kind: {kind}")

    for record in records:
        groups[key(record)].append(record)

    result = []
    for certificate, members in sorted(groups.items()):
        item: dict[str, object] = {
            "key_sha256": certificate,
            "count": len(members),
            "members": sorted(str(member["label"]) for member in members),
        }
        representative = members[0]
        if kind == "ht_gram_basin":
            item["unordered_row_column_gram_certificate_sha256s"] = (
                representative[
                    "unordered_row_column_gram_certificate_sha256s"
                ]
            )
        result.append(item)
    return result


def parse_matrix_argument(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError(
            "--matrix requires LABEL=PATH (for example reference=matrix.txt)"
        )
    label, path = value.split("=", 1)
    if not label.strip() or not path:
        raise argparse.ArgumentTypeError("--matrix requires nonempty LABEL=PATH")
    return label.strip(), Path(path)


def expand_inputs(
    explicit: Sequence[tuple[str, Path]],
    patterns: Sequence[str],
) -> list[tuple[str, Path]]:
    entries: list[tuple[str, Path]] = []
    for label, path in explicit:
        candidate = path if path.is_absolute() else Path.cwd() / path
        entries.append((label, candidate.resolve()))
    for pattern in patterns:
        matches = sorted(
            Path(match).resolve()
            for match in glob_module.glob(pattern, recursive=True)
        )
        if not matches:
            raise ValueError(f"glob matched no paths: {pattern}")
        for path in matches:
            if path.is_file():
                entries.append((display_path(path), path))
    if not entries:
        raise ValueError("provide repeated --matrix LABEL=PATH and/or --glob")

    seen_labels: dict[str, Path] = {}
    seen_paths: dict[Path, str] = {}
    result: list[tuple[str, Path]] = []
    for label, path in entries:
        if not path.is_file():
            raise ValueError(f"{label}: missing matrix file: {path}")
        if label in seen_labels:
            if seen_labels[label] == path:
                continue
            raise ValueError(
                f"duplicate label {label!r} names both "
                f"{seen_labels[label]} and {path}"
            )
        if path in seen_paths:
            continue
        seen_labels[label] = path
        seen_paths[path] = label
        result.append((label, path))
    return result


def apply_hadamard_operations(
    matrix: Matrix,
    row_permutation: Sequence[int],
    column_permutation: Sequence[int],
    row_signs: Sequence[int],
    column_signs: Sequence[int],
) -> Matrix:
    return [
        [
            row_signs[row] * column_signs[column]
            * matrix[row_permutation[row]][column_permutation[column]]
            for column in range(ORDER)
        ]
        for row in range(ORDER)
    ]


def synthetic_normalized_gram(default: int = -1) -> Gram:
    return [
        [ORDER if row == column else default for column in range(ORDER)]
        for row in range(ORDER)
    ]


def run_self_test() -> dict[str, object]:
    """Exercise invariants and edge-color identity deterministically."""

    if getattr(pynauty, "__version__", None) != PINNED_PYNAUTY:
        raise RuntimeError(f"requires pynauty=={PINNED_PYNAUTY}")
    base = [
        [-1 if row == column else 1 for column in range(ORDER)]
        for row in range(ORDER)
    ]
    expected = (ORDER - 2) * (2 ** (ORDER - 1))
    if determinant(base) != expected:
        raise AssertionError("Bareiss failed the J-2I determinant identity")

    row_permutation = list(range(ORDER))
    row_permutation[0], row_permutation[7] = 7, 0
    row_permutation[3], row_permutation[19] = 19, 3
    column_permutation = list(reversed(range(ORDER)))
    row_signs = [-1 if index % 3 == 0 else 1 for index in range(ORDER)]
    column_signs = [-1 if index % 5 in (1, 4) else 1 for index in range(ORDER)]
    transformed = apply_hadamard_operations(
        base,
        row_permutation,
        column_permutation,
        row_signs,
        column_signs,
    )
    base_row = certificate_record(gram(base))
    base_column = certificate_record(gram(transpose(base)))
    transformed_row = certificate_record(gram(transformed))
    transformed_column = certificate_record(gram(transpose(transformed)))
    if base_row["certificate_sha256"] != transformed_row["certificate_sha256"]:
        raise AssertionError("row signed-permutation invariance failed")
    if (
        base_column["certificate_sha256"]
        != transformed_column["certificate_sha256"]
    ):
        raise AssertionError("column signed-permutation invariance failed")

    column_only = apply_hadamard_operations(
        base,
        list(range(ORDER)),
        column_permutation,
        [1] * ORDER,
        column_signs,
    )
    if (
        base_row["certificate_sha256"]
        != certificate_record(gram(column_only))["certificate_sha256"]
    ):
        raise AssertionError("column operations changed the row Gram certificate")

    base_gram = gram(base)
    anchor_certificates = {
        certificate_record(base_gram, anchor=anchor)["certificate_sha256"]
        for anchor in range(ORDER)
    }
    if len(anchor_certificates) != 1:
        raise AssertionError("Gram normalization depends on the chosen anchor")

    base_basin = basin_payload(
        base_row["certificate_bytes"], base_column["certificate_bytes"]
    )
    transposed_basin = basin_payload(
        base_column["certificate_bytes"], base_row["certificate_bytes"]
    )
    if base_basin != transposed_basin:
        raise AssertionError("HT Gram-basin key is not transpose invariant")

    uniform_minus = certificate_record(synthetic_normalized_gram(-1))
    uniform_three = certificate_record(synthetic_normalized_gram(3))
    if (
        uniform_minus["certificate_sha256"]
        == uniform_three["certificate_sha256"]
    ):
        raise AssertionError("fixed palette failed to distinguish uniform labels")

    first = synthetic_normalized_gram(-1)
    second = synthetic_normalized_gram(-1)
    triangle = ((0, 1), (1, 2), (0, 2))
    matching = ((3, 4), (5, 6), (7, 8))
    for left, right in triangle:
        first[left][right] = first[right][left] = 3
        second[left][right] = second[right][left] = 7
    for left, right in matching:
        first[left][right] = first[right][left] = 7
        second[left][right] = second[right][left] = 3
    first_record = certificate_record(first)
    second_record = certificate_record(second)
    if first_record["label_histogram"] != second_record["label_histogram"]:
        raise AssertionError("equal-cardinality color-swap fixture is malformed")
    if first_record["certificate_sha256"] == second_record["certificate_sha256"]:
        raise AssertionError(
            "fixed ordered palette failed an asymmetric equal-cardinality swap"
        )

    one_change = synthetic_normalized_gram(-1)
    one_change[0][1] = one_change[1][0] = 3
    if (
        uniform_minus["certificate_sha256"]
        == certificate_record(one_change)["certificate_sha256"]
    ):
        raise AssertionError("one edge-label change did not alter the certificate")

    return {
        "status": "passed",
        "tests": [
            "Bareiss J-2I identity",
            "row sign/permutation invariance",
            "column sign/permutation invariance",
            "column-operation row-certificate invariance",
            "all 23 normalization anchors invariant",
            "transpose-invariant unordered basin key",
            "fixed edge-color identity with absent palette cells",
            "asymmetric equal-cardinality edge-color swap",
            "one-label-change separation",
            "descriptor explicitly marked noncanonical",
        ],
        "pynauty_version": getattr(pynauty, "__version__", "unknown"),
    }


def build_report(entries: Sequence[tuple[str, Path]]) -> dict[str, object]:
    records = [classify_matrix(label, path) for label, path in entries]
    row_groups = group_records(records, "row_gram")
    column_groups = group_records(records, "column_gram")
    basin_groups = group_records(records, "ht_gram_basin")
    return {
        "schema_version": 1,
        "method": (
            "exact Bareiss score; mod-4 normalized signed row/column Grams; "
            "fixed-palette colored-incidence pynauty certificates"
        ),
        "order": ORDER,
        "pynauty_version": getattr(pynauty, "__version__", "unknown"),
        "input_matrix_count": len(records),
        "row_gram_group_count": len(row_groups),
        "column_gram_group_count": len(column_groups),
        "ht_gram_basin_group_count": len(basin_groups),
        "records": records,
        "groups": {
            "row_gram": row_groups,
            "column_gram": column_groups,
            "ht_gram_basin": basin_groups,
        },
        "claim_boundary": {
            "scope": (
                "Class counts and novelty are exact only for the explicitly "
                "supplied local corpus."
            ),
            "different_keys": (
                "Different unordered row/column Gram-basin keys prove the "
                "matrices are not H/transpose equivalent."
            ),
            "equal_keys": (
                "Equal Gram-basin keys do not prove that two sign-matrix "
                "factors are H/transpose equivalent."
            ),
            "excluded_claims": (
                "This report does not classify all order-23 matrices and does "
                "not establish literature novelty, global optimality, or a "
                "world record."
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix",
        action="append",
        default=[],
        type=parse_matrix_argument,
        metavar="LABEL=PATH",
        help="add one labeled 23x23 sign matrix (repeatable)",
    )
    parser.add_argument(
        "--glob",
        action="append",
        default=[],
        help="add matrix paths matching a glob (repeatable; supports **)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="write JSON to a fresh path instead of stdout",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="run deterministic invariance and color-identity tests",
    )
    arguments = parser.parse_args()
    if getattr(pynauty, "__version__", None) != PINNED_PYNAUTY:
        parser.error(f"requires pynauty=={PINNED_PYNAUTY}")
    try:
        if arguments.self_test:
            if arguments.matrix or arguments.glob:
                parser.error("--self-test cannot be combined with matrix inputs")
            result = run_self_test()
        else:
            entries = expand_inputs(arguments.matrix, arguments.glob)
            result = build_report(entries)
    except (OSError, ValueError, ArithmeticError, RuntimeError) as error:
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
