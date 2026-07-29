#!/usr/bin/env python3
"""Exact-audit order-22 targets against a frozen set of known H-classes.

For a known H-class, the certificate of any one pivot-dephased incidence
graph must equal one of the pivot certificates of its representative.
Precomputing all known pivot certificates therefore reduces classification
of each candidate from 22^2 graph canonizations to one, without weakening the
test. Candidates absent from the lookup receive the full all-pivot
H-certificate and are reported as novel.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Any

import pynauty

from h_equivalence_audit import determinant, h_certificate, transpose


ORDER = 22
EXPECTED_ABS_DETERMINANT = 409_600_000_000_000


def sha256_bytes(contents: bytes) -> str:
    return hashlib.sha256(contents).hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_matrix(path: Path) -> list[list[int]]:
    try:
        matrix = [
            [int(token) for token in line.split()]
            for line in path.read_text(encoding="ascii").splitlines()
            if line.strip()
        ]
    except (OSError, UnicodeDecodeError, ValueError) as error:
        raise ValueError(f"cannot parse {path}: {error}") from error
    if len(matrix) != ORDER or any(
        len(row) != ORDER or any(value not in (-1, 1) for value in row)
        for row in matrix
    ):
        raise ValueError(f"{path}: expected exactly {ORDER}x{ORDER} signs")
    return matrix


def pivot_certificate(
    matrix: list[list[int]],
    pivot_row: int,
    pivot_column: int,
) -> bytes:
    side = ORDER - 1
    kept_rows = [row for row in range(ORDER) if row != pivot_row]
    kept_columns = [
        column for column in range(ORDER) if column != pivot_column
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
    graph = pynauty.Graph(
        number_of_vertices=2 * side,
        directed=False,
        adjacency_dict=adjacency,
        vertex_coloring=[
            set(range(side)),
            set(range(side, 2 * side)),
        ],
    )
    return pynauty.certificate(graph)


def certificate_sha256(matrix: list[list[int]]) -> str:
    return sha256_bytes(h_certificate(matrix))


def atomic_write(path: Path, contents: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
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
    parser.add_argument(
        "--input-dir",
        action="append",
        required=True,
        type=Path,
    )
    parser.add_argument("--known-report", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()

    input_dirs = [
        path.expanduser().absolute() for path in arguments.input_dir
    ]
    paths: list[tuple[str, Path]] = []
    for input_dir in input_dirs:
        matches = sorted(input_dir.glob("target-*.matrix.txt"))
        paths.extend((input_dir.name, path) for path in matches)
    if not paths:
        parser.error("input directories contain no target matrices")

    known_report_path = arguments.known_report.expanduser().absolute()
    known_report = json.loads(
        known_report_path.read_text(encoding="utf-8")
    )
    known_records = {
        str(record["h_certificate_sha256"]): record
        for record in known_report["classes"]
    }
    if len(known_records) != int(known_report["h_class_count"]):
        raise ValueError("known report has duplicate H certificates")

    pivot_lookup: dict[bytes, str] = {}
    pivot_lookup_multiplicities: Counter[str] = Counter()
    for direct, record in known_records.items():
        matrix = read_matrix(Path(str(record["path"])))
        checked = determinant(matrix)
        if abs(checked) != EXPECTED_ABS_DETERMINANT:
            raise ArithmeticError("known representative is off target")
        full = certificate_sha256(matrix)
        if full != direct:
            raise ArithmeticError("known representative certificate changed")
        for pivot_row in range(ORDER):
            for pivot_column in range(ORDER):
                certificate = pivot_certificate(
                    matrix,
                    pivot_row,
                    pivot_column,
                )
                previous = pivot_lookup.setdefault(certificate, direct)
                if previous != direct:
                    raise ArithmeticError(
                        "one pivot certificate maps to two H-classes"
                    )
                pivot_lookup_multiplicities[direct] += 1

    grouped: dict[str, dict[str, Any]] = {}
    determinant_signs: Counter[str] = Counter()
    arm_counts: Counter[str] = Counter()
    named_payload = hashlib.sha256()
    novel_records: list[dict[str, Any]] = []
    for index, (arm, path) in enumerate(paths, start=1):
        contents = path.read_bytes()
        named_payload.update(arm.encode("utf-8"))
        named_payload.update(b"/")
        named_payload.update(path.name.encode("ascii"))
        named_payload.update(b"\0")
        named_payload.update(contents)
        matrix = read_matrix(path)
        checked = determinant(matrix)
        if abs(checked) != EXPECTED_ABS_DETERMINANT:
            raise ArithmeticError(
                f"{path}: exact determinant is off target: {checked}"
            )
        determinant_signs["positive" if checked > 0 else "negative"] += 1
        arm_counts[arm] += 1

        one_pivot = pivot_certificate(matrix, 0, 0)
        direct = pivot_lookup.get(one_pivot)
        novel = direct is None
        if novel:
            direct = certificate_sha256(matrix)
            if direct in known_records:
                raise ArithmeticError(
                    "known H-class was absent from pivot lookup"
                )
            transposed = certificate_sha256(transpose(matrix))
            ht = min(direct, transposed)
            novel_records.append(
                {
                    "path": str(path),
                    "raw_sha256": sha256_bytes(contents),
                    "determinant": checked,
                    "h_certificate_sha256": direct,
                    "transpose_h_certificate_sha256": transposed,
                    "ht_certificate_sha256": ht,
                }
            )
        else:
            record = known_records[direct]
            transposed = str(record["transpose_h_certificate_sha256"])
            ht = str(record["ht_certificate_sha256"])

        group = grouped.setdefault(
            direct,
            {
                "h_certificate_sha256": direct,
                "transpose_h_certificate_sha256": transposed,
                "ht_certificate_sha256": ht,
                "known_h_class": not novel,
                "matrix_count": 0,
                "arm_counts": defaultdict(int),
                "representative": {
                    "path": str(path),
                    "raw_sha256": sha256_bytes(contents),
                    "determinant": checked,
                },
            },
        )
        if bool(group["known_h_class"]) == novel:
            raise ArithmeticError("H-class novelty status changed")
        group["matrix_count"] = int(group["matrix_count"]) + 1
        group["arm_counts"][arm] += 1
        if index % 1000 == 0:
            print(
                f"checked={index}/{len(paths)} "
                f"classes={len(grouped)} novel={len(novel_records)}",
                flush=True,
            )

    classes = []
    for direct in sorted(grouped):
        group = grouped[direct]
        group["arm_counts"] = dict(sorted(group["arm_counts"].items()))
        classes.append(group)
    report = {
        "schema_version": 1,
        "engine": "order22-known-class-pivot-audit-v1",
        "claim": (
            "Every input determinant was checked exactly. One pivot graph "
            "certificate per input was matched against every pivot graph "
            "certificate of each frozen known H-class. An equal pivot graph "
            "is a complete witness of H-equivalence."
        ),
        "order": ORDER,
        "target_absolute_determinant": EXPECTED_ABS_DETERMINANT,
        "input_dirs": [str(path) for path in input_dirs],
        "input_matrix_count": len(paths),
        "exact_determinants_checked": len(paths),
        "all_exact_determinants_on_target": True,
        "determinant_sign_counts": dict(sorted(determinant_signs.items())),
        "input_arm_counts": dict(sorted(arm_counts.items())),
        "named_payload_sha256": named_payload.hexdigest(),
        "known_report": {
            "path": str(known_report_path),
            "sha256": sha256(known_report_path),
            "h_class_count": len(known_records),
            "ht_class_count": int(known_report["ht_class_count"]),
        },
        "known_pivot_certificates_generated": ORDER
        * ORDER
        * len(known_records),
        "distinct_known_pivot_certificates": len(pivot_lookup),
        "known_h_classes_observed": sum(
            group["known_h_class"] is True for group in classes
        ),
        "observed_h_class_count": len(classes),
        "observed_ht_class_count": len(
            {
                str(group["ht_certificate_sha256"])
                for group in classes
            }
        ),
        "novel_h_class_count": sum(
            group["known_h_class"] is False for group in classes
        ),
        "novel_matrices": novel_records,
        "dependencies": {
            "pynauty": getattr(pynauty, "__version__", "unknown"),
        },
        "source": {
            "path": str(Path(__file__).absolute()),
            "sha256": sha256(Path(__file__).absolute()),
        },
        "classes": classes,
    }
    output = arguments.output.expanduser().absolute()
    atomic_write(
        output,
        (json.dumps(report, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    print(
        f"matrices={report['input_matrix_count']} "
        f"h_classes={report['observed_h_class_count']} "
        f"ht_classes={report['observed_ht_class_count']} "
        f"novel_h_classes={report['novel_h_class_count']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
