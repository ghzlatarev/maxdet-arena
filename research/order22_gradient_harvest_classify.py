#!/usr/bin/env python3
"""Exact-classify target matrices emitted by order22_gradient_harvest."""

from __future__ import annotations

import argparse
from collections import defaultdict
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


def certificate_sha256(matrix: list[list[int]]) -> str:
    return sha256_bytes(h_certificate(matrix))


def trace_gram_squared(matrix: list[list[int]]) -> int:
    gram = [
        [
            sum(
                matrix[inner][row] * matrix[inner][column]
                for inner in range(ORDER)
            )
            for column in range(ORDER)
        ]
        for row in range(ORDER)
    ]
    return sum(
        gram[row][column] * gram[column][row]
        for row in range(ORDER)
        for column in range(ORDER)
    )


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
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--known-report", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()

    input_dir = arguments.input_dir.expanduser().absolute()
    paths = sorted(input_dir.glob("target-*.matrix.txt"))
    if not paths:
        parser.error("input directory contains no target matrices")
    known_h_classes: set[str] = set()
    known_report_sha256 = None
    if arguments.known_report:
        known_path = arguments.known_report.expanduser().absolute()
        known = json.loads(known_path.read_text(encoding="utf-8"))
        known_h_classes = {
            str(record["h_certificate_sha256"])
            for record in known["classes"]
        }
        known_report_sha256 = sha256(known_path)

    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for path in paths:
        matrix = read_matrix(path)
        checked = determinant(matrix)
        if abs(checked) != EXPECTED_ABS_DETERMINANT:
            raise ArithmeticError(
                f"{path}: exact determinant is off target: {checked}"
            )
        direct = certificate_sha256(matrix)
        transposed = certificate_sha256(transpose(matrix))
        grouped[direct].append(
            {
                "path": str(path),
                "raw_sha256": sha256(path),
                "determinant": checked,
                "transpose_h_certificate_sha256": transposed,
                "ht_certificate_sha256": min(direct, transposed),
                "trace_gram_squared": trace_gram_squared(matrix),
            }
        )

    class_records = []
    for direct, members in sorted(grouped.items()):
        members.sort(key=lambda item: str(item["path"]))
        invariants = {
            (
                str(member["transpose_h_certificate_sha256"]),
                str(member["ht_certificate_sha256"]),
                int(member["trace_gram_squared"]),
            )
            for member in members
        }
        if len(invariants) != 1:
            raise ArithmeticError("H-class invariant mismatch")
        class_records.append(
            {
                "h_certificate_sha256": direct,
                "transpose_h_certificate_sha256": members[0][
                    "transpose_h_certificate_sha256"
                ],
                "ht_certificate_sha256": members[0][
                    "ht_certificate_sha256"
                ],
                "trace_gram_squared": members[0]["trace_gram_squared"],
                "known_h_class": direct in known_h_classes,
                "raw_representative_count": len(members),
                "representative": members[0],
            }
        )

    report = {
        "schema_version": 1,
        "engine": "order22-gradient-harvest-classification-v1",
        "order": ORDER,
        "target_absolute_determinant": EXPECTED_ABS_DETERMINANT,
        "input_dir": str(input_dir),
        "input_matrix_count": len(paths),
        "exact_determinants_checked": len(paths),
        "all_exact_determinants_on_target": True,
        "h_class_count": len(class_records),
        "ht_class_count": len(
            {
                str(record["ht_certificate_sha256"])
                for record in class_records
            }
        ),
        "known_h_class_count": sum(
            record["known_h_class"] is True for record in class_records
        ),
        "novel_h_class_count": sum(
            record["known_h_class"] is False for record in class_records
        ),
        "known_report": (
            {
                "path": str(arguments.known_report.expanduser().absolute()),
                "sha256": known_report_sha256,
            }
            if arguments.known_report
            else None
        ),
        "dependencies": {
            "pynauty": getattr(pynauty, "__version__", "unknown"),
        },
        "source": {
            "path": str(Path(__file__).absolute()),
            "sha256": sha256(Path(__file__).absolute()),
        },
        "classes": class_records,
    }
    output = arguments.output.expanduser().absolute()
    atomic_write(
        output,
        (json.dumps(report, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    print(
        f"matrices={report['input_matrix_count']} "
        f"h_classes={report['h_class_count']} "
        f"ht_classes={report['ht_class_count']} "
        f"novel_h_classes={report['novel_h_class_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
