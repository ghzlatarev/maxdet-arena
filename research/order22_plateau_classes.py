#!/usr/bin/env python3
"""Enumerate one-flip-connected H-classes on the order-22 maxdet plateau.

For each discovered H-equivalence class, one representative is sufficient:
signed row/column permutations map entry flips bijectively and preserve
absolute determinant.  Exact symbolic inverse entries identify every
single-entry flip that remains on the target determinant plateau.  Each new
neighbor is independently checked by Bareiss elimination and classified by
the complete pivoted-dephased pynauty certificate.

The completed result closes only the plateau components reached from the
supplied seeds; it does not prove that no disconnected components exist.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
import hashlib
import json
import os
import tempfile
import time
from pathlib import Path
from typing import Any

import pynauty
import sympy

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


def matrix_bytes(matrix: list[list[int]]) -> bytes:
    return "".join(" ".join(map(str, row)) + "\n" for row in matrix).encode(
        "ascii"
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


def certificate_sha256(matrix: list[list[int]]) -> str:
    return sha256_bytes(h_certificate(matrix))


def exact_plateau_neighbors(
    matrix: list[list[int]],
    signed_determinant: int,
) -> list[tuple[int, int, list[list[int]], int]]:
    """Return all one-entry flips preserving the target determinant magnitude."""

    symbolic = sympy.Matrix(matrix)
    symbolic_determinant = int(symbolic.det(method="domain-ge"))
    if symbolic_determinant != signed_determinant:
        raise ArithmeticError("SymPy and Bareiss determinant checks disagree")
    inverse = symbolic.inv(method="DM")
    if symbolic * inverse != sympy.eye(ORDER):
        raise ArithmeticError("exact symbolic inverse identity failed")

    neighbors = []
    for row in range(ORDER):
        for column in range(ORDER):
            ratio = (
                1
                - 2
                * matrix[row][column]
                * inverse[column, row]
            )
            proposed = sympy.Integer(signed_determinant) * ratio
            if proposed.q != 1:
                raise ArithmeticError("entry-flip determinant was not integral")
            proposed_determinant = int(proposed)
            if abs(proposed_determinant) != EXPECTED_ABS_DETERMINANT:
                continue
            neighbor = [source_row[:] for source_row in matrix]
            neighbor[row][column] = -neighbor[row][column]
            checked = determinant(neighbor)
            if checked != proposed_determinant:
                raise ArithmeticError("entry-flip determinant lemma check failed")
            neighbors.append((row, column, neighbor, checked))
    return neighbors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", action="append", required=True, type=Path)
    parser.add_argument("--include-seed-transposes", action="store_true")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--max-classes",
        type=int,
        help="diagnostic stop; omit to run each seeded component to closure",
    )
    arguments = parser.parse_args()
    if arguments.max_classes is not None and arguments.max_classes <= 0:
        parser.error("--max-classes must be positive")

    output_dir = arguments.output_dir.expanduser().absolute()
    if output_dir.exists() and any(output_dir.iterdir()):
        parser.error("--output-dir must be absent or empty")
    output_dir.mkdir(parents=True, exist_ok=True)
    classes_dir = output_dir / "classes"
    classes_dir.mkdir()
    report_path = output_dir / "report.json"

    started = time.monotonic()
    seed_records: list[dict[str, Any]] = []
    initial: list[tuple[list[list[int]], dict[str, Any]]] = []
    for raw_path in arguments.seed:
        path = raw_path.expanduser().absolute()
        matrix = read_matrix(path)
        checked = determinant(matrix)
        if abs(checked) != EXPECTED_ABS_DETERMINANT:
            raise ValueError(f"{path}: unexpected determinant {checked}")
        record = {
            "path": str(path),
            "sha256": sha256(path),
            "determinant": checked,
            "orientation": "direct",
        }
        seed_records.append(record)
        initial.append((matrix, record))
        if arguments.include_seed_transposes:
            transposed = transpose(matrix)
            transposed_record = {
                "path": str(path),
                "sha256": sha256(path),
                "determinant": determinant(transposed),
                "orientation": "transpose",
            }
            seed_records.append(transposed_record)
            initial.append((transposed, transposed_record))

    classes: dict[str, dict[str, Any]] = {}
    matrices: dict[str, list[list[int]]] = {}
    queue: deque[str] = deque()

    def retain(
        matrix: list[list[int]],
        signed_determinant: int,
        discovery: dict[str, Any],
    ) -> tuple[str, bool]:
        direct = certificate_sha256(matrix)
        if direct in classes:
            return direct, False
        transposed = certificate_sha256(transpose(matrix))
        index = len(classes)
        path = classes_dir / f"class-{index:03d}-{direct[:12]}.matrix.txt"
        contents = matrix_bytes(matrix)
        atomic_write(path, contents)
        record = {
            "index": index,
            "path": str(path),
            "raw_sha256": sha256_bytes(contents),
            "determinant": signed_determinant,
            "h_certificate_sha256": direct,
            "transpose_h_certificate_sha256": transposed,
            "ht_certificate_sha256": min(direct, transposed),
            "discovery": discovery,
            "explored": False,
        }
        classes[direct] = record
        matrices[direct] = matrix
        queue.append(direct)
        return direct, True

    for matrix, seed_record in initial:
        retain(
            matrix,
            determinant(matrix),
            {"kind": "seed", **seed_record},
        )

    adjacency: list[dict[str, Any]] = []
    neutral_labeled_neighbors = 0
    stopped_at_limit = False
    while queue:
        if (
            arguments.max_classes is not None
            and len(classes) >= arguments.max_classes
        ):
            stopped_at_limit = True
            break
        source_certificate = queue.popleft()
        source_record = classes[source_certificate]
        source_matrix = matrices[source_certificate]
        signed_determinant = determinant(source_matrix)
        neighbors = exact_plateau_neighbors(source_matrix, signed_determinant)
        neutral_labeled_neighbors += len(neighbors)
        neighbor_counts: Counter[str] = Counter()
        discoveries = 0
        for row, column, neighbor, neighbor_determinant in neighbors:
            target_certificate, novel = retain(
                neighbor,
                neighbor_determinant,
                {
                    "kind": "one-entry-flip",
                    "parent_h_certificate_sha256": source_certificate,
                    "row_one_based": row + 1,
                    "column_one_based": column + 1,
                },
            )
            neighbor_counts[target_certificate] += 1
            discoveries += int(novel)
        source_record["explored"] = True
        source_record["neutral_entry_flip_count"] = len(neighbors)
        adjacency.append(
            {
                "source_h_certificate_sha256": source_certificate,
                "neutral_entry_flip_count": len(neighbors),
                "novel_class_discoveries": discoveries,
                "neighbor_h_class_multiplicities": dict(
                    sorted(neighbor_counts.items())
                ),
            }
        )
        print(
            f"explored={len(adjacency)} discovered={len(classes)} "
            f"queue={len(queue)} neutral={len(neighbors)} "
            f"elapsed={time.monotonic() - started:.3f}s",
            flush=True,
        )

    ordered_classes = sorted(classes.values(), key=lambda item: item["index"])
    report = {
        "schema_version": 1,
        "engine": "order22-maxdet-plateau-h-class-bfs-v1",
        "claim": (
            "Every one-entry maxdet neighbor of one representative from each "
            "reached H-class was enumerated exactly. H-equivalence preserves "
            "the neighbor-class multiset."
        ),
        "claim_boundary": (
            "Closure applies only to components reached from the bound seeds; "
            "disconnected maxdet components may exist."
        ),
        "order": ORDER,
        "target_absolute_determinant": EXPECTED_ABS_DETERMINANT,
        "complete_seeded_component_closure": not queue and not stopped_at_limit,
        "stopped_at_class_limit": stopped_at_limit,
        "seed_transposes_included": arguments.include_seed_transposes,
        "seeds": seed_records,
        "dependencies": {
            "pynauty": getattr(pynauty, "__version__", "unknown"),
            "sympy": getattr(sympy, "__version__", "unknown"),
        },
        "source": {
            "path": str(Path(__file__).absolute()),
            "sha256": sha256(Path(__file__).absolute()),
        },
        "h_class_count": len(classes),
        "ht_class_count": len(
            {str(record["ht_certificate_sha256"]) for record in ordered_classes}
        ),
        "classes_explored": sum(
            record["explored"] is True for record in ordered_classes
        ),
        "entry_flip_determinants_checked": len(adjacency) * ORDER * ORDER,
        "neutral_labeled_neighbors": neutral_labeled_neighbors,
        "classes": ordered_classes,
        "adjacency": adjacency,
        "elapsed_seconds": time.monotonic() - started,
    }
    atomic_write(
        report_path,
        (json.dumps(report, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    print(
        f"complete={report['complete_seeded_component_closure']} "
        f"h_classes={report['h_class_count']} "
        f"ht_classes={report['ht_class_count']} "
        f"elapsed={report['elapsed_seconds']:.3f}s",
        flush=True,
    )
    return 0 if report["complete_seeded_component_closure"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
