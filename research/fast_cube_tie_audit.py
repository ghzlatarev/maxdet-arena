#!/usr/bin/env python3
"""Rerun tie-containing exact cubes and harvest bounded frontier masks."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from maxdet.exact import bareiss_determinant

from fast_cube_batch import (
    FRONTIER,
    atomic_json,
    atomic_write,
    read_matrix,
)


def matrix_bytes(matrix: list[list[int]]) -> bytes:
    return (
        "".join(" ".join(str(value) for value in row) + "\n" for row in matrix)
    ).encode("ascii")


def apply_mask(
    start: tuple[tuple[int, ...], ...],
    support_one_based: list[list[int]],
    mask: int,
) -> list[list[int]]:
    candidate = [list(row) for row in start]
    for index, (row, column) in enumerate(support_one_based):
        if mask & (1 << index):
            candidate[row - 1][column - 1] *= -1
    return candidate


def parse_verification(output: str) -> dict[str, str]:
    fields = {}
    for key, label in (
        ("determinant", "determinant"),
        ("matrix_sha256", "matrix sha256"),
        ("normalized_sha256", "normalized sha256"),
        ("receipt_sha256", "receipt sha256"),
    ):
        match = re.search(
            rf"^{re.escape(label)}:\s*(\S+)\s*$",
            output,
            flags=re.MULTILINE,
        )
        if not match:
            raise RuntimeError(f"arena verification omitted {label}")
        fields[key] = match.group(1)
    return fields


def audit_report(
    planned_cubes: int,
    completed_cubes: list[dict[str, Any]],
    ties_by_hash: dict[str, dict[str, Any]],
    *,
    complete: bool,
    wall_seconds: float,
) -> dict[str, Any]:
    captured_masks = sum(
        record["captured_nonzero_frontier_masks"]
        for record in completed_cubes
    )
    reported_masks = sum(
        record["reported_nonzero_frontier_ties"]
        for record in completed_cubes
    )
    known_union_ties = sum(
        bool(record["known_neutral_union"])
        for record in ties_by_hash.values()
    )
    return {
        "captured_nonzero_frontier_masks": captured_masks,
        "complete": complete,
        "completed_tie_cubes": len(completed_cubes),
        "engine": "fast-principal-minor-entry-cube-v1",
        "frontier_floor": str(FRONTIER),
        "method": "bounded-fast-cube-frontier-tie-audit-v1",
        "planned_tie_cubes": planned_cubes,
        "reported_nonzero_frontier_ties": reported_masks,
        "schema_version": 1,
        "tie_cubes": completed_cubes,
        "tie_masks_truncated": any(
            record["tie_masks_truncated"] for record in completed_cubes
        ),
        "total_rerun_engine_seconds": round(
            sum(
                record["engine_elapsed_seconds"]
                for record in completed_cubes
            ),
            6,
        ),
        "unique_ties_beyond_known_neutral_union":
            len(ties_by_hash) - known_union_ties,
        "unique_ties_known_neutral_union": known_union_ties,
        "unique_raw_frontier_matrices": len(ties_by_hash),
        "unique_ties": [
            ties_by_hash[digest] for digest in sorted(ties_by_hash)
        ],
        "wall_seconds": round(wall_seconds, 6),
    }


def run_audit(arguments: argparse.Namespace) -> int:
    root = REPOSITORY_ROOT
    batch_directory = (root / arguments.batch_dir).resolve()
    binary = (root / arguments.binary).resolve()
    manifest = json.loads((batch_directory / "manifest.json").read_text())
    aggregate = json.loads(
        (batch_directory / "aggregate-report.json").read_text()
    )
    if aggregate.get("complete") is not True:
        raise ValueError("batch aggregate must be complete before tie audit")
    manifest_runs = {record["id"]: record for record in manifest["runs"]}
    tie_runs = [
        record for record in aggregate["runs"] if int(record["best_ties"]) > 1
    ]
    audit_directory = batch_directory / "tie-audit"
    rerun_root = audit_directory / "reruns"
    tie_root = audit_directory / "ties"
    report_path = audit_directory / "report.json"
    known_start_hashes: dict[str, list[str]] = {}
    for record in manifest["starts"]:
        known_start_hashes.setdefault(record["raw_sha256"], []).append(
            record["label"]
        )
    known_frontier_report = json.loads(
        (root / arguments.known_frontier_report).read_text()
    )
    known_neutral_union_hashes = {
        record["raw_sha256"]
        for record in known_frontier_report["artifacts"]
        if record.get("kind") == "tie"
    }

    completed_cubes: list[dict[str, Any]] = []
    ties_by_hash: dict[str, dict[str, Any]] = {}
    started = time.monotonic()
    for summary in tie_runs:
        plan = manifest_runs[summary["id"]]
        run_directory = rerun_root / summary["id"]
        run_directory.mkdir(parents=True, exist_ok=True)
        output_path = run_directory / "best.matrix.txt"
        tie_output_path = run_directory / "first-tie.matrix.txt"
        log_path = run_directory / "search.jsonl"
        engine_report_path = run_directory / "report.json"
        engine_report: dict[str, Any] | None = None
        if engine_report_path.exists():
            possible = json.loads(engine_report_path.read_text())
            if possible.get("complete") is True:
                engine_report = possible
        if engine_report is None:
            command = [
                str(binary),
                "--start",
                str(root / plan["start_path"]),
                "--coordinates",
                str(
                    batch_directory
                    / plan["id"]
                    / "support.coords.txt"
                ),
                "--output",
                str(output_path),
                "--tie-output",
                str(tie_output_path),
                "--log",
                str(log_path),
                "--report",
                str(engine_report_path),
            ]
            result = subprocess.run(
                command,
                cwd=root,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            atomic_write(
                run_directory / "stdout.txt",
                result.stdout.encode("utf-8"),
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"tie rerun failed for {plan['id']}: {result.stdout}"
                )
            engine_report = json.loads(engine_report_path.read_text())

        if int(engine_report["best_absolute_determinant"]) != FRONTIER:
            raise RuntimeError("tie rerun changed the exact cube maximum")
        expected_nonzero = int(summary["best_ties"]) - 1
        reported_nonzero = int(engine_report["frontier_nonzero_ties"])
        if reported_nonzero != expected_nonzero:
            raise RuntimeError(
                f"tie count changed for {plan['id']}: "
                f"{reported_nonzero} != {expected_nonzero}"
            )
        masks = [
            int(mask)
            for mask in engine_report["frontier_tie_masks_decimal"]
        ]
        start_matrix = read_matrix(root / plan["start_path"])
        captured_hashes: list[str] = []
        for mask in masks:
            candidate = apply_mask(start_matrix, plan["support"], mask)
            determinant = bareiss_determinant(candidate)
            if abs(determinant) != FRONTIER:
                raise RuntimeError(
                    f"captured mask failed Bareiss: {plan['id']}/{mask}"
                )
            encoded = matrix_bytes(candidate)
            raw_hash = hashlib.sha256(encoded).hexdigest()
            captured_hashes.append(raw_hash)
            artifact_path = tie_root / f"{raw_hash}.matrix.txt"
            if artifact_path.exists():
                if artifact_path.read_bytes() != encoded:
                    raise RuntimeError("SHA-256 artifact collision")
            else:
                atomic_write(artifact_path, encoded)
            if raw_hash not in ties_by_hash:
                verification = subprocess.run(
                    ["./arena", "verify", str(artifact_path)],
                    cwd=root,
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                atomic_write(
                    tie_root / f"{raw_hash}.arena-verify.txt",
                    verification.stdout.encode("utf-8"),
                )
                if verification.returncode != 0:
                    raise RuntimeError(
                        f"arena rejected captured tie {raw_hash}"
                    )
                verification_fields = parse_verification(
                    verification.stdout
                )
                if verification_fields["matrix_sha256"] != raw_hash:
                    raise RuntimeError("arena raw hash disagrees")
                ties_by_hash[raw_hash] = {
                    "absolute_determinant": str(FRONTIER),
                    "artifact": str(
                        artifact_path.relative_to(root)
                    ),
                    "matches_batch_start_labels":
                        known_start_hashes.get(raw_hash, []),
                    "known_neutral_union":
                        raw_hash in known_neutral_union_hashes,
                    "normalized_sha256": verification_fields[
                        "normalized_sha256"
                    ],
                    "raw_sha256": raw_hash,
                    "receipt_sha256": verification_fields[
                        "receipt_sha256"
                    ],
                    "sources": [],
                }
            ties_by_hash[raw_hash]["sources"].append(
                {"cube_id": plan["id"], "mask_decimal": str(mask)}
            )

        first_tie_hash = hashlib.sha256(
            tie_output_path.read_bytes()
        ).hexdigest()
        if not captured_hashes or first_tie_hash != captured_hashes[0]:
            raise RuntimeError(
                f"first tie artifact mismatch for {plan['id']}"
            )
        cube_record = {
            "captured_nonzero_frontier_masks": len(masks),
            "cube_id": plan["id"],
            "engine_elapsed_seconds": engine_report["elapsed_seconds"],
            "first_tie_raw_sha256": first_tie_hash,
            "h_class_side": plan["h_class_side"],
            "reported_nonzero_frontier_ties": reported_nonzero,
            "strategy": plan["strategy"],
            "tie_masks_truncated": bool(
                engine_report["frontier_tie_masks_truncated"]
            ),
            "unique_raw_hashes_in_cube": len(set(captured_hashes)),
        }
        completed_cubes.append(cube_record)
        partial = audit_report(
            len(tie_runs),
            completed_cubes,
            ties_by_hash,
            complete=False,
            wall_seconds=time.monotonic() - started,
        )
        atomic_json(report_path, partial)
        print(
            f"tie-audit={len(completed_cubes)}/{len(tie_runs)} "
            f"cube={plan['id']} masks={len(masks)} "
            f"unique={len(ties_by_hash)}",
            flush=True,
        )

    final = audit_report(
        len(tie_runs),
        completed_cubes,
        ties_by_hash,
        complete=True,
        wall_seconds=time.monotonic() - started,
    )
    atomic_json(report_path, final)
    print(
        json.dumps(
            {
                "captured_masks": final[
                    "captured_nonzero_frontier_masks"
                ],
                "complete": True,
                "tie_cubes": len(tie_runs),
                "unique_raw_frontier_matrices": len(ties_by_hash),
                "wall_seconds": final["wall_seconds"],
            },
            sort_keys=True,
        ),
        flush=True,
    )
    return 0


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-dir", type=Path, required=True)
    parser.add_argument(
        "--binary",
        type=Path,
        default=Path("build/research/fast_principal_cube"),
    )
    parser.add_argument(
        "--known-frontier-report",
        type=Path,
        default=Path(
            "runs/direct-search/neutral-cycle/"
            "two-cycle-union-29952/report.json"
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    try:
        raise SystemExit(run_audit(parse_arguments()))
    except (OSError, ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
