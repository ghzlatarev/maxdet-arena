#!/usr/bin/env python3
"""Build and verify the retrospective radius-four/path-relink manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from maxdet import bareiss_determinant  # noqa: E402


RADIUS_DIRECTORY = ROOT / "runs/direct-search/radius4"
PATH_DIRECTORY = ROOT / "runs/direct-search/path-relink"
EXPECTED_RADIUS_COMBINATIONS = math.comb(23 * 23, 4)
RADIUS_SCREENING_SEMANTICS = (
    "all_combinations_floating_score_exact_margin_gate_"
    "no_rounding_certificate"
)

RADIUS_STARTS = {
    "class51": ROOT
    / "runs/direct-search/best-below/frontier-class51-28752.matrix.txt",
    "frontier": ROOT
    / "runs/direct-search/best-below/frontier-class14-late-28751.matrix.txt",
    "historical": ROOT
    / "runs/direct-search/reference-data/orrick-pre-april2003.matrix.txt",
    "reference": ROOT / "references/orrick-et-al-2003/matrix.txt",
}

FULL_RADIUS_CAMPAIGNS = {
    "class51": [f"class51-shard{index}" for index in range(4)],
    "frontier": [
        f"frontier-shard{index}-{29660 + index}" for index in range(4)
    ],
    "historical": [f"historical-shard{index}" for index in range(4)],
    "reference": [f"reference-shard{index}" for index in range(4)],
}

ALIGNMENT_METADATA = (
    "frontier-class14-aligned-to-near.json",
    "frontier-class14-aligned-to-class42-near.json",
)


class ManifestError(RuntimeError):
    """A retained artifact failed a consistency check."""


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def file_artifact(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    return {
        "path": relative(path),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size_bytes": len(payload),
    }


def read_matrix(path: Path) -> list[list[int]]:
    rows = [
        [int(token) for token in line.split()]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(rows) != 23 or any(len(row) != 23 for row in rows):
        raise ManifestError(f"{relative(path)} is not 23x23")
    if any(value not in (-1, 1) for row in rows for value in row):
        raise ManifestError(f"{relative(path)} contains a non-sign entry")
    return rows


def sign_bits_hex(matrix: list[list[int]]) -> str:
    bits = "".join("1" if value == 1 else "0" for row in matrix for value in row)
    bits += "0" * ((4 - len(bits) % 4) % 4)
    return "".join(f"{int(bits[index:index + 4], 2):x}" for index in range(0, len(bits), 4))


def matrix_artifact(path: Path) -> dict[str, Any]:
    matrix = read_matrix(path)
    return {
        **file_artifact(path),
        "absolute_determinant": str(abs(bareiss_determinant(matrix))),
        "row_major_sign_bits_hex": sign_bits_hex(matrix),
    }


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_bytes())
    if not isinstance(value, dict):
        raise ManifestError(f"{relative(path)} root is not an object")
    return value


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    records = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        if not line:
            continue
        value = json.loads(line)
        if not isinstance(value, dict):
            raise ManifestError(
                f"{relative(path)}:{line_number} is not an object"
            )
        records.append(value)
    if not records:
        raise ManifestError(f"{relative(path)} is empty")
    return records


def radius_start_key(stem: str) -> str:
    for key in ("class51", "historical", "reference", "frontier"):
        if stem.startswith(key):
            return key
    if stem in {"benchmark", "sanitizer", "smoke"}:
        return "frontier"
    raise ManifestError(f"no retained start mapping for radius run {stem}")


def radius_run(snapshot_path: Path) -> dict[str, Any]:
    suffix = ".snapshot.json"
    stem = snapshot_path.name[: -len(suffix)]
    log_path = RADIUS_DIRECTORY / f"{stem}.jsonl"
    output_path = RADIUS_DIRECTORY / f"{stem}.matrix.txt"
    research_path = RADIUS_DIRECTORY / f"{stem}-research.matrix.txt"
    start_path = RADIUS_STARTS[radius_start_key(stem)]

    snapshot = load_json(snapshot_path)
    log_records = load_jsonl(log_path)
    if log_records[-1] != snapshot:
        raise ManifestError(f"{stem}: final log record differs from snapshot")
    if snapshot.get("event") != "finished":
        raise ManifestError(f"{stem}: snapshot is not a finished event")
    if snapshot.get("screening_semantics") != RADIUS_SCREENING_SEMANTICS:
        raise ManifestError(f"{stem}: screening semantics changed")
    if snapshot.get("screened_combinations") != snapshot.get(
        "expected_combinations"
    ):
        raise ManifestError(f"{stem}: incomplete enumeration count")

    start = matrix_artifact(start_path)
    output = matrix_artifact(output_path)
    research = matrix_artifact(research_path)
    if output["absolute_determinant"] != snapshot.get("best_score"):
        raise ManifestError(f"{stem}: output score differs from snapshot")
    if research["absolute_determinant"] != snapshot.get(
        "research_best_score"
    ):
        raise ManifestError(f"{stem}: research score differs from snapshot")
    if snapshot.get("promotions") != 0:
        raise ManifestError(
            f"{stem}: retrospective start binding assumes zero promotions"
        )
    if output["sha256"] != start["sha256"]:
        raise ManifestError(
            f"{stem}: zero-promotion output bytes differ from intended start"
        )

    return {
        "run_id": stem,
        "semantic_start": start,
        "artifacts": {
            "snapshot": file_artifact(snapshot_path),
            "jsonl": file_artifact(log_path),
            "output": output,
            "research_output": research,
        },
        "evidence": {
            "snapshot_equals_final_log_record": True,
            "output_bytes_equal_intended_start": True,
            "promotions": snapshot["promotions"],
        },
        "telemetry": {
            key: snapshot[key]
            for key in (
                "complete",
                "termination",
                "seed",
                "shard_count",
                "shard_index",
                "expected_combinations",
                "screened_combinations",
                "top_pool",
                "calibration_samples",
                "calibration_maximum_absolute_error",
                "exact_margin",
                "near_gate_candidates",
                "near_gate_exact_checks",
                "final_pool_exact_checks",
                "best_score",
                "research_best_score",
                "elapsed_seconds",
            )
        },
    }


def full_radius_campaigns(
    runs_by_id: dict[str, dict[str, Any]]
) -> list[dict[str, Any]]:
    result = []
    for campaign_id, run_ids in FULL_RADIUS_CAMPAIGNS.items():
        runs = [runs_by_id[run_id] for run_id in run_ids]
        shard_counts = {run["telemetry"]["shard_count"] for run in runs}
        shard_indices = sorted(run["telemetry"]["shard_index"] for run in runs)
        expected = sum(
            run["telemetry"]["expected_combinations"] for run in runs
        )
        screened = sum(
            run["telemetry"]["screened_combinations"] for run in runs
        )
        starts = {run["semantic_start"]["sha256"] for run in runs}
        if (
            shard_counts != {4}
            or shard_indices != [0, 1, 2, 3]
            or expected != EXPECTED_RADIUS_COMBINATIONS
            or screened != EXPECTED_RADIUS_COMBINATIONS
            or len(starts) != 1
            or not all(run["telemetry"]["complete"] for run in runs)
        ):
            raise ManifestError(f"{campaign_id}: invalid four-shard coverage")
        result.append(
            {
                "campaign_id": campaign_id,
                "run_ids": run_ids,
                "semantic_start_sha256": next(iter(starts)),
                "expected_combinations": expected,
                "screened_combinations": screened,
                "complete_floating_screen": True,
                "exact_audit": False,
            }
        )
    return result


def path_run(snapshot_path: Path) -> dict[str, Any]:
    suffix = ".snapshot.json"
    stem = snapshot_path.name[: -len(suffix)]
    log_path = PATH_DIRECTORY / f"{stem}.jsonl"
    output_path = PATH_DIRECTORY / f"{stem}.matrix.txt"
    snapshot = load_json(snapshot_path)
    records = load_jsonl(log_path)
    start = records[0]
    if records[-1] != snapshot:
        raise ManifestError(f"{stem}: final log record differs from snapshot")
    if start.get("event") != "start" or snapshot.get("event") != "finished":
        raise ManifestError(f"{stem}: malformed start/final events")

    elites = []
    for elite in start.get("elites", []):
        if not isinstance(elite, dict) or not isinstance(elite.get("path"), str):
            raise ManifestError(f"{stem}: malformed elite record")
        artifact = matrix_artifact(ROOT / elite["path"])
        if artifact["absolute_determinant"] != elite.get(
            "absolute_determinant"
        ):
            raise ManifestError(f"{stem}: elite determinant changed")
        elites.append(artifact)
    if len(elites) != start.get("elite_count"):
        raise ManifestError(f"{stem}: elite count differs from start record")

    output = matrix_artifact(output_path)
    if output["absolute_determinant"] != snapshot.get("best_score"):
        raise ManifestError(f"{stem}: output score differs from final snapshot")
    return {
        "run_id": stem,
        "semantics": "pruned_heuristic_not_exhaustive",
        "elite_inputs": elites,
        "artifacts": {
            "snapshot": file_artifact(snapshot_path),
            "jsonl": file_artifact(log_path),
            "output": output,
        },
        "evidence": {"snapshot_equals_final_log_record": True},
        "parameters": {
            key: start.get(key)
            for key in (
                "seed",
                "beam_width",
                "exact_pool",
                "random_fraction",
                "bidirectional",
                "elite_count",
            )
        },
        "telemetry": {
            key: snapshot.get(key)
            for key in (
                "complete",
                "termination",
                "generated",
                "unique_candidates",
                "exact_checks",
                "singular_candidates",
                "retained_states",
                "completed_depths",
                "completed_pairs",
                "promotions",
                "best_score",
                "elapsed_seconds",
            )
        },
    }


def replay_alignment(
    source: list[list[int]], metadata: dict[str, Any]
) -> list[list[int]]:
    row_permutation = metadata["row_permutation_zero_based"]
    row_signs = metadata["row_signs"]
    column_permutation = metadata["column_permutation_zero_based"]
    column_signs = metadata["column_signs"]
    return [
        [
            row_signs[row]
            * column_signs[column]
            * source[row_permutation[row]][column_permutation[column]]
            for column in range(23)
        ]
        for row in range(23)
    ]


def alignment_record(metadata_path: Path) -> dict[str, Any]:
    metadata = load_json(metadata_path)
    target_path = ROOT / metadata["target"]
    source_path = ROOT / metadata["source"]
    output_path = ROOT / metadata["output"]
    target_matrix = read_matrix(target_path)
    source_matrix = read_matrix(source_path)
    output_matrix = read_matrix(output_path)
    if replay_alignment(source_matrix, metadata) != output_matrix:
        raise ManifestError(f"{metadata_path.name}: transform replay failed")
    hamming = sum(
        left != right
        for target_row, output_row in zip(target_matrix, output_matrix)
        for left, right in zip(target_row, output_row)
    )
    if hamming != metadata["hamming_distance"]:
        raise ManifestError(f"{metadata_path.name}: Hamming distance changed")
    return {
        "record_id": metadata_path.stem,
        "artifacts": {
            "metadata": file_artifact(metadata_path),
            "target": matrix_artifact(target_path),
            "source": matrix_artifact(source_path),
            "output": matrix_artifact(output_path),
        },
        "evidence": {
            "signed_permutation_replay": True,
            "hamming_distance": hamming,
            "determinant_preserved": (
                abs(bareiss_determinant(source_matrix))
                == abs(bareiss_determinant(output_matrix))
            ),
        },
    }


def build_manifest() -> dict[str, Any]:
    radius_runs = [
        radius_run(path)
        for path in sorted(RADIUS_DIRECTORY.glob("*.snapshot.json"))
    ]
    radius_by_id = {run["run_id"]: run for run in radius_runs}
    path_runs = [
        path_run(path)
        for path in sorted(PATH_DIRECTORY.glob("*.snapshot.json"))
    ]
    alignments = [
        alignment_record(PATH_DIRECTORY / name) for name in ALIGNMENT_METADATA
    ]
    return {
        "schema_version": 1,
        "manifest_id": "direct-search-radius4-path-20260728",
        "challenge_id": "maxdet-23-v1",
        "claim_boundary": (
            "This is a retrospective digest binding of currently retained "
            "files. Pre-hardening radius snapshots did not record their start "
            "hashes, and path logs recorded elite paths and exact scores but "
            "not hashes. Matching current bytes, determinants, zero-promotion "
            "radius outputs, and log/snapshot records strongly reconcile the "
            "artifacts but cannot prove that an input file was unchanged "
            "between historical execution and creation of this manifest."
        ),
        "radius4_semantics": {
            "classification": "complete_floating_screen_not_exact_audit",
            "screening_semantics": RADIUS_SCREENING_SEMANTICS,
            "exact_promotion": "Bareiss",
            "rounding_certificate": False,
        },
        "path_relink_semantics": {
            "classification": "pruned_heuristic_not_exhaustive",
            "exact_promotion": "Bareiss",
            "numerically_singular_states_may_be_pruned": True,
        },
        "full_radius4_campaigns": full_radius_campaigns(radius_by_id),
        "radius4_runs": radius_runs,
        "path_relink_runs": path_runs,
        "alignment_records": alignments,
    }


def canonical_bytes(value: dict[str, Any]) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n"
    ).encode("ascii")


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(payload)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--output", type=Path)
    mode.add_argument("--check", type=Path)
    arguments = parser.parse_args()
    payload = canonical_bytes(build_manifest())
    digest = hashlib.sha256(payload).hexdigest()
    if arguments.output is not None:
        output = arguments.output.resolve()
        atomic_write(output, payload)
        atomic_write(
            output.with_suffix(output.suffix + ".sha256"),
            f"{digest}  {output.name}\n".encode("ascii"),
        )
        print(f"{digest}  {output}")
        return 0

    assert arguments.check is not None
    actual = arguments.check.read_bytes()
    if actual != payload:
        raise ManifestError(
            f"{arguments.check} differs from a fresh validated manifest"
        )
    sidecar = arguments.check.with_suffix(arguments.check.suffix + ".sha256")
    expected_sidecar = f"{digest}  {arguments.check.name}\n".encode("ascii")
    if sidecar.read_bytes() != expected_sidecar:
        raise ManifestError(f"{sidecar} does not bind the manifest bytes")
    print(f"VALID {digest}  {arguments.check}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
