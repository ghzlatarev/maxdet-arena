#!/usr/bin/env python3
"""Apply explicit post-run semantics/provenance to completed cube reports.

This is intentionally a one-shot research-artifact migration.  It preserves
the original LNPS aggregate byte-for-byte before correcting ambiguous field
names, reconstructs the exact pre-correction driver source by inverting the
metadata-only patch, and adds a provenance sidecar to the already completed
H2 32-cube without changing that result.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from research.fast_cube_batch import (
    atomic_json,
    atomic_write,
    sha256_bytes,
)

LNPS = REPOSITORY_ROOT / (
    "runs/direct-search/fast-principal-cube/"
    "lnps-beam-h012-20260728"
)
H2_PARTITION = REPOSITORY_ROOT / (
    "runs/direct-search/fast-principal-cube/"
    "partition32-h2-bridge-20260728"
)
LNPS_SOURCE = REPOSITORY_ROOT / "research/fast_cube_lnps.py"
PARTITION_SOURCE = REPOSITORY_ROOT / "research/fast_cube_partition32.py"
ENGINE_SOURCE = REPOSITORY_ROOT / "research/fast_principal_cube.cpp"
ENGINE_BINARY = (
    REPOSITORY_ROOT / "build/research/fast_principal_cube_lnps"
)

EXPECTED_LNPS_PRE_CORRECTION_REPORT_SHA256 = (
    "8ac48ab4e0217f3942654e7129f0cec2980cde27f1740ac11f7744f5e0655104"
)
EXPECTED_H2_REPORT_SHA256 = (
    "4f746359c520f13383e8f4ad985d952bcf509358e7dc939482a60ed2a4e8b633"
)
H2_HISTORICAL_DRIVER_SHA256 = (
    "8d1462d2b702c9cb28828fdcf59ee45fa513e5c2f5c871320484ac98e8488400"
)


def replace_once(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise RuntimeError(
            f"historical reconstruction expected one match, got "
            f"{text.count(old)}"
        )
    return text.replace(old, new, 1)


def reconstruct_lnps_run_source(current: str) -> str:
    historical = current
    historical = replace_once(
        historical,
        """    total_assignment_visits = sum(
        int(report["assignments"]) for report in cube_reports
    )
    return {
        "affine_cubes_may_overlap": True,
        "assignment_visits_are_unique": False,
""",
        """    return {
""",
    )
    historical = replace_once(
        historical,
        """        "total_assignment_visits": total_assignment_visits,
""",
        """        "total_assignments": sum(
            int(report["assignments"]) for report in cube_reports
        ),
""",
    )
    historical = replace_once(
        historical,
        """        "wall_seconds": round(wall_seconds, 6),
        "claim_boundary": [
            "Assignment counts are visits across overlapping affine cubes, "
            "not unique matrices.",
            "Reserved fingerprints include planned cubes that may remain "
            "unrun at a soft time cap.",
            "Sign-normalized identity is not H/HT equivalence.",
            "Only an independently arena-verified strict score increase is "
            "a promotion.",
        ],
""",
        """        "wall_seconds": round(wall_seconds, 6),
""",
    )
    historical = replace_once(
        historical,
        """            "planned_or_reserved_affine_fingerprints_total": len(
                seen_fingerprints
            ),
""",
        """            "unique_affine_fingerprints_total": len(seen_fingerprints),
""",
    )
    historical = replace_once(
        historical,
        """    final["affine_fingerprints_prior"] = len(prior_fingerprints)
    final["affine_fingerprints_evaluated"] = len(
        {
            str(report["fingerprint_sha256"])
            for report in cube_reports
        }
    )
    final["affine_fingerprints_planned_or_reserved_total"] = len(
        seen_fingerprints
    )
    final["planned_cubes"] = sum(
        int(report["planned_cubes"]) for report in generation_reports
    )
    final["unrun_planned_cubes"] = (
        int(final["planned_cubes"]) - len(cube_reports)
    )
""",
        """    final["affine_fingerprints_prior"] = len(prior_fingerprints)
    final["affine_fingerprints_total_seen"] = len(seen_fingerprints)
""",
    )
    historical = replace_once(
        historical,
        """                "total_assignment_visits": final[
                    "total_assignment_visits"
                ],
""",
        """                "total_assignments": final["total_assignments"],
""",
    )
    return historical


def correct_lnps() -> dict[str, str | int]:
    report_path = LNPS / "aggregate-report.json"
    original = report_path.read_bytes()
    original_hash = sha256_bytes(original)
    if original_hash != EXPECTED_LNPS_PRE_CORRECTION_REPORT_SHA256:
        raise RuntimeError(
            f"unexpected LNPS pre-correction report hash: {original_hash}"
        )
    preserved = LNPS / "aggregate-report.pre-metadata-correction.json"
    if preserved.exists():
        raise FileExistsError(preserved)
    atomic_write(preserved, original)

    current_source = LNPS_SOURCE.read_text()
    historical_source = reconstruct_lnps_run_source(current_source)
    historical_bytes = historical_source.encode("utf-8")
    historical_hash = sha256_bytes(historical_bytes)
    historical_path = LNPS / "provenance/fast_cube_lnps.run-source.py"
    atomic_write(historical_path, historical_bytes)

    report = json.loads(original)
    visits = int(report.pop("total_assignments"))
    reserved_total = int(report.pop("affine_fingerprints_total_seen"))
    evaluated_fingerprints = {
        str(record["fingerprint_sha256"])
        for generation in LNPS.glob("generation-*/cubes/*/lnps-summary.json")
        for record in [json.loads(generation.read_text())]
    }
    planned = sum(
        int(generation["planned_cubes"])
        for generation in report["generation_reports"]
    )
    if (
        visits != 14_763_950_080
        or len(evaluated_fingerprints) != 110
        or planned != 126
        or reserved_total != 326
    ):
        raise RuntimeError("LNPS correction invariants changed")
    for generation in report["generation_reports"]:
        if "unique_affine_fingerprints_total" in generation:
            generation[
                "planned_or_reserved_affine_fingerprints_total"
            ] = generation.pop("unique_affine_fingerprints_total")
    report.update(
        {
            "affine_cubes_may_overlap": True,
            "affine_fingerprints_evaluated": len(evaluated_fingerprints),
            "affine_fingerprints_planned_or_reserved_total": reserved_total,
            "assignment_visits_are_unique": False,
            "claim_boundary": [
                "The 14,763,950,080 count is assignment-visits across "
                "overlapping affine cubes, not unique matrices.",
                "Exactly 110 of 126 planned cubes were evaluated; 16 were "
                "reserved in generation manifests but unrun at the soft cap.",
                "The prior count of 326 is 200 prior fingerprints plus 126 "
                "new planned/reserved fingerprints, not evaluated cubes.",
                "Sign-normalized identity is not H/HT equivalence.",
            ],
            "planned_cubes": planned,
            "post_run_metadata_correction": {
                "correction_only": True,
                "date": "2026-07-28",
                "evaluated_cube_results_changed": False,
                "original_report": str(preserved.relative_to(REPOSITORY_ROOT)),
                "original_report_raw_sha256": original_hash,
                "reason":
                    "rename ambiguous assignment/fingerprint fields and "
                    "record source freshness",
            },
            "provenance": {
                "engine_binary": str(
                    ENGINE_BINARY.relative_to(REPOSITORY_ROOT)
                ),
                "engine_binary_raw_sha256": sha256_bytes(
                    ENGINE_BINARY.read_bytes()
                ),
                "engine_source": str(
                    ENGINE_SOURCE.relative_to(REPOSITORY_ROOT)
                ),
                "engine_source_raw_sha256": sha256_bytes(
                    ENGINE_SOURCE.read_bytes()
                ),
                "run_driver_reconstructed": True,
                "run_driver_snapshot": str(
                    historical_path.relative_to(REPOSITORY_ROOT)
                ),
                "run_driver_snapshot_raw_sha256": historical_hash,
                "semantics_corrected_driver": str(
                    LNPS_SOURCE.relative_to(REPOSITORY_ROOT)
                ),
                "semantics_corrected_driver_raw_sha256": sha256_bytes(
                    LNPS_SOURCE.read_bytes()
                ),
            },
            "total_assignment_visits": visits,
            "unrun_planned_cubes": planned - len(evaluated_fingerprints),
        }
    )
    atomic_json(report_path, report)
    return {
        "corrected_report_raw_sha256": sha256_bytes(
            report_path.read_bytes()
        ),
        "historical_driver_raw_sha256": historical_hash,
        "original_report_raw_sha256": original_hash,
    }


def add_h2_sidecar() -> dict[str, str]:
    report_path = H2_PARTITION / "aggregate-report.json"
    report_hash = sha256_bytes(report_path.read_bytes())
    if report_hash != EXPECTED_H2_REPORT_SHA256:
        raise RuntimeError(f"unexpected H2 report hash: {report_hash}")
    sidecar_path = H2_PARTITION / "provenance.json"
    if sidecar_path.exists():
        raise FileExistsError(sidecar_path)
    sidecar = {
        "aggregate_report": str(report_path.relative_to(REPOSITORY_ROOT)),
        "aggregate_report_raw_sha256": report_hash,
        "driver_at_run_raw_sha256": H2_HISTORICAL_DRIVER_SHA256,
        "driver_note":
            "The driver was generalized after this run; the historical hash "
            "was captured immediately before that metadata-only/generalizing "
            "edit.",
        "driver_now": str(PARTITION_SOURCE.relative_to(REPOSITORY_ROOT)),
        "driver_now_raw_sha256": sha256_bytes(
            PARTITION_SOURCE.read_bytes()
        ),
        "engine_binary": str(ENGINE_BINARY.relative_to(REPOSITORY_ROOT)),
        "engine_binary_raw_sha256": sha256_bytes(
            ENGINE_BINARY.read_bytes()
        ),
        "engine_source": str(ENGINE_SOURCE.relative_to(REPOSITORY_ROOT)),
        "engine_source_raw_sha256": sha256_bytes(
            ENGINE_SOURCE.read_bytes()
        ),
        "manifest_raw_sha256": sha256_bytes(
            (H2_PARTITION / "manifest.json").read_bytes()
        ),
        "result_changed": False,
        "schema_version": 1,
    }
    atomic_json(sidecar_path, sidecar)
    return {
        "sidecar": str(sidecar_path.relative_to(REPOSITORY_ROOT)),
        "sidecar_raw_sha256": sha256_bytes(sidecar_path.read_bytes()),
    }


def main() -> int:
    result = {
        "h2_partition": add_h2_sidecar(),
        "lnps": correct_lnps(),
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
