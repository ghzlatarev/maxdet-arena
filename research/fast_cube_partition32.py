#!/usr/bin/env python3
"""Exhaust one 32-entry bridge cube as 32 disjoint 27-bit leaves.

The first 27 free entries are a calibrated bridge plus transverse entries.
Five exact pair-rescue entries form the outer prefix.  If a fixed-prefix leaf
base is singular, it is deterministically rerooted inside its 27-cube;
reported masks are XOR-mapped back to the one global 32-bit assignment space.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from maxdet.exact import bareiss_determinant
from research.fast_cube_batch import (
    DIMENSION,
    FRONTIER,
    ORDER,
    Features,
    atomic_json,
    atomic_write,
    canonical_matrix_bytes,
    compute_features,
    parse_max_rss,
    read_matrix,
    sha256_bytes,
)
from research.fast_cube_lnps import (
    DEFAULT_H2,
    DEFAULT_H2_BRIDGE,
    DEFAULT_REFERENCE,
    SearchState,
    affine_fingerprint,
    apply_mask,
    make_root_state,
    matrix_identity,
    parse_coordinate_file,
    popcount,
    support_bytes,
)

Coordinate = tuple[int, int]
Matrix = tuple[tuple[int, ...], ...]
FULL_DIMENSION = 32
OUTER_DIMENSION = FULL_DIMENSION - DIMENSION
LEAF_COUNT = 1 << OUTER_DIMENSION


@dataclass(frozen=True)
class GlobalCandidate:
    absolute_determinant: int
    signed_determinant: int
    global_mask: int
    leaf_id: str
    engine_mask: int
    matrix: Matrix
    source: str


def parse_arena(output: str) -> dict[str, str]:
    fields: dict[str, str] = {}
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
            raise RuntimeError(f"arena output omitted {label}")
        fields[key] = match.group(1)
    return fields


def build_feature_state(
    label: str, path: Path, root: Path, lineage: str
) -> SearchState:
    return make_root_state(lineage, path, root)


def choose_pair_rescues(
    base_support: Sequence[Coordinate],
    start_features: Features,
    endpoint_features: Features,
) -> tuple[tuple[Coordinate, ...], list[dict[str, Any]]]:
    selected = set(base_support)
    rescues: list[Coordinate] = []
    evidence: list[dict[str, Any]] = []
    rankings = (
        ("start_pair_synergy", start_features.pair_synergy_rank),
        ("endpoint_pair_synergy", endpoint_features.pair_synergy_rank),
        ("start_pair_score", start_features.pair_score_rank),
        ("endpoint_pair_score", endpoint_features.pair_score_rank),
    )
    maximum = max(len(ranking) for _, ranking in rankings)
    for rank in range(maximum):
        for source, ranking in rankings:
            if rank >= len(ranking):
                continue
            first, second = ranking[rank]
            first_selected = first in selected
            second_selected = second in selected
            if first_selected == second_selected:
                continue
            rescue = second if first_selected else first
            partner = first if first_selected else second
            if rescue in selected:
                continue
            selected.add(rescue)
            rescues.append(rescue)
            evidence.append(
                {
                    "coordinate": [rescue[0] + 1, rescue[1] + 1],
                    "partner_in_existing_support": [
                        partner[0] + 1,
                        partner[1] + 1,
                    ],
                    "start_exact_one_flip_absolute_determinant": str(
                        start_features.single_scores[rescue]
                    ),
                    "rank_within_exact_pair_order": rank + 1,
                    "endpoint_exact_one_flip_absolute_determinant": str(
                        endpoint_features.single_scores[rescue]
                    ),
                    "source": source,
                }
            )
            if len(rescues) == OUTER_DIMENSION:
                return tuple(rescues), evidence
    raise RuntimeError("could not select five pair-rescue coordinates")


def reroot_offset(
    fixed_leaf_base: Matrix,
    inner_support: Sequence[Coordinate],
) -> tuple[int, Matrix, int]:
    determinant = bareiss_determinant(fixed_leaf_base)
    if determinant != 0:
        return 0, fixed_leaf_base, determinant
    # Low-mask order is deterministic and recorded.  A single/pair reroot is
    # overwhelmingly likely; the bounded fallback remains exact.
    offsets = [1 << index for index in range(DIMENSION)]
    offsets.extend(
        (1 << first) | (1 << second)
        for first in range(DIMENSION)
        for second in range(first + 1, DIMENSION)
    )
    for offset in offsets:
        candidate = apply_mask(fixed_leaf_base, inner_support, offset)
        candidate_determinant = bareiss_determinant(candidate)
        if candidate_determinant != 0:
            return offset, candidate, candidate_determinant
    raise RuntimeError("could not find nonsingular reroot within leaf")


def global_mask(
    outer_mask: int, engine_mask: int, reroot_xor: int
) -> int:
    inner_original = engine_mask ^ reroot_xor
    return inner_original | (outer_mask << DIMENSION)


def candidate_from_matrix(
    matrix: Matrix,
    signed_determinant: int,
    mask: int,
    leaf_id: str,
    engine_mask: int,
    source: str,
) -> GlobalCandidate:
    return GlobalCandidate(
        absolute_determinant=abs(signed_determinant),
        signed_determinant=signed_determinant,
        global_mask=mask,
        leaf_id=leaf_id,
        engine_mask=engine_mask,
        matrix=matrix,
        source=source,
    )


def verify_promotion(
    matrix: Matrix,
    output_path: Path,
    root: Path,
) -> dict[str, Any]:
    determinant = bareiss_determinant(matrix)
    if abs(determinant) <= FRONTIER:
        raise RuntimeError("promotion verifier received a non-promotion")
    atomic_write(output_path, canonical_matrix_bytes(matrix))
    verification = subprocess.run(
        ["./arena", "verify", str(output_path)],
        cwd=root,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    atomic_write(
        output_path.with_suffix(".arena-verify.txt"),
        verification.stdout.encode("utf-8"),
    )
    if verification.returncode != 0:
        raise RuntimeError("strict promotion failed arena verification")
    return {
        "arena": parse_arena(verification.stdout),
        "signed_determinant": str(determinant),
    }


def aggregate_payload(
    *,
    complete: bool,
    reason: str,
    full_fingerprint: str,
    full_support: Sequence[Coordinate],
    leaf_reports: Sequence[dict[str, Any]],
    wall_seconds: float,
    rescue_evidence: Sequence[dict[str, Any]],
    global_top: Sequence[dict[str, Any]] = (),
    global_ties: Sequence[dict[str, Any]] = (),
    promotion: dict[str, Any] | None = None,
    provenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    best = max(
        (
            int(report["best_absolute_determinant"])
            for report in leaf_reports
        ),
        default=FRONTIER,
    )
    return {
        "best_absolute_determinant": str(best),
        "complete": complete,
        "completed_leaves": len(leaf_reports),
        "dimension": FULL_DIMENSION,
        "engine_leaf_dimension": DIMENSION,
        "frontier": str(FRONTIER),
        "frontier_gain": str(best - FRONTIER),
        "full_affine_fingerprint_sha256": full_fingerprint,
        "full_support": [
            [row + 1, column + 1] for row, column in full_support
        ],
        "global_ties": list(global_ties),
        "global_top_k": list(global_top),
        "leaf_partition": {
            "disjoint": True,
            "fixed_outer_bits": OUTER_DIMENSION,
            "leaf_assignments": 1 << DIMENSION,
            "leaf_count": LEAF_COUNT,
            "mask_mapping":
                "inner_global=engine_mask XOR reroot_xor; "
                "global=inner_global OR (outer_mask<<27)",
            "total_assignments": 1 << FULL_DIMENSION,
            "union_complete": (
                complete and len(leaf_reports) == LEAF_COUNT
            ),
        },
        "leaf_reports": list(leaf_reports),
        "method": "32-cube-as-32-disjoint-rerootable-27-cube-leaves-v1",
        "pair_rescue_evidence": list(rescue_evidence),
        "planned_leaves": LEAF_COUNT,
        "promotion": promotion,
        "provenance": provenance,
        "reason": reason,
        "schema_version": 1,
        "total_assignments": sum(
            int(report["assignments"]) for report in leaf_reports
        ),
        "total_engine_seconds": round(
            sum(float(report["elapsed_seconds"]) for report in leaf_reports),
            6,
        ),
        "total_frontier_tie_masks": sum(
            int(report["frontier_tie_count"])
            for report in leaf_reports
        ),
        "total_rerooted_leaves": sum(
            int(report["reroot_xor_mask_decimal"]) != 0
            for report in leaf_reports
        ),
        "total_zero_pivot_corrections": sum(
            int(report["zero_pivot_corrections"])
            for report in leaf_reports
        ),
        "wall_seconds": round(wall_seconds, 6),
    }


def run(arguments: argparse.Namespace) -> int:
    root = REPOSITORY_ROOT
    binary = (root / arguments.binary).resolve()
    if not binary.is_file():
        raise FileNotFoundError(f"missing evaluator binary: {binary}")
    output_directory = (root / arguments.output_dir).resolve()
    if output_directory.exists():
        raise FileExistsError(f"output exists: {output_directory}")
    output_directory.mkdir(parents=True)

    start = build_feature_state(
        arguments.start_label, arguments.start, root, arguments.start_label
    )
    endpoint_path = (root / arguments.endpoint).resolve()
    endpoint_matrix = read_matrix(endpoint_path)
    endpoint_determinant = bareiss_determinant(endpoint_matrix)
    calibrated = parse_coordinate_file(root / arguments.calibrated_support)
    if len(calibrated) != DIMENSION:
        raise ValueError("calibrated support must have 27 coordinates")
    if arguments.outer_support is not None:
        rescues = parse_coordinate_file(root / arguments.outer_support)
        if len(rescues) != OUTER_DIMENSION:
            raise ValueError("explicit outer support must have 5 coordinates")
        full_support = calibrated + rescues
        if (
            len(full_support) != FULL_DIMENSION
            or len(set(full_support)) != FULL_DIMENSION
        ):
            raise ValueError(
                "explicit inner and outer supports are not 32 unique entries"
            )
        rescue_evidence = [
            {
                "coordinate": [row + 1, column + 1],
                "source": "explicit_exact_connector",
            }
            for row, column in rescues
        ]
        outer_support_selection = "explicit_exact_connector"
    else:
        endpoint = build_feature_state(
            arguments.endpoint_label,
            arguments.endpoint,
            root,
            arguments.endpoint_label,
        )
        print(f"feature-build center={arguments.start_label}", flush=True)
        start_features = compute_features(start)  # type: ignore[arg-type]
        print(
            f"feature-build center={arguments.endpoint_label}", flush=True
        )
        endpoint_features = compute_features(
            endpoint
        )  # type: ignore[arg-type]
        rescues, rescue_evidence = choose_pair_rescues(
            calibrated, start_features, endpoint_features
        )
        full_support = calibrated + rescues
        outer_support_selection = "derived_exact_pair_rescues"
    if len(full_support) != FULL_DIMENSION or len(set(full_support)) != 32:
        raise RuntimeError("full support is not 32 unique coordinates")
    bridge_mask = (1 << arguments.bridge_size) - 1
    if apply_mask(
        start.matrix,
        full_support[: arguments.bridge_size],
        bridge_mask,
    ) != endpoint_matrix:
        raise ValueError(
            "leading full-support entries do not exactly reach endpoint"
        )
    endpoint_control = {
        "absolute_determinant": str(abs(endpoint_determinant)),
        "bridge_mask_decimal": str(bridge_mask),
        "bridge_size": arguments.bridge_size,
        "exact_matrix_reconstruction": True,
        "parsed_matrix_sha256": sha256_bytes(
            canonical_matrix_bytes(endpoint_matrix)
        ),
        "signed_determinant": str(endpoint_determinant),
    }
    full_fingerprint = affine_fingerprint(start.matrix, full_support)
    driver_path = Path(__file__).resolve()
    engine_source_path = root / "research/fast_principal_cube.cpp"
    runtime_dependency_paths = (
        root / "research/fast_cube_batch.py",
        root / "research/fast_cube_lnps.py",
        root / "maxdet/exact.py",
        root / "arena",
    )
    provenance = {
        "driver": str(driver_path.relative_to(root)),
        "driver_raw_sha256": sha256_bytes(driver_path.read_bytes()),
        "engine_binary": str(binary.relative_to(root)),
        "engine_binary_raw_sha256": sha256_bytes(binary.read_bytes()),
        "engine_source": str(engine_source_path.relative_to(root)),
        "engine_source_raw_sha256": sha256_bytes(
            engine_source_path.read_bytes()
        ),
        "runtime_dependencies": [
            {
                "path": str(path.relative_to(root)),
                "raw_sha256": sha256_bytes(path.read_bytes()),
            }
            for path in runtime_dependency_paths
        ],
    }
    atomic_write(
        output_directory / "support32.coords.txt",
        support_bytes(full_support),
    )
    atomic_json(
        output_directory / "manifest.json",
        {
            "full_affine_fingerprint_sha256": full_fingerprint,
            "full_support": [
                [row + 1, column + 1] for row, column in full_support
            ],
            "inner_support_size": DIMENSION,
            "outer_support_size": OUTER_DIMENSION,
            "outer_support_selection": outer_support_selection,
            "pair_rescue_evidence": rescue_evidence,
            "planned_leaf_ids": [
                f"leaf-{outer_mask:02d}-outer-{outer_mask:05b}"
                for outer_mask in range(LEAF_COUNT)
            ],
            "provenance": provenance,
            "schema_version": 1,
            "bridge_mask_decimal": str(bridge_mask),
            "bridge_size": arguments.bridge_size,
            "endpoint": str(arguments.endpoint),
            "endpoint_control": endpoint_control,
            "endpoint_label": arguments.endpoint_label,
            "endpoint_raw_sha256": sha256_bytes(
                endpoint_path.read_bytes()
            ),
            "start": str(arguments.start),
            "start_label": arguments.start_label,
            "start_raw_sha256": start.raw_sha256,
        },
    )

    leaf_reports: list[dict[str, Any]] = []
    global_candidates: list[GlobalCandidate] = []
    tie_candidates: dict[int, GlobalCandidate] = {}
    leaf_fingerprints: set[str] = set()
    aggregate_path = output_directory / "aggregate-report.json"
    started = time.monotonic()
    reason = "complete"
    promotion: dict[str, Any] | None = None

    for outer_mask in range(LEAF_COUNT):
        if time.monotonic() - started >= arguments.maximum_seconds:
            reason = "time_limit"
            break
        leaf_id = f"leaf-{outer_mask:02d}-outer-{outer_mask:05b}"
        leaf_directory = output_directory / "leaves" / leaf_id
        leaf_directory.mkdir(parents=True)
        fixed_base = apply_mask(start.matrix, rescues, outer_mask)
        reroot_xor, leaf_start, leaf_start_det = reroot_offset(
            fixed_base, calibrated
        )
        leaf_fingerprint = affine_fingerprint(leaf_start, calibrated)
        if leaf_fingerprint in leaf_fingerprints:
            raise RuntimeError("duplicate affine leaf fingerprint")
        leaf_fingerprints.add(leaf_fingerprint)
        atomic_write(
            leaf_directory / "start.matrix.txt",
            canonical_matrix_bytes(leaf_start),
        )
        atomic_write(
            leaf_directory / "support.coords.txt",
            support_bytes(calibrated),
        )

        reroot_global_mask = global_mask(
            outer_mask, 0, reroot_xor
        )
        if abs(leaf_start_det) > FRONTIER:
            promotion = verify_promotion(
                leaf_start,
                leaf_directory / "promotion.matrix.txt",
                root,
            )
            promotion.update(
                {
                    "global_mask_decimal": str(reroot_global_mask),
                    "leaf_id": leaf_id,
                    "source": "leaf_start",
                }
            )
            reason = "strict_promotion"
            break
        if reroot_global_mask != 0 and abs(leaf_start_det) < FRONTIER:
            global_candidates.append(
                candidate_from_matrix(
                    leaf_start,
                    leaf_start_det,
                    reroot_global_mask,
                    leaf_id,
                    0,
                    "leaf_start",
                )
            )
        if abs(leaf_start_det) == FRONTIER:
            tie_candidates[reroot_global_mask] = candidate_from_matrix(
                leaf_start,
                leaf_start_det,
                reroot_global_mask,
                leaf_id,
                0,
                "leaf_start",
            )

        command = [
            str(binary),
            "--start",
            str(leaf_directory / "start.matrix.txt"),
            "--coordinates",
            str(leaf_directory / "support.coords.txt"),
            "--output",
            str(leaf_directory / "best.matrix.txt"),
            "--log",
            str(leaf_directory / "search.jsonl"),
            "--report",
            str(leaf_directory / "report.json"),
            "--top-k",
            str(arguments.top_k),
        ]
        leaf_wall_started = time.monotonic()
        result = subprocess.run(
            command,
            cwd=root,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        leaf_wrapper_wall_seconds = time.monotonic() - leaf_wall_started
        atomic_write(
            leaf_directory / "stdout.txt", result.stdout.encode("utf-8")
        )
        atomic_write(
            leaf_directory / "timing.txt", result.stderr.encode("utf-8")
        )
        if result.returncode != 0:
            atomic_json(
                leaf_directory / "failure.json",
                {
                    "command": command,
                    "returncode": result.returncode,
                    "stderr": result.stderr,
                    "stdout": result.stdout,
                },
            )
            raise RuntimeError(f"leaf engine failed: {leaf_id}")
        report = json.loads(
            (leaf_directory / "report.json").read_text()
        )
        if report.get("complete") is not True:
            raise RuntimeError(f"incomplete leaf report: {leaf_id}")
        if report["start_parsed_matrix_sha256"] != sha256_bytes(
            canonical_matrix_bytes(leaf_start)
        ):
            raise RuntimeError("leaf start provenance mismatch")

        if int(report["best_absolute_determinant"]) > FRONTIER:
            best_matrix = read_matrix(leaf_directory / "best.matrix.txt")
            best_engine_mask = int(report["best_mask_decimal"])
            best_global_mask = global_mask(
                outer_mask, best_engine_mask, reroot_xor
            )
            promotion = verify_promotion(
                best_matrix,
                leaf_directory / "promotion.matrix.txt",
                root,
            )
            promotion.update(
                {
                    "global_mask_decimal": str(best_global_mask),
                    "leaf_id": leaf_id,
                    "source": "engine_best",
                }
            )
            reason = "strict_promotion"

        for record in report["top_k_candidates"]:
            engine_mask = int(record["mask_decimal"])
            mapped = global_mask(outer_mask, engine_mask, reroot_xor)
            matrix = apply_mask(leaf_start, calibrated, engine_mask)
            global_candidates.append(
                candidate_from_matrix(
                    matrix,
                    int(record["signed_determinant"]),
                    mapped,
                    leaf_id,
                    engine_mask,
                    "engine_top_k",
                )
            )

        for mask_text in report["frontier_tie_masks_decimal"]:
            engine_mask = int(mask_text)
            mapped = global_mask(outer_mask, engine_mask, reroot_xor)
            matrix = apply_mask(leaf_start, calibrated, engine_mask)
            determinant = bareiss_determinant(matrix)
            if abs(determinant) != FRONTIER:
                raise RuntimeError("mapped leaf tie failed exact verification")
            tie_candidates[mapped] = candidate_from_matrix(
                matrix,
                determinant,
                mapped,
                leaf_id,
                engine_mask,
                "engine_frontier_tie",
            )

        leaf_summary = {
            "assignments": int(report["assignments"]),
            "best_absolute_determinant": report[
                "best_absolute_determinant"
            ],
            "best_engine_mask_decimal": report["best_mask_decimal"],
            "best_global_mask_decimal": str(
                global_mask(
                    outer_mask,
                    int(report["best_mask_decimal"]),
                    reroot_xor,
                )
            ),
            "elapsed_seconds": float(report["elapsed_seconds"]),
            "frontier_tie_count": int(report["frontier_nonzero_ties"])
                + int(abs(leaf_start_det) == FRONTIER),
            "frontier_tie_masks_truncated": bool(
                report["frontier_tie_masks_truncated"]
            ),
            "leaf_affine_fingerprint_sha256": leaf_fingerprint,
            "leaf_id": leaf_id,
            "maximum_resident_set_size": parse_max_rss(result.stderr),
            "outer_mask_binary": f"{outer_mask:05b}",
            "outer_mask_decimal": str(outer_mask),
            "reroot_xor_mask_decimal": str(reroot_xor),
            "rerooted": reroot_xor != 0,
            "start_signed_determinant": str(leaf_start_det),
            "wrapper_wall_seconds": round(
                leaf_wrapper_wall_seconds, 6
            ),
            "zero_pivot_corrections": int(
                report["zero_pivot_corrections"]
            ),
        }
        leaf_reports.append(leaf_summary)
        atomic_json(leaf_directory / "partition-summary.json", leaf_summary)
        atomic_json(
            aggregate_path,
            aggregate_payload(
                complete=False,
                reason="running",
                full_fingerprint=full_fingerprint,
                full_support=full_support,
                leaf_reports=leaf_reports,
                wall_seconds=time.monotonic() - started,
                rescue_evidence=rescue_evidence,
                promotion=promotion,
                provenance=provenance,
            ),
        )
        print(
            f"leaf={len(leaf_reports)}/{LEAF_COUNT} id={leaf_id} "
            f"reroot={reroot_xor} "
            f"best={report['best_absolute_determinant']} "
            f"seconds={report['elapsed_seconds']:.3f}",
            flush=True,
        )
        if reason == "strict_promotion":
            break

    complete = reason == "complete" and len(leaf_reports) == LEAF_COUNT
    global_top_records: list[dict[str, Any]] = []
    global_tie_records: list[dict[str, Any]] = []
    if complete:
        # Thirty-two candidates per leaf suffice for the exact global top 32:
        # any omitted candidate has at least 32 no-worse candidates in its own
        # leaf.  Leaf mask zero was inserted separately when globally nonzero.
        global_candidates.sort(
            key=lambda item: (
                -item.absolute_determinant,
                item.global_mask,
                item.leaf_id,
                item.engine_mask,
            )
        )
        seen_masks: set[int] = set()
        top: list[GlobalCandidate] = []
        for candidate in global_candidates:
            if candidate.global_mask == 0:
                continue
            if candidate.global_mask in seen_masks:
                raise RuntimeError("global assignment appeared twice")
            seen_masks.add(candidate.global_mask)
            top.append(candidate)
            if len(top) == arguments.top_k:
                break
        if len(top) != arguments.top_k:
            raise RuntimeError("insufficient global top candidates")
        top_directory = output_directory / "global-top-k"
        for rank, candidate in enumerate(top, 1):
            exact = bareiss_determinant(candidate.matrix)
            if (
                exact != candidate.signed_determinant
                or abs(exact) != candidate.absolute_determinant
            ):
                raise RuntimeError("global top candidate failed Bareiss")
            raw, normalized, _ = matrix_identity(candidate.matrix)
            artifact = (
                top_directory
                / f"rank-{rank:02d}-mask-{candidate.global_mask}.matrix.txt"
            )
            atomic_write(artifact, canonical_matrix_bytes(candidate.matrix))
            global_top_records.append(
                {
                    "absolute_determinant": str(
                        candidate.absolute_determinant
                    ),
                    "artifact": str(artifact.relative_to(root)),
                    "engine_mask_decimal": str(candidate.engine_mask),
                    "global_hamming_weight": popcount(
                        candidate.global_mask
                    ),
                    "global_mask_decimal": str(candidate.global_mask),
                    "leaf_id": candidate.leaf_id,
                    "normalized_sha256": normalized,
                    "rank": rank,
                    "raw_sha256": raw,
                    "signed_determinant": str(
                        candidate.signed_determinant
                    ),
                    "source": candidate.source,
                }
            )

        tie_directory = output_directory / "global-frontier-ties"
        for mapped, candidate in sorted(tie_candidates.items()):
            exact = bareiss_determinant(candidate.matrix)
            if abs(exact) != FRONTIER:
                raise RuntimeError("global tie failed final Bareiss")
            raw, normalized, _ = matrix_identity(candidate.matrix)
            artifact = tie_directory / f"mask-{mapped}.matrix.txt"
            atomic_write(artifact, canonical_matrix_bytes(candidate.matrix))
            verification = subprocess.run(
                ["./arena", "verify", str(artifact)],
                cwd=root,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            atomic_write(
                artifact.with_suffix(".arena-verify.txt"),
                verification.stdout.encode("utf-8"),
            )
            if verification.returncode != 0:
                raise RuntimeError("arena rejected global frontier tie")
            global_tie_records.append(
                {
                    "artifact": str(artifact.relative_to(root)),
                    "global_mask_decimal": str(mapped),
                    "leaf_id": candidate.leaf_id,
                    "normalized_sha256": normalized,
                    "raw_sha256": raw,
                    "signed_determinant": str(exact),
                    "source": candidate.source,
                }
            )

    final = aggregate_payload(
        complete=complete,
        reason=reason,
        full_fingerprint=full_fingerprint,
        full_support=full_support,
        leaf_reports=leaf_reports,
        wall_seconds=time.monotonic() - started,
        rescue_evidence=rescue_evidence,
        global_top=global_top_records,
        global_ties=global_tie_records,
        promotion=promotion,
        provenance=provenance,
    )
    final["leaf_affine_fingerprints_unique"] = len(leaf_fingerprints)
    final["global_top_k_requested"] = arguments.top_k
    final["endpoint_control"] = endpoint_control
    final["outer_support_selection"] = outer_support_selection
    final["tie_mask_truncation_leaf_count"] = sum(
        bool(report["frontier_tie_masks_truncated"])
        for report in leaf_reports
    )
    atomic_json(aggregate_path, final)
    print(
        json.dumps(
            {
                "best_absolute_determinant": final[
                    "best_absolute_determinant"
                ],
                "complete": complete,
                "completed_leaves": len(leaf_reports),
                "global_frontier_ties": len(global_tie_records),
                "reason": reason,
                "total_assignments": final["total_assignments"],
                "wall_seconds": final["wall_seconds"],
            },
            sort_keys=True,
        ),
        flush=True,
    )
    return 0


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--binary",
        type=Path,
        default=Path("build/research/fast_principal_cube_lnps"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--start", "--h2", dest="start", type=Path, default=DEFAULT_H2
    )
    parser.add_argument(
        "--calibrated-support",
        "--h2-bridge",
        dest="calibrated_support",
        type=Path,
        default=DEFAULT_H2_BRIDGE,
    )
    parser.add_argument(
        "--outer-support",
        type=Path,
        help=(
            "optional explicit 5-coordinate outer support; with this option "
            "--bridge-size may be 28..32 and the complete support must "
            "reconstruct the endpoint"
        ),
    )
    parser.add_argument(
        "--endpoint",
        "--reference",
        dest="endpoint",
        type=Path,
        default=DEFAULT_REFERENCE,
    )
    parser.add_argument("--start-label", default="H2-QUBO")
    parser.add_argument("--endpoint-label", default="reference")
    parser.add_argument("--bridge-size", type=int, default=12)
    parser.add_argument("--top-k", type=int, default=32)
    parser.add_argument("--maximum-seconds", type=float, default=150.0)
    arguments = parser.parse_args()
    if not 1 <= arguments.top_k <= 256:
        parser.error("--top-k must lie within 1..256")
    maximum_bridge_size = (
        FULL_DIMENSION
        if arguments.outer_support is not None
        else DIMENSION
    )
    if not 1 <= arguments.bridge_size <= maximum_bridge_size:
        parser.error(
            f"--bridge-size must lie within 1..{maximum_bridge_size}"
        )
    if arguments.maximum_seconds <= 0:
        parser.error("--maximum-seconds must be positive")
    return arguments


if __name__ == "__main__":
    try:
        raise SystemExit(run(parse_arguments()))
    except (OSError, ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
