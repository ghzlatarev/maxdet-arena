#!/usr/bin/env python3
"""Audit and exhaust the two retained low-dimensional QD elite bridges.

This is a deliberately narrow, immutable campaign:

* elite 018 -> elite 022 is the exact raw 28-entry bridge;
* elite 018 -> elite 024 is a raw 22-entry bridge, expanded to one
  deterministic 27-entry search cube with five exact-score-ranked transverse
  coordinates.

Before a plan is written, all three affine objects (raw 28, raw 22, and
repaired 27) are compared exactly in the dephased 22x22 GF(2) space against
the 312 retained principal-cube searches, the 49 logical cubes in QD campaign
seed 36002, and the aligned elite004->elite021 connector.  The evaluator then
exhausts the two search cubes as three disjoint 27-bit tasks.  Every retained
best, top-K record, frontier tie, and endpoint control is replayed with exact
Bareiss arithmetic; bests and endpoints are also passed through ``./arena``.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import json
import re
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Iterable, Sequence

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from maxdet.exact import bareiss_determinant
import research.affine_gf2_audit as affine_audit
import research.fast_cube_connector_provenance as connector_provenance
import research.qd_cube_campaign as qd_campaign
from research.fast_cube_batch import (
    DIMENSION,
    FRONTIER,
    ORDER,
    atomic_json,
    atomic_write,
    canonical_matrix_bytes,
    compute_features,
    read_matrix,
    sha256_bytes,
)
from research.fast_cube_lnps import (
    affine_fingerprint,
    apply_mask,
    matrix_identity,
    support_bytes,
)

# ``pair_endpoint_rank`` belongs to the QD planner.  Keeping this explicit
# avoids silently changing the repair ranking if another feature helper grows
# a similarly named function.
pair_endpoint_rank = qd_campaign.pair_endpoint_rank
FeatureCenter = qd_campaign.FeatureCenter

Coordinate = tuple[int, int]
Matrix = tuple[tuple[int, ...], ...]

SCHEMA_VERSION = 1
METHOD = "remaining-qd-bridges-exact-campaign-v1"
TOP_K = 32
ASSIGNMENTS_PER_TASK = 1 << DIMENSION
EXPECTED_PRIOR_COUNTS = {
    "batch200": 200,
    "connector004021": 1,
    "deep32": 2,
    "lnps": 110,
    "qd36002": 49,
}
EXPECTED_PRIOR_TOTAL = sum(EXPECTED_PRIOR_COUNTS.values())
INVENTORY_TAG = b"maxdet-remaining-qd-bridge-prior-inventory-v1\0"
COMPARISON_TAG = b"maxdet-remaining-qd-bridge-comparisons-v1\0"
ARTIFACT_INVENTORY_TAG = (
    b"maxdet-remaining-qd-bridge-artifact-inventory-v1\0"
)

ELITE018 = Path(
    "runs/map-elites-pilot-20260728-seed35011/archive/"
    "elite-018-q650109375-n1-d6-b8-g2.matrix.txt"
)
ELITE022 = Path(
    "runs/map-elites-pilot-20260728-seed35011/archive/"
    "elite-022-q646078125-n1-d7-b8-g2.matrix.txt"
)
ELITE024 = Path(
    "runs/map-elites-pilot-20260728-seed35011/archive/"
    "elite-024-q637265625-n1-d5-b8-g2.matrix.txt"
)
QD_PRIOR_DIRECTORY = Path("runs/qd-selected-cubes-20260728-seed36002")
CONNECTOR_RUN_DIRECTORY = Path(
    "runs/direct-search/fast-principal-cube/"
    "qd-aligned-connector32-elite004-elite021-20260728"
)
CONNECTOR_INPUT_DIRECTORY = Path(
    "runs/direct-search/fast-principal-cube/"
    "qd-aligned-connector32-elite004-elite021-20260728-inputs"
)
CONNECTOR_TARGET = Path(
    "runs/map-elites-pilot-20260728-seed35011/archive/"
    "elite-004-q662671875-n0-d9-b8-g2.matrix.txt"
)
CONNECTOR_SOURCE = Path(
    "runs/map-elites-pilot-20260728-seed35011/archive/"
    "elite-021-q650109375-n2-d7-b8-g2.matrix.txt"
)
DEFAULT_BINARY = Path("build/research/fast_principal_cube_lnps")
DEFAULT_OUTPUT = Path("runs/qd-remaining-bridges-20260728")


@dataclass(frozen=True)
class LogicalCube:
    identifier: str
    endpoint_identifier: str
    start_path: Path
    endpoint_path: Path
    support: tuple[Coordinate, ...]
    raw_bridge_size: int
    mode: str

    @property
    def dimension(self) -> int:
        return len(self.support)

    @property
    def endpoint_logical_mask(self) -> int:
        return (1 << self.raw_bridge_size) - 1


def resolve_repository_path(path: Path) -> Path:
    resolved = (
        path.resolve()
        if path.is_absolute()
        else (REPOSITORY_ROOT / path).resolve()
    )
    resolved.relative_to(REPOSITORY_ROOT)
    return resolved


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPOSITORY_ROOT))


def read_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return value


def file_sha256(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def file_record(path: Path) -> dict[str, str | int]:
    return {
        "path": relative(path),
        "raw_sha256": file_sha256(path),
        "size_bytes": path.stat().st_size,
    }


def json_certificate(tag: bytes, value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return sha256_bytes(tag + encoded)


def parse_arena(output: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for key, label in (
        ("determinant", "determinant"),
        ("score", "score |det|"),
        ("matrix_sha256", "matrix sha256"),
        ("normalized_sha256", "normalized sha256"),
        ("receipt_sha256", "receipt sha256"),
    ):
        match = re.search(
            rf"^{re.escape(label)}:\s*(\S+)\s*$",
            output,
            flags=re.MULTILINE,
        )
        if match is None:
            raise RuntimeError(f"arena output omitted {label}")
        fields[key] = match.group(1)
    return fields


def exact_arena_verify(matrix_path: Path, output_path: Path) -> dict[str, str]:
    result = subprocess.run(
        ["./arena", "verify", str(matrix_path)],
        cwd=REPOSITORY_ROOT,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    atomic_write(output_path, result.stdout.encode("utf-8"))
    if result.returncode != 0:
        raise RuntimeError(f"arena rejected {matrix_path}")
    parsed = parse_arena(result.stdout)
    matrix = read_matrix(matrix_path)
    determinant = bareiss_determinant(matrix)
    if (
        int(parsed["determinant"]) != determinant
        or int(parsed["score"]) != abs(determinant)
        or parsed["matrix_sha256"] != file_sha256(matrix_path)
    ):
        raise RuntimeError(f"arena result disagrees with Bareiss: {matrix_path}")
    return parsed


def differing_coordinates(first: Matrix, second: Matrix) -> tuple[Coordinate, ...]:
    if len(first) != len(second):
        raise ValueError("bridge endpoints have different orders")
    return tuple(
        (row, column)
        for row in range(len(first))
        for column in range(len(first))
        if first[row][column] != second[row][column]
    )


def core_coordinate(coordinate: Coordinate) -> bool:
    row, column = coordinate
    return 1 <= row < ORDER and 1 <= column < ORDER


def one_based(support: Sequence[Coordinate]) -> tuple[Coordinate, ...]:
    return tuple((row + 1, column + 1) for row, column in support)


def select_transverse_repairs(
    start: Matrix,
    start_determinant: int,
    raw_bridge: Sequence[Coordinate],
) -> tuple[Coordinate, ...]:
    """Select five exact-pair-ranked, row/column-disjoint core entries."""

    features = compute_features(
        FeatureCenter("elite018", start, start_determinant)
    )  # type: ignore[arg-type]
    ranked = pair_endpoint_rank(features.pair_score_rank)
    bridge = set(raw_bridge)
    selected: list[Coordinate] = []
    selected_rows: set[int] = set()
    selected_columns: set[int] = set()
    for coordinate in ranked:
        row, column = coordinate
        if (
            coordinate in bridge
            or not core_coordinate(coordinate)
            or row in selected_rows
            or column in selected_columns
        ):
            continue
        selected.append(coordinate)
        selected_rows.add(row)
        selected_columns.add(column)
        if len(selected) == 5:
            break
    if len(selected) != 5:
        raise RuntimeError("could not select five transverse repair entries")
    return tuple(selected)


def derive_logical_cubes() -> tuple[
    tuple[LogicalCube, ...], tuple[Coordinate, ...]
]:
    start_path = resolve_repository_path(ELITE018)
    endpoint22_path = resolve_repository_path(ELITE022)
    endpoint24_path = resolve_repository_path(ELITE024)
    start = read_matrix(start_path)
    endpoint22 = read_matrix(endpoint22_path)
    endpoint24 = read_matrix(endpoint24_path)
    start_det = bareiss_determinant(start)
    raw28 = differing_coordinates(start, endpoint22)
    raw22 = differing_coordinates(start, endpoint24)
    if len(raw28) != 28 or len(raw22) != 22:
        raise RuntimeError(
            f"unexpected raw bridge sizes: {len(raw28)}, {len(raw22)}"
        )
    if not all(core_coordinate(coordinate) for coordinate in (*raw28, *raw22)):
        raise RuntimeError("raw QD bridge crossed the dephased boundary")
    repairs = select_transverse_repairs(start, start_det, raw22)
    repaired27 = raw22 + repairs
    if len(repaired27) != DIMENSION or len(set(repaired27)) != DIMENSION:
        raise RuntimeError("repaired bridge is not 27 unique coordinates")
    expected_repairs = (
        (13, 6),
        (21, 2),
        (2, 8),
        (3, 4),
        (7, 18),
    )
    if repairs != expected_repairs:
        raise RuntimeError(
            "deterministic repair selection changed; audit the feature order"
        )
    if apply_mask(start, raw28, (1 << 28) - 1) != endpoint22:
        raise RuntimeError("raw 28 bridge missed elite022")
    if apply_mask(start, raw22, (1 << 22) - 1) != endpoint24:
        raise RuntimeError("raw 22 bridge missed elite024")
    if apply_mask(start, repaired27, (1 << 22) - 1) != endpoint24:
        raise RuntimeError("repaired cube does not retain the raw endpoint mask")
    return (
        (
            LogicalCube(
                identifier="elite018-to-elite022-raw-d28",
                endpoint_identifier="elite022",
                start_path=start_path,
                endpoint_path=endpoint22_path,
                support=raw28,
                raw_bridge_size=28,
                mode="exact-raw-bridge",
            ),
            LogicalCube(
                identifier="elite018-to-elite024-repaired-d27",
                endpoint_identifier="elite024",
                start_path=start_path,
                endpoint_path=endpoint24_path,
                support=repaired27,
                raw_bridge_size=22,
                mode="raw-d22-plus-deterministic-transverse5",
            ),
        ),
        raw22,
    )


def require_unique_prior_inventory(
    cubes: Sequence[affine_audit.Cube],
) -> None:
    identifiers: set[tuple[str, str]] = set()
    affine_certificates: dict[str, tuple[str, str]] = {}
    for cube in cubes:
        identifier = (cube.family, cube.label)
        if identifier in identifiers:
            raise RuntimeError(f"duplicate prior identifier: {identifier}")
        identifiers.add(identifier)
        previous = affine_certificates.get(cube.affine_certificate_sha256)
        if previous is not None:
            raise RuntimeError(
                "duplicate prior affine cube: "
                f"{previous} and {identifier}"
            )
        affine_certificates[cube.affine_certificate_sha256] = identifier


def prior_inventory_certificate(
    cubes: Sequence[affine_audit.Cube],
) -> str:
    records = [
        {
            "affine_certificate_sha256":
                cube.affine_certificate_sha256,
            "family": cube.family,
            "label": cube.label,
            "rank": cube.rank,
            "start_sha256": cube.start_sha256,
            "support_set_certificate_sha256":
                cube.support_set_certificate_sha256,
        }
        for cube in sorted(cubes, key=lambda item: (item.family, item.label))
    ]
    return json_certificate(INVENTORY_TAG, records)


def validate_qd_prior() -> tuple[
    list[affine_audit.Cube], dict[str, Any]
]:
    directory = resolve_repository_path(QD_PRIOR_DIRECTORY)
    plan, plan_hash = qd_campaign.load_plan(directory)
    aggregate_path = directory / "aggregate-report.json"
    aggregate = read_object(aggregate_path)
    scheduled = {
        str(task["id"]): task
        for task in plan["tasks"]
        if task["disposition"] == "scheduled"
    }
    if (
        aggregate.get("complete") is not True
        or aggregate.get("plan_raw_sha256") != plan_hash
        or int(aggregate["counts"]["evaluated_task_count"]) != len(scheduled)
        or int(aggregate["counts"]["unrun_task_count"]) != 0
        or int(aggregate["counts"]["evaluated_assignment_visits"])
            != int(plan["counts"]["planned_assignment_visits"])
    ):
        raise RuntimeError("QD seed36002 aggregate is not complete")
    complete_paths = sorted(
        (directory / "shards").glob(
            "shard-*-of-*/tasks/*/complete.json"
        )
    )
    completed_ids: set[str] = set()
    completed_evidence: list[dict[str, str]] = []
    for complete_path in complete_paths:
        task_id = complete_path.parent.name
        if task_id not in scheduled or task_id in completed_ids:
            raise RuntimeError("QD prior has unknown/duplicate completed task")
        qd_campaign.validate_completed_task(
            complete_path, scheduled[task_id], plan, plan_hash
        )
        completed_ids.add(task_id)
        completed_evidence.append(
            {
                "path": relative(complete_path),
                "raw_sha256": file_sha256(complete_path),
            }
        )
    if completed_ids != set(scheduled):
        raise RuntimeError("QD prior does not have every scheduled task")

    centers = {
        str(record["id"]): record
        for record in plan["qd_frontier_elites"]
    }
    cubes: list[affine_audit.Cube] = []
    for logical in plan["logical_cubes"]:
        center = centers[str(logical["center_id"])]
        start_path = resolve_repository_path(Path(str(center["path"])))
        if file_sha256(start_path) != center["raw_sha256"]:
            raise RuntimeError("QD logical center changed")
        support = tuple(
            (int(row), int(column))
            for row, column in logical["support"]
        )
        cubes.append(
            affine_audit.make_cube(
                family="qd36002",
                label=str(logical["id"]),
                start_path=start_path,
                support=support,
            )
        )
    if len(cubes) != EXPECTED_PRIOR_COUNTS["qd36002"]:
        raise RuntimeError("QD prior logical-cube count changed")
    return cubes, {
        "aggregate": file_record(aggregate_path),
        "completed_task_count": len(completed_evidence),
        "completed_task_evidence_certificate_sha256": json_certificate(
            b"maxdet-qd36002-complete-evidence-v1\0",
            completed_evidence,
        ),
        "logical_cube_count": len(cubes),
        "plan": file_record(directory / "plan.json"),
        "plan_raw_sha256": plan_hash,
        "plan_sidecar": file_record(directory / "plan.sha256"),
    }


def validate_connector_prior() -> tuple[
    affine_audit.Cube, dict[str, Any]
]:
    run_directory = resolve_repository_path(CONNECTOR_RUN_DIRECTORY)
    input_directory = resolve_repository_path(CONNECTOR_INPUT_DIRECTORY)
    target = resolve_repository_path(CONNECTOR_TARGET)
    source = resolve_repository_path(CONNECTOR_SOURCE)
    alignment = input_directory / "alignment.json"
    verification = connector_provenance.verify_sidecar(
        SimpleNamespace(
            run_dir=run_directory,
            input_dir=input_directory,
            target=target,
            source=source,
            alignment=alignment,
        )
    )
    support_zero = connector_provenance.parse_coordinate_file(
        input_directory / "support32.coords.txt"
    )
    cube = affine_audit.make_cube(
        family="connector004021",
        label="qd-aligned-connector32-elite004-elite021",
        start_path=target,
        support=one_based(support_zero),
    )
    return cube, {
        "alignment": file_record(alignment),
        "input_support": file_record(
            input_directory / "support32.coords.txt"
        ),
        "provenance": file_record(run_directory / "provenance.json"),
        "provenance_sidecar": file_record(
            run_directory / "provenance.sha256"
        ),
        "verification": verification,
    }


def load_full_prior_corpus() -> tuple[
    list[affine_audit.Cube], dict[str, Any]
]:
    batch, batch_sources = affine_audit.load_batch_cubes(
        REPOSITORY_ROOT,
        affine_audit.DEFAULT_BATCH_DIRECTORY,
        EXPECTED_PRIOR_COUNTS["batch200"],
    )
    lnps, lnps_sources = affine_audit.load_lnps_cubes(
        REPOSITORY_ROOT,
        affine_audit.DEFAULT_LNPS_DIRECTORY,
        EXPECTED_PRIOR_COUNTS["lnps"],
    )
    deep, deep_sources = affine_audit.load_deep32_cubes(
        REPOSITORY_ROOT,
        affine_audit.DEFAULT_DEEP32_DIRECTORIES,
        EXPECTED_PRIOR_COUNTS["deep32"],
    )
    qd, qd_sources = validate_qd_prior()
    connector, connector_sources = validate_connector_prior()
    priors = [*batch, *lnps, *deep, *qd, connector]
    observed_counts = Counter(cube.family for cube in priors)
    if dict(sorted(observed_counts.items())) != EXPECTED_PRIOR_COUNTS:
        raise RuntimeError(
            f"prior family counts changed: {dict(observed_counts)}"
        )
    if len(priors) != EXPECTED_PRIOR_TOTAL:
        raise RuntimeError("prior corpus total changed")
    require_unique_prior_inventory(priors)
    return priors, {
        "batch200": batch_sources,
        "connector004021": connector_sources,
        "deep32": deep_sources,
        "lnps": lnps_sources,
        "qd36002": qd_sources,
    }


def comparison_certificate(records: Sequence[dict[str, Any]]) -> str:
    return json_certificate(COMPARISON_TAG, list(records))


def candidate_comparisons(
    candidate: affine_audit.Cube,
    priors: Sequence[affine_audit.Cube],
) -> tuple[list[dict[str, Any]], int]:
    records: list[dict[str, Any]] = []
    overlap = 0
    for prior in sorted(priors, key=lambda item: (item.family, item.label)):
        intersection = affine_audit.intersect_affine(candidate, prior)
        record = affine_audit.comparison_record(prior, intersection)
        records.append(record)
        overlap += intersection.state_count
    return records, overlap


def intersection_payload(
    first: affine_audit.Cube,
    second: affine_audit.Cube,
) -> dict[str, Any]:
    result = affine_audit.intersect_affine(first, second)
    return {
        "first": first.label,
        "first_contained_in_second": result.candidate_in_prior,
        "intersection_affine_certificate_sha256":
            result.affine_certificate_sha256,
        "intersection_dimension": result.dimension,
        "intersection_state_count": str(result.state_count),
        "intersects": result.intersects,
        "second": second.label,
        "second_contained_in_first": result.prior_in_candidate,
        "union_rank": result.union_rank,
    }


def build_novelty_audit(
    logical_cubes: Sequence[LogicalCube],
    raw22: Sequence[Coordinate],
) -> dict[str, Any]:
    priors, prior_sources = load_full_prior_corpus()
    start_path = logical_cubes[0].start_path
    raw28_cube = affine_audit.make_cube(
        family="candidate",
        label="elite018-to-elite022-raw-d28",
        start_path=start_path,
        support=one_based(logical_cubes[0].support),
    )
    raw22_cube = affine_audit.make_cube(
        family="candidate",
        label="elite018-to-elite024-raw-d22",
        start_path=start_path,
        support=one_based(raw22),
    )
    repaired27_cube = affine_audit.make_cube(
        family="candidate",
        label="elite018-to-elite024-repaired-d27",
        start_path=start_path,
        support=one_based(logical_cubes[1].support),
    )
    candidates = (raw28_cube, raw22_cube, repaired27_cube)
    candidate_records: list[dict[str, Any]] = []
    for candidate in candidates:
        comparisons, pairwise_overlap = candidate_comparisons(
            candidate, priors
        )
        nonempty = [record for record in comparisons if record["intersects"]]
        if pairwise_overlap != 0 or nonempty:
            raise RuntimeError(
                f"{candidate.label} intersects the pinned prior corpus"
            )
        candidate_records.append(
            {
                "affine_cube": affine_audit.cube_payload(
                    candidate, REPOSITORY_ROOT
                ),
                "candidate_new_dephased_state_count_relative_prior":
                    str(1 << candidate.rank),
                "comparisons": comparisons,
                "comparisons_certificate_sha256":
                    comparison_certificate(comparisons),
                "empty_intersection_count": len(comparisons),
                "nonempty_intersection_count": 0,
                "pairwise_intersection_state_sum": "0",
            }
        )

    raw28_raw22 = intersection_payload(raw28_cube, raw22_cube)
    raw28_repaired = intersection_payload(raw28_cube, repaired27_cube)
    raw22_repaired = intersection_payload(raw22_cube, repaired27_cube)
    if (
        int(raw28_repaired["intersection_state_count"]) != 8
        or raw28_repaired["intersection_dimension"] != 3
        or int(raw22_repaired["intersection_state_count"]) != 1 << 22
        or raw22_repaired["first_contained_in_second"] is not True
    ):
        raise RuntimeError("new-candidate overlap invariants changed")
    executed_total = (1 << 28) + (1 << 27)
    executed_unique = executed_total - 8
    prior_inventory = [
        affine_audit.cube_payload(cube, REPOSITORY_ROOT)
        for cube in sorted(priors, key=lambda item: (item.family, item.label))
    ]
    return {
        "candidate_audits": candidate_records,
        "claim_boundary": [
            (
                "Zero overlap is exact only in the canonical dephased "
                "22x22 GF(2) space and only against the 362 pinned logical "
                "cubes inventoried here."
            ),
            (
                "The audit does not quotient row/column permutations or "
                "transpose."
            ),
            (
                "The raw d22 bridge is not counted as a separate evaluator "
                "search; all 2^22 raw states are contained in the repaired "
                "d27 cube."
            ),
            (
                "Assignment visits count the two executed cubes separately; "
                "their exact mutual overlap is eight states."
            ),
        ],
        "coordinate_indexing": "one_based",
        "dephased_space": {
            "bit_count": affine_audit.CORE_BITS,
            "core_order": affine_audit.CORE_ORDER,
            "field": "GF(2)",
        },
        "executed_cube_union": {
            "assignment_visits": str(executed_total),
            "mutual_overlap_state_count": "8",
            "new_dephased_states_relative_prior_and_each_other":
                str(executed_unique),
            "pair": [
                raw28_cube.label,
                repaired27_cube.label,
            ],
        },
        "method": "exact-dephased-gf2-affine-full-corpus-audit-v1",
        "new_candidate_intersections": [
            raw28_raw22,
            raw28_repaired,
            raw22_repaired,
        ],
        "prior_corpus": {
            "count": len(priors),
            "counts_by_family": EXPECTED_PRIOR_COUNTS,
            "dephased_affine_certificates_unique": len(
                {cube.affine_certificate_sha256 for cube in priors}
            ),
            "inventory": prior_inventory,
            "inventory_certificate_sha256":
                prior_inventory_certificate(priors),
            "sources": prior_sources,
        },
        "provenance": {
            "affine_auditor": file_record(
                resolve_repository_path(
                    Path("research/affine_gf2_audit.py")
                )
            ),
            "campaign_driver": file_record(Path(__file__).resolve()),
            "connector_validator": file_record(
                resolve_repository_path(
                    Path(
                        "research/"
                        "fast_cube_connector_provenance.py"
                    )
                )
            ),
            "qd_validator": file_record(
                resolve_repository_path(
                    Path("research/qd_cube_campaign.py")
                )
            ),
        },
        "schema_version": SCHEMA_VERSION,
    }


def reroot_offset(
    fixed: Matrix, support: Sequence[Coordinate]
) -> tuple[int, Matrix, int]:
    determinant = bareiss_determinant(fixed)
    if determinant != 0:
        return 0, fixed, determinant
    offsets = [1 << index for index in range(len(support))]
    offsets.extend(
        (1 << first) | (1 << second)
        for first in range(len(support))
        for second in range(first + 1, len(support))
    )
    for offset in offsets:
        candidate = apply_mask(fixed, support, offset)
        determinant = bareiss_determinant(candidate)
        if determinant != 0:
            return offset, candidate, determinant
    raise RuntimeError("could not find a nonsingular task reroot")


def logical_record(cube: LogicalCube) -> dict[str, Any]:
    start = read_matrix(cube.start_path)
    endpoint = read_matrix(cube.endpoint_path)
    endpoint_mask = cube.endpoint_logical_mask
    if apply_mask(start, cube.support, endpoint_mask) != endpoint:
        raise RuntimeError("logical endpoint mask failed reconstruction")
    return {
        "dimension": cube.dimension,
        "endpoint": relative(cube.endpoint_path),
        "endpoint_absolute_determinant": str(
            abs(bareiss_determinant(endpoint))
        ),
        "endpoint_identifier": cube.endpoint_identifier,
        "endpoint_logical_mask_decimal": str(endpoint_mask),
        "endpoint_raw_sha256": file_sha256(cube.endpoint_path),
        "full_affine_fingerprint_sha256":
            affine_fingerprint(start, cube.support),
        "id": cube.identifier,
        "mode": cube.mode,
        "partition": {
            "disjoint": True,
            "inner_dimension": DIMENSION,
            "leaf_count": 1 << (cube.dimension - DIMENSION),
            "outer_dimension": cube.dimension - DIMENSION,
            "total_assignments": 1 << cube.dimension,
        },
        "raw_bridge_size": cube.raw_bridge_size,
        "start": relative(cube.start_path),
        "start_raw_sha256": file_sha256(cube.start_path),
        "support": [
            [row + 1, column + 1] for row, column in cube.support
        ],
        "support_set_sha256": sha256_bytes(support_bytes(cube.support)),
    }


def build_tasks(
    logical_cubes: Sequence[LogicalCube],
) -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    fingerprints: set[str] = set()
    scheduled_index = 0
    for cube in logical_cubes:
        start = read_matrix(cube.start_path)
        inner = cube.support[:DIMENSION]
        outer = cube.support[DIMENSION:]
        if len(inner) != DIMENSION:
            raise RuntimeError("evaluator task does not have 27 inner bits")
        endpoint_inner = cube.endpoint_logical_mask & (
            (1 << DIMENSION) - 1
        )
        endpoint_outer = cube.endpoint_logical_mask >> DIMENSION
        for outer_mask in range(1 << len(outer)):
            fixed = apply_mask(start, outer, outer_mask)
            reroot, task_start, signed = reroot_offset(fixed, inner)
            fingerprint = affine_fingerprint(task_start, inner)
            if fingerprint in fingerprints:
                raise RuntimeError("duplicate evaluator task cube")
            fingerprints.add(fingerprint)
            endpoint_expected = outer_mask == endpoint_outer
            task_id = (
                f"{cube.identifier}-leaf-{outer_mask:02d}"
                if outer
                else cube.identifier
            )
            tasks.append(
                {
                    "affine_fingerprint_sha256": fingerprint,
                    "assignment_visits": ASSIGNMENTS_PER_TASK,
                    "dimension": DIMENSION,
                    "endpoint_expected": endpoint_expected,
                    "expected_endpoint_engine_mask_decimal": (
                        str(endpoint_inner ^ reroot)
                        if endpoint_expected
                        else None
                    ),
                    "id": task_id,
                    "logical_cube_id": cube.identifier,
                    "outer_mask_decimal": str(outer_mask),
                    "reroot_xor_mask_decimal": str(reroot),
                    "scheduled_index": scheduled_index,
                    "start_parsed_matrix_sha256": sha256_bytes(
                        canonical_matrix_bytes(task_start)
                    ),
                    "start_signed_determinant": str(signed),
                }
            )
            scheduled_index += 1
    return tasks


def plan_campaign(output_directory: Path, binary: Path) -> dict[str, Any]:
    output_directory = resolve_repository_path(output_directory)
    binary = resolve_repository_path(binary)
    if output_directory.exists():
        raise FileExistsError(f"output exists: {output_directory}")
    if not binary.is_file():
        raise FileNotFoundError(f"missing evaluator: {binary}")

    logical_cubes, raw22 = derive_logical_cubes()
    novelty = build_novelty_audit(logical_cubes, raw22)
    logical_records = [
        logical_record(cube) for cube in logical_cubes
    ]
    tasks = build_tasks(logical_cubes)
    if len(tasks) != 3:
        raise RuntimeError("expected exactly three evaluator tasks")

    output_directory.mkdir(parents=True)
    inputs = output_directory / "inputs"
    inputs.mkdir()
    atomic_write(
        inputs / "elite018-elite022-raw28.coords.txt",
        support_bytes(logical_cubes[0].support),
    )
    atomic_write(
        inputs / "elite018-elite024-raw22.coords.txt",
        support_bytes(raw22),
    )
    atomic_write(
        inputs / "elite018-elite024-repair5.coords.txt",
        support_bytes(logical_cubes[1].support[22:]),
    )
    atomic_write(
        inputs / "elite018-elite024-repaired27.coords.txt",
        support_bytes(logical_cubes[1].support),
    )
    novelty_path = output_directory / "novelty-audit.json"
    atomic_json(novelty_path, novelty)

    dependency_paths = (
        Path(__file__).resolve(),
        resolve_repository_path(Path("research/affine_gf2_audit.py")),
        resolve_repository_path(Path("research/fast_cube_batch.py")),
        resolve_repository_path(Path("research/fast_cube_lnps.py")),
        resolve_repository_path(
            Path("research/fast_cube_connector_provenance.py")
        ),
        resolve_repository_path(Path("research/qd_cube_campaign.py")),
        resolve_repository_path(Path("maxdet/exact.py")),
        resolve_repository_path(Path("arena")),
    )
    engine_source = resolve_repository_path(
        Path("research/fast_principal_cube.cpp")
    )
    plan = {
        "claim_boundary": novelty["claim_boundary"],
        "counts": {
            "assignment_visits": sum(
                int(task["assignment_visits"]) for task in tasks
            ),
            "evaluator_task_count": len(tasks),
            "logical_search_cube_count": len(logical_cubes),
            "raw_bridge_control_count": 2,
        },
        "engine": "fast-principal-minor-entry-cube-v1",
        "frontier": str(FRONTIER),
        "logical_cubes": logical_records,
        "method": METHOD,
        "novelty_audit": file_record(novelty_path),
        "provenance": {
            "engine_binary": relative(binary),
            "engine_binary_raw_sha256": file_sha256(binary),
            "engine_source": relative(engine_source),
            "engine_source_raw_sha256": file_sha256(engine_source),
            "runtime_dependencies": [
                file_record(path) for path in dependency_paths
            ],
        },
        "raw_d22_control": {
            "dimension": 22,
            "endpoint": relative(logical_cubes[1].endpoint_path),
            "executed_as_separate_cube": False,
            "state_count": str(1 << 22),
            "subset_of": logical_cubes[1].identifier,
            "support": [
                [row + 1, column + 1] for row, column in raw22
            ],
        },
        "schema_version": SCHEMA_VERSION,
        "tasks": tasks,
        "top_k": TOP_K,
    }
    plan_path = output_directory / "plan.json"
    atomic_json(plan_path, plan)
    plan_hash = file_sha256(plan_path)
    atomic_write(
        output_directory / "plan.sha256",
        f"{plan_hash}  plan.json\n".encode("ascii"),
    )
    return plan


def load_plan(output_directory: Path) -> tuple[dict[str, Any], str]:
    output_directory = resolve_repository_path(output_directory)
    plan_path = output_directory / "plan.json"
    sidecar_path = output_directory / "plan.sha256"
    expected = sidecar_path.read_text(encoding="ascii").split()[0]
    actual = file_sha256(plan_path)
    if expected != actual:
        raise RuntimeError("immutable plan hash mismatch")
    plan = read_object(plan_path)
    if (
        plan.get("method") != METHOD
        or int(plan.get("top_k", 0)) != TOP_K
        or int(plan["counts"]["evaluator_task_count"]) != 3
    ):
        raise RuntimeError("unexpected campaign plan")
    provenance = plan["provenance"]
    for path_key, hash_key in (
        ("engine_binary", "engine_binary_raw_sha256"),
        ("engine_source", "engine_source_raw_sha256"),
    ):
        path = resolve_repository_path(Path(provenance[path_key]))
        if file_sha256(path) != provenance[hash_key]:
            raise RuntimeError(f"pinned provenance changed: {path_key}")
    for record in provenance["runtime_dependencies"]:
        path = resolve_repository_path(Path(record["path"]))
        if file_sha256(path) != record["raw_sha256"]:
            raise RuntimeError(
                f"pinned runtime dependency changed: {record['path']}"
            )
    novelty_path = resolve_repository_path(
        Path(plan["novelty_audit"]["path"])
    )
    if file_sha256(novelty_path) != plan["novelty_audit"]["raw_sha256"]:
        raise RuntimeError("novelty audit changed after planning")
    return plan, actual


def logical_by_id(plan: dict[str, Any]) -> dict[str, dict[str, Any]]:
    records = {
        str(record["id"]): record for record in plan["logical_cubes"]
    }
    if len(records) != len(plan["logical_cubes"]):
        raise RuntimeError("duplicate logical cube ID")
    return records


def reconstruct_task(
    plan: dict[str, Any], task: dict[str, Any]
) -> tuple[Matrix, tuple[Coordinate, ...], dict[str, Any]]:
    logical = logical_by_id(plan)[str(task["logical_cube_id"])]
    start_path = resolve_repository_path(Path(logical["start"]))
    if file_sha256(start_path) != logical["start_raw_sha256"]:
        raise RuntimeError("logical start changed")
    start = read_matrix(start_path)
    support = tuple(
        (int(row) - 1, int(column) - 1)
        for row, column in logical["support"]
    )
    inner = support[:DIMENSION]
    outer = support[DIMENSION:]
    fixed = apply_mask(start, outer, int(task["outer_mask_decimal"]))
    reroot = int(task["reroot_xor_mask_decimal"])
    task_start = apply_mask(fixed, inner, reroot)
    determinant = bareiss_determinant(task_start)
    if determinant == 0 or str(determinant) != task["start_signed_determinant"]:
        raise RuntimeError("task start determinant changed")
    if sha256_bytes(canonical_matrix_bytes(task_start)) != task[
        "start_parsed_matrix_sha256"
    ]:
        raise RuntimeError("task start hash mismatch")
    if affine_fingerprint(task_start, inner) != task[
        "affine_fingerprint_sha256"
    ]:
        raise RuntimeError("task affine fingerprint mismatch")
    return task_start, inner, logical


def resolve_report_artifact(path_value: str, attempt: Path) -> Path:
    path = Path(path_value)
    resolved = path.resolve() if path.is_absolute() else (
        REPOSITORY_ROOT / path
    ).resolve()
    resolved.relative_to(attempt.resolve())
    if not resolved.is_file():
        raise RuntimeError(f"missing report artifact: {resolved}")
    return resolved


def replay_engine_result(
    plan: dict[str, Any],
    task: dict[str, Any],
    attempt: Path,
) -> dict[str, Any]:
    start, support, logical = reconstruct_task(plan, task)
    report_path = attempt / "report.json"
    best_path = attempt / "best.matrix.txt"
    report = read_object(report_path)
    if (
        report.get("engine") != plan["engine"]
        or report.get("complete") is not True
        or report.get("all_assignments_bound_checked") is not True
        or int(report.get("dimension", 0)) != DIMENSION
        or int(report.get("assignments", 0)) != ASSIGNMENTS_PER_TASK
        or report.get("start_parsed_matrix_sha256")
            != task["start_parsed_matrix_sha256"]
        or report.get("coordinate_file_raw_sha256")
            != sha256_bytes(support_bytes(support))
    ):
        raise RuntimeError(f"invalid engine report: {task['id']}")
    report_support = tuple(
        (int(row) - 1, int(column) - 1)
        for row, column in report["support"]
    )
    if report_support != support:
        raise RuntimeError("engine report support changed")

    best_mask = int(report["best_mask_decimal"])
    if not 0 <= best_mask < 1 << DIMENSION:
        raise RuntimeError("engine best mask is outside the task")
    reconstructed_best = apply_mask(start, support, best_mask)
    artifact_best = read_matrix(best_path)
    if reconstructed_best != artifact_best:
        raise RuntimeError("best artifact does not match its mask")
    best_exact = bareiss_determinant(reconstructed_best)
    best_hash = sha256_bytes(canonical_matrix_bytes(reconstructed_best))
    if (
        best_exact != int(report["best_signed_determinant"])
        or abs(best_exact) != int(report["best_absolute_determinant"])
        or report["best_matrix_sha256"] != best_hash
        or report["output_raw_sha256"] != best_hash
    ):
        raise RuntimeError("engine best failed exact replay")

    top_records = report.get("top_k_candidates", [])
    if (
        not isinstance(top_records, list)
        or len(top_records) != int(report["top_k_captured"])
        or int(report["top_k_requested"]) != TOP_K
    ):
        raise RuntimeError("engine top-K accounting mismatch")
    previous_key: tuple[int, int] | None = None
    for expected_rank, record in enumerate(top_records, 1):
        mask = int(record["mask_decimal"])
        if not 0 < mask < 1 << DIMENSION:
            raise RuntimeError("top-K mask is outside/nonzero gate")
        matrix = apply_mask(start, support, mask)
        exact = bareiss_determinant(matrix)
        artifact = resolve_report_artifact(str(record["artifact"]), attempt)
        artifact_matrix = read_matrix(artifact)
        key = (-abs(exact), mask)
        if (
            exact != int(record["signed_determinant"])
            or abs(exact) != int(record["absolute_determinant"])
            or int(record["rank"]) != expected_rank
            or artifact_matrix != matrix
            or file_sha256(artifact) != record["artifact_raw_sha256"]
            or previous_key is not None
            and key < previous_key
        ):
            raise RuntimeError("top-K record failed exact ordered replay")
        previous_key = key

    tie_masks = [int(value) for value in report["frontier_tie_masks_decimal"]]
    if len(tie_masks) != len(set(tie_masks)):
        raise RuntimeError("engine returned duplicate frontier tie masks")
    for mask in tie_masks:
        if not 0 < mask < 1 << DIMENSION:
            raise RuntimeError("frontier tie mask is outside/nonzero gate")
        if abs(bareiss_determinant(apply_mask(start, support, mask))) != FRONTIER:
            raise RuntimeError("frontier tie failed exact replay")
    if (
        bool(report["frontier_tie_masks_truncated"])
        and len(tie_masks) >= int(report["frontier_nonzero_ties"])
    ):
        raise RuntimeError("tie truncation metadata is inconsistent")
    if (
        not bool(report["frontier_tie_masks_truncated"])
        and len(tie_masks) != int(report["frontier_nonzero_ties"])
    ):
        raise RuntimeError("complete tie list count is inconsistent")

    endpoint_confirmed: bool | None = None
    if task["endpoint_expected"]:
        endpoint_mask = int(task["expected_endpoint_engine_mask_decimal"])
        endpoint = apply_mask(start, support, endpoint_mask)
        endpoint_path = resolve_repository_path(Path(logical["endpoint"]))
        endpoint_artifact = read_matrix(endpoint_path)
        endpoint_confirmed = (
            endpoint == endpoint_artifact
            and file_sha256(endpoint_path) == logical["endpoint_raw_sha256"]
            and abs(bareiss_determinant(endpoint))
                == int(logical["endpoint_absolute_determinant"])
        )
        if not endpoint_confirmed:
            raise RuntimeError("endpoint control failed")

    logical_best_mask = (
        (best_mask ^ int(task["reroot_xor_mask_decimal"]))
        | (int(task["outer_mask_decimal"]) << DIMENSION)
    )
    return {
        "best_absolute_determinant": str(abs(best_exact)),
        "best_engine_mask_decimal": str(best_mask),
        "best_logical_mask_decimal": str(logical_best_mask),
        "best_matrix_raw_sha256": best_hash,
        "best_signed_determinant": str(best_exact),
        "endpoint_control_confirmed": endpoint_confirmed,
        "frontier_nonzero_tie_count": int(
            report["frontier_nonzero_ties"]
        ),
        "frontier_tie_masks": tie_masks,
        "frontier_tie_masks_truncated": bool(
            report["frontier_tie_masks_truncated"]
        ),
        "replayed_best_records": 1,
        "replayed_returned_tie_masks": len(tie_masks),
        "replayed_top_k_records": len(top_records),
    }


def next_attempt(task_directory: Path) -> Path:
    maximum = 0
    if task_directory.exists():
        for child in task_directory.iterdir():
            match = re.fullmatch(r"attempt-(\d{3})", child.name)
            if match is not None:
                maximum = max(maximum, int(match.group(1)))
    attempt = task_directory / f"attempt-{maximum + 1:03d}"
    attempt.mkdir(parents=True, exist_ok=False)
    return attempt


def validate_arena_artifact(
    matrix_path: Path,
    output_path: Path,
) -> dict[str, str]:
    parsed = parse_arena(output_path.read_text(encoding="utf-8"))
    matrix = read_matrix(matrix_path)
    determinant = bareiss_determinant(matrix)
    if (
        not output_path.read_text(encoding="utf-8").startswith("VERIFIED\n")
        or int(parsed["determinant"]) != determinant
        or int(parsed["score"]) != abs(determinant)
        or parsed["matrix_sha256"] != file_sha256(matrix_path)
    ):
        raise RuntimeError("stored arena verification disagrees")
    return parsed


def validate_complete(
    complete_path: Path,
    task: dict[str, Any],
    plan: dict[str, Any],
    plan_hash: str,
) -> dict[str, Any]:
    complete = read_object(complete_path)
    if (
        complete.get("complete") is not True
        or complete.get("plan_raw_sha256") != plan_hash
        or complete.get("task_id") != task["id"]
        or int(complete.get("assignment_visits", 0))
            != ASSIGNMENTS_PER_TASK
    ):
        raise RuntimeError(f"invalid task completion marker: {complete_path}")
    attempt = resolve_repository_path(
        Path(complete["fresh_attempt_directory"])
    )
    attempt.relative_to(complete_path.parent.resolve())
    for field, filename in (
        ("engine_report_raw_sha256", "report.json"),
        ("best_matrix_raw_sha256", "best.matrix.txt"),
        ("best_arena_verification_raw_sha256", "best.arena-verify.txt"),
    ):
        if file_sha256(attempt / filename) != complete[field]:
            raise RuntimeError(f"completed artifact changed: {filename}")
    replay = replay_engine_result(plan, task, attempt)
    for key in (
        "best_absolute_determinant",
        "best_engine_mask_decimal",
        "best_logical_mask_decimal",
        "best_matrix_raw_sha256",
        "best_signed_determinant",
        "endpoint_control_confirmed",
        "frontier_nonzero_tie_count",
        "frontier_tie_masks_truncated",
        "replayed_best_records",
        "replayed_returned_tie_masks",
        "replayed_top_k_records",
    ):
        if complete.get(key) != replay[key]:
            raise RuntimeError(f"completion/replay mismatch: {key}")
    best_arena = validate_arena_artifact(
        attempt / "best.matrix.txt",
        attempt / "best.arena-verify.txt",
    )
    if best_arena != complete["best_arena"]:
        raise RuntimeError("stored best arena fields changed")
    if task["endpoint_expected"]:
        endpoint_matrix = attempt / "endpoint-control.matrix.txt"
        endpoint_arena_path = attempt / "endpoint-control.arena-verify.txt"
        if (
            file_sha256(endpoint_matrix)
                != complete["endpoint_control_matrix_raw_sha256"]
            or file_sha256(endpoint_arena_path)
                != complete[
                    "endpoint_control_arena_verification_raw_sha256"
                ]
        ):
            raise RuntimeError("endpoint control artifact changed")
        endpoint_arena = validate_arena_artifact(
            endpoint_matrix, endpoint_arena_path
        )
        if endpoint_arena != complete["endpoint_control_arena"]:
            raise RuntimeError("endpoint arena fields changed")
    elif complete.get("endpoint_control_arena") is not None:
        raise RuntimeError("unexpected endpoint arena record")
    tie_records = complete["frontier_tie_arena_verifications"]
    if len(tie_records) != len(replay["frontier_tie_masks"]):
        raise RuntimeError("frontier tie verification count changed")
    for record, mask in zip(tie_records, replay["frontier_tie_masks"]):
        if int(record["mask_decimal"]) != mask:
            raise RuntimeError("frontier tie arena mask changed")
        matrix_path = resolve_repository_path(Path(record["matrix"]))
        arena_path = resolve_repository_path(Path(record["arena_verification"]))
        if (
            file_sha256(matrix_path) != record["matrix_raw_sha256"]
            or file_sha256(arena_path)
                != record["arena_verification_raw_sha256"]
            or validate_arena_artifact(matrix_path, arena_path)
                != record["arena"]
        ):
            raise RuntimeError("frontier tie arena artifact changed")
    return complete


def write_progress(
    output_directory: Path,
    plan_hash: str,
    completed: Sequence[dict[str, Any]],
    started: float,
) -> None:
    atomic_json(
        output_directory / "run-progress.json",
        {
            "assignment_visits_completed": sum(
                int(record["assignment_visits"]) for record in completed
            ),
            "complete": len(completed) == 3,
            "completed_task_count": len(completed),
            "elapsed_seconds": round(time.monotonic() - started, 6),
            "plan_raw_sha256": plan_hash,
            "schema_version": SCHEMA_VERSION,
        },
    )


def run_campaign(output_directory: Path, resume: bool) -> dict[str, Any]:
    output_directory = resolve_repository_path(output_directory)
    plan, plan_hash = load_plan(output_directory)
    binary = resolve_repository_path(
        Path(plan["provenance"]["engine_binary"])
    )
    tasks_root = output_directory / "tasks"
    if tasks_root.exists() and not resume:
        raise FileExistsError("tasks exist; pass --resume")
    tasks_root.mkdir(exist_ok=True)
    completed: list[dict[str, Any]] = []
    started = time.monotonic()

    for task in plan["tasks"]:
        task_directory = tasks_root / str(task["id"])
        complete_path = task_directory / "complete.json"
        if complete_path.is_file():
            completed.append(
                validate_complete(
                    complete_path, task, plan, plan_hash
                )
            )
            continue
        attempt = next_attempt(task_directory)
        start, support, logical = reconstruct_task(plan, task)
        start_path = attempt / "start.matrix.txt"
        support_path = attempt / "support.coords.txt"
        best_path = attempt / "best.matrix.txt"
        report_path = attempt / "report.json"
        atomic_write(start_path, canonical_matrix_bytes(start))
        atomic_write(support_path, support_bytes(support))
        command = [
            str(binary),
            "--start",
            str(start_path),
            "--coordinates",
            str(support_path),
            "--output",
            str(best_path),
            "--tie-output",
            str(attempt / "frontier-tie.matrix.txt"),
            "--log",
            str(attempt / "search.jsonl"),
            "--report",
            str(report_path),
            "--top-k",
            str(TOP_K),
            "--top-k-output-dir",
            str(attempt / "top-k"),
        ]
        result = subprocess.run(
            command,
            cwd=REPOSITORY_ROOT,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        atomic_write(
            attempt / "stdout.txt", result.stdout.encode("utf-8")
        )
        atomic_write(
            attempt / "stderr.txt", result.stderr.encode("utf-8")
        )
        if result.returncode != 0:
            atomic_json(
                attempt / "failure.json",
                {
                    "command": command,
                    "returncode": result.returncode,
                    "stderr": result.stderr,
                    "stdout": result.stdout,
                },
            )
            raise RuntimeError(f"engine failed: {task['id']}")

        replay = replay_engine_result(plan, task, attempt)
        best_arena_path = attempt / "best.arena-verify.txt"
        best_arena = exact_arena_verify(best_path, best_arena_path)
        if int(best_arena["score"]) != int(
            replay["best_absolute_determinant"]
        ):
            raise RuntimeError("best arena score disagrees with replay")

        endpoint_arena: dict[str, str] | None = None
        endpoint_matrix_hash: str | None = None
        endpoint_arena_hash: str | None = None
        if task["endpoint_expected"]:
            endpoint_mask = int(
                task["expected_endpoint_engine_mask_decimal"]
            )
            endpoint = apply_mask(start, support, endpoint_mask)
            endpoint_matrix_path = (
                attempt / "endpoint-control.matrix.txt"
            )
            endpoint_arena_path = (
                attempt / "endpoint-control.arena-verify.txt"
            )
            atomic_write(
                endpoint_matrix_path, canonical_matrix_bytes(endpoint)
            )
            endpoint_arena = exact_arena_verify(
                endpoint_matrix_path, endpoint_arena_path
            )
            endpoint_matrix_hash = file_sha256(endpoint_matrix_path)
            endpoint_arena_hash = file_sha256(endpoint_arena_path)
            if int(endpoint_arena["score"]) != int(
                logical["endpoint_absolute_determinant"]
            ):
                raise RuntimeError("endpoint arena score changed")

        tie_verifications: list[dict[str, Any]] = []
        for mask in replay["frontier_tie_masks"]:
            tie = apply_mask(start, support, mask)
            tie_directory = attempt / "frontier-tie-replay"
            matrix_path = tie_directory / f"mask-{mask}.matrix.txt"
            arena_path = tie_directory / f"mask-{mask}.arena-verify.txt"
            atomic_write(matrix_path, canonical_matrix_bytes(tie))
            arena = exact_arena_verify(matrix_path, arena_path)
            if int(arena["score"]) != FRONTIER:
                raise RuntimeError("frontier tie arena score changed")
            tie_verifications.append(
                {
                    "arena": arena,
                    "arena_verification": relative(arena_path),
                    "arena_verification_raw_sha256":
                        file_sha256(arena_path),
                    "mask_decimal": str(mask),
                    "matrix": relative(matrix_path),
                    "matrix_raw_sha256": file_sha256(matrix_path),
                }
            )

        complete = {
            "affine_fingerprint_sha256":
                task["affine_fingerprint_sha256"],
            "assignment_visits": ASSIGNMENTS_PER_TASK,
            "best_absolute_determinant":
                replay["best_absolute_determinant"],
            "best_arena": best_arena,
            "best_arena_verification_raw_sha256":
                file_sha256(best_arena_path),
            "best_engine_mask_decimal":
                replay["best_engine_mask_decimal"],
            "best_logical_mask_decimal":
                replay["best_logical_mask_decimal"],
            "best_matrix_raw_sha256":
                replay["best_matrix_raw_sha256"],
            "best_signed_determinant":
                replay["best_signed_determinant"],
            "complete": True,
            "endpoint_control_arena": endpoint_arena,
            "endpoint_control_arena_verification_raw_sha256":
                endpoint_arena_hash,
            "endpoint_control_confirmed":
                replay["endpoint_control_confirmed"],
            "endpoint_control_matrix_raw_sha256": endpoint_matrix_hash,
            "engine_report_raw_sha256": file_sha256(report_path),
            "fresh_attempt_directory": relative(attempt),
            "frontier_nonzero_tie_count":
                replay["frontier_nonzero_tie_count"],
            "frontier_tie_arena_verifications": tie_verifications,
            "frontier_tie_masks_truncated":
                replay["frontier_tie_masks_truncated"],
            "logical_cube_id": logical["id"],
            "plan_raw_sha256": plan_hash,
            "replayed_best_records": replay["replayed_best_records"],
            "replayed_returned_tie_masks":
                replay["replayed_returned_tie_masks"],
            "replayed_top_k_records":
                replay["replayed_top_k_records"],
            "schema_version": SCHEMA_VERSION,
            "task_id": task["id"],
        }
        atomic_json(complete_path, complete)
        completed.append(complete)
        write_progress(
            output_directory, plan_hash, completed, started
        )
        print(
            f"task={task['id']} best="
            f"{complete['best_absolute_determinant']} "
            f"ties={complete['frontier_nonzero_tie_count']}",
            flush=True,
        )

    write_progress(output_directory, plan_hash, completed, started)
    return read_object(output_directory / "run-progress.json")


def collect_completed(
    output_directory: Path,
    plan: dict[str, Any],
    plan_hash: str,
) -> list[dict[str, Any]]:
    completed: list[dict[str, Any]] = []
    expected_paths: set[Path] = set()
    for task in plan["tasks"]:
        path = output_directory / "tasks" / str(task["id"]) / "complete.json"
        expected_paths.add(path.resolve())
        if not path.is_file():
            raise RuntimeError(f"task is incomplete: {task['id']}")
        completed.append(
            validate_complete(path, task, plan, plan_hash)
        )
    observed = {
        path.resolve()
        for path in (output_directory / "tasks").glob("*/complete.json")
    }
    if observed != expected_paths:
        raise RuntimeError("unexpected/missing task completion marker")
    return completed


def build_summary(
    plan: dict[str, Any],
    plan_hash: str,
    completed: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    logical_results: list[dict[str, Any]] = []
    for logical in plan["logical_cubes"]:
        records = [
            record
            for record in completed
            if record["logical_cube_id"] == logical["id"]
        ]
        logical_results.append(
            {
                "assignment_visits": str(
                    sum(
                        int(record["assignment_visits"])
                        for record in records
                    )
                ),
                "best_absolute_determinant": str(
                    max(
                        int(record["best_absolute_determinant"])
                        for record in records
                    )
                ),
                "dimension": logical["dimension"],
                "endpoint_absolute_determinant":
                    logical["endpoint_absolute_determinant"],
                "endpoint_controls_confirmed": sum(
                    record["endpoint_control_confirmed"] is True
                    for record in records
                ),
                "id": logical["id"],
                "task_count": len(records),
            }
        )
    best = max(
        int(record["best_absolute_determinant"])
        for record in completed
    )
    total_visits = sum(
        int(record["assignment_visits"]) for record in completed
    )
    novelty = read_object(
        resolve_repository_path(Path(plan["novelty_audit"]["path"]))
    )
    return {
        "arena_verification": {
            "best_artifacts_verified": len(completed),
            "endpoint_controls_verified": sum(
                record["endpoint_control_arena"] is not None
                for record in completed
            ),
            "frontier_ties_verified": sum(
                len(record["frontier_tie_arena_verifications"])
                for record in completed
            ),
        },
        "assignment_visits": str(total_visits),
        "assignment_visits_are_pairwise_unique": False,
        "best_absolute_determinant": str(best),
        "complete": len(completed) == len(plan["tasks"]),
        "executed_cube_exact_mutual_overlap_state_count":
            novelty["executed_cube_union"]["mutual_overlap_state_count"],
        "frontier": str(FRONTIER),
        "frontier_gain": str(best - FRONTIER),
        "logical_results": logical_results,
        "method": METHOD,
        "new_dephased_states_relative_pinned_prior_and_each_other":
            novelty["executed_cube_union"][
                "new_dephased_states_relative_prior_and_each_other"
            ],
        "novelty_audit_raw_sha256":
            plan["novelty_audit"]["raw_sha256"],
        "plan_raw_sha256": plan_hash,
        "promotions": [
            {
                "best_absolute_determinant":
                    record["best_absolute_determinant"],
                "task_id": record["task_id"],
            }
            for record in completed
            if int(record["best_absolute_determinant"]) > FRONTIER
        ],
        "raw_d22_control": {
            "executed_as_separate_cube": False,
            "state_count": str(1 << 22),
            "states_contained_in_repaired_d27": str(1 << 22),
        },
        "replay_audit": {
            "replayed_best_records": sum(
                int(record["replayed_best_records"])
                for record in completed
            ),
            "replayed_returned_tie_masks": sum(
                int(record["replayed_returned_tie_masks"])
                for record in completed
            ),
            "replayed_top_k_records": sum(
                int(record["replayed_top_k_records"])
                for record in completed
            ),
            "tie_mask_truncation_task_count": sum(
                record["frontier_tie_masks_truncated"] is True
                for record in completed
            ),
        },
        "schema_version": SCHEMA_VERSION,
        "task_count": len(completed),
    }


def artifact_inventory(
    output_directory: Path,
) -> tuple[list[dict[str, str | int]], str]:
    excluded = {"provenance.json", "provenance.sha256"}
    records = [
        file_record(path)
        for path in sorted(output_directory.rglob("*"))
        if path.is_file() and path.name not in excluded
    ]
    return records, json_certificate(ARTIFACT_INVENTORY_TAG, records)


def write_provenance(
    output_directory: Path,
    plan: dict[str, Any],
    plan_hash: str,
) -> dict[str, Any]:
    records, digest = artifact_inventory(output_directory)
    provenance = {
        "aggregate_report": file_record(
            output_directory / "aggregate-report.json"
        ),
        "artifact_inventory": {
            "count": len(records),
            "sha256": digest,
        },
        "immutable_plan_raw_sha256": plan_hash,
        "method": METHOD,
        "novelty_audit": plan["novelty_audit"],
        "result_changed": False,
        "run_directory": relative(output_directory),
        "schema_version": SCHEMA_VERSION,
        "tool": file_record(Path(__file__).resolve()),
    }
    path = output_directory / "provenance.json"
    atomic_json(path, provenance)
    digest_path = output_directory / "provenance.sha256"
    atomic_write(
        digest_path,
        f"{file_sha256(path)}  provenance.json\n".encode("ascii"),
    )
    return provenance


def summarize_campaign(output_directory: Path) -> dict[str, Any]:
    output_directory = resolve_repository_path(output_directory)
    plan, plan_hash = load_plan(output_directory)
    completed = collect_completed(
        output_directory, plan, plan_hash
    )
    summary = build_summary(plan, plan_hash, completed)
    atomic_json(output_directory / "aggregate-report.json", summary)
    write_provenance(output_directory, plan, plan_hash)
    return summary


def verify_campaign(output_directory: Path) -> dict[str, Any]:
    output_directory = resolve_repository_path(output_directory)
    plan, plan_hash = load_plan(output_directory)
    completed = collect_completed(
        output_directory, plan, plan_hash
    )
    expected_summary = build_summary(plan, plan_hash, completed)
    aggregate_path = output_directory / "aggregate-report.json"
    if read_object(aggregate_path) != expected_summary:
        raise RuntimeError("aggregate report does not reproduce")

    logical_cubes, raw22 = derive_logical_cubes()
    reproduced_novelty = build_novelty_audit(logical_cubes, raw22)
    novelty_path = output_directory / "novelty-audit.json"
    if read_object(novelty_path) != reproduced_novelty:
        raise RuntimeError("novelty audit does not reproduce")

    provenance_path = output_directory / "provenance.json"
    provenance_hash_path = output_directory / "provenance.sha256"
    expected_hash = provenance_hash_path.read_text(
        encoding="ascii"
    ).split()[0]
    if file_sha256(provenance_path) != expected_hash:
        raise RuntimeError("provenance sidecar hash mismatch")
    provenance = read_object(provenance_path)
    records, digest = artifact_inventory(output_directory)
    if (
        int(provenance["artifact_inventory"]["count"]) != len(records)
        or provenance["artifact_inventory"]["sha256"] != digest
        or provenance["immutable_plan_raw_sha256"] != plan_hash
        or provenance["aggregate_report"]["raw_sha256"]
            != file_sha256(aggregate_path)
    ):
        raise RuntimeError("artifact provenance inventory changed")
    return {
        "aggregate_report_raw_sha256": file_sha256(aggregate_path),
        "artifact_inventory_count": len(records),
        "artifact_inventory_sha256": digest,
        "assignment_visits": expected_summary["assignment_visits"],
        "best_absolute_determinant":
            expected_summary["best_absolute_determinant"],
        "novelty_audit_raw_sha256": file_sha256(novelty_path),
        "plan_raw_sha256": plan_hash,
        "provenance_raw_sha256": file_sha256(provenance_path),
        "verified": True,
    }


def self_test() -> dict[str, Any]:
    logical_cubes, raw22 = derive_logical_cubes()
    if (
        [cube.dimension for cube in logical_cubes] != [28, 27]
        or len(raw22) != 22
        or logical_cubes[1].support[:22] != raw22
    ):
        raise AssertionError("bridge derivation self-test failed")
    for cube in logical_cubes:
        inner_mask = cube.endpoint_logical_mask & (
            (1 << DIMENSION) - 1
        )
        outer_mask = cube.endpoint_logical_mask >> DIMENSION
        reconstructed = inner_mask | (outer_mask << DIMENSION)
        if reconstructed != cube.endpoint_logical_mask:
            raise AssertionError("partition mask round-trip failed")
    overlap = len(
        set(logical_cubes[0].support)
        & set(logical_cubes[1].support)
    )
    if overlap != 3:
        raise AssertionError("executed support intersection changed")
    return {
        "checks": [
            "raw bridge sizes and endpoint reconstruction",
            "deterministic transverse repair coordinates",
            "27+outer partition mask round-trip",
            "executed support intersection dimension",
        ],
        "passed": True,
        "schema_version": SCHEMA_VERSION,
        "tool": file_record(Path(__file__).resolve()),
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("self-test")
    plan = subparsers.add_parser("plan")
    plan.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    plan.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    run = subparsers.add_parser("run")
    run.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    run.add_argument("--resume", action="store_true")
    summarize = subparsers.add_parser("summarize")
    summarize.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT
    )
    verify = subparsers.add_parser("verify")
    verify.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    if arguments.command == "self-test":
        result = self_test()
    elif arguments.command == "plan":
        result = plan_campaign(arguments.output_dir, arguments.binary)
    elif arguments.command == "run":
        result = run_campaign(arguments.output_dir, arguments.resume)
    elif arguments.command == "summarize":
        result = summarize_campaign(arguments.output_dir)
    elif arguments.command == "verify":
        result = verify_campaign(arguments.output_dir)
    else:  # pragma: no cover
        raise AssertionError(arguments.command)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
