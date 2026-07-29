#!/usr/bin/env python3
"""Plan and run a resumable, shardable QD-selected exact-cube campaign.

The planner selects distance- and orientation-diverse exact-frontier matrices
exported by the corrected MAP-Elites run.  It derives four 27-entry supports
per selected center and adds exact frontier-to-frontier bridge controls of
dimension 28..32.  Every bridge is partitioned into disjoint 27-bit leaves.

The C++ fast-principal-minor evaluator remains the arithmetic authority.  This
driver pins all source/binary/input hashes, deduplicates affine cubes, gives
planned/evaluated/unrun counters distinct names, and never reuses an engine
output path.  A resumed task gets a fresh numbered attempt directory.
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
from typing import Any, Iterable, Sequence

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
    read_matrix,
    sha256_bytes,
)
from research.fast_cube_lnps import (
    affine_fingerprint,
    apply_mask,
    matrix_identity,
    popcount,
    support_bytes,
)

Coordinate = tuple[int, int]
Matrix = tuple[tuple[int, ...], ...]
ASSIGNMENTS_PER_TASK = 1 << DIMENSION
DIRECT_MODES = (
    "exact_single_gain",
    "exact_pair_score",
    "exact_pair_synergy",
    "balanced_single_synergy",
)
PREFERRED_QD_RANKS = (3, 6, 11, 16, 17)
DEPRIORITIZED_QD_RANKS = frozenset((0, 8, 13, 14, 15))
DEFAULT_QD_ARCHIVE = Path(
    "runs/map-elites-pilot-20260728-seed35011/archive"
)
DEFAULT_QD_SUMMARY = Path(
    "runs/map-elites-pilot-20260728-seed35011/summary.json"
)
DEFAULT_PRIOR_ROOT = Path(
    "runs/direct-search/fast-principal-cube"
)
DEFAULT_BINARY = Path(
    "build/research/fast_principal_cube_lnps"
)
ELITE_PATTERN = re.compile(
    r"^elite-(?P<rank>\d+)-q(?P<quotient>\d+)"
    r"-n(?P<nearest>\d+)-d(?P<distance>\d+)"
    r"-b(?P<balance>\d+)-g(?P<gram>\d+)\.matrix\.txt$"
)


@dataclass(frozen=True)
class Elite:
    id: str
    archive_rank: int
    path: Path
    matrix: Matrix
    signed_determinant: int
    raw_sha256: str
    normalized_sha256: str
    normalized_bits: int
    nearest_seed: int
    distance_bin: int
    balance_bin: int
    gram_bin: int


@dataclass(frozen=True)
class FeatureCenter:
    """The minimal interface expected by compute_features()."""

    label: str
    matrix: Matrix
    determinant: int


def repository_relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPOSITORY_ROOT))


def resolve_repository_path(path: Path) -> Path:
    resolved = (
        path.resolve()
        if path.is_absolute()
        else (REPOSITORY_ROOT / path).resolve()
    )
    resolved.relative_to(REPOSITORY_ROOT)
    return resolved


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return value


def parse_arena(output: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for key, label in (
        ("determinant", "determinant"),
        ("score", "score |det|"),
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
        result[key] = match.group(1)
    return result


def core_coordinate(coordinate: Coordinate) -> bool:
    """QD evolves the dephased 22x22 core, never row/column zero."""

    row, column = coordinate
    return 1 <= row < ORDER and 1 <= column < ORDER


def assert_core_support(support: Sequence[Coordinate], label: str) -> None:
    if len(support) != len(set(support)):
        raise ValueError(f"{label} contains duplicate coordinates")
    if not all(core_coordinate(coordinate) for coordinate in support):
        raise ValueError(
            f"{label} crosses the fixed dephased boundary"
        )


def gf2_direction_rank(support: Sequence[Coordinate]) -> int:
    """Rank entry-flip directions in the dephased 22x22 core."""

    vectors = [
        1 << ((row - 1) * (ORDER - 1) + column - 1)
        for row, column in support
    ]
    pivots: dict[int, int] = {}
    for vector in vectors:
        reduced = vector
        while reduced:
            pivot = reduced.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = reduced
                break
            reduced ^= pivots[pivot]
    return len(pivots)


def load_frontier_elites(
    archive_directory: Path,
) -> list[Elite]:
    elites: list[Elite] = []
    normalized_seen: set[str] = set()
    for path in sorted(archive_directory.glob("*.matrix.txt")):
        match = ELITE_PATTERN.match(path.name)
        if match is None:
            raise ValueError(f"unrecognized QD archive filename: {path}")
        matrix = read_matrix(path)
        determinant = bareiss_determinant(matrix)
        if abs(determinant) != FRONTIER:
            continue
        # Correct QD archive outputs are dephased.  Enforcing this prevents
        # compute_features() from silently selecting fixed boundary entries.
        if any(matrix[0][column] != 1 for column in range(ORDER)):
            raise ValueError(f"QD elite first row is not dephased: {path}")
        if any(matrix[row][0] != 1 for row in range(ORDER)):
            raise ValueError(f"QD elite first column is not dephased: {path}")
        raw, normalized, bits = matrix_identity(matrix)
        if normalized in normalized_seen:
            continue
        normalized_seen.add(normalized)
        fields = {
            key: int(value)
            for key, value in match.groupdict().items()
        }
        elites.append(
            Elite(
                id=f"qd-{fields['rank']:03d}-{normalized[:12]}",
                archive_rank=fields["rank"],
                path=path,
                matrix=matrix,
                signed_determinant=determinant,
                raw_sha256=raw,
                normalized_sha256=normalized,
                normalized_bits=bits,
                nearest_seed=fields["nearest"],
                distance_bin=fields["distance"],
                balance_bin=fields["balance"],
                gram_bin=fields["gram"],
            )
        )
    if not elites:
        raise ValueError("QD archive contains no unique frontier elites")
    return elites


def hamming(first: Elite, second: Elite) -> int:
    return popcount(first.normalized_bits ^ second.normalized_bits)


def select_centers(
    elites: Sequence[Elite], requested: int
) -> tuple[list[Elite], dict[str, Any]]:
    if requested <= 0 or requested > len(elites):
        raise ValueError(
            f"center count must lie within 1..{len(elites)}"
        )
    selected: list[Elite] = []
    # Begin with the independently audited clean HT-diverse roots.  These
    # avoid replaying the known bridge endpoints when enough alternatives
    # exist.
    by_rank = {elite.archive_rank: elite for elite in elites}
    for rank in PREFERRED_QD_RANKS:
        elite = by_rank.get(rank)
        if elite is not None and len(selected) < requested:
            selected.append(elite)
    # Then cover any QD seed lineage not represented above.
    for nearest_seed in sorted({elite.nearest_seed for elite in elites}):
        if any(
            elite.nearest_seed == nearest_seed for elite in selected
        ):
            continue
        candidates = [
            elite
            for elite in elites
            if (
                elite.nearest_seed == nearest_seed
                and elite.archive_rank not in DEPRIORITIZED_QD_RANKS
            )
        ]
        if not candidates:
            candidates = [
                elite
                for elite in elites
                if elite.nearest_seed == nearest_seed
            ]
        choice = min(
            candidates,
            key=lambda elite: (
                -elite.distance_bin,
                elite.normalized_sha256,
            ),
        )
        if choice not in selected and len(selected) < requested:
            selected.append(choice)
    while len(selected) < requested:
        remaining = [
            elite for elite in elites if elite not in selected
        ]
        choice = min(
            remaining,
            key=lambda elite: (
                elite.archive_rank in DEPRIORITIZED_QD_RANKS,
                -min(hamming(elite, chosen) for chosen in selected),
                -elite.distance_bin,
                elite.nearest_seed,
                elite.normalized_sha256,
            ),
        )
        selected.append(choice)
    minimum_distance = (
        min(
            hamming(first, second)
            for index, first in enumerate(selected)
            for second in selected[index + 1 :]
        )
        if len(selected) > 1
        else 0
    )
    return selected, {
        "algorithm":
            "audited-roots-then-seed-coverage-and-max-min-hamming-v1",
        "available_unique_frontier_elites": len(elites),
        "minimum_selected_pairwise_normalized_hamming_distance":
            minimum_distance,
        "requested_centers": requested,
        "selected_nearest_seed_counts": {
            str(nearest): sum(
                elite.nearest_seed == nearest for elite in selected
            )
            for nearest in sorted(
                {elite.nearest_seed for elite in selected}
            )
        },
    }


def pair_endpoint_rank(
    pairs: Iterable[tuple[Coordinate, Coordinate]],
) -> list[Coordinate]:
    result: list[Coordinate] = []
    seen: set[Coordinate] = set()
    for first, second in pairs:
        for coordinate in (first, second):
            if core_coordinate(coordinate) and coordinate not in seen:
                seen.add(coordinate)
                result.append(coordinate)
    return result


def fill_support(
    ranked: Sequence[Coordinate],
    *,
    dimension: int = DIMENSION,
    row_cap: int = 2,
    column_cap: int = 2,
) -> tuple[Coordinate, ...]:
    selected: list[Coordinate] = []
    selected_set: set[Coordinate] = set()

    def fill(cap_rows: int | None, cap_columns: int | None) -> None:
        row_counts = [0] * ORDER
        column_counts = [0] * ORDER
        for row, column in selected:
            row_counts[row] += 1
            column_counts[column] += 1
        for coordinate in ranked:
            if len(selected) == dimension:
                return
            if (
                coordinate in selected_set
                or not core_coordinate(coordinate)
            ):
                continue
            row, column = coordinate
            if (
                cap_rows is not None
                and row_counts[row] >= cap_rows
            ):
                continue
            if (
                cap_columns is not None
                and column_counts[column] >= cap_columns
            ):
                continue
            selected.append(coordinate)
            selected_set.add(coordinate)
            row_counts[row] += 1
            column_counts[column] += 1

    fill(row_cap, column_cap)
    if len(selected) < dimension:
        fill(3, 3)
    if len(selected) < dimension:
        fill(None, None)
    if len(selected) != dimension:
        raise ValueError(
            f"could select only {len(selected)} of {dimension} coordinates"
        )
    support = tuple(selected)
    assert_core_support(support, "derived support")
    return support


def derive_supports(
    elite: Elite, features: Features
) -> dict[str, tuple[Coordinate, ...]]:
    singles = [
        coordinate
        for coordinate in features.single_rank
        if core_coordinate(coordinate)
    ]
    pair_score = pair_endpoint_rank(features.pair_score_rank)
    pair_synergy = pair_endpoint_rank(features.pair_synergy_rank)
    supports: dict[str, tuple[Coordinate, ...]] = {
        "exact_single_gain": fill_support(singles),
        "exact_pair_score": fill_support(pair_score),
        "exact_pair_synergy": fill_support(pair_synergy),
    }
    usage = {
        coordinate: sum(
            coordinate in support for support in supports.values()
        )
        for coordinate in (
            (row, column)
            for row in range(1, ORDER)
            for column in range(1, ORDER)
        )
    }
    single_rank = {
        coordinate: index for index, coordinate in enumerate(singles)
    }
    synergy_rank = {
        coordinate: index
        for index, coordinate in enumerate(pair_synergy)
    }
    hybrid_rank = sorted(
        usage,
        key=lambda coordinate: (
            usage[coordinate],
            min(
                single_rank.get(coordinate, ORDER * ORDER),
                synergy_rank.get(coordinate, ORDER * ORDER),
            ),
            single_rank.get(coordinate, ORDER * ORDER)
            + synergy_rank.get(coordinate, ORDER * ORDER),
            coordinate,
        ),
    )
    supports["balanced_single_synergy"] = fill_support(hybrid_rank)
    if tuple(supports) != DIRECT_MODES:
        raise RuntimeError("direct support mode order changed")
    for mode, support in supports.items():
        if len(support) != DIMENSION:
            raise RuntimeError(f"{elite.id}/{mode} is not 27-dimensional")
        assert_core_support(support, f"{elite.id}/{mode}")
        if gf2_direction_rank(support) != DIMENSION:
            raise RuntimeError(f"{elite.id}/{mode} failed GF(2) rank gate")
    return supports


def differing_coordinates(
    first: Elite, second: Elite
) -> tuple[Coordinate, ...]:
    differences = tuple(
        (row, column)
        for row in range(ORDER)
        for column in range(ORDER)
        if first.matrix[row][column] != second.matrix[row][column]
    )
    assert_core_support(
        differences, f"bridge {first.id}<->{second.id}"
    )
    if apply_mask(
        first.matrix, differences, (1 << len(differences)) - 1
    ) != second.matrix:
        raise RuntimeError("bridge coordinate application missed endpoint")
    return differences


def select_bridge_controls(
    elites: Sequence[Elite],
    requested: int,
    minimum_dimension: int,
    maximum_dimension: int,
) -> list[tuple[Elite, Elite, tuple[Coordinate, ...]]]:
    candidates: list[
        tuple[int, str, str, Elite, Elite, tuple[Coordinate, ...]]
    ] = []
    for first_index, first in enumerate(elites):
        for second in elites[first_index + 1 :]:
            support = differing_coordinates(first, second)
            dimension = len(support)
            if minimum_dimension <= dimension <= maximum_dimension:
                candidates.append(
                    (
                        dimension,
                        first.normalized_sha256,
                        second.normalized_sha256,
                        first,
                        second,
                        support,
                    )
                )
    candidates.sort(key=lambda item: item[:3])
    selected: list[
        tuple[Elite, Elite, tuple[Coordinate, ...]]
    ] = []
    used_endpoints: set[str] = set()
    used_supports: set[tuple[Coordinate, ...]] = set()
    for _, _, _, first, second, support in candidates:
        if (
            first.normalized_sha256 in used_endpoints
            or second.normalized_sha256 in used_endpoints
            or support in used_supports
        ):
            continue
        selected.append((first, second, support))
        used_endpoints.add(first.normalized_sha256)
        used_endpoints.add(second.normalized_sha256)
        used_supports.add(support)
        if len(selected) == requested:
            break
    if len(selected) != requested:
        raise ValueError(
            f"found only {len(selected)} disjoint bridge controls in "
            f"dimensions {minimum_dimension}..{maximum_dimension}"
        )
    return selected


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
        rerooted = apply_mask(fixed, support, offset)
        determinant = bareiss_determinant(rerooted)
        if determinant != 0:
            return offset, rerooted, determinant
    raise RuntimeError("could not find a nonsingular one/two-bit reroot")


def load_prior_fingerprints(
    prior_root: Path,
) -> dict[str, dict[str, Any]]:
    """Load fingerprints backed by completed evaluator artifacts only."""

    fingerprints: dict[str, dict[str, Any]] = {}
    if not prior_root.exists():
        return fingerprints

    def retain(
        fingerprint: str,
        evidence_paths: Sequence[Path],
    ) -> None:
        evidence = [
            {
                "path": repository_relative(path),
                "raw_sha256": sha256_bytes(path.read_bytes()),
            }
            for path in evidence_paths
        ]
        fingerprints[fingerprint] = {
            "evidence": evidence,
            "primary": evidence[0]["path"],
        }

    for path in prior_root.rglob("lnps-summary.json"):
        payload = read_json(path)
        fingerprint = payload.get("fingerprint_sha256")
        report_path = path.parent / "report.json"
        if not report_path.is_file():
            continue
        report = read_json(report_path)
        if (
            fingerprint
            and int(payload.get("assignments", 0))
                == ASSIGNMENTS_PER_TASK
            and report.get("complete") is True
            and int(report.get("assignments", 0))
                == ASSIGNMENTS_PER_TASK
        ):
            retain(str(fingerprint), (path, report_path))
    for path in prior_root.rglob("partition-summary.json"):
        payload = read_json(path)
        fingerprint = payload.get(
            "leaf_affine_fingerprint_sha256"
        )
        report_path = path.parent / "report.json"
        if not report_path.is_file():
            continue
        report = read_json(report_path)
        if (
            fingerprint
            and int(payload.get("assignments", 0))
                == ASSIGNMENTS_PER_TASK
            and report.get("complete") is True
            and int(report.get("assignments", 0))
                == ASSIGNMENTS_PER_TASK
        ):
            retain(str(fingerprint), (path, report_path))
    for path in prior_root.rglob("manifest.json"):
        payload = read_json(path)
        if (
            payload.get("method")
            != "deterministic-mixed-support-cube-batch-v1"
        ):
            continue
        for record in payload.get("runs", []):
            run_id = str(record.get("id", ""))
            report_path = path.parent / run_id / "report.json"
            if not report_path.is_file():
                continue
            report = read_json(report_path)
            if (
                report.get("complete") is True
                and int(report.get("assignments", 0))
                    == ASSIGNMENTS_PER_TASK
            ):
                fingerprint = str(record["fingerprint_sha256"])
                retain(fingerprint, (path, report_path))
    return fingerprints


def elite_record(elite: Elite) -> dict[str, Any]:
    return {
        "absolute_determinant": str(abs(elite.signed_determinant)),
        "archive_rank": elite.archive_rank,
        "balance_bin": elite.balance_bin,
        "distance_bin": elite.distance_bin,
        "gram_bin": elite.gram_bin,
        "id": elite.id,
        "nearest_seed": elite.nearest_seed,
        "normalized_sha256": elite.normalized_sha256,
        "path": repository_relative(elite.path),
        "raw_sha256": elite.raw_sha256,
        "signed_determinant": str(elite.signed_determinant),
    }


def build_plan(arguments: argparse.Namespace) -> dict[str, Any]:
    output_directory = resolve_repository_path(arguments.output_dir)
    if output_directory.exists():
        raise FileExistsError(
            f"refusing to overwrite campaign: {output_directory}"
        )
    binary = resolve_repository_path(arguments.binary)
    engine_source = resolve_repository_path(
        Path("research/fast_principal_cube.cpp")
    )
    qd_archive = resolve_repository_path(arguments.qd_archive)
    qd_summary = resolve_repository_path(arguments.qd_summary)
    prior_root = resolve_repository_path(arguments.prior_root)
    if not binary.is_file():
        raise FileNotFoundError(f"missing evaluator binary: {binary}")
    if not qd_summary.is_file():
        raise FileNotFoundError(f"missing QD summary: {qd_summary}")
    summary = read_json(qd_summary)
    if (
        summary.get("complete") is not True
        or int(summary.get("promotions", -1)) != 0
    ):
        raise ValueError("QD summary is not the retained complete pilot")

    elites = load_frontier_elites(qd_archive)
    centers, selection = select_centers(elites, arguments.center_count)
    bridges = select_bridge_controls(
        elites,
        arguments.bridge_controls,
        arguments.bridge_min_dimension,
        arguments.bridge_max_dimension,
    )
    prior_fingerprints = load_prior_fingerprints(prior_root)

    logical_cubes: list[dict[str, Any]] = []
    for center in centers:
        print(f"feature-build center={center.id}", flush=True)
        features = compute_features(
            FeatureCenter(
                center.id, center.matrix, center.signed_determinant
            )
        )  # type: ignore[arg-type]
        for mode, support in derive_supports(center, features).items():
            logical_cubes.append(
                {
                    "center_id": center.id,
                    "dimension": DIMENSION,
                    "endpoint_id": None,
                    "full_affine_fingerprint_sha256":
                        affine_fingerprint(center.matrix, support),
                    "id": f"direct-{center.id}-{mode}",
                    "kind": "direct",
                    "mode": mode,
                    "gf2_direction_rank": gf2_direction_rank(support),
                    "support": [
                        [row + 1, column + 1]
                        for row, column in support
                    ],
                    "support_sha256": sha256_bytes(
                        support_bytes(support)
                    ),
                }
            )
    for index, (first, second, support) in enumerate(bridges):
        if not DIMENSION < len(support) <= 32:
            raise RuntimeError("bridge planner emitted unsupported dimension")
        if (
            abs(bareiss_determinant(first.matrix)) != FRONTIER
            or abs(bareiss_determinant(second.matrix)) != FRONTIER
        ):
            raise RuntimeError("bridge parent failed Bareiss preflight")
        logical_cubes.append(
            {
                "center_id": first.id,
                "dimension": len(support),
                "endpoint_id": second.id,
                "endpoint_exact_match": True,
                "endpoint_signed_determinant":
                    str(second.signed_determinant),
                "full_affine_fingerprint_sha256":
                    affine_fingerprint(first.matrix, support),
                "id": (
                    f"bridge-{index:02d}-d{len(support)}-"
                    f"{first.id}-{second.id}"
                ),
                "kind": "bridge",
                "mode": "exact-frontier-bridge-control",
                "gf2_direction_rank": gf2_direction_rank(support),
                "support": [
                    [row + 1, column + 1]
                    for row, column in support
                ],
                "support_sha256": sha256_bytes(
                    support_bytes(support)
                ),
            }
        )

    elite_by_id = {elite.id: elite for elite in elites}
    candidate_tasks: list[dict[str, Any]] = []
    campaign_fingerprints: dict[str, str] = {}
    scheduled_index = 0
    for logical in logical_cubes:
        center = elite_by_id[str(logical["center_id"])]
        support = tuple(
            (int(row) - 1, int(column) - 1)
            for row, column in logical["support"]
        )
        assert_core_support(support, str(logical["id"]))
        if gf2_direction_rank(support) != len(support):
            raise RuntimeError(
                f"{logical['id']} failed the dephased GF(2) rank gate"
            )
        inner = support[:DIMENSION]
        outer = support[DIMENSION:]
        leaf_count = 1 << len(outer)
        logical["candidate_task_count"] = leaf_count
        logical["partition"] = {
            "disjoint": True,
            "fixed_outer_bits": len(outer),
            "inner_dimension": len(inner),
            "leaf_count": leaf_count,
            "total_logical_assignments": 1 << len(support),
        }
        if len(inner) != DIMENSION:
            raise RuntimeError("all evaluator tasks must be 27-bit")
        for outer_mask in range(leaf_count):
            fixed = apply_mask(center.matrix, outer, outer_mask)
            reroot, task_start, signed = reroot_offset(fixed, inner)
            fingerprint = affine_fingerprint(task_start, inner)
            task_id = (
                f"{logical['id']}-leaf-{outer_mask:02d}"
                if leaf_count > 1
                else str(logical["id"])
            )
            disposition = "scheduled"
            duplicate_of: str | None = None
            if fingerprint in prior_fingerprints:
                disposition = "duplicate_prior"
                duplicate_of = str(
                    prior_fingerprints[fingerprint]["primary"]
                )
            elif fingerprint in campaign_fingerprints:
                disposition = "duplicate_campaign"
                duplicate_of = campaign_fingerprints[fingerprint]
            else:
                campaign_fingerprints[fingerprint] = task_id
            endpoint_expected = (
                logical["kind"] == "bridge"
                and outer_mask == leaf_count - 1
            )
            task = {
                "affine_fingerprint_sha256": fingerprint,
                "assignment_visits": ASSIGNMENTS_PER_TASK,
                "dimension": DIMENSION,
                "disposition": disposition,
                "duplicate_of": duplicate_of,
                "endpoint_expected": endpoint_expected,
                "expected_endpoint_engine_mask_decimal": (
                    str(((1 << DIMENSION) - 1) ^ reroot)
                    if endpoint_expected
                    else None
                ),
                "id": task_id,
                "kind": logical["kind"],
                "logical_cube_id": logical["id"],
                "outer_mask_decimal": str(outer_mask),
                "reroot_xor_mask_decimal": str(reroot),
                "scheduled_index": (
                    scheduled_index
                    if disposition == "scheduled"
                    else None
                ),
                "start_parsed_matrix_sha256": sha256_bytes(
                    canonical_matrix_bytes(task_start)
                ),
                "start_signed_determinant": str(signed),
            }
            candidate_tasks.append(task)
            if disposition == "scheduled":
                scheduled_index += 1

    direct_logical = sum(
        logical["kind"] == "direct" for logical in logical_cubes
    )
    bridge_logical = sum(
        logical["kind"] == "bridge" for logical in logical_cubes
    )
    scheduled = [
        task
        for task in candidate_tasks
        if task["disposition"] == "scheduled"
    ]
    counts = {
        "bridge_logical_cube_count": bridge_logical,
        "candidate_evaluator_task_count": len(candidate_tasks),
        "direct_logical_cube_count": direct_logical,
        "duplicate_campaign_task_count": sum(
            task["disposition"] == "duplicate_campaign"
            for task in candidate_tasks
        ),
        "duplicate_prior_task_count": sum(
            task["disposition"] == "duplicate_prior"
            for task in candidate_tasks
        ),
        "logical_cube_count": len(logical_cubes),
        "planned_assignment_visits": sum(
            int(task["assignment_visits"]) for task in scheduled
        ),
        "planned_bridge_assignment_visits": sum(
            int(task["assignment_visits"])
            for task in scheduled
            if task["kind"] == "bridge"
        ),
        "planned_direct_assignment_visits": sum(
            int(task["assignment_visits"])
            for task in scheduled
            if task["kind"] == "direct"
        ),
        "scheduled_unique_evaluator_task_count": len(scheduled),
    }
    driver = Path(__file__).resolve()
    runtime_dependencies = [
        resolve_repository_path(path)
        for path in (
            Path("research/fast_cube_batch.py"),
            Path("research/fast_cube_lnps.py"),
            Path("maxdet/exact.py"),
            Path("arena"),
        )
    ]
    prior_ledger = [
        {
            "affine_fingerprint_sha256": fingerprint,
            **record,
        }
        for fingerprint, record in sorted(prior_fingerprints.items())
    ]
    prior_ledger_bytes = (
        json.dumps(
            prior_ledger,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("utf-8")
    direct_overlap_ledger: list[dict[str, Any]] = []
    direct_cubes = [
        logical
        for logical in logical_cubes
        if logical["kind"] == "direct"
    ]
    for first_index, first in enumerate(direct_cubes):
        first_support = {
            (int(row), int(column))
            for row, column in first["support"]
        }
        for second in direct_cubes[first_index + 1 :]:
            if first["center_id"] != second["center_id"]:
                continue
            second_support = {
                (int(row), int(column))
                for row, column in second["support"]
            }
            intersection = len(first_support & second_support)
            direct_overlap_ledger.append(
                {
                    "first_logical_cube_id": first["id"],
                    "intersection_dimension": intersection,
                    "pairwise_assignment_overlap": 1 << intersection,
                    "second_logical_cube_id": second["id"],
                }
            )
    plan = {
        "claim_boundary": [
            "Assignment visits count exhaustive evaluator visits and are "
            "not claimed to be globally unique matrices.",
            "Affine fingerprints deduplicate identical fixed/free entry "
            "cubes, not H/HT-equivalence classes.",
            "The direct-support overlap ledger is pairwise evidence, not an "
            "inclusion-exclusion or unique-assignment count.",
            "Only an independently arena-verified strict score increase is "
            "a promotion.",
        ],
        "counts": counts,
        "direct_support_modes": list(DIRECT_MODES),
        "direct_support_pairwise_overlap_ledger":
            direct_overlap_ledger,
        "engine": "fast-principal-minor-entry-cube-v1",
        "frontier": str(FRONTIER),
        "logical_cubes": logical_cubes,
        "method": "qd-selected-shardable-exact-cubes-v1",
        "prior_fingerprint_count": len(prior_fingerprints),
        "prior_fingerprint_ledger": prior_ledger,
        "prior_fingerprint_ledger_raw_sha256": sha256_bytes(
            prior_ledger_bytes
        ),
        "provenance": {
            "driver": repository_relative(driver),
            "driver_raw_sha256": sha256_bytes(driver.read_bytes()),
            "engine_binary": repository_relative(binary),
            "engine_binary_raw_sha256": sha256_bytes(
                binary.read_bytes()
            ),
            "engine_source": repository_relative(engine_source),
            "engine_source_raw_sha256": sha256_bytes(
                engine_source.read_bytes()
            ),
            "prior_root": repository_relative(prior_root),
            "qd_archive": repository_relative(qd_archive),
            "qd_summary": repository_relative(qd_summary),
            "qd_summary_raw_sha256": sha256_bytes(
                qd_summary.read_bytes()
            ),
            "runtime_dependencies": [
                {
                    "path": repository_relative(path),
                    "raw_sha256": sha256_bytes(path.read_bytes()),
                }
                for path in runtime_dependencies
            ],
        },
        "qd_frontier_elites": [
            elite_record(elite) for elite in elites
        ],
        "schema_version": 1,
        "selected_center_ids": [elite.id for elite in centers],
        "selection": selection,
        "tasks": candidate_tasks,
    }
    output_directory.mkdir(parents=True)
    plan_path = output_directory / "plan.json"
    atomic_json(plan_path, plan)
    plan_hash = sha256_bytes(plan_path.read_bytes())
    atomic_write(
        output_directory / "plan.sha256",
        f"{plan_hash}  plan.json\n".encode("ascii"),
    )
    return plan


def load_plan(output_directory: Path) -> tuple[dict[str, Any], str]:
    plan_path = output_directory / "plan.json"
    sidecar_path = output_directory / "plan.sha256"
    if not plan_path.is_file() or not sidecar_path.is_file():
        raise FileNotFoundError("campaign plan or hash sidecar is missing")
    expected = sidecar_path.read_text().split()[0]
    actual = sha256_bytes(plan_path.read_bytes())
    if expected != actual:
        raise RuntimeError("campaign plan SHA-256 mismatch")
    plan = read_json(plan_path)
    provenance = plan["provenance"]
    for path_key, hash_key in (
        ("driver", "driver_raw_sha256"),
        ("engine_binary", "engine_binary_raw_sha256"),
        ("engine_source", "engine_source_raw_sha256"),
        ("qd_summary", "qd_summary_raw_sha256"),
    ):
        path = resolve_repository_path(Path(provenance[path_key]))
        if sha256_bytes(path.read_bytes()) != provenance[hash_key]:
            raise RuntimeError(f"pinned provenance changed: {path_key}")
    for record in provenance["runtime_dependencies"]:
        path = resolve_repository_path(Path(record["path"]))
        if sha256_bytes(path.read_bytes()) != record["raw_sha256"]:
            raise RuntimeError(
                f"pinned runtime dependency changed: {record['path']}"
            )
    ledger_bytes = (
        json.dumps(
            plan["prior_fingerprint_ledger"],
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("utf-8")
    if sha256_bytes(ledger_bytes) != plan[
        "prior_fingerprint_ledger_raw_sha256"
    ]:
        raise RuntimeError("prior fingerprint ledger digest mismatch")
    for record in plan["prior_fingerprint_ledger"]:
        for evidence in record["evidence"]:
            path = resolve_repository_path(Path(evidence["path"]))
            if sha256_bytes(path.read_bytes()) != evidence["raw_sha256"]:
                raise RuntimeError(
                    "prior fingerprint evidence changed: "
                    f"{evidence['path']}"
                )
    return plan, actual


def reconstruct_task(
    plan: dict[str, Any], task: dict[str, Any]
) -> tuple[Matrix, tuple[Coordinate, ...], dict[str, Any]]:
    centers = {
        str(record["id"]): record
        for record in plan["qd_frontier_elites"]
    }
    logicals = {
        str(record["id"]): record
        for record in plan["logical_cubes"]
    }
    logical = logicals[str(task["logical_cube_id"])]
    center_record = centers[str(logical["center_id"])]
    center_path = resolve_repository_path(Path(center_record["path"]))
    center = read_matrix(center_path)
    if sha256_bytes(canonical_matrix_bytes(center)) != center_record[
        "raw_sha256"
    ]:
        raise RuntimeError("center matrix changed after planning")
    if abs(bareiss_determinant(center)) != FRONTIER:
        raise RuntimeError("center failed Bareiss at task reconstruction")
    support = tuple(
        (int(row) - 1, int(column) - 1)
        for row, column in logical["support"]
    )
    assert_core_support(support, str(logical["id"]))
    inner = support[:DIMENSION]
    outer = support[DIMENSION:]
    fixed = apply_mask(center, outer, int(task["outer_mask_decimal"]))
    reroot = int(task["reroot_xor_mask_decimal"])
    start = apply_mask(fixed, inner, reroot)
    if bareiss_determinant(start) == 0:
        raise RuntimeError("planned task start became singular")
    if sha256_bytes(canonical_matrix_bytes(start)) != task[
        "start_parsed_matrix_sha256"
    ]:
        raise RuntimeError("planned task start hash mismatch")
    if affine_fingerprint(start, inner) != task[
        "affine_fingerprint_sha256"
    ]:
        raise RuntimeError("planned task affine fingerprint mismatch")
    return start, inner, logical


def replay_engine_result(
    plan: dict[str, Any],
    task: dict[str, Any],
    report_path: Path,
    best_path: Path,
) -> dict[str, Any]:
    """Independently Bareiss-replay every retained engine record."""

    start, support, logical = reconstruct_task(plan, task)
    report = read_json(report_path)
    if (
        report.get("complete") is not True
        or report.get("all_assignments_bound_checked") is not True
        or int(report.get("dimension", 0)) != DIMENSION
        or int(report.get("assignments", 0))
            != ASSIGNMENTS_PER_TASK
        or report.get("start_parsed_matrix_sha256")
            != task["start_parsed_matrix_sha256"]
        or report.get("coordinate_file_raw_sha256")
            != sha256_bytes(support_bytes(support))
    ):
        raise RuntimeError(
            f"engine report/provenance mismatch: {task['id']}"
        )

    best_mask = int(report["best_mask_decimal"])
    if not 0 <= best_mask < (1 << DIMENSION):
        raise RuntimeError("engine best mask lies outside its cube")
    reconstructed_best = apply_mask(start, support, best_mask)
    artifact_best = read_matrix(best_path)
    if reconstructed_best != artifact_best:
        raise RuntimeError("engine best artifact does not match best mask")
    best_exact = bareiss_determinant(reconstructed_best)
    best_raw_sha256 = sha256_bytes(
        canonical_matrix_bytes(reconstructed_best)
    )
    if (
        best_exact != int(report["best_signed_determinant"])
        or abs(best_exact) != int(report["best_absolute_determinant"])
        or report["best_matrix_sha256"] != best_raw_sha256
        or report["output_raw_sha256"] != best_raw_sha256
    ):
        raise RuntimeError("engine best failed independent Bareiss replay")

    top_replayed = 0
    for record in report.get("top_k_candidates", []):
        mask = int(record["mask_decimal"])
        if not 0 < mask < (1 << DIMENSION):
            raise RuntimeError("top-k mask lies outside its cube")
        exact = bareiss_determinant(apply_mask(start, support, mask))
        if (
            exact != int(record["signed_determinant"])
            or abs(exact) != int(record["absolute_determinant"])
        ):
            raise RuntimeError("top-k record failed Bareiss replay")
        top_replayed += 1
    if top_replayed != int(report.get("top_k_captured", 0)):
        raise RuntimeError("top-k replay count disagrees with report")

    returned_tie_masks = [
        int(mask)
        for mask in report["frontier_tie_masks_decimal"]
    ]
    for mask in returned_tie_masks:
        if not 0 < mask < (1 << DIMENSION):
            raise RuntimeError("returned tie mask lies outside its cube")
        if abs(bareiss_determinant(
            apply_mask(start, support, mask)
        )) != FRONTIER:
            raise RuntimeError("returned tie failed Bareiss replay")

    start_exact = bareiss_determinant(start)
    start_is_frontier = abs(start_exact) == FRONTIER
    engine_nonzero_ties = int(report["frontier_nonzero_ties"])
    engine_frontier_assignments = (
        engine_nonzero_ties + int(start_is_frontier)
    )
    logical_nonzero_frontier_assignments = (
        engine_frontier_assignments
        - int(int(task["outer_mask_decimal"]) == 0)
    )

    endpoint_confirmed: bool | None = None
    endpoint_listed: bool | None = None
    if task["endpoint_expected"]:
        expected_mask = int(
            task["expected_endpoint_engine_mask_decimal"]
        )
        expected = apply_mask(start, support, expected_mask)
        centers = {
            str(record["id"]): record
            for record in plan["qd_frontier_elites"]
        }
        endpoint_record = centers[str(logical["endpoint_id"])]
        endpoint = read_matrix(
            resolve_repository_path(Path(endpoint_record["path"]))
        )
        endpoint_exact = bareiss_determinant(expected)
        endpoint_confirmed = (
            expected == endpoint
            and endpoint_exact
                == int(endpoint_record["signed_determinant"])
            and abs(endpoint_exact) == FRONTIER
        )
        if not endpoint_confirmed:
            raise RuntimeError(
                f"bridge endpoint control failed: {task['id']}"
            )
        endpoint_listed = expected_mask in set(returned_tie_masks)

    return {
        "best_absolute_determinant": str(abs(best_exact)),
        "best_matrix_raw_sha256": best_raw_sha256,
        "best_signed_determinant": str(best_exact),
        "endpoint_control_confirmed": endpoint_confirmed,
        "endpoint_listed_in_returned_tie_masks": endpoint_listed,
        "engine_frontier_assignment_occurrences":
            engine_frontier_assignments,
        "engine_nonzero_frontier_tie_count": engine_nonzero_ties,
        "engine_reported_tie_masks_truncated": bool(
            report["frontier_tie_masks_truncated"]
        ),
        "logical_nonzero_frontier_assignment_occurrences":
            logical_nonzero_frontier_assignments,
        "replayed_best_records": 1,
        "replayed_returned_tie_masks": len(returned_tie_masks),
        "replayed_top_k_records": top_replayed,
    }


def task_attempt_directory(task_directory: Path) -> Path:
    maximum = 0
    if task_directory.exists():
        for child in task_directory.iterdir():
            match = re.fullmatch(r"attempt-(\d{3})", child.name)
            if match is not None:
                maximum = max(maximum, int(match.group(1)))
    attempt = task_directory / f"attempt-{maximum + 1:03d}"
    attempt.mkdir(parents=True, exist_ok=False)
    return attempt


def validate_completed_task(
    path: Path,
    task: dict[str, Any],
    plan: dict[str, Any],
    plan_hash: str,
) -> dict[str, Any]:
    payload = read_json(path)
    if (
        payload.get("complete") is not True
        or payload.get("plan_raw_sha256") != plan_hash
        or payload.get("task_id") != task["id"]
        or payload.get("affine_fingerprint_sha256")
            != task["affine_fingerprint_sha256"]
        or int(payload.get("assignment_visits", 0))
            != ASSIGNMENTS_PER_TASK
    ):
        raise RuntimeError(f"invalid completed task marker: {path}")
    report_path = resolve_repository_path(
        Path(payload["engine_report"])
    )
    if sha256_bytes(report_path.read_bytes()) != payload[
        "engine_report_raw_sha256"
    ]:
        raise RuntimeError(f"completed task report changed: {path}")
    best_path = resolve_repository_path(
        Path(payload["best_matrix_artifact"])
    )
    if sha256_bytes(best_path.read_bytes()) != payload[
        "best_matrix_artifact_raw_sha256"
    ]:
        raise RuntimeError(f"completed best artifact changed: {path}")
    replay = replay_engine_result(
        plan, task, report_path, best_path
    )
    report = read_json(report_path)
    best_engine_mask = int(report["best_mask_decimal"])
    best_logical_mask = (
        (
            best_engine_mask
            ^ int(task["reroot_xor_mask_decimal"])
        )
        | (int(task["outer_mask_decimal"]) << DIMENSION)
    )
    if (
        payload.get("best_engine_mask_decimal")
            != str(best_engine_mask)
        or payload.get("best_logical_mask_decimal")
            != str(best_logical_mask)
        or payload.get("logical_cube_id")
            != task["logical_cube_id"]
    ):
        raise RuntimeError(
            f"completed marker mask/logical-cube mismatch: {path}"
        )
    for key in (
        "best_absolute_determinant",
        "best_matrix_raw_sha256",
        "best_signed_determinant",
        "endpoint_control_confirmed",
        "endpoint_listed_in_returned_tie_masks",
        "engine_frontier_assignment_occurrences",
        "engine_nonzero_frontier_tie_count",
        "engine_reported_tie_masks_truncated",
        "logical_nonzero_frontier_assignment_occurrences",
        "replayed_best_records",
        "replayed_returned_tie_masks",
        "replayed_top_k_records",
    ):
        if payload.get(key) != replay[key]:
            raise RuntimeError(
                f"completed marker disagrees with report: {path} ({key})"
            )
    if payload.get("kind") != task["kind"]:
        raise RuntimeError(f"completed marker kind mismatch: {path}")
    promotion = payload.get("promotion")
    if int(replay["best_absolute_determinant"]) > FRONTIER:
        if not isinstance(promotion, dict):
            raise RuntimeError("strict promotion marker is missing")
        verification_path = resolve_repository_path(
            Path(promotion["arena_verification_artifact"])
        )
        if sha256_bytes(verification_path.read_bytes()) != promotion[
            "arena_verification_artifact_raw_sha256"
        ]:
            raise RuntimeError("arena verification artifact changed")
        parsed = parse_arena(verification_path.read_text())
        for key in (
            "determinant",
            "score",
            "normalized_sha256",
            "receipt_sha256",
        ):
            if promotion.get(key) != parsed[key]:
                raise RuntimeError(
                    f"promotion marker disagrees with arena output: {key}"
                )
        if int(parsed["score"]) != int(
            replay["best_absolute_determinant"]
        ):
            raise RuntimeError("promotion score disagrees with engine")
        if promotion.get(
            "best_matrix_artifact_raw_sha256"
        ) != replay["best_matrix_raw_sha256"]:
            raise RuntimeError("promotion artifact hash disagrees with best")
    elif promotion is not None:
        raise RuntimeError("non-promotion task has promotion metadata")
    return payload


def write_shard_report(
    path: Path,
    plan_hash: str,
    shard_index: int,
    shard_count: int,
    selected_tasks: Sequence[dict[str, Any]],
    completed_before: int,
    completed_now: Sequence[dict[str, Any]],
    started: float,
) -> None:
    all_completed = completed_before + len(completed_now)
    atomic_json(
        path,
        {
            "assignment_visits_completed_this_invocation": sum(
                int(record["assignment_visits"])
                for record in completed_now
            ),
            "complete": all_completed == len(selected_tasks),
            "completed_before_this_invocation": completed_before,
            "completed_this_invocation": len(completed_now),
            "elapsed_seconds": round(time.monotonic() - started, 6),
            "plan_raw_sha256": plan_hash,
            "schema_version": 1,
            "selected_assignment_visits": sum(
                int(task["assignment_visits"])
                for task in selected_tasks
            ),
            "selected_task_count": len(selected_tasks),
            "shard_count": shard_count,
            "shard_index": shard_index,
            "unrun_selected_task_count":
                len(selected_tasks) - all_completed,
        },
    )


def run_shard(arguments: argparse.Namespace) -> dict[str, Any]:
    output_directory = resolve_repository_path(arguments.output_dir)
    plan, plan_hash = load_plan(output_directory)
    if not 0 <= arguments.shard_index < arguments.shard_count:
        raise ValueError("shard index must lie within 0..shard-count-1")
    tasks = [
        task
        for task in plan["tasks"]
        if task["disposition"] == "scheduled"
        and int(task["scheduled_index"]) % arguments.shard_count
            == arguments.shard_index
    ]
    shard_directory = (
        output_directory
        / "shards"
        / (
            f"shard-{arguments.shard_index:03d}-of-"
            f"{arguments.shard_count:03d}"
        )
    )
    shard_manifest = shard_directory / "shard-manifest.json"
    if shard_directory.exists():
        if not arguments.resume:
            raise FileExistsError(
                f"shard exists; pass --resume: {shard_directory}"
            )
        existing = read_json(shard_manifest)
        if (
            existing.get("plan_raw_sha256") != plan_hash
            or int(existing.get("shard_count", -1))
                != arguments.shard_count
            or int(existing.get("shard_index", -1))
                != arguments.shard_index
        ):
            raise RuntimeError("existing shard manifest mismatch")
    else:
        shard_directory.mkdir(parents=True)
        atomic_json(
            shard_manifest,
            {
                "fresh_output_enforcement": True,
                "plan_raw_sha256": plan_hash,
                "schema_version": 1,
                "selected_task_ids": [task["id"] for task in tasks],
                "shard_count": arguments.shard_count,
                "shard_index": arguments.shard_index,
            },
        )

    binary = resolve_repository_path(
        Path(plan["provenance"]["engine_binary"])
    )
    completed_before = 0
    completed_now: list[dict[str, Any]] = []
    started = time.monotonic()
    attempted_now = 0
    for task in tasks:
        task_directory = shard_directory / "tasks" / str(task["id"])
        complete_path = task_directory / "complete.json"
        if complete_path.is_file():
            validate_completed_task(
                complete_path, task, plan, plan_hash
            )
            completed_before += 1
            continue
        if (
            arguments.max_tasks is not None
            and attempted_now >= arguments.max_tasks
        ):
            break
        attempt = task_attempt_directory(task_directory)
        start_matrix, support, logical = reconstruct_task(
            plan, task
        )
        start_path = attempt / "start.matrix.txt"
        support_path = attempt / "support.coords.txt"
        best_path = attempt / "best.matrix.txt"
        tie_path = attempt / "frontier-tie.matrix.txt"
        log_path = attempt / "search.jsonl"
        report_path = attempt / "report.json"
        top_directory = attempt / "top-k"
        stdout_path = attempt / "stdout.txt"
        stderr_path = attempt / "stderr.txt"
        atomic_write(start_path, canonical_matrix_bytes(start_matrix))
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
            str(tie_path),
            "--log",
            str(log_path),
            "--report",
            str(report_path),
            "--top-k",
            str(arguments.top_k),
            "--top-k-output-dir",
            str(top_directory),
        ]
        result = subprocess.run(
            command,
            cwd=REPOSITORY_ROOT,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        atomic_write(stdout_path, result.stdout.encode("utf-8"))
        atomic_write(stderr_path, result.stderr.encode("utf-8"))
        attempted_now += 1
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
            raise RuntimeError(f"engine failed for task {task['id']}")
        report = read_json(report_path)
        replay = replay_engine_result(
            plan, task, report_path, best_path
        )

        promotion: dict[str, Any] | None = None
        best_score = int(replay["best_absolute_determinant"])
        if best_score > FRONTIER:
            verification = subprocess.run(
                ["./arena", "verify", str(best_path)],
                cwd=REPOSITORY_ROOT,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            atomic_write(
                attempt / "arena-verify.txt",
                verification.stdout.encode("utf-8"),
            )
            if verification.returncode != 0:
                raise RuntimeError("promotion failed arena verification")
            verification_path = attempt / "arena-verify.txt"
            promotion = {
                **parse_arena(verification.stdout),
                "arena_verification_artifact":
                    repository_relative(verification_path),
                "arena_verification_artifact_raw_sha256":
                    sha256_bytes(verification_path.read_bytes()),
                "best_matrix_artifact_raw_sha256":
                    replay["best_matrix_raw_sha256"],
            }

        engine_mask = int(report["best_mask_decimal"])
        logical_mask = (
            (engine_mask ^ int(task["reroot_xor_mask_decimal"]))
            | (
                int(task["outer_mask_decimal"])
                << DIMENSION
            )
        )
        completed = {
            "affine_fingerprint_sha256":
                task["affine_fingerprint_sha256"],
            "assignment_visits": ASSIGNMENTS_PER_TASK,
            "best_absolute_determinant": str(best_score),
            "best_engine_mask_decimal": str(engine_mask),
            "best_logical_mask_decimal": str(logical_mask),
            "best_matrix_artifact": repository_relative(best_path),
            "best_matrix_artifact_raw_sha256":
                replay["best_matrix_raw_sha256"],
            "best_matrix_raw_sha256":
                replay["best_matrix_raw_sha256"],
            "best_signed_determinant":
                replay["best_signed_determinant"],
            "complete": True,
            "endpoint_control_confirmed":
                replay["endpoint_control_confirmed"],
            "endpoint_listed_in_returned_tie_masks":
                replay["endpoint_listed_in_returned_tie_masks"],
            "engine_frontier_assignment_occurrences":
                replay["engine_frontier_assignment_occurrences"],
            "engine_nonzero_frontier_tie_count":
                replay["engine_nonzero_frontier_tie_count"],
            "engine_report": repository_relative(report_path),
            "engine_report_raw_sha256": sha256_bytes(
                report_path.read_bytes()
            ),
            "engine_reported_tie_masks_truncated":
                replay["engine_reported_tie_masks_truncated"],
            "fresh_attempt_directory": repository_relative(attempt),
            "kind": logical["kind"],
            "logical_nonzero_frontier_assignment_occurrences":
                replay[
                    "logical_nonzero_frontier_assignment_occurrences"
                ],
            "logical_cube_id": logical["id"],
            "plan_raw_sha256": plan_hash,
            "promotion": promotion,
            "replayed_best_records":
                replay["replayed_best_records"],
            "replayed_returned_tie_masks":
                replay["replayed_returned_tie_masks"],
            "replayed_top_k_records":
                replay["replayed_top_k_records"],
            "schema_version": 1,
            "task_id": task["id"],
        }
        atomic_json(complete_path, completed)
        completed_now.append(completed)
        print(
            f"task={task['id']} kind={logical['kind']} "
            f"best={best_score} engine_nonzero_ties="
            f"{report['frontier_nonzero_ties']} "
            f"seconds={report['elapsed_seconds']}",
            flush=True,
        )
        write_shard_report(
            shard_directory / "shard-report.json",
            plan_hash,
            arguments.shard_index,
            arguments.shard_count,
            tasks,
            completed_before,
            completed_now,
            started,
        )
        if promotion is not None:
            break

    write_shard_report(
        shard_directory / "shard-report.json",
        plan_hash,
        arguments.shard_index,
        arguments.shard_count,
        tasks,
        completed_before,
        completed_now,
        started,
    )
    return read_json(shard_directory / "shard-report.json")


def summarize_campaign(arguments: argparse.Namespace) -> dict[str, Any]:
    output_directory = resolve_repository_path(arguments.output_dir)
    plan, plan_hash = load_plan(output_directory)
    scheduled = {
        str(task["id"]): task
        for task in plan["tasks"]
        if task["disposition"] == "scheduled"
    }
    completed: dict[str, dict[str, Any]] = {}
    for path in sorted(
        (output_directory / "shards").glob(
            "shard-*-of-*/tasks/*/complete.json"
        )
    ):
        task_id = path.parent.name
        if task_id not in scheduled:
            raise RuntimeError(f"unexpected completed task: {path}")
        if task_id in completed:
            raise RuntimeError(f"task completed in two shards: {task_id}")
        completed[task_id] = validate_completed_task(
            path, scheduled[task_id], plan, plan_hash
        )
    evaluated = list(completed.values())
    planned_counts = plan["counts"]
    evaluated_direct = [
        record for record in evaluated if record["kind"] == "direct"
    ]
    evaluated_bridge = [
        record for record in evaluated if record["kind"] == "bridge"
    ]
    best = max(
        (
            int(record["best_absolute_determinant"])
            for record in evaluated
        ),
        default=FRONTIER,
    )
    expected_controls = sum(
        task["endpoint_expected"] for task in scheduled.values()
    )
    confirmed_controls = sum(
        record["endpoint_control_confirmed"] is True
        for record in evaluated_bridge
    )
    report = {
        "assignment_visits_are_unique": False,
        "best_absolute_determinant": str(best),
        "complete": len(evaluated) == len(scheduled),
        "counts": {
            "candidate_evaluator_task_count":
                planned_counts["candidate_evaluator_task_count"],
            "duplicate_campaign_task_count":
                planned_counts["duplicate_campaign_task_count"],
            "duplicate_prior_task_count":
                planned_counts["duplicate_prior_task_count"],
            "duplicate_prior_tasks_skipped_not_reevaluated":
                planned_counts["duplicate_prior_task_count"],
            "evaluated_assignment_visits": sum(
                int(record["assignment_visits"])
                for record in evaluated
            ),
            "evaluated_bridge_assignment_visits": sum(
                int(record["assignment_visits"])
                for record in evaluated_bridge
            ),
            "evaluated_bridge_task_count": len(evaluated_bridge),
            "evaluated_direct_assignment_visits": sum(
                int(record["assignment_visits"])
                for record in evaluated_direct
            ),
            "evaluated_direct_task_count": len(evaluated_direct),
            "evaluated_task_count": len(evaluated),
            "logical_bridge_cube_count":
                planned_counts["bridge_logical_cube_count"],
            "logical_cube_count":
                planned_counts["logical_cube_count"],
            "logical_direct_cube_count":
                planned_counts["direct_logical_cube_count"],
            "planned_assignment_visits":
                planned_counts["planned_assignment_visits"],
            "planned_bridge_assignment_visits":
                planned_counts["planned_bridge_assignment_visits"],
            "planned_direct_assignment_visits":
                planned_counts["planned_direct_assignment_visits"],
            "scheduled_unique_evaluator_task_count": len(scheduled),
            "unrun_assignment_visits": (
                int(planned_counts["planned_assignment_visits"])
                - sum(
                    int(record["assignment_visits"])
                    for record in evaluated
                )
            ),
            "unrun_task_count": len(scheduled) - len(evaluated),
        },
        "endpoint_controls": {
            "confirmed": confirmed_controls,
            "expected": expected_controls,
        },
        "frontier": str(FRONTIER),
        "frontier_gain": str(best - FRONTIER),
        "method": plan["method"],
        "plan_raw_sha256": plan_hash,
        "promotions": [
            {
                "promotion": record["promotion"],
                "task_id": record["task_id"],
            }
            for record in evaluated
            if record["promotion"] is not None
        ],
        "replay_audit": {
            "replayed_best_records": sum(
                int(record["replayed_best_records"])
                for record in evaluated
            ),
            "replayed_returned_tie_masks": sum(
                int(record["replayed_returned_tie_masks"])
                for record in evaluated
            ),
            "replayed_top_k_records": sum(
                int(record["replayed_top_k_records"])
                for record in evaluated
            ),
        },
        "schema_version": 1,
        "tie_occurrence_accounting": {
            "engine_frontier_assignment_occurrences": sum(
                int(record["engine_frontier_assignment_occurrences"])
                for record in evaluated
            ),
            "engine_nonzero_frontier_tie_count": sum(
                int(record["engine_nonzero_frontier_tie_count"])
                for record in evaluated
            ),
            "logical_nonzero_frontier_assignment_occurrences": sum(
                int(
                    record[
                        "logical_nonzero_frontier_assignment_occurrences"
                    ]
                )
                for record in evaluated
            ),
            "reported_tie_masks_truncated_task_count": sum(
                record["engine_reported_tie_masks_truncated"] is True
                for record in evaluated
            ),
            "semantics": (
                "Occurrence counts sum per-cube assignments and are not "
                "unique matrices; logical counts apply reroot/outer-mask "
                "mapping and exclude logical mask zero."
            ),
        },
    }
    atomic_json(output_directory / "aggregate-report.json", report)
    return report


def self_test() -> None:
    matrix = tuple(
        tuple(
            1 if row == 0 or column == 0 else (
                1 if (row + 2 * column) % 3 else -1
            )
            for column in range(ORDER)
        )
        for row in range(ORDER)
    )
    ranked = [
        (row, column)
        for row in range(1, ORDER)
        for column in range(1, ORDER)
    ]
    support = fill_support(ranked)
    assert len(support) == DIMENSION
    assert all(core_coordinate(coordinate) for coordinate in support)
    assert gf2_direction_rank(support) == DIMENSION
    rerooted = apply_mask(matrix, support, 0b10101)
    assert affine_fingerprint(matrix, support) == affine_fingerprint(
        rerooted, support
    )
    outside = next(
        coordinate
        for coordinate in ranked
        if coordinate not in set(support)
    )
    assert affine_fingerprint(
        matrix, support
    ) != affine_fingerprint(
        apply_mask(matrix, (outside,), 1), support
    )
    print(
        json.dumps(
            {
                "affine_reroot_invariant": True,
                "core_boundary_excluded": True,
                "fresh_support_unique": True,
            },
            sort_keys=True,
        )
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    plan = subparsers.add_parser("plan")
    plan.add_argument("--output-dir", type=Path, required=True)
    plan.add_argument(
        "--binary", type=Path, default=DEFAULT_BINARY
    )
    plan.add_argument(
        "--qd-archive", type=Path, default=DEFAULT_QD_ARCHIVE
    )
    plan.add_argument(
        "--qd-summary", type=Path, default=DEFAULT_QD_SUMMARY
    )
    plan.add_argument(
        "--prior-root", type=Path, default=DEFAULT_PRIOR_ROOT
    )
    plan.add_argument("--center-count", type=int, default=12)
    plan.add_argument("--bridge-controls", type=int, default=1)
    plan.add_argument(
        "--bridge-min-dimension", type=int, default=28
    )
    plan.add_argument(
        "--bridge-max-dimension", type=int, default=32
    )

    run = subparsers.add_parser("run")
    run.add_argument("--output-dir", type=Path, required=True)
    run.add_argument("--shard-count", type=int, required=True)
    run.add_argument("--shard-index", type=int, required=True)
    run.add_argument("--resume", action="store_true")
    run.add_argument("--top-k", type=int, default=32)
    run.add_argument("--max-tasks", type=int)

    summarize = subparsers.add_parser("summarize")
    summarize.add_argument("--output-dir", type=Path, required=True)
    subparsers.add_parser("self-test")

    arguments = parser.parse_args()
    if arguments.command == "plan":
        if not 28 <= arguments.bridge_min_dimension:
            parser.error("--bridge-min-dimension must be at least 28")
        if not (
            arguments.bridge_min_dimension
            <= arguments.bridge_max_dimension
            <= 32
        ):
            parser.error("bridge dimension range must end within 28..32")
        if arguments.bridge_controls <= 0:
            parser.error("--bridge-controls must be positive")
    if arguments.command == "run":
        if arguments.shard_count <= 0:
            parser.error("--shard-count must be positive")
        if not 1 <= arguments.top_k <= 256:
            parser.error("--top-k must lie within 1..256")
        if (
            arguments.max_tasks is not None
            and arguments.max_tasks <= 0
        ):
            parser.error("--max-tasks must be positive")
    return arguments


def main() -> int:
    arguments = parse_arguments()
    if arguments.command == "self-test":
        self_test()
        return 0
    if arguments.command == "plan":
        plan = build_plan(arguments)
        print(json.dumps(plan["counts"], sort_keys=True))
        return 0
    if arguments.command == "run":
        report = run_shard(arguments)
        print(json.dumps(report, sort_keys=True))
        return 0
    if arguments.command == "summarize":
        report = summarize_campaign(arguments)
        print(json.dumps(report, sort_keys=True))
        return 0
    raise RuntimeError("unreachable command")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
