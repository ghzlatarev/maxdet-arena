#!/usr/bin/env python3
"""Deterministic large-neighborhood beam search over exact 27-entry cubes.

The C++ fast-principal-minor evaluator is the arithmetic authority for each
complete cube.  This driver:

* starts from one verified frontier representative in each H0/H1/H2 lineage;
* calibrates the exact H0<->H1 and reference<->QUBO 12-flip bridges;
* captures the best nonzero, strictly subfrontier masks by Hamming weight;
* keeps a sign-normalized, lineage-balanced, Hamming-diverse beam;
* recenters fresh 27-entry cubes with controlled 6/9/12/15 overlaps; and
* independently Bareiss-verifies every state admitted to the beam.

It deliberately does not modify the trusted verifier or the GOMEA search.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from maxdet.exact import bareiss_determinant, normalize_signs
from research.fast_cube_batch import (
    DIMENSION,
    FRONTIER,
    ORDER,
    TWO22,
    Features,
    atomic_json,
    atomic_write,
    canonical_matrix_bytes,
    compute_features,
    parse_max_rss,
    read_matrix,
    sha256_bytes,
)

Coordinate = tuple[int, int]
Matrix = tuple[tuple[int, ...], ...]

DEFAULT_H0 = Path(
    "runs/direct-search/h24-deletion-remaining-classes/elites/"
    "class-9-r10-c24.matrix.txt"
)
DEFAULT_H1 = Path(
    "runs/direct-search/multiflip/"
    "class9-depth12-tie-replay-29831/tie.matrix.txt"
)
DEFAULT_H2 = Path(
    "runs/qubo-trust-pilot-20260728-seed31003/"
    "best-proposal.matrix.txt"
)
DEFAULT_CYCLE_REPORT = Path(
    "runs/direct-search/neutral-cycle/"
    "class9-six-generator-29920/report.json"
)
DEFAULT_H2_BRIDGE = Path(
    "research/h2_reference_qubo_bridge27.coords.txt"
)
DEFAULT_REFERENCE = Path("references/orrick-et-al-2003/matrix.txt")
DEFAULT_PRIOR_MANIFEST = Path(
    "runs/direct-search/fast-principal-cube/"
    "batch200-mixed-h012-20260728/manifest.json"
)

QUOTIENT_BAND_WIDTH = 16_593_750
DEFAULT_LINEAGES = ("H0", "H1", "H2")


@dataclass(frozen=True)
class SearchState:
    """A concrete search center and its reproducible lineage metadata."""

    id: str
    label: str
    lineage: str
    generation: int
    path: Path
    matrix: Matrix
    determinant: int
    raw_sha256: str
    normalized_sha256: str
    normalized_bits: int
    parent_id: str | None
    source_cube_id: str | None
    source_mask: int
    active_flips: tuple[Coordinate, ...]
    incoming_support: tuple[Coordinate, ...]
    ancestor_normalized_hashes: frozenset[str]


@dataclass(frozen=True)
class Candidate:
    """A subfrontier state proposed by an exhaustively scanned cube."""

    lineage: str
    generation: int
    parent: SearchState
    source_cube_id: str
    source_mask: int
    source_weight: int
    matrix: Matrix
    determinant: int
    absolute_determinant: int
    raw_sha256: str
    normalized_sha256: str
    normalized_bits: int
    active_flips: tuple[Coordinate, ...]
    incoming_support: tuple[Coordinate, ...]


@dataclass(frozen=True)
class CubePlan:
    id: str
    generation: int
    parent: SearchState
    mode: str
    overlap: int | None
    support: tuple[Coordinate, ...]
    support_sha256: str
    fingerprint_sha256: str
    prior_fingerprint_duplicate: bool
    required_bridge_mask: int | None


def matrix_identity(matrix: Matrix) -> tuple[str, str, int]:
    raw = canonical_matrix_bytes(matrix)
    normalized = tuple(
        tuple(value for value in row) for row in normalize_signs(matrix)
    )
    normalized_raw = canonical_matrix_bytes(normalized)
    bits = 0
    index = 0
    for row in normalized:
        for value in row:
            if value == 1:
                bits |= 1 << index
            index += 1
    return sha256_bytes(raw), sha256_bytes(normalized_raw), bits


def popcount(value: int) -> int:
    """Return the population count on Python versions predating int.bit_count."""

    return bin(value).count("1")


def make_root_state(
    lineage: str, relative_path: Path, root: Path
) -> SearchState:
    path = root / relative_path
    matrix = read_matrix(path)
    determinant = bareiss_determinant(matrix)
    if abs(determinant) != FRONTIER:
        raise ValueError(f"{lineage} root is not frontier: {relative_path}")
    raw_hash, normalized_hash, normalized_bits = matrix_identity(matrix)
    return SearchState(
        id=f"root-{lineage.lower()}",
        label=f"{lineage}-root",
        lineage=lineage,
        generation=0,
        path=relative_path,
        matrix=matrix,
        determinant=determinant,
        raw_sha256=raw_hash,
        normalized_sha256=normalized_hash,
        normalized_bits=normalized_bits,
        parent_id=None,
        source_cube_id=None,
        source_mask=0,
        active_flips=tuple(),
        incoming_support=tuple(),
        ancestor_normalized_hashes=frozenset({normalized_hash}),
    )


def apply_mask(
    matrix: Matrix, support: Sequence[Coordinate], mask: int
) -> Matrix:
    work = [list(row) for row in matrix]
    for index, (row, column) in enumerate(support):
        if (mask >> index) & 1:
            work[row][column] *= -1
    return tuple(tuple(row) for row in work)


def parse_coordinate_file(path: Path) -> tuple[Coordinate, ...]:
    coordinates: list[Coordinate] = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        stripped = line.partition("#")[0].strip()
        if not stripped:
            continue
        tokens = stripped.split()
        if len(tokens) != 2:
            raise ValueError(f"bad coordinate line {path}:{line_number}")
        row, column = (int(token) for token in tokens)
        coordinate = (row - 1, column - 1)
        if not (
            0 <= coordinate[0] < ORDER
            and 0 <= coordinate[1] < ORDER
        ):
            raise ValueError(f"coordinate out of range {path}:{line_number}")
        coordinates.append(coordinate)
    if len(coordinates) != len(set(coordinates)):
        raise ValueError(f"duplicate coordinates in {path}")
    return tuple(coordinates)


def support_bytes(support: Sequence[Coordinate]) -> bytes:
    return "".join(
        f"{row + 1} {column + 1}\n" for row, column in support
    ).encode("ascii")


def affine_fingerprint(
    matrix: Matrix, support: Sequence[Coordinate]
) -> str:
    free = set(support)
    ternary = bytearray()
    for row in range(ORDER):
        for column in range(ORDER):
            if (row, column) in free:
                ternary.append(2)
            else:
                ternary.append(1 if matrix[row][column] == 1 else 0)
    return sha256_bytes(b"maxdet-entry-cube-v1\0" + bytes(ternary))


def load_bridges(
    root: Path,
    roots: dict[str, SearchState],
    cycle_report_path: Path,
    h2_bridge_path: Path,
    reference_path: Path,
) -> tuple[
    dict[str, tuple[Coordinate, ...]],
    tuple[Coordinate, ...],
    dict[str, Any],
]:
    cycle_report = json.loads((root / cycle_report_path).read_text())
    a0 = tuple(
        (int(row) - 1, int(column) - 1)
        for row, column in cycle_report["generators"]["A0"]
    )
    if len(a0) != 12 or len(set(a0)) != 12:
        raise ValueError("A0 bridge must contain 12 unique coordinates")

    h0_to_h1 = apply_mask(roots["H0"].matrix, a0, (1 << 12) - 1)
    h1_to_h0 = apply_mask(roots["H1"].matrix, a0, (1 << 12) - 1)
    if h0_to_h1 != roots["H1"].matrix or h1_to_h0 != roots["H0"].matrix:
        raise ValueError("A0 coordinates do not exactly connect H0 and H1")
    if (
        abs(bareiss_determinant(h0_to_h1)) != FRONTIER
        or abs(bareiss_determinant(h1_to_h0)) != FRONTIER
    ):
        raise ValueError("H0/H1 bridge endpoints are not exact frontier ties")

    h2_coordinates = parse_coordinate_file(root / h2_bridge_path)
    if len(h2_coordinates) != DIMENSION:
        raise ValueError("H2 bridge file must contain 12+15 coordinates")
    h2_bridge = h2_coordinates[:12]
    h2_transverse = h2_coordinates[12:]
    reference = read_matrix(root / reference_path)
    reference_determinant = bareiss_determinant(reference)
    if abs(reference_determinant) != FRONTIER:
        raise ValueError("Orrick reference is not an exact frontier matrix")
    h2_to_reference = apply_mask(
        roots["H2"].matrix, h2_bridge, (1 << 12) - 1
    )
    reference_to_h2 = apply_mask(
        reference, h2_bridge, (1 << 12) - 1
    )
    if h2_to_reference != reference or reference_to_h2 != roots["H2"].matrix:
        raise ValueError(
            "12 listed H2 coordinates do not exactly connect QUBO/reference"
        )
    if (
        abs(bareiss_determinant(h2_to_reference)) != FRONTIER
        or abs(bareiss_determinant(reference_to_h2)) != FRONTIER
    ):
        raise ValueError("H2/reference bridge endpoints are not exact ties")

    calibration = {
        "complete": True,
        "frontier": str(FRONTIER),
        "h0_h1": {
            "bridge_size": len(a0),
            "forward_exact_matrix_match": True,
            "forward_signed_determinant": str(
                bareiss_determinant(h0_to_h1)
            ),
            "inverse_exact_matrix_match": True,
            "inverse_signed_determinant": str(
                bareiss_determinant(h1_to_h0)
            ),
            "mask_decimal": str((1 << 12) - 1),
        },
        "h2_reference": {
            "bridge_size": len(h2_bridge),
            "forward_exact_matrix_match": True,
            "forward_signed_determinant": str(
                bareiss_determinant(h2_to_reference)
            ),
            "inverse_exact_matrix_match": True,
            "inverse_signed_determinant": str(
                bareiss_determinant(reference_to_h2)
            ),
            "mask_decimal": str((1 << 12) - 1),
            "reference_path": str(reference_path),
            "transverse_size": len(h2_transverse),
        },
    }
    return (
        {"H0": a0, "H1": a0, "H2": h2_bridge},
        h2_transverse,
        calibration,
    )


def deterministic_rng(
    seed: int,
    *parts: str | int,
) -> random.Random:
    digest = hashlib.sha256(
        "\0".join((str(seed), *(str(part) for part in parts))).encode(
            "ascii"
        )
    ).digest()
    return random.Random(int.from_bytes(digest[:8], "big"))


def jittered_coordinates(
    ranked: Sequence[Coordinate],
    rng: random.Random,
    *,
    pool_size: int = 240,
    noise: float = 0.18,
) -> list[Coordinate]:
    pool = ranked[: min(pool_size, len(ranked))]
    decorated = [
        (index + noise * len(pool) * rng.random(), coordinate)
        for index, coordinate in enumerate(pool)
    ]
    decorated.sort()
    return [coordinate for _, coordinate in decorated] + list(
        ranked[len(pool) :]
    )


def pair_endpoint_rank(
    pairs: Sequence[tuple[Coordinate, Coordinate]],
) -> list[Coordinate]:
    result: list[Coordinate] = []
    seen: set[Coordinate] = set()
    for first, second in pairs:
        for coordinate in (first, second):
            if coordinate not in seen:
                seen.add(coordinate)
                result.append(coordinate)
    return result


def fill_to_dimension(
    initial: Iterable[Coordinate],
    ranked: Sequence[Coordinate],
    *,
    excluded: Iterable[Coordinate] = (),
    row_cap: int | None = 2,
    column_cap: int | None = 2,
) -> tuple[Coordinate, ...]:
    selected = list(dict.fromkeys(initial))
    selected_set = set(selected)
    excluded_set = set(excluded)
    if selected_set & excluded_set:
        raise ValueError("initial support intersects excluded coordinates")
    if len(selected) > DIMENSION:
        raise ValueError("initial support is larger than cube dimension")

    def add_with_caps(
        cap_rows: int | None, cap_columns: int | None
    ) -> None:
        row_counts = [0] * ORDER
        column_counts = [0] * ORDER
        for row, column in selected:
            row_counts[row] += 1
            column_counts[column] += 1
        for coordinate in ranked:
            if len(selected) == DIMENSION:
                return
            if coordinate in selected_set or coordinate in excluded_set:
                continue
            row, column = coordinate
            if cap_rows is not None and row_counts[row] >= cap_rows:
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

    add_with_caps(row_cap, column_cap)
    if len(selected) < DIMENSION:
        add_with_caps(3, 3)
    if len(selected) < DIMENSION:
        add_with_caps(None, None)
    if len(selected) != DIMENSION:
        raise ValueError(f"could fill only {len(selected)} support entries")
    return tuple(selected)


def root_support(
    state: SearchState,
    bridge: tuple[Coordinate, ...],
    h2_transverse: tuple[Coordinate, ...],
    features: Features,
    variant: int,
    attempt: int,
    seed: int,
) -> tuple[Coordinate, ...]:
    rng = deterministic_rng(
        seed, "root", state.raw_sha256, variant, attempt
    )
    singles = jittered_coordinates(features.single_rank, rng)
    pair_score = jittered_coordinates(
        pair_endpoint_rank(features.pair_score_rank), rng
    )
    pair_synergy = jittered_coordinates(
        pair_endpoint_rank(features.pair_synergy_rank), rng
    )

    if variant == 0:
        if state.lineage == "H2" and attempt == 0:
            support = bridge + h2_transverse
        else:
            # Preserve bridge coordinates as mask bits 0..11.  H0 and H1 use
            # distinct deterministic transverse orderings, so their affine
            # cubes are not redundant even though they share A0.
            support = fill_to_dimension(
                bridge,
                singles,
                excluded=(),
                row_cap=None,
                column_cap=None,
            )
        if support[:12] != bridge:
            raise ValueError("required root bridge lost mask bits 0..11")
        return support
    if variant == 1:
        return fill_to_dimension((), singles)
    if variant == 2:
        return fill_to_dimension((), pair_score)
    if variant == 3:
        return fill_to_dimension((), pair_synergy)

    overlap = 6 if variant == 4 else 12
    bridge_rank = sorted(
        bridge,
        key=lambda coordinate: (
            -features.single_scores[coordinate],
            coordinate,
        ),
    )
    anchors = bridge_rank[:overlap]
    excluded = set(bridge) - set(anchors)
    ranked = pair_synergy if variant == 5 else pair_score
    return fill_to_dimension(
        anchors,
        ranked + singles,
        excluded=excluded,
        row_cap=None,
        column_cap=None,
    )


def child_support(
    state: SearchState,
    features: Features,
    generation: int,
    variant: int,
    attempt: int,
    seed: int,
) -> tuple[tuple[Coordinate, ...], int]:
    if not state.incoming_support or not state.active_flips:
        raise ValueError("nonroot center lacks incoming cube metadata")
    overlaps_by_generation = {
        1: (6, 9, 12),
        2: (9, 12, 15),
    }
    overlaps = overlaps_by_generation.get(generation, (6, 12, 15))
    overlap = overlaps[variant % len(overlaps)]
    rng = deterministic_rng(
        seed,
        "child",
        state.raw_sha256,
        generation,
        variant,
        attempt,
    )

    active_rank = sorted(
        state.active_flips,
        key=lambda coordinate: (
            -features.single_scores[coordinate],
            coordinate,
        ),
    )
    # Omitting an actually active transition coordinate makes the immediate
    # parent unreachable from this child cube.
    forbidden = active_rank[
        (variant + attempt) % len(active_rank)
    ]
    previous_available = [
        coordinate
        for coordinate in state.incoming_support
        if coordinate != forbidden
    ]
    previous_available.sort(
        key=lambda coordinate: (
            0 if coordinate in state.active_flips else 1,
            -features.single_scores[coordinate],
            coordinate,
        )
    )
    anchors = previous_available[:overlap]
    if len(anchors) != overlap:
        raise ValueError("could not construct requested support overlap")
    excluded = set(state.incoming_support) - set(anchors)
    excluded.add(forbidden)

    if variant % 3 == 0:
        primary = list(features.single_rank)
    elif variant % 3 == 1:
        primary = pair_endpoint_rank(features.pair_synergy_rank)
    else:
        primary = pair_endpoint_rank(features.pair_score_rank)
    primary = jittered_coordinates(primary, rng)
    singles = jittered_coordinates(features.single_rank, rng)
    support = fill_to_dimension(
        anchors,
        primary + singles,
        excluded=excluded,
        row_cap=None,
        column_cap=None,
    )
    actual_overlap = len(set(support) & set(state.incoming_support))
    if actual_overlap != overlap:
        raise ValueError(
            f"support overlap {actual_overlap}, expected {overlap}"
        )
    if DIMENSION - actual_overlap < 12:
        raise ValueError("child support has fewer than 12 fresh entries")
    if forbidden in support:
        raise ValueError("child cube can reach its immediate parent")
    return support, overlap


def load_prior_fingerprints(root: Path, manifest_path: Path) -> set[str]:
    path = root / manifest_path
    if not path.exists():
        return set()
    manifest = json.loads(path.read_text())
    return {
        str(record["fingerprint_sha256"])
        for record in manifest.get("runs", [])
        if record.get("fingerprint_sha256")
    }


def plan_generation(
    centers: Sequence[SearchState],
    bridges: dict[str, tuple[Coordinate, ...]],
    h2_transverse: tuple[Coordinate, ...],
    generation: int,
    root_cube_count: int,
    child_cube_count: int,
    seed: int,
    seen_fingerprints: set[str],
    prior_fingerprints: set[str],
) -> list[CubePlan]:
    plans: list[CubePlan] = []
    for center in sorted(centers, key=lambda state: state.id):
        print(
            f"feature-build generation={generation} center={center.id}",
            flush=True,
        )
        features = compute_features(center)  # type: ignore[arg-type]
        count = root_cube_count if generation == 0 else child_cube_count
        for variant in range(count):
            required_bridge = generation == 0 and variant == 0
            for attempt in range(200):
                if generation == 0:
                    support = root_support(
                        center,
                        bridges[center.lineage],
                        h2_transverse,
                        features,
                        variant,
                        attempt,
                        seed,
                    )
                    overlap = (
                        12 if variant in (0, 5) else (
                            6 if variant == 4 else None
                        )
                    )
                else:
                    support, overlap = child_support(
                        center,
                        features,
                        generation,
                        variant,
                        attempt,
                        seed,
                    )
                if len(support) != DIMENSION or len(set(support)) != DIMENSION:
                    raise ValueError("planned support is not 27 unique entries")
                fingerprint = affine_fingerprint(center.matrix, support)
                # Required bridge cubes may deliberately repeat an earlier
                # calibration.  All exploratory cubes are strictly deduped.
                if fingerprint not in seen_fingerprints or required_bridge:
                    if (
                        required_bridge
                        and fingerprint in seen_fingerprints
                        and attempt == 0
                    ):
                        break
                    if fingerprint not in seen_fingerprints:
                        break
                if required_bridge:
                    # A distinct transverse set keeps bridge mask 4095 exact.
                    continue
            else:
                raise RuntimeError(
                    f"could not find unique support for {center.id}/{variant}"
                )
            prior_duplicate = fingerprint in prior_fingerprints
            seen_fingerprints.add(fingerprint)
            encoded = support_bytes(support)
            mode = (
                (
                    "bridge_calibration",
                    "exact_single_gain",
                    "exact_pair_score",
                    "exact_pair_synergy",
                    "bridge_overlap_6",
                    "destroy_repair_overlap_12",
                )[variant]
                if generation == 0
                else (
                    "overlap_single_gain",
                    "overlap_pair_synergy",
                    "overlap_pair_score",
                )[variant % 3]
            )
            plans.append(
                CubePlan(
                    id=(
                        f"g{generation}-"
                        f"{center.lineage.lower()}-{center.id}-c{variant}"
                    ),
                    generation=generation,
                    parent=center,
                    mode=mode,
                    overlap=overlap,
                    support=support,
                    support_sha256=sha256_bytes(encoded),
                    fingerprint_sha256=fingerprint,
                    prior_fingerprint_duplicate=prior_duplicate,
                    required_bridge_mask=(
                        (1 << 12) - 1 if required_bridge else None
                    ),
                )
            )
    return plans


def state_distance(first: Candidate, second: Candidate) -> int:
    return popcount(first.normalized_bits ^ second.normalized_bits)


def quotient_band(score: int) -> int:
    deficit = FRONTIER - score
    if deficit <= 0:
        return 0
    width = TWO22 * QUOTIENT_BAND_WIDTH
    return (deficit + width - 1) // width


def candidate_key(candidate: Candidate) -> tuple[int, int, str, str, int]:
    return (
        quotient_band(candidate.absolute_determinant),
        -candidate.absolute_determinant,
        candidate.normalized_sha256,
        candidate.source_cube_id,
        candidate.source_mask,
    )


def select_diverse_beam(
    candidates: Sequence[Candidate],
    lineage: str,
    count: int,
    globally_seen_hashes: set[str],
) -> tuple[list[Candidate], dict[str, int | float | str]]:
    lineage_candidates = [
        candidate
        for candidate in candidates
        if candidate.lineage == lineage
        and candidate.normalized_sha256 not in globally_seen_hashes
        and candidate.normalized_sha256
        not in candidate.parent.ancestor_normalized_hashes
    ]
    best_by_hash: dict[str, Candidate] = {}
    for candidate in lineage_candidates:
        incumbent = best_by_hash.get(candidate.normalized_sha256)
        if incumbent is None or candidate_key(candidate) < candidate_key(
            incumbent
        ):
            best_by_hash[candidate.normalized_sha256] = candidate
    unique = sorted(best_by_hash.values(), key=candidate_key)

    floor_fraction = 0.90
    eligible = [
        candidate
        for candidate in unique
        if candidate.absolute_determinant >= int(FRONTIER * floor_fraction)
    ]
    if len(eligible) < count:
        floor_fraction = 0.85
        eligible = [
            candidate
            for candidate in unique
            if candidate.absolute_determinant
            >= int(FRONTIER * floor_fraction)
        ]
    if len(eligible) < count:
        floor_fraction = 0.0
        eligible = unique
    if len(eligible) < count:
        raise RuntimeError(
            f"only {len(eligible)} unique {lineage} candidates for beam {count}"
        )

    selected: list[Candidate] = []
    used_cubes: set[str] = set()
    while len(selected) < count:
        distinct_cube = [
            candidate
            for candidate in eligible
            if candidate not in selected
            and candidate.source_cube_id not in used_cubes
        ]
        available = distinct_cube or [
            candidate
            for candidate in eligible
            if candidate not in selected
        ]
        if not selected:
            choice = min(available, key=candidate_key)
        else:
            choice = min(
                available,
                key=lambda candidate: (
                    -min(
                        state_distance(candidate, chosen)
                        for chosen in selected
                    ),
                    quotient_band(candidate.absolute_determinant),
                    -candidate.absolute_determinant,
                    candidate.normalized_sha256,
                    candidate.source_cube_id,
                    candidate.source_mask,
                ),
            )
        selected.append(choice)
        used_cubes.add(choice.source_cube_id)

    minimum_distance = (
        min(
            state_distance(first, second)
            for index, first in enumerate(selected)
            for second in selected[index + 1 :]
        )
        if len(selected) > 1
        else 0
    )
    return selected, {
        "eligible_candidates": len(eligible),
        "floor_fraction": floor_fraction,
        "minimum_pairwise_sign_normalized_hamming_distance":
            minimum_distance,
        "raw_candidates": len(lineage_candidates),
        "unique_sign_normalized_candidates": len(unique),
    }


def candidate_from_record(
    plan: CubePlan,
    weight: int,
    record: dict[str, Any],
) -> Candidate:
    mask = int(record["mask_decimal"])
    if mask == 0 or popcount(mask) != weight:
        raise ValueError("engine top-per-weight mask metadata mismatch")
    matrix = apply_mask(plan.parent.matrix, plan.support, mask)
    raw_hash, normalized_hash, normalized_bits = matrix_identity(matrix)
    active = tuple(
        coordinate
        for index, coordinate in enumerate(plan.support)
        if (mask >> index) & 1
    )
    return Candidate(
        lineage=plan.parent.lineage,
        generation=plan.generation + 1,
        parent=plan.parent,
        source_cube_id=plan.id,
        source_mask=mask,
        source_weight=weight,
        matrix=matrix,
        determinant=int(record["signed_determinant"]),
        absolute_determinant=int(record["absolute_determinant"]),
        raw_sha256=raw_hash,
        normalized_sha256=normalized_hash,
        normalized_bits=normalized_bits,
        active_flips=active,
        incoming_support=plan.support,
    )


def verify_and_materialize_beam(
    selected: Sequence[Candidate],
    output_directory: Path,
    generation: int,
    globally_seen_hashes: set[str],
) -> list[SearchState]:
    beam_directory = (
        output_directory / f"generation-{generation:02d}" / "beam"
    )
    beam: list[SearchState] = []
    for index, candidate in enumerate(
        sorted(
            selected,
            key=lambda item: (
                item.lineage,
                candidate_key(item),
            ),
        )
    ):
        exact = bareiss_determinant(candidate.matrix)
        if (
            exact != candidate.determinant
            or abs(exact) != candidate.absolute_determinant
        ):
            raise RuntimeError(
                "selected beam candidate failed exact Bareiss verification"
            )
        if candidate.normalized_sha256 in globally_seen_hashes:
            raise RuntimeError("beam materialization encountered duplicate state")
        state_id = (
            f"g{generation + 1}-{candidate.lineage.lower()}-"
            f"{index:02d}-{candidate.normalized_sha256[:12]}"
        )
        relative_path = (
            Path(f"generation-{generation:02d}")
            / "beam"
            / f"{state_id}.matrix.txt"
        )
        atomic_write(
            output_directory / relative_path,
            canonical_matrix_bytes(candidate.matrix),
        )
        ancestors = set(candidate.parent.ancestor_normalized_hashes)
        ancestors.add(candidate.parent.normalized_sha256)
        ancestors.add(candidate.normalized_sha256)
        state = SearchState(
            id=state_id,
            label=state_id,
            lineage=candidate.lineage,
            generation=generation + 1,
            path=relative_path,
            matrix=candidate.matrix,
            determinant=exact,
            raw_sha256=candidate.raw_sha256,
            normalized_sha256=candidate.normalized_sha256,
            normalized_bits=candidate.normalized_bits,
            parent_id=candidate.parent.id,
            source_cube_id=candidate.source_cube_id,
            source_mask=candidate.source_mask,
            active_flips=candidate.active_flips,
            incoming_support=candidate.incoming_support,
            ancestor_normalized_hashes=frozenset(ancestors),
        )
        beam.append(state)
        globally_seen_hashes.add(candidate.normalized_sha256)
    return beam


def state_record(state: SearchState) -> dict[str, Any]:
    return {
        "absolute_determinant": str(abs(state.determinant)),
        "active_flip_count": len(state.active_flips),
        "generation": state.generation,
        "id": state.id,
        "incoming_support_size": len(state.incoming_support),
        "lineage_h_class": state.lineage,
        "normalized_sha256": state.normalized_sha256,
        "parent_id": state.parent_id,
        "path": str(state.path),
        "raw_sha256": state.raw_sha256,
        "signed_determinant": str(state.determinant),
        "source_cube_id": state.source_cube_id,
        "source_mask_decimal": str(state.source_mask),
    }


def plan_record(plan: CubePlan) -> dict[str, Any]:
    return {
        "fingerprint_sha256": plan.fingerprint_sha256,
        "generation": plan.generation,
        "id": plan.id,
        "lineage_h_class": plan.parent.lineage,
        "mode": plan.mode,
        "overlap_with_parent_support": plan.overlap,
        "parent_id": plan.parent.id,
        "prior_fingerprint_duplicate": plan.prior_fingerprint_duplicate,
        "required_bridge_mask_decimal": (
            str(plan.required_bridge_mask)
            if plan.required_bridge_mask is not None
            else None
        ),
        "support": [
            [row + 1, column + 1] for row, column in plan.support
        ],
        "support_sha256": plan.support_sha256,
    }


def run_cube(
    plan: CubePlan,
    binary: Path,
    root: Path,
    output_directory: Path,
    top_per_weight: int,
    minimum_weight: int,
    maximum_weight: int,
    known_frontier_hashes: set[str],
) -> tuple[list[Candidate], dict[str, Any], bool]:
    generation_directory = (
        output_directory / f"generation-{plan.generation:02d}"
    )
    run_directory = generation_directory / "cubes" / plan.id
    run_directory.mkdir(parents=True, exist_ok=False)
    support_path = run_directory / "support.coords.txt"
    start_path = run_directory / "start.matrix.txt"
    output_path = run_directory / "best.matrix.txt"
    tie_path = run_directory / "frontier-tie.matrix.txt"
    report_path = run_directory / "report.json"
    log_path = run_directory / "search.jsonl"
    stdout_path = run_directory / "stdout.txt"
    timing_path = run_directory / "timing.txt"
    atomic_write(support_path, support_bytes(plan.support))
    atomic_write(start_path, canonical_matrix_bytes(plan.parent.matrix))

    command = [
        "/usr/bin/time",
        "-lp",
        str(binary),
        "--start",
        str(start_path),
        "--coordinates",
        str(support_path),
        "--output",
        str(output_path),
        "--tie-output",
        str(tie_path),
        "--log",
        str(log_path),
        "--report",
        str(report_path),
        "--top-k-per-weight",
        str(top_per_weight),
    ]
    started = time.monotonic()
    result = subprocess.run(
        command,
        cwd=root,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    wall_seconds = time.monotonic() - started
    atomic_write(stdout_path, result.stdout.encode("utf-8"))
    atomic_write(timing_path, result.stderr.encode("utf-8"))
    if result.returncode != 0:
        atomic_json(
            run_directory / "failure.json",
            {
                "command": command,
                "returncode": result.returncode,
                "stderr": result.stderr,
                "stdout": result.stdout,
            },
        )
        raise RuntimeError(f"cube engine failed: {plan.id}")
    report = json.loads(report_path.read_text())
    if report.get("complete") is not True:
        raise RuntimeError(f"incomplete cube engine report: {plan.id}")
    if (
        report["coordinate_file_raw_sha256"] != plan.support_sha256
        or report["start_parsed_matrix_sha256"]
        != plan.parent.raw_sha256
        or int(report["dimension"]) != DIMENSION
    ):
        raise RuntimeError(f"cube provenance mismatch: {plan.id}")

    strict_promotion = int(report["best_absolute_determinant"]) > FRONTIER
    arena_verification: dict[str, Any] | None = None
    if strict_promotion:
        exact = bareiss_determinant(read_matrix(output_path))
        if abs(exact) != int(report["best_absolute_determinant"]):
            raise RuntimeError("strict promotion failed exact determinant check")
        verification = subprocess.run(
            ["./arena", "verify", str(output_path)],
            cwd=root,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        atomic_write(
            run_directory / "arena-verify.txt",
            verification.stdout.encode("utf-8"),
        )
        arena_verification = {
            "returncode": verification.returncode,
            "stdout": verification.stdout,
        }
        if verification.returncode != 0:
            raise RuntimeError(
                "strict promotion failed independent arena verification"
            )

    bridge_confirmed: bool | None = None
    if plan.required_bridge_mask is not None:
        bridge_matrix = apply_mask(
            plan.parent.matrix, plan.support, plan.required_bridge_mask
        )
        bridge_determinant = bareiss_determinant(bridge_matrix)
        bridge_confirmed = abs(bridge_determinant) == FRONTIER
        if not bridge_confirmed:
            raise RuntimeError(f"required bridge failed in {plan.id}")

    new_tie: dict[str, Any] | None = None
    if tie_path.exists():
        tie_matrix = read_matrix(tie_path)
        tie_determinant = bareiss_determinant(tie_matrix)
        if abs(tie_determinant) != FRONTIER:
            raise RuntimeError("captured tie is not exact frontier")
        tie_raw, tie_normalized, _ = matrix_identity(tie_matrix)
        unseen = tie_normalized not in known_frontier_hashes
        if unseen:
            known_frontier_hashes.add(tie_normalized)
            preserved = (
                output_directory
                / "frontier-ties"
                / f"{tie_normalized}.matrix.txt"
            )
            atomic_write(preserved, canonical_matrix_bytes(tie_matrix))
        new_tie = {
            "normalized_sha256": tie_normalized,
            "preserved": unseen,
            "raw_sha256": tie_raw,
            "signed_determinant": str(tie_determinant),
        }

    candidates: list[Candidate] = []
    for group in report.get("top_k_per_weight_candidates", []):
        weight = int(group["hamming_weight"])
        if not minimum_weight <= weight <= maximum_weight:
            continue
        for record in group["candidates"]:
            candidate = candidate_from_record(plan, weight, record)
            if (
                candidate.absolute_determinant >= FRONTIER
                or candidate.source_mask == 0
            ):
                raise RuntimeError("invalid top-per-weight candidate")
            candidates.append(candidate)

    summary = {
        "arena_verification": arena_verification,
        "assignments": int(report["assignments"]),
        "best_absolute_determinant": report[
            "best_absolute_determinant"
        ],
        "best_mask_decimal": report["best_mask_decimal"],
        "best_ties": int(report["best_ties"]),
        "bridge_mask_confirmed": bridge_confirmed,
        "candidate_count": len(candidates),
        "elapsed_seconds": float(report["elapsed_seconds"]),
        "fingerprint_sha256": plan.fingerprint_sha256,
        "frontier_nonzero_ties": int(report["frontier_nonzero_ties"]),
        "id": plan.id,
        "lineage_h_class": plan.parent.lineage,
        "maximum_resident_set_size": parse_max_rss(result.stderr),
        "mode": plan.mode,
        "new_frontier_tie": new_tie,
        "overlap_with_parent_support": plan.overlap,
        "parent_id": plan.parent.id,
        "principal_minor_seconds": float(
            report["principal_minor_seconds"]
        ),
        "prior_fingerprint_duplicate": plan.prior_fingerprint_duplicate,
        "scan_seconds": float(report["scan_seconds"]),
        "strict_promotion": strict_promotion,
        "support_sha256": plan.support_sha256,
        "top_k_per_weight_requested": top_per_weight,
        "wall_seconds": round(wall_seconds, 6),
        "zero_pivot_corrections": int(
            report["zero_pivot_corrections"]
        ),
    }
    atomic_json(run_directory / "lnps-summary.json", summary)
    return candidates, summary, strict_promotion


def aggregate_report(
    arguments: argparse.Namespace,
    roots: Sequence[SearchState],
    calibration: dict[str, Any],
    generation_reports: Sequence[dict[str, Any]],
    cube_reports: Sequence[dict[str, Any]],
    *,
    complete: bool,
    reason: str,
    wall_seconds: float,
) -> dict[str, Any]:
    best = max(
        (
            int(report["best_absolute_determinant"])
            for report in cube_reports
        ),
        default=FRONTIER,
    )
    total_assignment_visits = sum(
        int(report["assignments"]) for report in cube_reports
    )
    return {
        "affine_cubes_may_overlap": True,
        "assignment_visits_are_unique": False,
        "best_absolute_determinant": str(best),
        "bridge_calibration": calibration,
        "complete": complete,
        "completed_cubes": len(cube_reports),
        "dimension": DIMENSION,
        "engine": "fast-principal-minor-entry-cube-v1",
        "frontier": str(FRONTIER),
        "frontier_gain": str(best - FRONTIER),
        "generation_reports": list(generation_reports),
        "method": "deterministic-top-per-weight-beam-lnps-v1",
        "options": {
            "beam_per_lineage": arguments.beam_per_lineage,
            "child_cubes": arguments.child_cubes,
            "generations": arguments.generations,
            "maximum_seconds": arguments.maximum_seconds,
            "maximum_weight": arguments.maximum_weight,
            "minimum_weight": arguments.minimum_weight,
            "root_cubes": arguments.root_cubes,
            "seed": arguments.seed,
            "top_per_weight": arguments.top_per_weight,
        },
        "reason": reason,
        "roots": [state_record(state) for state in roots],
        "schema_version": 1,
        "total_assignment_visits": total_assignment_visits,
        "total_candidate_records": sum(
            int(report["candidate_count"]) for report in cube_reports
        ),
        "total_engine_seconds": round(
            sum(float(report["elapsed_seconds"]) for report in cube_reports),
            6,
        ),
        "total_frontier_nonzero_ties": sum(
            int(report["frontier_nonzero_ties"])
            for report in cube_reports
        ),
        "total_zero_pivot_corrections": sum(
            int(report["zero_pivot_corrections"])
            for report in cube_reports
        ),
        "wall_seconds": round(wall_seconds, 6),
        "claim_boundary": [
            "Assignment counts are visits across overlapping affine cubes, "
            "not unique matrices.",
            "Reserved fingerprints include planned cubes that may remain "
            "unrun at a soft time cap.",
            "Sign-normalized identity is not H/HT equivalence.",
            "Only an independently arena-verified strict score increase is "
            "a promotion.",
        ],
    }


def run_search(arguments: argparse.Namespace) -> int:
    root = REPOSITORY_ROOT
    binary = (root / arguments.binary).resolve()
    if not binary.is_file():
        raise FileNotFoundError(f"missing evaluator binary: {binary}")
    output_directory = (root / arguments.output_dir).resolve()
    if output_directory.exists():
        raise FileExistsError(
            f"output directory already exists: {output_directory}"
        )
    output_directory.mkdir(parents=True)

    roots = [
        make_root_state("H0", arguments.h0, root),
        make_root_state("H1", arguments.h1, root),
        make_root_state("H2", arguments.h2, root),
    ]
    roots_by_lineage = {state.lineage: state for state in roots}
    bridges, h2_transverse, calibration = load_bridges(
        root,
        roots_by_lineage,
        arguments.cycle_report,
        arguments.h2_bridge,
        arguments.reference,
    )
    atomic_json(output_directory / "bridge-calibration.json", calibration)

    prior_fingerprints = load_prior_fingerprints(
        root, arguments.prior_manifest
    )
    seen_fingerprints = set(prior_fingerprints)
    globally_seen_hashes = {
        state.normalized_sha256 for state in roots
    }
    known_frontier_hashes = set(globally_seen_hashes)
    expanded_hashes: set[str] = set()
    centers: list[SearchState] = roots
    cube_reports: list[dict[str, Any]] = []
    generation_reports: list[dict[str, Any]] = []
    aggregate_path = output_directory / "aggregate-report.json"
    started = time.monotonic()
    reason = "complete"

    for generation in range(arguments.generations):
        plans = plan_generation(
            centers,
            bridges,
            h2_transverse,
            generation,
            arguments.root_cubes,
            arguments.child_cubes,
            arguments.seed,
            seen_fingerprints,
            prior_fingerprints,
        )
        generation_directory = (
            output_directory / f"generation-{generation:02d}"
        )
        generation_directory.mkdir(parents=True, exist_ok=False)
        atomic_json(
            generation_directory / "manifest.json",
            {
                "centers": [state_record(state) for state in centers],
                "generation": generation,
                "planned_cubes": len(plans),
                "plans": [plan_record(plan) for plan in plans],
                "schema_version": 1,
            },
        )

        generation_candidates: list[Candidate] = []
        completed_this_generation = 0
        strict_promotion = False
        for plan in plans:
            if time.monotonic() - started >= arguments.maximum_seconds:
                reason = "time_limit"
                break
            candidates, report, strict_promotion = run_cube(
                plan,
                binary,
                root,
                output_directory,
                arguments.top_per_weight,
                arguments.minimum_weight,
                arguments.maximum_weight,
                known_frontier_hashes,
            )
            generation_candidates.extend(candidates)
            cube_reports.append(report)
            completed_this_generation += 1
            print(
                f"cube={len(cube_reports)} generation={generation} "
                f"lineage={plan.parent.lineage} "
                f"best={report['best_absolute_determinant']} "
                f"candidates={len(candidates)} "
                f"seconds={report['elapsed_seconds']:.3f}",
                flush=True,
            )
            atomic_json(
                aggregate_path,
                aggregate_report(
                    arguments,
                    roots,
                    calibration,
                    generation_reports,
                    cube_reports,
                    complete=False,
                    reason="running",
                    wall_seconds=time.monotonic() - started,
                ),
            )
            if strict_promotion:
                reason = "strict_promotion"
                break
        if reason != "complete":
            generation_reports.append(
                {
                    "beam": [],
                    "candidate_records": len(generation_candidates),
                    "complete": False,
                    "completed_cubes": completed_this_generation,
                    "generation": generation,
                    "planned_cubes": len(plans),
                    "reason": reason,
                }
            )
            break

        selected_all: list[Candidate] = []
        selection_seen_hashes = set(globally_seen_hashes)
        diversity: dict[str, Any] = {}
        for lineage in DEFAULT_LINEAGES:
            selected, metrics = select_diverse_beam(
                generation_candidates,
                lineage,
                arguments.beam_per_lineage,
                selection_seen_hashes,
            )
            selected_all.extend(selected)
            selection_seen_hashes.update(
                candidate.normalized_sha256 for candidate in selected
            )
            diversity[lineage] = metrics
        next_centers = verify_and_materialize_beam(
            selected_all,
            output_directory,
            generation,
            globally_seen_hashes,
        )
        expanded_hashes.update(
            state.normalized_sha256 for state in centers
        )
        if any(
            state.normalized_sha256 in expanded_hashes
            for state in next_centers
        ):
            raise RuntimeError("beam selected an already expanded center")
        generation_report = {
            "beam": [state_record(state) for state in next_centers],
            "candidate_records": len(generation_candidates),
            "complete": True,
            "completed_cubes": completed_this_generation,
            "diversity": diversity,
            "generation": generation,
            "planned_cubes": len(plans),
            "planned_or_reserved_affine_fingerprints_total": len(
                seen_fingerprints
            ),
        }
        generation_reports.append(generation_report)
        atomic_json(
            generation_directory / "generation-report.json",
            generation_report,
        )
        centers = next_centers

    complete = (
        reason == "complete"
        and len(generation_reports) == arguments.generations
        and all(report["complete"] for report in generation_reports)
    )
    final = aggregate_report(
        arguments,
        roots,
        calibration,
        generation_reports,
        cube_reports,
        complete=complete,
        reason=reason,
        wall_seconds=time.monotonic() - started,
    )
    final["affine_fingerprints_prior"] = len(prior_fingerprints)
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
    final["frontier_sign_normalized_hashes_seen"] = len(
        known_frontier_hashes
    )
    atomic_json(aggregate_path, final)
    print(
        json.dumps(
            {
                "best_absolute_determinant": final[
                    "best_absolute_determinant"
                ],
                "complete": complete,
                "completed_cubes": final["completed_cubes"],
                "reason": reason,
                "total_assignment_visits": final[
                    "total_assignment_visits"
                ],
                "wall_seconds": final["wall_seconds"],
            },
            sort_keys=True,
        ),
        flush=True,
    )
    return 0


def self_test() -> None:
    matrix = tuple(
        tuple(1 if (row + column) % 3 else -1 for column in range(ORDER))
        for row in range(ORDER)
    )
    support = tuple(
        [(0, column) for column in range(ORDER)]
        + [(1, column) for column in range(DIMENSION - ORDER)]
    )
    reroot = apply_mask(matrix, support, 0b1010101)
    assert affine_fingerprint(matrix, support) == affine_fingerprint(
        reroot, support
    )
    outside = next(
        (row, column)
        for row in range(ORDER)
        for column in range(ORDER)
        if (row, column) not in set(support)
    )
    outside_matrix = apply_mask(matrix, (outside,), 1)
    assert affine_fingerprint(matrix, support) != affine_fingerprint(
        outside_matrix, support
    )

    sign_equivalent = [list(row) for row in matrix]
    for column in range(ORDER):
        sign_equivalent[3][column] *= -1
    for row in range(ORDER):
        sign_equivalent[row][7] *= -1
    _, normalized_first, bits_first = matrix_identity(matrix)
    _, normalized_second, bits_second = matrix_identity(
        tuple(tuple(row) for row in sign_equivalent)
    )
    assert normalized_first == normalized_second
    assert bits_first == bits_second

    filled = fill_to_dimension(
        support[:6],
        tuple(
            (row, column)
            for row in range(ORDER)
            for column in range(ORDER)
        ),
        excluded=support[6:12],
    )
    assert len(filled) == DIMENSION
    assert len(set(filled)) == DIMENSION
    assert set(filled).isdisjoint(set(support[6:12]))
    assert quotient_band(FRONTIER - TWO22 * QUOTIENT_BAND_WIDTH) == 1
    assert quotient_band(
        FRONTIER - TWO22 * (QUOTIENT_BAND_WIDTH + 1)
    ) == 2
    print(
        json.dumps(
            {
                "affine_reroot_invariant": True,
                "deterministic_support_fill": True,
                "quotient_band_boundaries": True,
                "sign_normalization_invariant": True,
            },
            sort_keys=True,
        )
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument(
        "--binary",
        type=Path,
        default=Path("build/research/fast_principal_cube_lnps"),
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--h0", type=Path, default=DEFAULT_H0)
    parser.add_argument("--h1", type=Path, default=DEFAULT_H1)
    parser.add_argument("--h2", type=Path, default=DEFAULT_H2)
    parser.add_argument(
        "--cycle-report", type=Path, default=DEFAULT_CYCLE_REPORT
    )
    parser.add_argument(
        "--h2-bridge", type=Path, default=DEFAULT_H2_BRIDGE
    )
    parser.add_argument(
        "--reference", type=Path, default=DEFAULT_REFERENCE
    )
    parser.add_argument(
        "--prior-manifest", type=Path, default=DEFAULT_PRIOR_MANIFEST
    )
    parser.add_argument("--generations", type=int, default=3)
    parser.add_argument("--root-cubes", type=int, default=6)
    parser.add_argument("--child-cubes", type=int, default=3)
    parser.add_argument("--beam-per-lineage", type=int, default=6)
    parser.add_argument("--top-per-weight", type=int, default=2)
    parser.add_argument("--minimum-weight", type=int, default=4)
    parser.add_argument("--maximum-weight", type=int, default=23)
    parser.add_argument("--maximum-seconds", type=float, default=315.0)
    parser.add_argument("--seed", type=int, default=31_280)
    arguments = parser.parse_args()
    if arguments.self_test:
        return arguments
    if arguments.output_dir is None:
        parser.error("--output-dir is required unless --self-test is used")
    if arguments.root_cubes != 6:
        parser.error("--root-cubes must be 6 for the fixed root schedule")
    for name in (
        "generations",
        "root_cubes",
        "child_cubes",
        "beam_per_lineage",
        "top_per_weight",
    ):
        if getattr(arguments, name) <= 0:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    if not 1 <= arguments.minimum_weight <= arguments.maximum_weight <= 27:
        parser.error("weight range must lie within 1..27")
    if arguments.maximum_seconds <= 0:
        parser.error("--maximum-seconds must be positive")
    return arguments


if __name__ == "__main__":
    try:
        parsed = parse_arguments()
        if parsed.self_test:
            self_test()
            raise SystemExit(0)
        raise SystemExit(run_search(parsed))
    except (OSError, ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
