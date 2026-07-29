#!/usr/bin/env python3
"""Deterministic batch driver for exact 27-bit principal-minor cubes.

The C++ evaluator remains the arithmetic authority.  This driver only chooses
supports, invokes it, preserves per-cube timing/evidence, and atomically
updates an aggregate report.  Four verified frontier starts cover both sides
of the two recovered neutral circuits, and a fifth start covers the newly
classified QUBO frontier H-class.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from maxdet.exact import bareiss_determinant

PRIME = 4_294_967_291
TWO22 = 1 << 22
QUOTIENT_BOUND = 1_089_457_290
FRONTIER = 2_779_447_296_000_000
ORDER = 23
DIMENSION = 27
STRATEGIES = (
    "top_single_gain",
    "pair_interaction_linkage",
    "neutral_switch_transverse",
    "row_column_balanced_random",
    "destroyed_repaired_overlap",
)

Coordinate = tuple[int, int]


@dataclass(frozen=True)
class Start:
    label: str
    h_class_side: str
    path: Path
    matrix: tuple[tuple[int, ...], ...]
    raw_sha256: str
    parsed_sha256: str
    determinant: int
    generators: tuple[tuple[Coordinate, ...], ...]
    known_union: frozenset[Coordinate]


@dataclass(frozen=True)
class Features:
    single_rank: tuple[Coordinate, ...]
    single_scores: dict[Coordinate, int]
    pair_score_rank: tuple[tuple[Coordinate, Coordinate], ...]
    pair_synergy_rank: tuple[tuple[Coordinate, Coordinate], ...]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_matrix_bytes(matrix: tuple[tuple[int, ...], ...]) -> bytes:
    return (
        "".join(" ".join(str(value) for value in row) + "\n" for row in matrix)
    ).encode("ascii")


def read_matrix(path: Path) -> tuple[tuple[int, ...], ...]:
    raw = path.read_bytes()
    rows = tuple(
        tuple(int(token) for token in line.split())
        for line in raw.decode("ascii").splitlines()
        if line.strip()
    )
    if (
        len(rows) != ORDER
        or any(len(row) != ORDER for row in rows)
        or any(value not in (-1, 1) for row in rows for value in row)
    ):
        raise ValueError(f"invalid 23x23 sign matrix: {path}")
    return rows


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / (
        f".{path.name}.fast-cube-batch-{os.getpid()}-{time.time_ns()}.tmp"
    )
    descriptor = os.open(
        temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600
    )
    installed = False
    try:
        offset = 0
        while offset < len(data):
            offset += os.write(descriptor, data[offset:])
        os.fsync(descriptor)
        os.close(descriptor)
        descriptor = -1
        os.replace(temporary, path)
        installed = True
        directory_descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if not installed:
            temporary.unlink(missing_ok=True)


def atomic_json(path: Path, value: Any) -> None:
    atomic_write(
        path,
        (
            json.dumps(value, indent=2, sort_keys=True, separators=(",", ": "))
            + "\n"
        ).encode("utf-8"),
    )


def matrix_mod_inverse(
    matrix: tuple[tuple[int, ...], ...],
) -> tuple[int, tuple[tuple[int, ...], ...]]:
    width = 2 * ORDER
    work = [
        [
            *(value % PRIME for value in row),
            *(1 if row_index == column else 0 for column in range(ORDER)),
        ]
        for row_index, row in enumerate(matrix)
    ]
    determinant = 1
    sign = 1
    for column in range(ORDER):
        pivot_row = next(
            (
                row
                for row in range(column, ORDER)
                if work[row][column] != 0
            ),
            None,
        )
        if pivot_row is None:
            raise ValueError("frontier start is singular modulo exact prime")
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        determinant = determinant * pivot % PRIME
        pivot_inverse = pow(pivot, PRIME - 2, PRIME)
        work[column] = [value * pivot_inverse % PRIME for value in work[column]]
        for row in range(ORDER):
            if row == column or work[row][column] == 0:
                continue
            multiplier = work[row][column]
            work[row] = [
                (work[row][inner] - multiplier * work[column][inner]) % PRIME
                for inner in range(width)
            ]
    determinant = determinant * sign % PRIME
    inverse = tuple(
        tuple(work[row][ORDER:]) for row in range(ORDER)
    )
    return determinant, inverse


def recover_absolute(residue: int) -> int:
    quotient_residue = (
        residue * pow(TWO22 % PRIME, PRIME - 2, PRIME) % PRIME
    )
    quotient = (
        quotient_residue
        if quotient_residue <= PRIME // 2
        else quotient_residue - PRIME
    )
    if abs(quotient) > QUOTIENT_BOUND:
        raise ValueError("modular candidate violates Hadamard quotient bound")
    return abs(quotient) * TWO22


def compute_features(start: Start) -> Features:
    determinant_mod, inverse = matrix_mod_inverse(start.matrix)
    if determinant_mod != start.determinant % PRIME:
        raise ValueError(f"modular determinant mismatch for {start.label}")
    coordinates = tuple(
        (row, column)
        for row in range(ORDER)
        for column in range(ORDER)
    )
    deltas = {
        coordinate: (
            -2 * start.matrix[coordinate[0]][coordinate[1]]
        )
        % PRIME
        for coordinate in coordinates
    }
    diagonals = {
        coordinate: (
            1
            + deltas[coordinate]
            * inverse[coordinate[1]][coordinate[0]]
        )
        % PRIME
        for coordinate in coordinates
    }
    single_scores = {
        coordinate: recover_absolute(
            determinant_mod * diagonals[coordinate] % PRIME
        )
        for coordinate in coordinates
    }
    single_rank = tuple(
        sorted(
            coordinates,
            key=lambda coordinate: (
                -single_scores[coordinate],
                coordinate,
            ),
        )
    )

    pair_records: list[
        tuple[int, int, Coordinate, Coordinate]
    ] = []
    for first_index, first in enumerate(coordinates):
        first_delta = deltas[first]
        for second in coordinates[first_index + 1 :]:
            cross_first_second = (
                first_delta * inverse[first[1]][second[0]] % PRIME
            )
            cross_second_first = (
                deltas[second] * inverse[second[1]][first[0]] % PRIME
            )
            ratio = (
                diagonals[first] * diagonals[second]
                - cross_first_second * cross_second_first
            ) % PRIME
            score = recover_absolute(determinant_mod * ratio % PRIME)
            synergy = score - max(
                single_scores[first], single_scores[second]
            )
            pair_records.append((score, synergy, first, second))
    pair_score_rank = tuple(
        (first, second)
        for _, _, first, second in sorted(
            pair_records,
            key=lambda record: (
                -record[0],
                -record[1],
                record[2],
                record[3],
            ),
        )
    )
    pair_synergy_rank = tuple(
        (first, second)
        for _, _, first, second in sorted(
            pair_records,
            key=lambda record: (
                -record[1],
                -record[0],
                record[2],
                record[3],
            ),
        )
    )
    return Features(
        single_rank=single_rank,
        single_scores=single_scores,
        pair_score_rank=pair_score_rank,
        pair_synergy_rank=pair_synergy_rank,
    )


def coordinate_counts(
    selected: set[Coordinate],
) -> tuple[list[int], list[int]]:
    row_counts = [0] * ORDER
    column_counts = [0] * ORDER
    for row, column in selected:
        row_counts[row] += 1
        column_counts[column] += 1
    return row_counts, column_counts


def fill_ranked(
    selected: set[Coordinate],
    ranked: tuple[Coordinate, ...] | list[Coordinate],
    count: int,
    *,
    row_cap: int | None,
    column_cap: int | None,
    excluded: frozenset[Coordinate] | set[Coordinate] = frozenset(),
) -> set[Coordinate]:
    row_counts, column_counts = coordinate_counts(selected)
    for coordinate in ranked:
        if len(selected) >= count:
            break
        if coordinate in selected or coordinate in excluded:
            continue
        row, column = coordinate
        if row_cap is not None and row_counts[row] >= row_cap:
            continue
        if column_cap is not None and column_counts[column] >= column_cap:
            continue
        selected.add(coordinate)
        row_counts[row] += 1
        column_counts[column] += 1
    if len(selected) < count and (row_cap is not None or column_cap is not None):
        return fill_ranked(
            selected,
            ranked,
            count,
            row_cap=None,
            column_cap=None,
            excluded=excluded,
        )
    if len(selected) != count:
        raise ValueError(f"could select only {len(selected)} of {count} entries")
    return selected


def jittered_rank(
    ranked: tuple[Coordinate, ...],
    rng: random.Random,
    pool_size: int,
    noise: float,
) -> list[Coordinate]:
    pool = ranked[: min(pool_size, len(ranked))]
    decorated = [
        (rank + noise * pool_size * rng.random(), coordinate)
        for rank, coordinate in enumerate(pool)
    ]
    decorated.sort()
    return [coordinate for _, coordinate in decorated]


def support_top_single(
    features: Features,
    rng: random.Random,
    variant: int,
) -> set[Coordinate]:
    cap = None if variant == 0 else (2 if variant % 2 else 3)
    ranked = jittered_rank(
        features.single_rank,
        rng,
        60 + 10 * variant,
        0.08 + 0.06 * variant,
    )
    return fill_ranked(
        set(),
        ranked + list(features.single_rank),
        DIMENSION,
        row_cap=cap,
        column_cap=cap,
    )


def support_pair_linkage(
    features: Features,
    rng: random.Random,
    variant: int,
) -> set[Coordinate]:
    edges = (
        features.pair_score_rank
        if variant % 2 == 0
        else features.pair_synergy_rank
    )
    pool_size = 800 + 300 * variant
    pool = list(edges[:pool_size])
    decorated = [
        (
            rank + (0.10 + 0.04 * variant) * pool_size * rng.random(),
            edge,
        )
        for rank, edge in enumerate(pool)
    ]
    decorated.sort()
    selected: set[Coordinate] = set()
    row_counts = [0] * ORDER
    column_counts = [0] * ORDER
    cap = 2 if variant % 3 else 3
    for _, edge in decorated:
        additions = [coordinate for coordinate in edge if coordinate not in selected]
        if not additions or len(selected) + len(additions) > DIMENSION:
            continue
        if any(
            row_counts[row] >= cap or column_counts[column] >= cap
            for row, column in additions
        ):
            continue
        for coordinate in additions:
            selected.add(coordinate)
            row_counts[coordinate[0]] += 1
            column_counts[coordinate[1]] += 1
        if len(selected) == DIMENSION:
            break
    return fill_ranked(
        selected,
        features.single_rank,
        DIMENSION,
        row_cap=cap,
        column_cap=cap,
    )


def support_neutral_transverse(
    start: Start,
    features: Features,
    rng: random.Random,
    variant: int,
) -> set[Coordinate]:
    if start.generators:
        first_index = variant % len(start.generators)
        selected = set(start.generators[first_index])
    else:
        # The new QUBO H-class has no recovered neutral switch yet.  Give it
        # a start-specific 12-entry structural anchor instead of transplanting
        # a coordinate mask from an unrelated class.
        anchor_rank = jittered_rank(
            features.single_rank, rng, 180, 0.12 + 0.03 * variant
        )
        selected = fill_ranked(
            set(),
            anchor_rank + list(features.single_rank),
            12,
            row_cap=1,
            column_cap=1,
        )
        first_index = 0
    if start.generators and variant >= 6:
        second_index = (
            first_index + 1 + (variant - 6)
        ) % len(start.generators)
        selected.update(start.generators[second_index])
    ranked = jittered_rank(
        features.single_rank,
        rng,
        180,
        0.15 + 0.04 * variant,
    )
    return fill_ranked(
        selected,
        ranked + list(features.single_rank),
        DIMENSION,
        row_cap=None,
        column_cap=None,
        excluded=start.known_union,
    )


def support_balanced_random(
    features: Features,
    rng: random.Random,
    variant: int,
) -> set[Coordinate]:
    del features, variant
    # A random perfect matching covers all 23 rows and columns once.  Four
    # additional disjoint row/column endpoints make exactly 27 coordinates
    # while keeping every line degree one or two.
    columns = list(range(ORDER))
    rng.shuffle(columns)
    selected = {(row, columns[row]) for row in range(ORDER)}
    extra_rows = rng.sample(range(ORDER), DIMENSION - ORDER)
    for _ in range(100):
        extra_columns = rng.sample(range(ORDER), DIMENSION - ORDER)
        rng.shuffle(extra_columns)
        extras = set(zip(extra_rows, extra_columns))
        if len(extras) == DIMENSION - ORDER and selected.isdisjoint(extras):
            selected.update(extras)
            return selected
    raise RuntimeError("could not construct balanced random support")


def support_destroyed_repaired(
    start: Start,
    features: Features,
    rng: random.Random,
    variant: int,
) -> set[Coordinate]:
    if start.generators:
        generator = start.generators[variant % len(start.generators)]
        choice_offset = variant // len(start.generators)
    else:
        anchor_rank = jittered_rank(
            features.single_rank, rng, 220, 0.18 + 0.04 * variant
        )
        generator = tuple(
            sorted(
                fill_ranked(
                    set(),
                    anchor_rank + list(features.single_rank),
                    12,
                    row_cap=1,
                    column_cap=1,
                )
            )
        )
        choice_offset = variant
    generator_set = set(generator)
    selected = set(generator)
    edge_rank = (
        features.pair_synergy_rank
        if variant % 2
        else features.pair_score_rank
    )
    incident: dict[Coordinate, list[Coordinate]] = {
        coordinate: [] for coordinate in generator
    }
    for first, second in edge_rank:
        if first in generator_set and second not in start.known_union:
            incident[first].append(second)
        if second in generator_set and first not in start.known_union:
            incident[second].append(first)
        if all(len(values) >= 12 for values in incident.values()):
            break
    generator_order = list(generator)
    rng.shuffle(generator_order)
    for coordinate in generator_order:
        candidates = incident[coordinate]
        if not candidates:
            continue
        window = candidates[: min(8, len(candidates))]
        alternate = window[
            (choice_offset + rng.randrange(len(window))) % len(window)
        ]
        selected.add(alternate)
        if len(selected) >= 24:
            break
    linked_endpoints: list[Coordinate] = []
    for first, second in edge_rank[:5000]:
        if first in selected and second not in start.known_union:
            linked_endpoints.append(second)
        if second in selected and first not in start.known_union:
            linked_endpoints.append(first)
    return fill_ranked(
        selected,
        linked_endpoints + list(features.single_rank),
        DIMENSION,
        row_cap=None,
        column_cap=None,
        excluded=start.known_union - generator_set,
    )


def canonical_support(support: set[Coordinate]) -> tuple[Coordinate, ...]:
    if len(support) != DIMENSION:
        raise ValueError(f"support has dimension {len(support)}, expected 27")
    return tuple(sorted(support))


def support_bytes(support: tuple[Coordinate, ...]) -> bytes:
    return "".join(
        f"{row + 1} {column + 1}\n" for row, column in support
    ).encode("ascii")


def affine_cube_fingerprint(
    start: Start, support: tuple[Coordinate, ...]
) -> str:
    free = set(support)
    ternary = bytearray()
    for row in range(ORDER):
        for column in range(ORDER):
            coordinate = (row, column)
            if coordinate in free:
                ternary.append(2)
            else:
                ternary.append(
                    1 if start.matrix[row][column] == 1 else 0
                )
    return sha256_bytes(b"maxdet-entry-cube-v1\0" + bytes(ternary))


def deterministic_rng(
    seed: int,
    start: Start,
    strategy: str,
    variant: int,
    attempt: int,
) -> random.Random:
    digest = hashlib.sha256(
        (
            f"{seed}\0{start.raw_sha256}\0{strategy}\0"
            f"{variant}\0{attempt}"
        ).encode("ascii")
    ).digest()
    return random.Random(int.from_bytes(digest[:8], "big"))


def choose_support(
    start: Start,
    features: Features,
    strategy: str,
    variant: int,
    attempt: int,
    seed: int,
) -> tuple[Coordinate, ...]:
    rng = deterministic_rng(seed, start, strategy, variant, attempt)
    effective_variant = variant + attempt
    if strategy == "top_single_gain":
        support = support_top_single(features, rng, effective_variant)
    elif strategy == "pair_interaction_linkage":
        support = support_pair_linkage(features, rng, effective_variant)
    elif strategy == "neutral_switch_transverse":
        support = support_neutral_transverse(
            start, features, rng, effective_variant
        )
    elif strategy == "row_column_balanced_random":
        support = support_balanced_random(
            features, rng, effective_variant
        )
    elif strategy == "destroyed_repaired_overlap":
        support = support_destroyed_repaired(
            start, features, rng, effective_variant
        )
    else:
        raise ValueError(f"unknown strategy: {strategy}")
    return canonical_support(support)


def load_starts(
    root: Path,
    first_report: Path,
    second_report: Path,
    third_start: Path,
) -> list[Start]:
    report_specs = []
    for report_path, generator_key, family_label in (
        (first_report, "generators", "cycle1"),
        (second_report, "selected_generators", "cycle2"),
    ):
        report = json.loads((root / report_path).read_text())
        generators_raw = report[generator_key]
        generators = tuple(
            tuple((row - 1, column - 1) for row, column in coordinates)
            for _, coordinates in sorted(generators_raw.items())
        )
        known_union = frozenset(
            coordinate
            for generator in generators
            for coordinate in generator
        )
        for input_key, side in (("base", "H0"), ("a0_endpoint", "H1")):
            input_record = report["inputs"][input_key]
            report_specs.append(
                (
                    f"{family_label}-{side.lower()}",
                    side,
                    Path(input_record["path"]),
                    input_record["raw_sha256"],
                    generators,
                    known_union,
                )
            )

    starts: list[Start] = []
    for (
        label,
        side,
        relative_path,
        expected_hash,
        generators,
        known_union,
    ) in report_specs:
        path = root / relative_path
        raw = path.read_bytes()
        raw_hash = sha256_bytes(raw)
        if raw_hash != expected_hash:
            raise ValueError(f"start hash mismatch: {relative_path}")
        matrix = read_matrix(path)
        determinant = bareiss_determinant(matrix)
        if abs(determinant) != FRONTIER:
            raise ValueError(f"start is not at frontier: {relative_path}")
        starts.append(
            Start(
                label=label,
                h_class_side=side,
                path=relative_path,
                matrix=matrix,
                raw_sha256=raw_hash,
                parsed_sha256=sha256_bytes(
                    canonical_matrix_bytes(matrix)
                ),
                determinant=determinant,
                generators=generators,
                known_union=known_union,
            )
        )
    third_path = root / third_start
    third_raw = third_path.read_bytes()
    third_matrix = read_matrix(third_path)
    third_determinant = bareiss_determinant(third_matrix)
    if abs(third_determinant) != FRONTIER:
        raise ValueError(
            f"third-class start is not at frontier: {third_start}"
        )
    starts.append(
        Start(
            label="qubo-h2",
            h_class_side="H2",
            path=third_start,
            matrix=third_matrix,
            raw_sha256=sha256_bytes(third_raw),
            parsed_sha256=sha256_bytes(
                canonical_matrix_bytes(third_matrix)
            ),
            determinant=third_determinant,
            generators=tuple(),
            known_union=frozenset(),
        )
    )
    return starts


def build_manifest(
    starts: list[Start],
    variants_per_start: int,
    seed: int,
) -> dict[str, Any]:
    runs: list[dict[str, Any]] = []
    seen_fingerprints: set[str] = set()
    seen_support_hashes: set[str] = set()
    for start in starts:
        print(f"feature-build start={start.label}", flush=True)
        features = compute_features(start)
        for strategy in STRATEGIES:
            for variant in range(variants_per_start):
                for attempt in range(100):
                    support = choose_support(
                        start,
                        features,
                        strategy,
                        variant,
                        attempt,
                        seed,
                    )
                    encoded = support_bytes(support)
                    support_hash = sha256_bytes(encoded)
                    fingerprint = affine_cube_fingerprint(start, support)
                    if fingerprint not in seen_fingerprints:
                        break
                else:
                    raise RuntimeError(
                        f"could not deduplicate {start.label}/{strategy}/"
                        f"{variant}"
                    )
                seen_fingerprints.add(fingerprint)
                seen_support_hashes.add(support_hash)
                runs.append(
                    {
                        "fingerprint_sha256": fingerprint,
                        "h_class_side": start.h_class_side,
                        "id": "",
                        "start_label": start.label,
                        "start_path": str(start.path),
                        "start_raw_sha256": start.raw_sha256,
                        "strategy": strategy,
                        "support": [
                            [row + 1, column + 1]
                            for row, column in support
                        ],
                        "support_sha256": support_hash,
                        "variant": variant,
                    }
                )
    start_order = {start.label: index for index, start in enumerate(starts)}
    strategy_order = {
        strategy: index for index, strategy in enumerate(STRATEGIES)
    }
    runs.sort(
        key=lambda record: (
            record["variant"],
            strategy_order[record["strategy"]],
            start_order[record["start_label"]],
        )
    )
    for ordered_index, record in enumerate(runs):
        record["id"] = (
            f"cube-{ordered_index:04d}-{record['start_label']}-"
            f"{record['strategy'].replace('_', '-')}"
        )
    expected = len(starts) * len(STRATEGIES) * variants_per_start
    if len(runs) != expected:
        raise RuntimeError("manifest run count mismatch")
    return {
        "dimension": DIMENSION,
        "engine": "fast-principal-minor-entry-cube-v1",
        "frontier_floor": str(FRONTIER),
        "fingerprint_method":
            "sha256-maxdet-entry-cube-v1-ternary-fixed-signs-and-free-cells",
        "method": "deterministic-mixed-support-cube-batch-v1",
        "planned_runs": len(runs),
        "schema_version": 1,
        "seed": seed,
        "starts": [
            {
                "h_class_side": start.h_class_side,
                "label": start.label,
                "path": str(start.path),
                "raw_sha256": start.raw_sha256,
            }
            for start in starts
        ],
        "strategies": list(STRATEGIES),
        "support_fingerprints_unique": len(seen_fingerprints),
        "support_set_hashes_unique": len(seen_support_hashes),
        "variants_per_start": variants_per_start,
        "runs": runs,
    }


def parse_max_rss(stderr: str) -> int | None:
    match = re.search(
        r"^\s*(\d+)\s+maximum resident set size\s*$",
        stderr,
        flags=re.MULTILINE,
    )
    return int(match.group(1)) if match else None


def aggregate_report(
    manifest: dict[str, Any],
    completed: list[dict[str, Any]],
    *,
    complete: bool,
    reason: str,
    wall_seconds: float,
) -> dict[str, Any]:
    strategy_counts = {strategy: 0 for strategy in STRATEGIES}
    side_counts = {
        side: 0
        for side in sorted(
            {start["h_class_side"] for start in manifest["starts"]}
        )
    }
    best_score = 0
    total_assignments = 0
    total_engine_seconds = 0.0
    total_principal_seconds = 0.0
    total_scan_seconds = 0.0
    total_ties = 0
    total_corrections = 0
    total_promotions = 0
    maximum_rss = 0
    for record in completed:
        strategy_counts[record["strategy"]] += 1
        side_counts[record["h_class_side"]] += 1
        best_score = max(best_score, int(record["best_absolute_determinant"]))
        total_assignments += int(record["assignments"])
        total_engine_seconds += float(record["elapsed_seconds"])
        total_principal_seconds += float(record["principal_minor_seconds"])
        total_scan_seconds += float(record["scan_seconds"])
        total_ties += int(record["best_ties"])
        total_corrections += int(record["zero_pivot_corrections"])
        total_promotions += int(record["promotions"])
        if record.get("maximum_resident_set_size") is not None:
            maximum_rss = max(
                maximum_rss, int(record["maximum_resident_set_size"])
            )
    return {
        "best_absolute_determinant": str(best_score),
        "complete": complete,
        "completed_runs": len(completed),
        "dimension": DIMENSION,
        "engine": manifest["engine"],
        "frontier_floor": str(FRONTIER),
        "frontier_gain": str(best_score - FRONTIER if completed else 0),
        "h_class_side_counts": side_counts,
        "maximum_resident_set_size": maximum_rss,
        "method": manifest["method"],
        "planned_runs": manifest["planned_runs"],
        "reason": reason,
        "runs": completed,
        "schema_version": 1,
        "seed": manifest["seed"],
        "strategy_counts": strategy_counts,
        "support_fingerprints_unique": manifest[
            "support_fingerprints_unique"
        ],
        "support_set_hashes_unique": manifest[
            "support_set_hashes_unique"
        ],
        "total_assignments": total_assignments,
        "total_best_ties": total_ties,
        "total_engine_seconds": round(total_engine_seconds, 6),
        "total_principal_minor_seconds": round(
            total_principal_seconds, 6
        ),
        "total_promotions": total_promotions,
        "total_scan_seconds": round(total_scan_seconds, 6),
        "total_zero_pivot_corrections": total_corrections,
        "wall_seconds": round(wall_seconds, 6),
    }


def run_batch(arguments: argparse.Namespace) -> int:
    root = REPOSITORY_ROOT
    binary = (root / arguments.binary).resolve()
    if not binary.is_file():
        raise FileNotFoundError(f"missing evaluator binary: {binary}")
    output_directory = (root / arguments.output_dir).resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    manifest_path = output_directory / "manifest.json"
    aggregate_path = output_directory / "aggregate-report.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if (
            manifest.get("seed") != arguments.seed
            or manifest.get("variants_per_start")
            != arguments.variants_per_start
        ):
            raise ValueError("existing manifest does not match batch options")
    else:
        starts = load_starts(
            root,
            arguments.first_report,
            arguments.second_report,
            arguments.third_start,
        )
        manifest = build_manifest(
            starts, arguments.variants_per_start, arguments.seed
        )
        atomic_json(manifest_path, manifest)

    completed_records: list[dict[str, Any]] = []
    started = time.monotonic()
    reason = "complete"
    for run_record in manifest["runs"]:
        run_directory = output_directory / run_record["id"]
        run_directory.mkdir(parents=True, exist_ok=True)
        support_path = run_directory / "support.coords.txt"
        expected_support = "".join(
            f"{row} {column}\n"
            for row, column in run_record["support"]
        ).encode("ascii")
        if (
            support_path.exists()
            and support_path.read_bytes() != expected_support
        ):
            raise ValueError(f"support mismatch in {run_record['id']}")
        if not support_path.exists():
            atomic_write(support_path, expected_support)
        report_path = run_directory / "report.json"
        output_path = run_directory / "best.matrix.txt"
        log_path = run_directory / "search.jsonl"
        timing_path = run_directory / "timing.txt"
        stdout_path = run_directory / "stdout.txt"

        engine_report: dict[str, Any] | None = None
        maximum_rss: int | None = None
        if report_path.exists():
            possible = json.loads(report_path.read_text())
            if possible.get("complete") is True:
                engine_report = possible
                if timing_path.exists():
                    maximum_rss = parse_max_rss(timing_path.read_text())
        if engine_report is None:
            elapsed = time.monotonic() - started
            if elapsed >= arguments.maximum_seconds:
                reason = "time_limit"
                break
            command = [
                "/usr/bin/time",
                "-lp",
                str(binary),
                "--start",
                str(root / run_record["start_path"]),
                "--coordinates",
                str(support_path),
                "--output",
                str(output_path),
                "--log",
                str(log_path),
                "--report",
                str(report_path),
            ]
            result = subprocess.run(
                command,
                cwd=root,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            atomic_write(stdout_path, result.stdout.encode("utf-8"))
            atomic_write(timing_path, result.stderr.encode("utf-8"))
            if result.returncode != 0:
                reason = "engine_failure"
                failure = {
                    "command": command,
                    "returncode": result.returncode,
                    "run": run_record,
                    "stderr": result.stderr,
                    "stdout": result.stdout,
                }
                atomic_json(run_directory / "failure.json", failure)
                aggregate = aggregate_report(
                    manifest,
                    completed_records,
                    complete=False,
                    reason=reason,
                    wall_seconds=time.monotonic() - started,
                )
                atomic_json(aggregate_path, aggregate)
                return 1
            engine_report = json.loads(report_path.read_text())
            if engine_report.get("complete") is not True:
                raise RuntimeError(
                    f"incomplete engine report: {run_record['id']}"
                )
            maximum_rss = parse_max_rss(result.stderr)

        if (
            engine_report["coordinate_file_raw_sha256"]
            != run_record["support_sha256"]
            or engine_report["start_raw_sha256"]
            != run_record["start_raw_sha256"]
            or int(engine_report["dimension"]) != DIMENSION
        ):
            raise ValueError(f"engine provenance mismatch: {run_record['id']}")
        summary = {
            "assignments": engine_report["assignments"],
            "best_absolute_determinant": engine_report[
                "best_absolute_determinant"
            ],
            "best_mask_decimal": engine_report["best_mask_decimal"],
            "best_matrix_sha256": engine_report["best_matrix_sha256"],
            "best_ties": engine_report["best_ties"],
            "elapsed_seconds": engine_report["elapsed_seconds"],
            "fingerprint_sha256": run_record["fingerprint_sha256"],
            "h_class_side": run_record["h_class_side"],
            "id": run_record["id"],
            "maximum_resident_set_size": maximum_rss,
            "principal_minor_seconds": engine_report[
                "principal_minor_seconds"
            ],
            "promotions": engine_report["promotions"],
            "scan_seconds": engine_report["scan_seconds"],
            "start_label": run_record["start_label"],
            "strategy": run_record["strategy"],
            "support_sha256": run_record["support_sha256"],
            "zero_pivot_corrections": engine_report[
                "zero_pivot_corrections"
            ],
        }
        completed_records.append(summary)
        aggregate = aggregate_report(
            manifest,
            completed_records,
            complete=False,
            reason="running",
            wall_seconds=time.monotonic() - started,
        )
        atomic_json(aggregate_path, aggregate)
        print(
            f"completed={len(completed_records)}/{manifest['planned_runs']} "
            f"strategy={summary['strategy']} start={summary['start_label']} "
            f"best={summary['best_absolute_determinant']} "
            f"seconds={summary['elapsed_seconds']}",
            flush=True,
        )

        if int(summary["best_absolute_determinant"]) > FRONTIER:
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
            if verification.returncode != 0:
                raise RuntimeError(
                    "strict promotion failed independent arena verification"
                )
            reason = "strict_promotion"
            break

    fully_complete = (
        reason == "complete"
        and len(completed_records) == manifest["planned_runs"]
    )
    final_report = aggregate_report(
        manifest,
        completed_records,
        complete=fully_complete,
        reason=reason,
        wall_seconds=time.monotonic() - started,
    )
    atomic_json(aggregate_path, final_report)
    print(
        json.dumps(
            {
                "best_absolute_determinant": final_report[
                    "best_absolute_determinant"
                ],
                "complete": fully_complete,
                "completed_runs": len(completed_records),
                "reason": reason,
                "total_assignments": final_report["total_assignments"],
                "wall_seconds": final_report["wall_seconds"],
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
        default=Path("build/research/fast_principal_cube"),
    )
    parser.add_argument(
        "--first-report",
        type=Path,
        default=Path(
            "runs/direct-search/neutral-cycle/"
            "class9-six-generator-29920/report.json"
        ),
    )
    parser.add_argument(
        "--second-report",
        type=Path,
        default=Path(
            "runs/direct-search/neutral-cycle/"
            "sphere32-six-generator-29943/report.json"
        ),
    )
    parser.add_argument(
        "--third-start",
        type=Path,
        default=Path(
            "runs/qubo-trust-pilot-20260728-seed31003/"
            "best-proposal.matrix.txt"
        ),
    )
    parser.add_argument("--maximum-seconds", type=float, default=600.0)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=30_280)
    parser.add_argument("--variants-per-start", type=int, default=8)
    arguments = parser.parse_args()
    if arguments.maximum_seconds <= 0:
        parser.error("--maximum-seconds must be positive")
    if arguments.variants_per_start <= 0:
        parser.error("--variants-per-start must be positive")
    return arguments


if __name__ == "__main__":
    try:
        raise SystemExit(run_batch(parse_arguments()))
    except (OSError, ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
