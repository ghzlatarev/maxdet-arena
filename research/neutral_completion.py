#!/usr/bin/env python3
"""Complete two observed neutral arm assignments into a six-generator cube.

Each observed 12-flip neighbor has three degree-2 and six degree-1 rows,
and the same pattern on columns. Its unused arm rectangles admit two
cyclic row assignments and two cyclic column assignments. There are two
ways to couple their orientations, so two observed A/B families give
four possible six-generator completions. This program exactly screens
all four, requires a unique completion whose frontier states form a
12-cycle, archives the selected cube's ties, and optionally merges that
cycle with an earlier neutral-graph report.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

from neutral_cycle import (
    Matrix,
    atomic_write,
    apply_mask,
    degree_profile,
    determinant,
    difference,
    mask_record,
    matrix_text,
    read_matrix,
    sha256,
    xor_masks,
)
from neutral_extension import connected_components, generator_rank


Coordinate = tuple[int, int]
Mask = frozenset[Coordinate]


def completion_modes(
    name: str, first: Mask
) -> tuple[dict[str, tuple[Mask, Mask]], dict[str, object]]:
    if len(first) != 12:
        raise ValueError(f"{name} observed assignment must have 12 flips")
    row_degrees, column_degrees = degree_profile(first)
    expected_degrees = [1] * 6 + [2] * 3
    if (
        sorted(row_degrees.values()) != expected_degrees
        or sorted(column_degrees.values()) != expected_degrees
    ):
        raise ValueError(
            f"{name} must have the (2x3,1x6) row/column profile"
        )

    double_rows = sorted(
        row for row, degree in row_degrees.items() if degree == 2
    )
    single_rows = sorted(
        row for row, degree in row_degrees.items() if degree == 1
    )
    double_columns = sorted(
        column
        for column, degree in column_degrees.items()
        if degree == 2
    )
    single_columns = sorted(
        column
        for column, degree in column_degrees.items()
        if degree == 1
    )
    column_pairs = [
        tuple(
            sorted(
                column
                for row, column in first
                if row == double_row
            )
        )
        for double_row in double_rows
    ]
    row_pairs = [
        tuple(
            sorted(
                row
                for row, column in first
                if column == double_column
            )
        )
        for double_column in double_columns
    ]
    if (
        len({value for pair in column_pairs for value in pair}) != 6
        or len({value for pair in row_pairs for value in pair}) != 6
    ):
        raise ValueError(f"{name} observed arms do not form three pairs")

    candidates: dict[tuple[int, int], Mask] = {}
    for row_shift in (1, 2):
        for column_shift in (1, 2):
            edges: set[Coordinate] = set()
            for index, row in enumerate(double_rows):
                edges.update(
                    (row, column)
                    for column in column_pairs[
                        (index + row_shift) % 3
                    ]
                )
            for index, column in enumerate(double_columns):
                edges.update(
                    (row, column)
                    for row in row_pairs[
                        (index + column_shift) % 3
                    ]
                )
            candidate = frozenset(edges)
            candidate_rows, candidate_columns = degree_profile(candidate)
            if (
                len(candidate) != 12
                or candidate_rows != row_degrees
                or candidate_columns != column_degrees
                or candidate & first
            ):
                raise ValueError(
                    f"{name} cyclic completion failed validation"
                )
            candidates[(row_shift, column_shift)] = candidate

    modes = {
        "parallel": (
            candidates[(1, 1)],
            candidates[(2, 2)],
        ),
        "crossed": (
            candidates[(1, 2)],
            candidates[(2, 1)],
        ),
    }
    for mode, (second, third) in modes.items():
        if second & third:
            raise ValueError(f"{name} {mode} completions overlap")
    metadata: dict[str, object] = {
        "double_rows_one_based": [row + 1 for row in double_rows],
        "single_rows_one_based": [row + 1 for row in single_rows],
        "double_columns_one_based": [
            column + 1 for column in double_columns
        ],
        "single_columns_one_based": [
            column + 1 for column in single_columns
        ],
        "observed_column_pairs_one_based": [
            [column + 1 for column in pair]
            for pair in column_pairs
        ],
        "observed_row_pairs_one_based": [
            [row + 1 for row in pair] for pair in row_pairs
        ],
    }
    return modes, metadata


def exact_cube(
    base: Matrix,
    generators: list[Mask],
    score_cache: dict[Mask, int],
    matrix_cache: dict[Mask, Matrix],
) -> tuple[list[dict[str, object]], list[int], list[tuple[int, str, int]]]:
    states: list[dict[str, object]] = []
    for state in range(1 << len(generators)):
        mask = xor_masks(
            [
                generator
                for index, generator in enumerate(generators)
                if state & (1 << index)
            ]
        )
        if mask not in score_cache:
            matrix = apply_mask(base, mask)
            matrix_cache[mask] = matrix
            score_cache[mask] = abs(determinant(matrix))
        states.append(
            {
                "state_hex": f"{state:02x}",
                "mask": mask,
                "hamming_from_base": len(mask),
                "absolute_determinant": score_cache[mask],
            }
        )

    base_score = abs(determinant(base))
    ties = [
        state
        for state in range(1 << len(generators))
        if states[state]["absolute_determinant"] == base_score
    ]
    tie_set = set(ties)
    edges: list[tuple[int, str, int]] = []
    names = ["A0", "A1", "A2", "B0", "B1", "B2"]
    for state in ties:
        for index, name in enumerate(names):
            neighbor = state ^ (1 << index)
            if neighbor in tie_set and state < neighbor:
                edges.append((state, name, neighbor))
    return states, ties, edges


def is_single_twelve_cycle(
    ties: list[int], edges: list[tuple[int, str, int]]
) -> bool:
    if len(ties) != 12 or len(edges) != 12:
        return False
    degrees = {state: 0 for state in ties}
    graph_edges: set[tuple[str, str]] = set()
    for first, _, second in edges:
        degrees[first] += 1
        degrees[second] += 1
        graph_edges.add((f"{first:02x}", f"{second:02x}"))
    components = connected_components(
        {f"{state:02x}" for state in ties}, graph_edges
    )
    return (
        all(degree == 2 for degree in degrees.values())
        and len(components) == 1
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--a0-endpoint", type=Path, required=True)
    parser.add_argument("--b0-endpoint", type=Path, required=True)
    parser.add_argument("--prior-report", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    base, base_raw = read_matrix(args.base)
    a0_endpoint, a0_raw = read_matrix(args.a0_endpoint)
    b0_endpoint, b0_raw = read_matrix(args.b0_endpoint)
    base_score = abs(determinant(base))
    if (
        abs(determinant(a0_endpoint)) != base_score
        or abs(determinant(b0_endpoint)) != base_score
    ):
        raise ValueError("observed endpoints must tie the base exactly")

    a0 = difference(base, a0_endpoint)
    b0 = difference(base, b0_endpoint)
    a_modes, a_metadata = completion_modes("A", a0)
    b_modes, b_metadata = completion_modes("B", b0)

    if args.output_dir.exists():
        raise FileExistsError(
            f"fresh output directory already exists: {args.output_dir}"
        )
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir()

    score_cache: dict[Mask, int] = {}
    matrix_cache: dict[Mask, Matrix] = {}
    completion_results: list[dict[str, object]] = []
    complete_cubes: list[
        tuple[
            str,
            str,
            list[Mask],
            list[dict[str, object]],
            list[int],
            list[tuple[int, str, int]],
        ]
    ] = []
    all_tie_masks: set[Mask] = set()
    all_promotion_masks: set[Mask] = set()
    for a_mode in ("parallel", "crossed"):
        for b_mode in ("parallel", "crossed"):
            generators = [
                a0,
                *a_modes[a_mode],
                b0,
                *b_modes[b_mode],
            ]
            if generator_rank(generators) != 6:
                raise ValueError(
                    f"{a_mode}/{b_mode} generators are dependent"
                )
            states, ties, edges = exact_cube(
                base, generators, score_cache, matrix_cache
            )
            for state in ties:
                all_tie_masks.add(states[state]["mask"])
            for state in states:
                if int(state["absolute_determinant"]) > base_score:
                    all_promotion_masks.add(state["mask"])
            full_cycle = is_single_twelve_cycle(ties, edges)
            completion_results.append(
                {
                    "a_mode": a_mode,
                    "b_mode": b_mode,
                    "selector_states": len(states),
                    "tie_states": len(ties),
                    "tie_state_hexes": [
                        f"{state:02x}" for state in ties
                    ],
                    "tie_edges": len(edges),
                    "single_twelve_cycle": full_cycle,
                    "best_absolute_determinant": str(
                        max(
                            int(state["absolute_determinant"])
                            for state in states
                        )
                    ),
                    "promotion_states": sum(
                        int(state["absolute_determinant"]) > base_score
                        for state in states
                    ),
                }
            )
            if full_cycle:
                complete_cubes.append(
                    (
                        a_mode,
                        b_mode,
                        generators,
                        states,
                        ties,
                        edges,
                    )
                )

    if len(complete_cubes) != 1:
        raise ValueError(
            "completion orientation was not uniquely determined "
            f"by a 12-cycle: found {len(complete_cubes)}"
        )
    (
        selected_a_mode,
        selected_b_mode,
        selected_generators,
        selected_states,
        selected_ties,
        selected_edges,
    ) = complete_cubes[0]
    names = ["A0", "A1", "A2", "B0", "B1", "B2"]

    local_tie_masks = {
        selected_states[state]["mask"] for state in selected_ties
    }
    local_hashes = {
        sha256(
            matrix_text(matrix_cache[mask]).encode("utf-8")
        ): mask
        for mask in local_tie_masks
    }
    local_edges: set[tuple[str, str, str]] = set()
    for first, generator, second in selected_edges:
        first_mask = selected_states[first]["mask"]
        second_mask = selected_states[second]["mask"]
        first_hash = sha256(
            matrix_text(matrix_cache[first_mask]).encode("utf-8")
        )
        second_hash = sha256(
            matrix_text(matrix_cache[second_mask]).encode("utf-8")
        )
        low, high = sorted((first_hash, second_hash))
        local_edges.add((low, generator, high))

    prior_raw: bytes | None = None
    prior_tie_hashes: set[str] = set()
    prior_edges: set[tuple[str, str, str]] = set()
    if args.prior_report is not None:
        prior_raw = args.prior_report.read_bytes()
        prior = json.loads(prior_raw.decode("utf-8"))
        if prior.get("complete") is not True:
            raise ValueError("prior report is not complete")
        prior_tie_hashes = {
            str(artifact["raw_sha256"])
            for artifact in prior.get("artifacts", [])
            if artifact.get("kind") == "tie"
        }
        for edge in prior.get("tie_edges", []):
            first = str(edge["first_raw_sha256"])
            second = str(edge["second_raw_sha256"])
            low, high = sorted((first, second))
            prior_edges.add((low, str(edge["generator"]), high))

    combined_hashes = prior_tie_hashes | set(local_hashes)
    combined_labeled_edges = prior_edges | {
        (first, f"S32/{generator}", second)
        for first, generator, second in local_edges
    }
    combined_graph_edges = {
        (first, second)
        for first, _, second in combined_labeled_edges
    }
    combined_components = connected_components(
        combined_hashes, combined_graph_edges
    )

    artifact_masks = sorted(
        all_tie_masks | all_promotion_masks,
        key=lambda mask: (
            0 if mask in all_tie_masks else 1,
            min(
                (
                    int(state["state_hex"], 16)
                    for state in selected_states
                    if state["mask"] == mask
                ),
                default=64,
            ),
            sha256(
                matrix_text(matrix_cache[mask]).encode("utf-8")
            ),
        ),
    )
    artifacts: list[dict[str, object]] = []
    for index, mask in enumerate(artifact_masks):
        score = score_cache[mask]
        kind = "promotion" if score > base_score else "tie"
        text = matrix_text(matrix_cache[mask])
        raw_hash = sha256(text.encode("utf-8"))
        filename = f"{kind}-{index:03d}-{raw_hash[:12]}.matrix.txt"
        atomic_write(args.output_dir / filename, text)
        artifacts.append(
            {
                "kind": kind,
                "path": filename,
                "raw_sha256": raw_hash,
                "absolute_determinant": str(score),
                "hamming_from_base": len(mask),
                "selected_state_hexes": [
                    state["state_hex"]
                    for state in selected_states
                    if state["mask"] == mask
                ],
            }
        )

    report = {
        "schema_version": 1,
        "event": "finished",
        "complete": True,
        "method": "neutral-arm-completion-cube-v1",
        "base_absolute_determinant": str(base_score),
        "completion_selector_states": 4 * 64,
        "completion_unique_matrices": len(score_cache),
        "exact_checks": len(score_cache),
        "selected_a_mode": selected_a_mode,
        "selected_b_mode": selected_b_mode,
        "selected_cube_states": len(selected_states),
        "selected_cube_ties": len(selected_ties),
        "selected_cube_promotions": sum(
            int(state["absolute_determinant"]) > base_score
            for state in selected_states
        ),
        "all_completion_unique_ties": len(all_tie_masks),
        "all_completion_unique_promotions": len(
            all_promotion_masks
        ),
        "completion_results": completion_results,
        "families": {
            "A": a_metadata,
            "B": b_metadata,
        },
        "selected_generators": {
            name: mask_record(generator)
            for name, generator in zip(names, selected_generators)
        },
        "selected_states": [
            {
                key: (
                    str(value)
                    if key == "absolute_determinant"
                    else value
                )
                for key, value in state.items()
                if key != "mask"
            }
            for state in selected_states
        ],
        "selected_tie_edges": [
            {
                "first_raw_sha256": first,
                "generator": generator,
                "second_raw_sha256": second,
            }
            for first, generator, second in sorted(local_edges)
        ],
        "prior_graph": (
            None
            if args.prior_report is None
            else {
                "path": str(args.prior_report),
                "raw_sha256": hashlib.sha256(prior_raw).hexdigest(),
                "tie_matrices": len(prior_tie_hashes),
                "tie_edges": len(prior_edges),
            }
        ),
        "graph_extension": {
            "local_tie_matrices": len(local_hashes),
            "local_tie_edges": len(local_edges),
            "intersection_with_prior": len(
                set(local_hashes) & prior_tie_hashes
            ),
            "new_ties_beyond_prior": len(
                set(local_hashes) - prior_tie_hashes
            ),
            "combined_tie_matrices": len(combined_hashes),
            "combined_tie_edges": len(combined_graph_edges),
            "combined_components": len(combined_components),
            "combined_component_sizes": [
                len(component) for component in combined_components
            ],
        },
        "inputs": {
            "base": {
                "path": str(args.base),
                "raw_sha256": sha256(base_raw),
            },
            "a0_endpoint": {
                "path": str(args.a0_endpoint),
                "raw_sha256": sha256(a0_raw),
            },
            "b0_endpoint": {
                "path": str(args.b0_endpoint),
                "raw_sha256": sha256(b0_raw),
            },
        },
        "artifacts": artifacts,
    }
    atomic_write(
        args.output_dir / "report.json",
        json.dumps(report, indent=2, sort_keys=True) + "\n",
    )
    print(
        json.dumps(
            {
                "all_completion_unique_promotions": len(
                    all_promotion_masks
                ),
                "combined_tie_matrices": len(combined_hashes),
                "exact_checks": len(score_cache),
                "new_ties_beyond_prior": len(
                    set(local_hashes) - prior_tie_hashes
                ),
                "selected_cube_ties": len(selected_ties),
                "selected_modes": [
                    selected_a_mode,
                    selected_b_mode,
                ],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
