#!/usr/bin/env python3
"""Exhaust the selector union of two six-generator neutral cycles.

The first cycle is expressed around the class-9 base. The second cycle
is expressed around a bridge endpoint S. Their twelve switch masks plus
the class-9-to-S bridge form thirteen GF(2) generators. This program
enumerates all selector states, deduplicates exact matrices, evaluates
every unique state with Bareiss elimination, archives every frontier tie
or promotion, and compares the resulting neutral graph with the two
known cycles joined by the bridge.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

from neutral_cycle import (
    Matrix,
    atomic_write,
    apply_mask,
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


def read_report(path: Path, method: str) -> tuple[dict[str, object], bytes]:
    raw = path.read_bytes()
    report = json.loads(raw.decode("utf-8"))
    if report.get("complete") is not True or report.get("method") != method:
        raise ValueError(f"{path} is not a complete {method} report")
    return report, raw


def parse_mask(record: object, name: str) -> Mask:
    if not isinstance(record, list):
        raise ValueError(f"missing generator {name}")
    result = frozenset(
        (int(row) - 1, int(column) - 1)
        for row, column in record
    )
    if len(result) != 12:
        raise ValueError(f"generator {name} must contain 12 flips")
    return result


def matrix_for_mask(base: Matrix, mask: Mask) -> tuple[Matrix, str, str]:
    matrix = apply_mask(base, mask)
    text = matrix_text(matrix)
    return matrix, text, sha256(text.encode("utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--first-cycle-report", type=Path, required=True)
    parser.add_argument("--second-cycle-report", type=Path, required=True)
    parser.add_argument("--bridge-endpoint", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    base, base_raw = read_matrix(args.base)
    bridge_endpoint, bridge_raw = read_matrix(args.bridge_endpoint)
    first_report, first_raw = read_report(
        args.first_cycle_report,
        "six-generator-neutral-compound-cube-v1",
    )
    second_report, second_raw = read_report(
        args.second_cycle_report,
        "neutral-arm-completion-cube-v1",
    )
    if (
        first_report.get("inputs", {})
        .get("base", {})
        .get("raw_sha256")
        != sha256(base_raw)
    ):
        raise ValueError("first cycle does not use the supplied base")
    if (
        second_report.get("inputs", {})
        .get("base", {})
        .get("raw_sha256")
        != sha256(bridge_raw)
    ):
        raise ValueError(
            "second cycle does not use the supplied bridge endpoint"
        )

    local_names = ["A0", "A1", "A2", "B0", "B1", "B2"]
    first_generators = [
        parse_mask(
            first_report.get("generators", {}).get(name),
            f"C0/{name}",
        )
        for name in local_names
    ]
    second_generators = [
        parse_mask(
            second_report.get("selected_generators", {}).get(name),
            f"C1/{name}",
        )
        for name in local_names
    ]
    bridge = difference(base, bridge_endpoint)
    if not bridge:
        raise ValueError("bridge generator must be nonempty")
    generator_names = (
        [f"C0/{name}" for name in local_names]
        + [f"C1/{name}" for name in local_names]
        + ["S"]
    )
    generators = first_generators + second_generators + [bridge]
    rank = generator_rank(generators)

    if args.output_dir.exists():
        raise FileExistsError(
            f"fresh output directory already exists: {args.output_dir}"
        )
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir()

    selector_count = 1 << len(generators)
    selector_masks: list[Mask] = []
    mask_aliases: dict[Mask, list[int]] = defaultdict(list)
    for selector in range(selector_count):
        mask = xor_masks(
            [
                generator
                for index, generator in enumerate(generators)
                if selector & (1 << index)
            ]
        )
        selector_masks.append(mask)
        mask_aliases[mask].append(selector)

    base_score = abs(determinant(base))
    mask_scores: dict[Mask, int] = {}
    for mask in mask_aliases:
        matrix = apply_mask(base, mask)
        mask_scores[mask] = abs(determinant(matrix))
    selector_scores = [
        mask_scores[mask] for mask in selector_masks
    ]
    tie_selectors = {
        selector
        for selector, score in enumerate(selector_scores)
        if score == base_score
    }
    promotion_selectors = {
        selector
        for selector, score in enumerate(selector_scores)
        if score > base_score
    }

    first_artifact_by_state = {
        int(str(artifact["state_hex"]), 16): str(
            artifact["raw_sha256"]
        )
        for artifact in first_report.get("tie_artifacts", [])
    }
    first_expected_states = set(first_artifact_by_state)
    second_artifact_by_state: dict[int, str] = {}
    for artifact in second_report.get("artifacts", []):
        if artifact.get("kind") != "tie":
            continue
        states = artifact.get("selected_state_hexes", [])
        if len(states) != 1:
            raise ValueError("second-cycle tie state is ambiguous")
        second_artifact_by_state[int(str(states[0]), 16)] = str(
            artifact["raw_sha256"]
        )
    second_expected_local_states = set(second_artifact_by_state)

    second_offset = 6
    bridge_bit = 1 << 12
    expected_first_states = set(first_expected_states)
    expected_second_states = {
        bridge_bit | (state << second_offset)
        for state in second_expected_local_states
    }
    expected_tie_states = expected_first_states | expected_second_states

    state_hashes: dict[int, str] = {}
    state_texts: dict[int, str] = {}
    for selector in sorted(tie_selectors | promotion_selectors):
        _, text, raw_hash = matrix_for_mask(
            base, selector_masks[selector]
        )
        state_hashes[selector] = raw_hash
        state_texts[selector] = text

    for state, expected_hash in first_artifact_by_state.items():
        if state_hashes.get(state) != expected_hash:
            raise ValueError(
                f"first-cycle state {state:02x} hash mismatch"
            )
    for local_state, expected_hash in second_artifact_by_state.items():
        state = bridge_bit | (local_state << second_offset)
        if state_hashes.get(state) != expected_hash:
            raise ValueError(
                f"second-cycle state {local_state:02x} hash mismatch"
            )

    tie_edges: set[tuple[int, str, int]] = set()
    for selector in tie_selectors:
        for index, name in enumerate(generator_names):
            neighbor = selector ^ (1 << index)
            if neighbor in tie_selectors and selector < neighbor:
                tie_edges.add((selector, name, neighbor))

    expected_edges: set[tuple[int, str, int]] = set()
    for edge in first_report.get("tie_edges", []):
        first = int(str(edge["first_state_hex"]), 16)
        second = int(str(edge["second_state_hex"]), 16)
        low, high = sorted((first, second))
        expected_edges.add(
            (low, f"C0/{edge['generator']}", high)
        )
    second_hash_to_local_state = {
        raw_hash: state
        for state, raw_hash in second_artifact_by_state.items()
    }
    for edge in second_report.get("selected_tie_edges", []):
        first_local = second_hash_to_local_state[
            str(edge["first_raw_sha256"])
        ]
        second_local = second_hash_to_local_state[
            str(edge["second_raw_sha256"])
        ]
        first = bridge_bit | (first_local << second_offset)
        second = bridge_bit | (second_local << second_offset)
        low, high = sorted((first, second))
        expected_edges.add(
            (low, f"C1/{edge['generator']}", high)
        )
    expected_edges.add((0, "S", bridge_bit))

    graph_nodes = {f"{state:04x}" for state in tie_selectors}
    graph_edges = {
        (f"{first:04x}", f"{second:04x}")
        for first, _, second in tie_edges
    }
    components = connected_components(graph_nodes, graph_edges)
    additional_states = tie_selectors - expected_tie_states
    missing_expected_states = expected_tie_states - tie_selectors
    additional_edges = tie_edges - expected_edges
    missing_expected_edges = expected_edges - tie_edges

    mixed_tie_states = {
        state
        for state in tie_selectors
        if (state & 0x03F) != 0
        and (state & 0xFC0) != 0
    }
    translated_first_ties = {
        state
        for state in tie_selectors
        if (state & 0x03F) != 0
        and (state & 0xFC0) == 0
        and (state & bridge_bit) != 0
    }
    unbridged_second_ties = {
        state
        for state in tie_selectors
        if (state & 0x03F) == 0
        and (state & 0xFC0) != 0
        and (state & bridge_bit) == 0
    }

    artifacts: list[dict[str, object]] = []
    artifact_states = sorted(
        tie_selectors | promotion_selectors,
        key=lambda state: (
            0 if state in tie_selectors else 1,
            state,
        ),
    )
    for index, state in enumerate(artifact_states):
        kind = "promotion" if state in promotion_selectors else "tie"
        raw_hash = state_hashes[state]
        filename = (
            f"{kind}-{index:03d}-{state:04x}-"
            f"{raw_hash[:12]}.matrix.txt"
        )
        atomic_write(args.output_dir / filename, state_texts[state])
        artifacts.append(
            {
                "kind": kind,
                "path": filename,
                "selector_state_hex": f"{state:04x}",
                "raw_sha256": raw_hash,
                "absolute_determinant": str(
                    selector_scores[state]
                ),
                "hamming_from_base": len(selector_masks[state]),
                "expected_prior_tie": state in expected_tie_states,
            }
        )

    overlap_records: list[dict[str, object]] = []
    for left, left_name in enumerate(generator_names):
        for right in range(left):
            overlap = generators[left] & generators[right]
            if overlap:
                overlap_records.append(
                    {
                        "first": generator_names[right],
                        "second": left_name,
                        "overlap": len(overlap),
                        "coordinates_one_based": mask_record(
                            frozenset(overlap)
                        ),
                    }
                )

    below_scores = [
        score for score in selector_scores if score < base_score
    ]
    states = [
        {
            "selector_state_hex": f"{state:04x}",
            "absolute_determinant": str(selector_scores[state]),
            "hamming_from_base": len(selector_masks[state]),
        }
        for state in range(selector_count)
    ]
    score_histogram = Counter(selector_scores)
    report = {
        "schema_version": 1,
        "event": "finished",
        "complete": True,
        "method": "two-cycle-thirteen-generator-union-v1",
        "base_absolute_determinant": str(base_score),
        "best_absolute_determinant": str(max(selector_scores)),
        "best_below_frontier": (
            None if not below_scores else str(max(below_scores))
        ),
        "generators": len(generators),
        "generator_gf2_rank": rank,
        "selector_states": selector_count,
        "unique_matrices": len(mask_aliases),
        "deduplicated_selector_states": (
            selector_count - len(mask_aliases)
        ),
        "exact_checks": len(mask_scores),
        "tie_matrices": len(tie_selectors),
        "promotion_matrices": len(promotion_selectors),
        "expected_tie_matrices": len(expected_tie_states),
        "additional_tie_matrices": len(additional_states),
        "missing_expected_tie_matrices": len(
            missing_expected_states
        ),
        "tie_edges": len(tie_edges),
        "expected_tie_edges": len(expected_edges),
        "additional_tie_edges": len(additional_edges),
        "missing_expected_tie_edges": len(
            missing_expected_edges
        ),
        "tie_components": len(components),
        "tie_component_sizes": [
            len(component) for component in components
        ],
        "mixed_cycle_ties": len(mixed_tie_states),
        "translated_first_cycle_ties": len(
            translated_first_ties
        ),
        "unbridged_second_cycle_ties": len(
            unbridged_second_ties
        ),
        "mixed_states_create_additional_links": bool(
            additional_states or additional_edges
        ),
        "inputs": {
            "base": {
                "path": str(args.base),
                "raw_sha256": sha256(base_raw),
            },
            "bridge_endpoint": {
                "path": str(args.bridge_endpoint),
                "raw_sha256": sha256(bridge_raw),
            },
            "first_cycle_report": {
                "path": str(args.first_cycle_report),
                "raw_sha256": hashlib.sha256(first_raw).hexdigest(),
            },
            "second_cycle_report": {
                "path": str(args.second_cycle_report),
                "raw_sha256": hashlib.sha256(second_raw).hexdigest(),
            },
        },
        "generator_records": {
            name: {
                "flips": len(generator),
                "coordinates_one_based": mask_record(generator),
            }
            for name, generator in zip(generator_names, generators)
        },
        "generator_overlaps": overlap_records,
        "states": states,
        "score_histogram": [
            {
                "absolute_determinant": str(score),
                "states": count,
            }
            for score, count in sorted(
                score_histogram.items(), reverse=True
            )
        ],
        "tie_edge_records": [
            {
                "first_state_hex": f"{first:04x}",
                "generator": generator,
                "second_state_hex": f"{second:04x}",
                "expected_prior_edge": (
                    first,
                    generator,
                    second,
                )
                in expected_edges,
            }
            for first, generator, second in sorted(tie_edges)
        ],
        "additional_tie_state_hexes": [
            f"{state:04x}" for state in sorted(additional_states)
        ],
        "missing_expected_tie_state_hexes": [
            f"{state:04x}"
            for state in sorted(missing_expected_states)
        ],
        "additional_tie_edge_records": [
            {
                "first_state_hex": f"{first:04x}",
                "generator": generator,
                "second_state_hex": f"{second:04x}",
            }
            for first, generator, second in sorted(additional_edges)
        ],
        "missing_expected_tie_edge_records": [
            {
                "first_state_hex": f"{first:04x}",
                "generator": generator,
                "second_state_hex": f"{second:04x}",
            }
            for first, generator, second in sorted(
                missing_expected_edges
            )
        ],
        "artifacts": artifacts,
    }
    atomic_write(
        args.output_dir / "report.json",
        json.dumps(report, indent=2, sort_keys=True) + "\n",
    )
    print(
        json.dumps(
            {
                "additional_tie_edges": len(additional_edges),
                "additional_tie_matrices": len(additional_states),
                "best_absolute_determinant": str(max(selector_scores)),
                "exact_checks": len(mask_scores),
                "generator_gf2_rank": rank,
                "promotion_matrices": len(promotion_selectors),
                "selector_states": selector_count,
                "tie_matrices": len(tie_selectors),
                "unique_matrices": len(mask_aliases),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
