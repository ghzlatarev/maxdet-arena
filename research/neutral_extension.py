#!/usr/bin/env python3
"""Extend a recovered neutral cycle by two arbitrary XOR generators.

The first six generators are loaded from a `neutral_cycle.py` report.
Generator S is the difference between the cycle base and an outside
frontier endpoint. Generator C is the difference between that endpoint
and its next frontier neighbor. All 2^8 selector states are materialized,
deduplicated by their exact flip mask, scored with Bareiss determinants,
and connected through one-generator neutral edges. Every unique tie or
strict promotion is written atomically to a fresh output directory.
"""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
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


def generator_rank(
    generators: list[frozenset[tuple[int, int]]],
) -> int:
    basis: dict[int, int] = {}
    rank = 0
    for generator in generators:
        vector = 0
        for row, column in generator:
            vector |= 1 << (int(row) * 23 + int(column))
        while vector:
            pivot = vector.bit_length() - 1
            if pivot in basis:
                vector ^= basis[pivot]
            else:
                basis[pivot] = vector
                rank += 1
                break
    return rank


def connected_components(
    nodes: set[str],
    edges: set[tuple[str, str]],
) -> list[list[str]]:
    adjacency: dict[str, set[str]] = {
        node: set() for node in nodes
    }
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    remaining = set(nodes)
    components: list[list[str]] = []
    while remaining:
        root = min(remaining)
        queue = deque([root])
        component: set[str] = set()
        while queue:
            node = queue.popleft()
            if node in component:
                continue
            component.add(node)
            queue.extend(adjacency[node] - component)
        remaining -= component
        components.append(sorted(component))
    components.sort(key=lambda component: (-len(component), component))
    return components


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--cycle-report", type=Path, required=True)
    parser.add_argument("--s-endpoint", type=Path, required=True)
    parser.add_argument("--c-endpoint", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    base, base_raw = read_matrix(args.base)
    s_endpoint, s_raw = read_matrix(args.s_endpoint)
    c_endpoint, c_raw = read_matrix(args.c_endpoint)
    cycle_report_raw = args.cycle_report.read_bytes()
    cycle_report = json.loads(cycle_report_raw.decode("utf-8"))
    if (
        cycle_report.get("complete") is not True
        or cycle_report.get("method")
        != "six-generator-neutral-compound-cube-v1"
    ):
        raise ValueError("cycle report is not a complete supported report")
    expected_base_hash = (
        cycle_report.get("inputs", {})
        .get("base", {})
        .get("raw_sha256")
    )
    if expected_base_hash != sha256(base_raw):
        raise ValueError("base does not match the cycle report")

    cycle_names = ["A0", "A1", "A2", "B0", "B1", "B2"]
    cycle_generators: list[frozenset[tuple[int, int]]] = []
    for name in cycle_names:
        record = cycle_report.get("generators", {}).get(name)
        if not isinstance(record, list):
            raise ValueError(f"cycle report is missing generator {name}")
        generator = frozenset(
            (int(row) - 1, int(column) - 1)
            for row, column in record
        )
        if len(generator) != 12:
            raise ValueError(f"cycle generator {name} is not 12 flips")
        cycle_generators.append(generator)

    s_generator = difference(base, s_endpoint)
    c_generator = difference(s_endpoint, c_endpoint)
    if not s_generator or not c_generator:
        raise ValueError("S and C generators must both be nonempty")
    generator_names = cycle_names + ["S", "C"]
    generators = cycle_generators + [s_generator, c_generator]

    if args.output_dir.exists():
        raise FileExistsError(
            f"fresh output directory already exists: {args.output_dir}"
        )
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir()

    selector_masks: dict[int, frozenset[tuple[int, int]]] = {}
    mask_aliases: dict[
        frozenset[tuple[int, int]], list[int]
    ] = defaultdict(list)
    for selector in range(1 << len(generators)):
        mask = xor_masks(
            [
                generators[index]
                for index in range(len(generators))
                if selector & (1 << index)
            ]
        )
        selector_masks[selector] = mask
        mask_aliases[mask].append(selector)

    base_score = abs(determinant(base))
    mask_matrices: dict[frozenset[tuple[int, int]], Matrix] = {}
    mask_scores: dict[frozenset[tuple[int, int]], int] = {}
    mask_hashes: dict[frozenset[tuple[int, int]], str] = {}
    for mask in mask_aliases:
        matrix = apply_mask(base, mask)
        text = matrix_text(matrix)
        mask_matrices[mask] = matrix
        mask_scores[mask] = abs(determinant(matrix))
        mask_hashes[mask] = sha256(text.encode("utf-8"))

    tie_masks = sorted(
        (
            mask
            for mask, score in mask_scores.items()
            if score == base_score
        ),
        key=lambda mask: (
            min(mask_aliases[mask]),
            mask_hashes[mask],
        ),
    )
    promotion_masks = sorted(
        (
            mask
            for mask, score in mask_scores.items()
            if score > base_score
        ),
        key=lambda mask: (
            -mask_scores[mask],
            min(mask_aliases[mask]),
            mask_hashes[mask],
        ),
    )
    tie_hashes = {mask_hashes[mask] for mask in tie_masks}
    cycle_hashes = {
        str(artifact["raw_sha256"])
        for artifact in cycle_report.get("tie_artifacts", [])
    }

    selector_ties = {
        selector
        for selector, mask in selector_masks.items()
        if mask_scores[mask] == base_score
    }
    edge_records: set[tuple[str, str, str]] = set()
    graph_edges: set[tuple[str, str]] = set()
    for selector in selector_ties:
        first_hash = mask_hashes[selector_masks[selector]]
        for index, name in enumerate(generator_names):
            neighbor = selector ^ (1 << index)
            if neighbor not in selector_ties:
                continue
            second_hash = mask_hashes[selector_masks[neighbor]]
            if first_hash == second_hash:
                continue
            low, high = sorted((first_hash, second_hash))
            graph_edges.add((low, high))
            edge_records.add((low, name, high))

    components = connected_components(tie_hashes, graph_edges)
    base_hash = sha256(matrix_text(base).encode("utf-8"))
    base_component = next(
        component for component in components if base_hash in component
    )
    s_hash = sha256(matrix_text(s_endpoint).encode("utf-8"))
    c_hash = sha256(matrix_text(c_endpoint).encode("utf-8"))

    artifacts: list[dict[str, object]] = []
    for kind, masks in (("tie", tie_masks), ("promotion", promotion_masks)):
        for index, mask in enumerate(masks):
            text = matrix_text(mask_matrices[mask])
            raw_hash = mask_hashes[mask]
            filename = (
                f"{kind}-{index:03d}-{raw_hash[:12]}.matrix.txt"
            )
            atomic_write(args.output_dir / filename, text)
            artifacts.append(
                {
                    "kind": kind,
                    "path": filename,
                    "raw_sha256": raw_hash,
                    "absolute_determinant": str(mask_scores[mask]),
                    "hamming_from_base": len(mask),
                    "selector_state_hexes": [
                        f"{selector:02x}"
                        for selector in mask_aliases[mask]
                    ],
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

    states: list[dict[str, object]] = []
    for selector in range(1 << len(generators)):
        mask = selector_masks[selector]
        states.append(
            {
                "selector_state_hex": f"{selector:02x}",
                "generators": [
                    generator_names[index]
                    for index in range(len(generators))
                    if selector & (1 << index)
                ],
                "matrix_raw_sha256": mask_hashes[mask],
                "hamming_from_base": len(mask),
                "absolute_determinant": str(mask_scores[mask]),
            }
        )

    report = {
        "schema_version": 1,
        "event": "finished",
        "complete": True,
        "method": "eight-generator-neutral-compound-cube-v1",
        "base_absolute_determinant": str(base_score),
        "best_absolute_determinant": str(max(mask_scores.values())),
        "selector_states": 1 << len(generators),
        "generator_gf2_rank": generator_rank(generators),
        "unique_matrices": len(mask_aliases),
        "deduplicated_selector_states": (
            (1 << len(generators)) - len(mask_aliases)
        ),
        "exact_checks": len(mask_aliases),
        "unique_tie_matrices": len(tie_masks),
        "tie_selector_states": len(selector_ties),
        "promotion_matrices": len(promotion_masks),
        "new_ties_beyond_cycle": len(tie_hashes - cycle_hashes),
        "tie_graph_edges": len(graph_edges),
        "tie_graph_components": len(components),
        "base_component_size": len(base_component),
        "outside_pair_connected_to_base": (
            s_hash in base_component and c_hash in base_component
        ),
        "inputs": {
            "base": {
                "path": str(args.base),
                "raw_sha256": sha256(base_raw),
            },
            "cycle_report": {
                "path": str(args.cycle_report),
                "raw_sha256": hashlib.sha256(
                    cycle_report_raw
                ).hexdigest(),
            },
            "s_endpoint": {
                "path": str(args.s_endpoint),
                "raw_sha256": sha256(s_raw),
            },
            "c_endpoint": {
                "path": str(args.c_endpoint),
                "raw_sha256": sha256(c_raw),
            },
        },
        "generators": {
            name: {
                "flips": len(generator),
                "coordinates_one_based": mask_record(generator),
            }
            for name, generator in zip(generator_names, generators)
        },
        "generator_overlaps": overlap_records,
        "states": states,
        "tie_edges": [
            {
                "first_raw_sha256": first,
                "generator": generator,
                "second_raw_sha256": second,
            }
            for first, generator, second in sorted(edge_records)
        ],
        "tie_components": components,
        "artifacts": artifacts,
    }
    atomic_write(
        args.output_dir / "report.json",
        json.dumps(report, indent=2, sort_keys=True) + "\n",
    )
    print(
        json.dumps(
            {
                "base_component_size": len(base_component),
                "best_absolute_determinant": str(
                    max(mask_scores.values())
                ),
                "exact_checks": len(mask_aliases),
                "new_ties_beyond_cycle": len(
                    tie_hashes - cycle_hashes
                ),
                "promotion_matrices": len(promotion_masks),
                "selector_states": 1 << len(generators),
                "tie_graph_components": len(components),
                "unique_matrices": len(mask_aliases),
                "unique_tie_matrices": len(tie_masks),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
