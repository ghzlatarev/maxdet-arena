#!/usr/bin/env python3
"""Summarize exact pairwise geometry of the ten local frontier portals.

The analyzer never mutates its inputs. It prints deterministic JSON to stdout
by default, or atomically installs it at ``--output``.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable


PORTALS = (
    "1e4b14334f15",
    "4072ead420f0",
    "9035bdf2a85b",
    "b584c923ea12",
    "b64c33090c9a",
    "db2cddf4b8f1",
    "de7642266b69",
    "df0b940533f8",
    "eb138a",
    "ff1b5d3735bd",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write(path: Path, value: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(value)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def portal_label(path: str) -> str:
    matches = [label for label in PORTALS if label in path]
    if len(matches) != 1:
        raise ValueError(f"cannot identify one portal in {path!r}")
    return matches[0]


def support_components(
    coordinates: Iterable[Iterable[int]],
) -> list[dict[str, Any]]:
    adjacency: dict[tuple[str, int], set[tuple[str, int]]] = {}
    for pair in coordinates:
        row, column = (int(value) for value in pair)
        row_vertex = ("row", row)
        column_vertex = ("column", column)
        adjacency.setdefault(row_vertex, set()).add(column_vertex)
        adjacency.setdefault(column_vertex, set()).add(row_vertex)

    seen: set[tuple[str, int]] = set()
    components: list[dict[str, Any]] = []
    for start in sorted(adjacency):
        if start in seen:
            continue
        queue = deque([start])
        seen.add(start)
        rows: set[int] = set()
        columns: set[int] = set()
        degree_sum = 0
        while queue:
            vertex = queue.popleft()
            degree_sum += len(adjacency[vertex])
            (rows if vertex[0] == "row" else columns).add(vertex[1])
            for neighbor in adjacency[vertex]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        edge_count = degree_sum // 2
        components.append(
            {
                "columns": len(columns),
                "complete_bipartite": edge_count
                == len(rows) * len(columns),
                "edges": edge_count,
                "rows": len(rows),
            }
        )
    return sorted(
        components,
        key=lambda item: (
            item["rows"],
            item["columns"],
            item["edges"],
        ),
    )


def support_shape(components: list[dict[str, Any]]) -> str:
    signatures = Counter(
        (item["rows"], item["columns"])
        for item in components
        if item["complete_bipartite"]
    )
    if len(signatures) == len(
        {
            (item["rows"], item["columns"])
            for item in components
        }
    ) and all(item["complete_bipartite"] for item in components):
        if signatures == {(2, 2): 8}:
            return "8*K2,2"
        if signatures == {(2, 4): 2, (4, 2): 2}:
            return "2*K2,4 + 2*K4,2"
        if signatures == {(4, 4): 2}:
            return "2*K4,4"
    return " + ".join(
        (
            f"K{item['rows']},{item['columns']}"
            if item["complete_bipartite"]
            else (
                f"B({item['rows']}r,{item['columns']}c,"
                f"{item['edges']}e)"
            )
        )
        for item in components
    )


def load_pairs(repo: Path) -> tuple[dict[tuple[str, str], dict[str, Any]], list[dict[str, str]]]:
    exploitation = (
        repo
        / "runs/direct-search/frontier-portal-exploitation-20260729"
    )
    aggregate_paths = sorted(
        exploitation.glob("exact-from-*/aggregate-report.json")
    )
    aggregate_paths.append(
        repo
        / "runs/direct-search/frontier-factor-class-expansion-20260728"
        / "eb138a-bridges/gram-aut-exhaustive/aggregate-report.json"
    )

    pairs: dict[tuple[str, str], dict[str, Any]] = {}
    inputs: list[dict[str, str]] = []
    for report_path in aggregate_paths:
        report = json.loads(report_path.read_text("utf-8"))
        if (
            report.get("complete") is not True
            or int(report.get("aut_g_elements_enumerated", 0)) != 442_368
        ):
            raise ValueError(f"incomplete alignment report: {report_path}")
        inputs.append(
            {
                "path": str(report_path.relative_to(repo)),
                "sha256": sha256(report_path),
            }
        )
        source = portal_label(report["source"])
        for alignment in report["alignments"]:
            target = portal_label(alignment["target"])
            key = tuple(sorted((source, target)))
            distance = int(alignment["hamming_distance"])
            metadata_path = repo / alignment["metadata"]
            metadata = json.loads(metadata_path.read_text("utf-8"))
            coordinates = metadata["differing_coordinates_one_based"]
            if len(coordinates) != distance:
                raise ValueError(f"support length mismatch for {key}")
            components = support_components(coordinates)
            witness = {
                "components": components,
                "metadata": str(metadata_path.relative_to(repo)),
                "metadata_sha256": sha256(metadata_path),
                "shape": support_shape(components),
            }
            if key not in pairs:
                pairs[key] = {
                    "distance": distance,
                    "witnesses": [witness],
                }
            else:
                if pairs[key]["distance"] != distance:
                    raise ValueError(f"asymmetric distance for {key}")
                if all(
                    item["metadata_sha256"] != witness["metadata_sha256"]
                    for item in pairs[key]["witnesses"]
                ):
                    pairs[key]["witnesses"].append(witness)

    expected = len(PORTALS) * (len(PORTALS) - 1) // 2
    if len(pairs) != expected:
        raise ValueError(f"found {len(pairs)} pairs, expected {expected}")
    return pairs, inputs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    pairs, inputs = load_pairs(repo)

    connector_neighbors = {label: [] for label in PORTALS}
    shape_pairs: dict[str, list[list[str]]] = {}
    shape_witness_counts: Counter[str] = Counter()
    connector_shape_pairs: dict[str, list[list[str]]] = {}
    connector_shape_witness_counts: Counter[str] = Counter()
    records = []
    for (first, second), record in sorted(pairs.items()):
        if record["distance"] <= 32:
            connector_neighbors[first].append(second)
            connector_neighbors[second].append(first)
        observed_shapes = sorted(
            {witness["shape"] for witness in record["witnesses"]}
        )
        for shape in observed_shapes:
            shape_pairs.setdefault(shape, []).append([first, second])
            if record["distance"] <= 32:
                connector_shape_pairs.setdefault(shape, []).append(
                    [first, second]
                )
        for witness in record["witnesses"]:
            shape_witness_counts[witness["shape"]] += 1
            if record["distance"] <= 32:
                connector_shape_witness_counts[witness["shape"]] += 1
        records.append(
            {
                "first": first,
                "second": second,
                "distance": record["distance"],
                "observed_minimizer_shapes": observed_shapes,
                "witnesses": sorted(
                    record["witnesses"],
                    key=lambda item: item["metadata"],
                ),
            }
        )

    distance_matrix = {
        first: {
            second: (
                0
                if first == second
                else pairs[tuple(sorted((first, second)))]["distance"]
            )
            for second in PORTALS
        }
        for first in PORTALS
    }
    distance_counts = Counter(record["distance"] for record in pairs.values())

    closure_path = (
        repo
        / "runs/direct-search/frontier-portal-exploitation-20260729"
        / "connector32-de764-1e4b/aggregate-report.json"
    )
    closure = json.loads(closure_path.read_text("utf-8"))
    closure_provenance_path = closure_path.parent / "provenance.json"
    closure_summary = {
        "best_strict_interior_absolute_determinant": closure[
            "global_top_k"
        ][0]["absolute_determinant"],
        "complete": closure["complete"],
        "frontier_gain": closure["frontier_gain"],
        "frontier_masks": [
            item["global_mask_decimal"] for item in closure["global_ties"]
        ],
        "pair": ["1e4b14334f15", "de7642266b69"],
        "report": str(closure_path.relative_to(repo)),
        "report_sha256": sha256(closure_path),
        "provenance": str(closure_provenance_path.relative_to(repo)),
        "provenance_sha256": sha256(closure_provenance_path),
        "total_assignments": closure["total_assignments"],
        "wall_seconds": closure["wall_seconds"],
    }

    result = {
        "aut_g_order": 442_368,
        "claim_boundary": (
            "Distances are exact over Aut(G) row actions and all signed "
            "column permutations, not over row actions outside Aut(G). "
            "Each retained alignment independently passed the arena verifier."
        ),
        "complete": True,
        "connector_graph_max_dimension": 32,
        "connector_graph": {
            "degrees": {
                label: len(connector_neighbors[label])
                for label in PORTALS
            },
            "edge_count": sum(
                len(neighbors)
                for neighbors in connector_neighbors.values()
            )
            // 2,
            "neighbors": {
                label: sorted(connector_neighbors[label])
                for label in PORTALS
            },
            "universal_hubs": [
                label
                for label in PORTALS
                if len(connector_neighbors[label]) == len(PORTALS) - 1
            ],
        },
        "connector_support_shape_counts": {
            shape: len(members)
            for shape, members in sorted(connector_shape_pairs.items())
        },
        "connector_support_shape_pairs": {
            shape: members
            for shape, members in sorted(connector_shape_pairs.items())
        },
        "connector_support_shape_witness_counts": {
            shape: count
            for shape, count in sorted(
                connector_shape_witness_counts.items()
            )
        },
        "distance_counts": {
            str(distance): count
            for distance, count in sorted(distance_counts.items())
        },
        "distance_matrix": distance_matrix,
        "exact_closed_connector": closure_summary,
        "input_reports": inputs,
        "pair_count": len(pairs),
        "pairs": records,
        "portal_count": len(PORTALS),
        "portals": list(PORTALS),
        "schema_version": 1,
        "support_shape_counts": {
            shape: len(members)
            for shape, members in sorted(shape_pairs.items())
        },
        "support_shape_pairs": {
            shape: members for shape, members in sorted(shape_pairs.items())
        },
        "support_shape_witness_counts": {
            shape: count
            for shape, count in sorted(shape_witness_counts.items())
        },
    }
    report_bytes = (
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    if arguments.output is None:
        print(report_bytes.decode("utf-8"), end="")
    else:
        atomic_write(arguments.output.resolve(), report_bytes)
        print(
            json.dumps(
                {
                    "output": str(arguments.output),
                    "sha256": hashlib.sha256(report_bytes).hexdigest(),
                },
                sort_keys=True,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
