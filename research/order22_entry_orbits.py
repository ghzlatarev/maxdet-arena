#!/usr/bin/env python3
"""Exact order-22 entry orbits under signed row/column automorphisms."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from collections import Counter
from pathlib import Path

try:
    import pynauty
except ImportError as error:
    raise SystemExit(
        "error: install pynauty==2.8.8.1 in an isolated environment"
    ) from error


ORDER = 22
PINNED_PYNAUTY = "2.8.8.1"
if getattr(pynauty, "__version__", None) != PINNED_PYNAUTY:
    raise SystemExit(f"error: this certificate requires pynauty=={PINNED_PYNAUTY}")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def matrix_slug(path: Path) -> str:
    suffix = ".matrix.txt"
    name = path.name
    slug = name[: -len(suffix)] if name.endswith(suffix) else path.stem
    if not slug or any(
        character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_"
        for character in slug
    ):
        raise ValueError(f"{path}: cannot derive safe deterministic export slug")
    return slug


def atomic_write_text(path: Path, contents: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            output.write(contents)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def read_matrix(path: Path) -> list[list[int]]:
    rows = [
        [int(token) for token in line.split()]
        for line in path.read_text(encoding="ascii").splitlines()
        if line.strip()
    ]
    if len(rows) != ORDER or any(
        len(row) != ORDER or any(value not in (-1, 1) for value in row)
        for row in rows
    ):
        raise ValueError(f"{path}: expected exactly {ORDER}x{ORDER} signs")
    return rows


def factor_graph(
    matrix: list[list[int]],
) -> tuple[pynauty.Graph, dict[str, object]]:
    """Build a colored signed double cover with explicit antipodal fibers."""
    row_port_base = 0
    column_port_base = 2 * ORDER
    row_fiber_base = 4 * ORDER
    column_fiber_base = 5 * ORDER
    vertex_count = 6 * ORDER
    adjacency: dict[int, list[int]] = {
        vertex: [] for vertex in range(vertex_count)
    }

    def row_port(row: int, negative: bool) -> int:
        return row_port_base + 2 * row + int(negative)

    def column_port(column: int, negative: bool) -> int:
        return column_port_base + 2 * column + int(negative)

    def add_edge(first: int, second: int) -> None:
        adjacency[first].append(second)
        adjacency[second].append(first)

    for row in range(ORDER):
        fiber = row_fiber_base + row
        add_edge(fiber, row_port(row, False))
        add_edge(fiber, row_port(row, True))
    for column in range(ORDER):
        fiber = column_fiber_base + column
        add_edge(fiber, column_port(column, False))
        add_edge(fiber, column_port(column, True))
    for row in range(ORDER):
        for column in range(ORDER):
            crossed = matrix[row][column] < 0
            add_edge(
                row_port(row, False),
                column_port(column, crossed),
            )
            add_edge(
                row_port(row, True),
                column_port(column, not crossed),
            )
    color_classes = [
        set(range(row_port_base, column_port_base)),
        set(range(column_port_base, row_fiber_base)),
        set(range(row_fiber_base, column_fiber_base)),
        set(range(column_fiber_base, vertex_count)),
    ]
    graph = pynauty.Graph(
        number_of_vertices=vertex_count,
        directed=False,
        adjacency_dict=adjacency,
        vertex_coloring=color_classes,
    )
    canonical = {
        "encoding": "signed-bipartite-double-cover-v1",
        "order": ORDER,
        "vertex_count": vertex_count,
        "color_classes": [sorted(color) for color in color_classes],
        "edges": [
            [first, second]
            for first in range(vertex_count)
            for second in sorted(adjacency[first])
            if first < second
        ],
    }
    return graph, canonical


def decode_generator(
    permutation: list[int],
    matrix: list[list[int]],
) -> dict[str, object]:
    column_port_base = 2 * ORDER
    row_fiber_base = 4 * ORDER
    column_fiber_base = 5 * ORDER
    row_permutation = [
        permutation[row_fiber_base + row] - row_fiber_base
        for row in range(ORDER)
    ]
    column_permutation = [
        permutation[column_fiber_base + column] - column_fiber_base
        for column in range(ORDER)
    ]
    if sorted(row_permutation) != list(range(ORDER)):
        raise ValueError("generator failed to permute row fibers")
    if sorted(column_permutation) != list(range(ORDER)):
        raise ValueError("generator failed to permute column fibers")
    row_switches: list[bool] = []
    column_switches: list[bool] = []
    for row, image_row in enumerate(row_permutation):
        plus_image = permutation[2 * row]
        minus_image = permutation[2 * row + 1]
        switched = plus_image == 2 * image_row + 1
        if plus_image != 2 * image_row + int(switched):
            raise ValueError("row plus-port escaped image fiber")
        if minus_image != 2 * image_row + int(not switched):
            raise ValueError("row antipodal pair was not preserved")
        row_switches.append(switched)
    for column, image_column in enumerate(column_permutation):
        plus_image = permutation[column_port_base + 2 * column]
        minus_image = permutation[column_port_base + 2 * column + 1]
        expected_plus = column_port_base + 2 * image_column
        switched = plus_image == expected_plus + 1
        if plus_image != expected_plus + int(switched):
            raise ValueError("column plus-port escaped image fiber")
        if minus_image != expected_plus + int(not switched):
            raise ValueError("column antipodal pair was not preserved")
        column_switches.append(switched)
    for row in range(ORDER):
        for column in range(ORDER):
            image_value = matrix[
                row_permutation[row]
            ][column_permutation[column]]
            if row_switches[row]:
                image_value = -image_value
            if column_switches[column]:
                image_value = -image_value
            if image_value != matrix[row][column]:
                raise ValueError("decoded generator does not stabilize factor")
    return {
        "row_permutation": row_permutation,
        "row_switch_mask": sum(
            switched << row for row, switched in enumerate(row_switches)
        ),
        "column_permutation": column_permutation,
        "column_switch_mask": sum(
            switched << column
            for column, switched in enumerate(column_switches)
        ),
    }


def act_cell(
    cell: tuple[int, int], generator: dict[str, object]
) -> tuple[int, int]:
    row_permutation = generator["row_permutation"]
    column_permutation = generator["column_permutation"]
    if not isinstance(row_permutation, list) or not isinstance(
        column_permutation, list
    ):
        raise TypeError("invalid decoded generator")
    row, column = cell
    return row_permutation[row], column_permutation[column]


def cell_orbits(
    generators: list[dict[str, object]],
) -> list[list[tuple[int, int]]]:
    unseen = {
        (row, column)
        for row in range(ORDER)
        for column in range(ORDER)
    }
    result: list[list[tuple[int, int]]] = []
    while unseen:
        start = min(unseen)
        orbit = {start}
        queue = [start]
        while queue:
            cell = queue.pop()
            for generator in generators:
                image = act_cell(cell, generator)
                if image not in orbit:
                    orbit.add(image)
                    queue.append(image)
        unseen.difference_update(orbit)
        result.append(sorted(orbit))
    result.sort(key=lambda orbit: (len(orbit), orbit[0]))
    if sum(map(len, result)) != ORDER * ORDER:
        raise ValueError("cell orbits do not cover all entries")
    return result


def unordered_cell_pair_orbits(
    generators: list[dict[str, object]],
) -> list[list[tuple[int, int]]]:
    """Radius-2 flip-set orbits, with each cell encoded row*ORDER+column."""
    cell_count = ORDER * ORDER
    unseen = {
        (first, second)
        for first in range(cell_count)
        for second in range(first + 1, cell_count)
    }
    result: list[list[tuple[int, int]]] = []

    def act_index(index: int, generator: dict[str, object]) -> int:
        row, column = divmod(index, ORDER)
        image_row, image_column = act_cell((row, column), generator)
        return image_row * ORDER + image_column

    while unseen:
        start = min(unseen)
        orbit = {start}
        queue = [start]
        while queue:
            first, second = queue.pop()
            for generator in generators:
                image = tuple(
                    sorted(
                        (
                            act_index(first, generator),
                            act_index(second, generator),
                        )
                    )
                )
                if image not in orbit:
                    orbit.add(image)
                    queue.append(image)
        unseen.difference_update(orbit)
        result.append(sorted(orbit))
    result.sort(key=lambda orbit: (len(orbit), orbit[0]))
    expected = cell_count * (cell_count - 1) // 2
    if sum(map(len, result)) != expected:
        raise ValueError("cell-pair orbits do not cover all radius-2 sets")
    return result


def analyze(path: Path, include_radius_two: bool) -> dict[str, object]:
    matrix = read_matrix(path)
    graph, canonical_graph = factor_graph(matrix)
    raw = pynauty.autgrp(graph)
    generators = [
        decode_generator(permutation, matrix) for permutation in raw[0]
    ]
    orbits = cell_orbits(generators)
    canonical_graph_bytes = (
        json.dumps(
            canonical_graph,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("ascii")
    orbit_partition = [
        [[row, column] for row, column in orbit] for orbit in orbits
    ]
    orbit_bytes = (
        json.dumps(
            orbit_partition,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("ascii")
    group_order = int(round(raw[1] * (10 ** raw[2])))
    result: dict[str, object] = {
        "order": ORDER,
        "source_matrix": str(path),
        "source_matrix_sha256": sha256(path),
        "graph_encoding": canonical_graph["encoding"],
        "graph_sha256": sha256_bytes(canonical_graph_bytes),
        "automorphism_group_order_on_double_cover": group_order,
        "effective_factor_automorphism_group_order": group_order // 2,
        "central_global_switch_kernel_order": 2,
        "generator_count": len(generators),
        "generators": generators,
        "entry_orbit_count": len(orbits),
        "entry_orbit_sizes": [len(orbit) for orbit in orbits],
        "entry_orbit_size_histogram": {
            str(size): count
            for size, count in sorted(Counter(map(len, orbits)).items())
        },
        "entry_orbit_representatives_zero_based": [
            list(orbit[0]) for orbit in orbits
        ],
        "entry_orbit_representatives_one_based": [
            [orbit[0][0] + 1, orbit[0][1] + 1] for orbit in orbits
        ],
        "entry_orbit_partition": orbit_partition,
        "entry_orbit_partition_sha256": sha256_bytes(orbit_bytes),
    }
    if include_radius_two:
        pair_orbits = unordered_cell_pair_orbits(generators)
        pair_partition = [
            [[first, second] for first, second in orbit]
            for orbit in pair_orbits
        ]
        pair_bytes = (
            json.dumps(
                pair_partition,
                sort_keys=True,
                separators=(",", ":"),
            )
            + "\n"
        ).encode("ascii")
        result.update(
            {
                "radius_two_orbit_count": len(pair_orbits),
                "radius_two_orbit_sizes": [
                    len(orbit) for orbit in pair_orbits
                ],
                "radius_two_orbit_size_histogram": {
                    str(size): count
                    for size, count in sorted(
                        Counter(map(len, pair_orbits)).items()
                    )
                },
                "radius_two_representatives_zero_based_flat": [
                    list(orbit[0]) for orbit in pair_orbits
                ],
                "radius_two_representatives_one_based": [
                    [
                        [
                            orbit[0][0] // ORDER + 1,
                            orbit[0][0] % ORDER + 1,
                        ],
                        [
                            orbit[0][1] // ORDER + 1,
                            orbit[0][1] % ORDER + 1,
                        ],
                    ]
                    for orbit in pair_orbits
                ],
                "radius_two_orbit_partition_sha256": sha256_bytes(pair_bytes),
            }
        )
    return result


def export_representatives(
    analysis: dict[str, object],
    source_path: Path,
    destination: Path,
) -> None:
    slug = matrix_slug(source_path)
    destination.mkdir(parents=True, exist_ok=True)
    radius_one_path = destination / f"{slug}-radius1-representatives.tsv"
    radius_one_representatives = analysis[
        "entry_orbit_representatives_one_based"
    ]
    if not isinstance(radius_one_representatives, list):
        raise TypeError("invalid radius-one representative list")
    radius_one_text = "".join(
        f"{row}\t{column}\n"
        for row, column in radius_one_representatives
    )
    atomic_write_text(radius_one_path, radius_one_text)
    exports: dict[str, object] = {
        "radius_one": {
            "format": "1-based rows: r c",
            "line_count": len(radius_one_representatives),
            "path": str(radius_one_path),
            "sha256": sha256(radius_one_path),
        }
    }
    radius_two_representatives = analysis.get(
        "radius_two_representatives_one_based"
    )
    if radius_two_representatives is not None:
        if not isinstance(radius_two_representatives, list):
            raise TypeError("invalid radius-two representative list")
        radius_two_path = (
            destination / f"{slug}-radius2-representatives.tsv"
        )
        radius_two_text = "".join(
            f"{first[0]}\t{first[1]}\t{second[0]}\t{second[1]}\n"
            for first, second in radius_two_representatives
        )
        atomic_write_text(radius_two_path, radius_two_text)
        exports["radius_two"] = {
            "format": "1-based rows: r1 c1 r2 c2",
            "line_count": len(radius_two_representatives),
            "path": str(radius_two_path),
            "sha256": sha256(radius_two_path),
        }
    analysis["representative_exports"] = exports


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", action="append", required=True, type=Path)
    parser.add_argument("--include-radius-two", action="store_true")
    parser.add_argument("--representative-dir", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    slugs = [matrix_slug(path) for path in args.matrix]
    if len(set(slugs)) != len(slugs):
        parser.error("--matrix paths must produce distinct export slugs")
    analyses = [
        analyze(path, args.include_radius_two) for path in args.matrix
    ]
    if args.representative_dir is not None:
        for analysis, source_path in zip(analyses, args.matrix):
            export_representatives(
                analysis,
                source_path,
                args.representative_dir,
            )
    engine_source = Path(__file__)
    report = {
        "schema_version": 1,
        "order": ORDER,
        "claim": (
            "Colored-graph automorphisms are exactly the side-preserving "
            "signed row/column stabilizer of each bound factor. Entry and "
            "unordered-entry-pair orbits are exact under that group."
        ),
        "claim_boundary": (
            "Transpose/side-swapping equivalences are deliberately excluded. "
            "The parser and graph construction fail closed unless every input "
            "is exactly a 22x22 sign matrix."
        ),
        "engine_source": str(engine_source),
        "engine_source_sha256": sha256(engine_source),
        "dependency": {
            "pynauty": getattr(pynauty, "__version__", "unknown"),
        },
        "analyses": analyses,
    }
    atomic_write_text(
        args.output,
        json.dumps(report, indent=2, sort_keys=True) + "\n",
    )
    for item in analyses:
        print(
            item["source_matrix"],
            "group",
            item["effective_factor_automorphism_group_order"],
            "entry-orbits",
            item["entry_orbit_count"],
            "sizes",
            item["entry_orbit_sizes"],
            "radius-two-orbits",
            item.get("radius_two_orbit_count"),
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
