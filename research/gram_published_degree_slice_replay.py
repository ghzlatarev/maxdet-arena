#!/usr/bin/env python3
"""Independent exact replay for gram_published_degree_slice.cpp output."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


ORDER = 23
FRONTIER = 2_779_447_296_000_000
EXPECTED_LABELED = 3_488_400
EXPECTED_ORBITS = 20
REQUIRED_DIVISOR = 1 << 22


class ReplayError(RuntimeError):
    """The retained report failed an independent exact check."""


def bareiss_determinant(matrix: list[list[int]]) -> int:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ReplayError("Bareiss input must be a nonempty square matrix")
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot_row = next(
            (
                row
                for row in range(column, size)
                if work[row][column] != 0
            ),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, size):
            for inner in range(column + 1, size):
                numerator = (
                    work[row][inner] * pivot
                    - work[row][column] * work[column][inner]
                )
                quotient, remainder = divmod(numerator, previous)
                if remainder:
                    raise ReplayError("Bareiss division was not exact")
                work[row][inner] = quotient
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def build_gram(added_edges: Any) -> list[list[int]]:
    if not isinstance(added_edges, list) or len(added_edges) != 12:
        raise ReplayError("each orbit must retain exactly 12 added edges")
    edges: set[tuple[int, int]] = set()
    for left in range(3):
        for right in range(left + 1, 3):
            edges.add((left, right))
    for block in range(5):
        offset = 3 + 4 * block
        for left in range(4):
            for right in range(left + 1, 4):
                edges.add((offset + left, offset + right))

    added: set[tuple[int, int]] = set()
    for raw_edge in added_edges:
        if (
            not isinstance(raw_edge, list)
            or len(raw_edge) != 2
            or any(type(vertex) is not int for vertex in raw_edge)
        ):
            raise ReplayError("malformed added edge")
        left, right = raw_edge[0] - 1, raw_edge[1] - 1
        if not (0 <= left < 3 <= right < ORDER):
            raise ReplayError(
                "added edge is not an absent K3-to-K4-side edge"
            )
        if (left, right) in added:
            raise ReplayError("duplicate added edge")
        added.add((left, right))
    edges |= added
    if len(edges) != 45:
        raise ReplayError("reconstructed defect graph must have 45 edges")

    gram = [
        [23 if row == column else -1 for column in range(ORDER)]
        for row in range(ORDER)
    ]
    for left, right in edges:
        gram[left][right] = 3
        gram[right][left] = 3

    added_degrees = [0] * ORDER
    for left, right in added:
        added_degrees[left] += 1
        added_degrees[right] += 1
    if added_degrees[:3] != [4, 4, 4]:
        raise ReplayError("K3 added degrees are not 4,4,4")
    if Counter(added_degrees[3:]) != Counter({0: 14, 2: 6}):
        raise ReplayError("K4-side added degrees are not 2^6,0^14")
    return gram


def positive_definite(gram: list[list[int]]) -> bool:
    return all(
        bareiss_determinant([row[:size] for row in gram[:size]]) > 0
        for size in range(1, ORDER + 1)
    )


def block_options() -> list[tuple[int, int, int]]:
    return [
        (first, second, third)
        for first in range(3)
        for second in range(3)
        for third in range(3)
        if first + second + third <= 4
    ]


def ordered_profiles(
    block: int = 0,
    remaining: tuple[int, int, int] = (2, 2, 2),
    prefix: tuple[tuple[int, int, int], ...] = (),
) -> Iterable[tuple[tuple[int, int, int], ...]]:
    if block == 5:
        if remaining == (0, 0, 0):
            yield prefix
        return
    for option in block_options():
        if all(option[color] <= remaining[color] for color in range(3)):
            yield from ordered_profiles(
                block + 1,
                tuple(
                    remaining[color] - option[color]
                    for color in range(3)
                ),
                prefix + (option,),
            )


def profile_key(profile: Iterable[tuple[int, int, int]]) -> int:
    blocks = tuple(profile)
    candidates = []
    for permutation in itertools.permutations(range(3)):
        transformed = sorted(
            tuple(block[color] for color in permutation)
            for block in blocks
        )
        key = 0
        for first, second, third in transformed:
            key = key * 125 + 25 * first + 5 * second + third
        candidates.append(key)
    return min(candidates)


def minimum_radius(profile: list[list[int]]) -> int:
    best_common_edges = 0
    for assigned_blocks in itertools.permutations(range(5), 3):
        common_edges = sum(
            min(sum(profile[assigned_blocks[color]]), 2)
            + min(profile[assigned_blocks[color]][color], 2)
            for color in range(3)
        )
        best_common_edges = max(best_common_edges, common_edges)
    return 12 - best_common_edges


def replay(report_path: Path) -> dict[str, Any]:
    report_bytes = report_path.read_bytes()
    report = json.loads(report_bytes)
    if not isinstance(report, dict):
        raise ReplayError("report root must be an object")
    family = report.get("family")
    orbits = report.get("orbits")
    if not isinstance(family, dict) or not isinstance(orbits, list):
        raise ReplayError("report family/orbits fields are malformed")
    if family.get("direct_labeled_count") != EXPECTED_LABELED:
        raise ReplayError("report labeled count mismatch")
    if len(orbits) != EXPECTED_ORBITS:
        raise ReplayError("report orbit count mismatch")
    if family.get("no_base_edge_deletions") is not True:
        raise ReplayError("report does not pin the no-deletion scope")

    generated = list(ordered_profiles())
    block_multisets = {tuple(sorted(profile)) for profile in generated}
    generated_keys = {profile_key(profile) for profile in generated}
    report_keys = {
        int(orbit["canonical_key"])
        for orbit in orbits
    }
    if len(block_multisets) != 61:
        raise ReplayError("expected 61 block multisets before S3 quotient")
    if generated_keys != report_keys or len(generated_keys) != EXPECTED_ORBITS:
        raise ReplayError("independent profile quotient disagrees with report")

    exact: list[dict[str, Any]] = []
    radius_histogram: Counter[int] = Counter()
    for expected_index, orbit in enumerate(orbits):
        if orbit.get("orbit_index") != expected_index:
            raise ReplayError("orbit indices are not contiguous")
        gram = build_gram(orbit.get("added_edges"))
        determinant = bareiss_determinant(gram)
        if str(determinant) != orbit.get("determinant"):
            raise ReplayError(
                f"orbit {expected_index} determinant mismatch"
            )
        root = math.isqrt(determinant) if determinant >= 0 else 0
        square = determinant >= 0 and root * root == determinant
        if orbit.get("is_square") is not square:
            raise ReplayError(f"orbit {expected_index} square mismatch")
        if square and orbit.get("square_root") != str(root):
            raise ReplayError(f"orbit {expected_index} root mismatch")
        divisible = square and root % REQUIRED_DIVISOR == 0
        if orbit.get("divisible_by_2_22") is not divisible:
            raise ReplayError(
                f"orbit {expected_index} divisibility mismatch"
            )
        pd = positive_definite(gram)
        if orbit.get("positive_definite") is not pd:
            raise ReplayError(
                f"orbit {expected_index} positive-definite mismatch"
            )
        profile = orbit.get("block_missing_color_counts")
        if (
            not isinstance(profile, list)
            or len(profile) != 5
            or minimum_radius(profile)
            != orbit.get("minimum_swap_radius_from_published")
        ):
            raise ReplayError(
                f"orbit {expected_index} minimum-radius mismatch"
            )
        radius_histogram[minimum_radius(profile)] += 1
        exact.append(
            {
                "orbit_index": expected_index,
                "determinant": str(determinant),
                "positive_definite": pd,
                "is_square": square,
                **({"square_root": str(root)} if square else {}),
            }
        )

    determinants = [int(item["determinant"]) for item in exact]
    square_roots = [
        int(item["square_root"]) for item in exact if item["is_square"]
    ]
    expected_statistics = {
        "exact_determinants": EXPECTED_ORBITS,
        "positive_determinants": sum(value > 0 for value in determinants),
        "positive_definite_orbits": sum(
            bool(item["positive_definite"]) for item in exact
        ),
        "above_frontier_determinants": sum(
            value > FRONTIER * FRONTIER for value in determinants
        ),
        "exact_squares": len(square_roots),
        "frontier_ties": square_roots.count(FRONTIER),
        "above_frontier_squares": sum(
            root > FRONTIER for root in square_roots
        ),
    }
    statistics = report.get("statistics")
    if not isinstance(statistics, dict) or any(
        statistics.get(key) != value
        for key, value in expected_statistics.items()
    ):
        raise ReplayError("report aggregate statistics mismatch")

    maximum = max(
        exact, key=lambda item: int(item["determinant"])
    )
    return {
        "schema_version": 1,
        "engine": "gram-published-degree-slice-independent-replay",
        "complete": True,
        "report": {
            "path": str(report_path),
            "sha256": hashlib.sha256(report_bytes).hexdigest(),
        },
        "checks": {
            "arbitrary_precision_bareiss_orbits": len(exact),
            "exact_sylvester_pd_orbits": sum(
                bool(item["positive_definite"]) for item in exact
            ),
            "ordered_block_profiles": len(generated),
            "block_multisets_before_s3": len(block_multisets),
            "base_automorphism_orbits": len(generated_keys),
            "radius_histogram": {
                str(radius): radius_histogram[radius]
                for radius in sorted(radius_histogram)
            },
            "scope_and_degree_constraints_checked": True,
        },
        "maximum": maximum,
        "square_roots": [str(root) for root in sorted(square_roots)],
        "statistics": expected_statistics,
        "termination": "completed",
    }


def atomic_write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(value, handle, sort_keys=True, separators=(",", ":"))
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    result = replay(arguments.report)
    atomic_write_json(arguments.output, result)
    print(
        "independent replay passed "
        f"orbits={result['statistics']['exact_determinants']} "
        f"squares={result['statistics']['exact_squares']} "
        f"output={arguments.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
