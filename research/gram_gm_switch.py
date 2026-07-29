#!/usr/bin/env python3
"""Exhaust and audit one-step Godsil--McKay switches of Gram defect graphs.

The source normalization is G = 24I - J + 4A.  For every Hasse-surviving
source hit, the companion enumerator exhausts every proper, even set X such
that:

* A[X] is regular; and
* every vertex outside X has 0, |X|/2, or |X| neighbors in X.

Sets with no |X|/2 outside vertex are omitted because their switch is a no-op.
Every emitted mate is independently reconstructed in Python.  Exact graph
isomorphism (NetworkX VF2) deduplicates mates, and every retained class is
checked with integer Bareiss determinants, Sylvester positive-definiteness,
and the scaled identity S G S^T = |X|^2 G', where S = |X| Q.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import networkx as nx


ORDER = 23
EDGE_COUNT = ORDER * (ORDER - 1) // 2
NORMALIZATION = "G=24I-J+4A"
FRONTIER_ROOT = 2_779_447_296_000_000
SOURCE_ENGINE = "gram-tabu"
ROUTE_ENGINE = "gram-gm-switch"
CHALLENGE_ID = "maxdet-23-v1"


class InputError(ValueError):
    """A bound input, checkpoint, or enumerator record is malformed."""


def sha256_bytes(contents: bytes) -> str:
    return hashlib.sha256(contents).hexdigest()


def sha256_path(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_json(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def atomic_write_json(path: Path, value: Any) -> None:
    contents = json.dumps(value, indent=2, sort_keys=True) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as output:
            output.write(contents)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def load_json_bytes(path: Path) -> tuple[dict[str, Any], bytes]:
    try:
        contents = path.read_bytes()
        value = json.loads(contents)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise InputError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise InputError(f"{path} must contain a JSON object")
    return value, contents


def parse_decimal(value: Any, field: str) -> int:
    if not isinstance(value, str) or not value or not value.isascii():
        raise InputError(f"{field} must be an ASCII decimal string")
    if not value.isdigit():
        raise InputError(f"{field} must contain decimal digits only")
    return int(value)


def edge_index(order: int, first: int, second: int) -> int:
    if not 0 <= first < second < order:
        raise InputError("edge endpoints are not in strict upper-triangle order")
    return first * (2 * order - first - 1) // 2 + second - first - 1


def popcount(value: int) -> int:
    return bin(value).count("1")


def edge_mask_from_pairs(
    order: int, pairs: Iterable[tuple[int, int]]
) -> int:
    result = 0
    for first, second in pairs:
        bit = 1 << edge_index(order, first, second)
        if result & bit:
            raise InputError("duplicate edge")
        result |= bit
    return result


def edge_pairs_from_mask(order: int, edge_mask: int) -> list[tuple[int, int]]:
    edge_count = order * (order - 1) // 2
    if type(edge_mask) is not int or edge_mask < 0 or edge_mask >= 1 << edge_count:
        raise InputError("edge mask is outside the graph")
    result: list[tuple[int, int]] = []
    bit = 0
    for first in range(order):
        for second in range(first + 1, order):
            if edge_mask >> bit & 1:
                result.append((first, second))
            bit += 1
    return result


def adjacency_masks(order: int, edge_mask: int) -> list[int]:
    result = [0] * order
    for first, second in edge_pairs_from_mask(order, edge_mask):
        result[first] |= 1 << second
        result[second] |= 1 << first
    return result


def parse_hit_edges(hit: dict[str, Any]) -> int:
    edges = hit.get("edges")
    edge_count = hit.get("edge_count")
    if not isinstance(edges, list):
        raise InputError("hit.edges must be a list")
    if type(edge_count) is not int or edge_count != len(edges):
        raise InputError("hit.edge_count disagrees with hit.edges")
    parsed: list[tuple[int, int]] = []
    for index, edge in enumerate(edges):
        if (
            not isinstance(edge, list)
            or len(edge) != 2
            or type(edge[0]) is not int
            or type(edge[1]) is not int
        ):
            raise InputError(f"hit.edges[{index}] must contain two integers")
        first, second = edge
        if not 1 <= first < second <= ORDER:
            raise InputError(
                f"hit.edges[{index}] must satisfy 1 <= first < second <= {ORDER}"
            )
        parsed.append((first - 1, second - 1))
    return edge_mask_from_pairs(ORDER, parsed)


def build_gram(order: int, edge_mask: int) -> list[list[int]]:
    adjacency = adjacency_masks(order, edge_mask)
    return [
        [
            order
            if row == column
            else 3
            if adjacency[row] >> column & 1
            else -1
            for column in range(order)
        ]
        for row in range(order)
    ]


def bareiss_determinant(matrix: list[list[int]]) -> int:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise InputError("determinant input must be a non-empty square matrix")
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot_row = next(
            (row for row in range(column, size) if work[row][column] != 0),
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
                    raise ArithmeticError("non-exact Bareiss division")
                work[row][inner] = quotient
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def exact_positive_definite(matrix: list[list[int]]) -> bool:
    return all(
        bareiss_determinant([row[:size] for row in matrix[:size]]) > 0
        for size in range(1, len(matrix) + 1)
    )


def matrix_multiply(
    left: list[list[int]], right: list[list[int]]
) -> list[list[int]]:
    if not left or not right or len(left[0]) != len(right):
        raise InputError("incompatible matrix product")
    return [
        [
            sum(left[row][inner] * right[inner][column] for inner in range(len(right)))
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def transpose(matrix: list[list[int]]) -> list[list[int]]:
    return [list(column) for column in zip(*matrix)]


def scaled_switch_matrix(order: int, switch_set: int) -> tuple[list[list[int]], int]:
    if switch_set <= 0 or switch_set >= 1 << order:
        raise InputError("switching set is empty or outside the graph")
    size = popcount(switch_set)
    result = [[0] * order for _ in range(order)]
    for row in range(order):
        for column in range(order):
            row_inside = switch_set >> row & 1
            column_inside = switch_set >> column & 1
            if row_inside and column_inside:
                result[row][column] = 2 - (size if row == column else 0)
            elif not row_inside and not column_inside and row == column:
                result[row][column] = size
    return result, size


def validate_switch_set(
    order: int, edge_mask: int, switch_set: int
) -> tuple[list[int], int]:
    if (
        type(switch_set) is not int
        or switch_set <= 0
        or switch_set >= 1 << order
    ):
        raise InputError("invalid switching set mask")
    size = popcount(switch_set)
    if size < 2 or size % 2 or size == order:
        raise InputError("switching set must be proper, even, and have size >= 2")
    adjacency = adjacency_masks(order, edge_mask)
    inside = [vertex for vertex in range(order) if switch_set >> vertex & 1]
    induced_degrees = {
        popcount(adjacency[vertex] & switch_set) for vertex in inside
    }
    if len(induced_degrees) != 1:
        raise InputError("switching set does not induce a regular graph")
    half = size // 2
    half_vertices: list[int] = []
    for vertex in range(order):
        if switch_set >> vertex & 1:
            continue
        count = popcount(adjacency[vertex] & switch_set)
        if count == half:
            half_vertices.append(vertex)
        elif count not in (0, size):
            raise InputError("outside vertex violates the GM neighbor condition")
    if not half_vertices:
        raise InputError("switch is a no-op: there is no half-neighbor vertex")
    return half_vertices, next(iter(induced_degrees))


def switch_edge_mask(order: int, edge_mask: int, switch_set: int) -> int:
    half_vertices, _ = validate_switch_set(order, edge_mask, switch_set)
    result = edge_mask
    for outside in half_vertices:
        remaining = switch_set
        while remaining:
            lowest = remaining & -remaining
            inside = lowest.bit_length() - 1
            first, second = sorted((inside, outside))
            result ^= 1 << edge_index(order, first, second)
            remaining ^= lowest
    return result


def graph_from_mask(order: int, edge_mask: int) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(order))
    graph.add_edges_from(edge_pairs_from_mask(order, edge_mask))
    return graph


def graph_bucket_key(graph: nx.Graph) -> tuple[Any, ...]:
    degrees = tuple(sorted(degree for _, degree in graph.degree()))
    wl_hash = nx.weisfeiler_lehman_graph_hash(
        graph, iterations=5, digest_size=16
    )
    return (graph.number_of_edges(), degrees, wl_hash)


def exact_isomorphic(left: nx.Graph, right: nx.Graph) -> bool:
    return nx.algorithms.isomorphism.GraphMatcher(left, right).is_isomorphic()


def parse_expected_hash(value: str, field: str) -> str:
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
        raise InputError(f"{field} must be 64 lowercase hex digits")
    return value


def validate_source_hit(hit: dict[str, Any], hit_index: int) -> dict[str, Any]:
    edge_mask = parse_hit_edges(hit)
    gram = build_gram(ORDER, edge_mask)
    determinant = bareiss_determinant(gram)
    claimed_determinant = parse_decimal(
        hit.get("determinant"), f"hits[{hit_index}].determinant"
    )
    if determinant != claimed_determinant:
        raise InputError(f"hits[{hit_index}] determinant mismatch")
    root = math.isqrt(determinant)
    if root * root != determinant:
        raise InputError(f"hits[{hit_index}] determinant is not a square")
    if parse_decimal(hit.get("square_root"), f"hits[{hit_index}].square_root") != root:
        raise InputError(f"hits[{hit_index}] square-root mismatch")
    positive = exact_positive_definite(gram)
    divisible = root % (1 << 22) == 0
    if hit.get("positive_definite") is not positive:
        raise InputError(f"hits[{hit_index}] positive-definite flag mismatch")
    if hit.get("divisible_by_2_22") is not divisible:
        raise InputError(f"hits[{hit_index}] divisibility flag mismatch")
    if hit.get("qualified") is not (positive and divisible and root > FRONTIER_ROOT):
        raise InputError(f"hits[{hit_index}] qualified flag mismatch")
    return {
        "hit_index": hit_index,
        "edge_mask": edge_mask,
        "edge_mask_hex": f"{edge_mask:064x}",
        "determinant": determinant,
        "square_root": root,
        "hit_identity_sha256": sha256_bytes(canonical_json(hit)),
    }


def validate_bound_sources(
    snapshot_path: Path,
    hasse_path: Path,
    expected_snapshot_sha256: str,
    expected_hasse_sha256: str,
) -> tuple[
    dict[str, Any],
    list[dict[str, Any]],
    dict[str, Any],
    bytes,
    bytes,
]:
    snapshot, snapshot_bytes = load_json_bytes(snapshot_path)
    hasse, hasse_bytes = load_json_bytes(hasse_path)
    expected_snapshot_sha256 = parse_expected_hash(
        expected_snapshot_sha256, "--expected-snapshot-sha256"
    )
    expected_hasse_sha256 = parse_expected_hash(
        expected_hasse_sha256, "--expected-hasse-sha256"
    )
    if sha256_bytes(snapshot_bytes) != expected_snapshot_sha256:
        raise InputError("source snapshot SHA-256 mismatch")
    if sha256_bytes(hasse_bytes) != expected_hasse_sha256:
        raise InputError("Hasse report SHA-256 mismatch")
    if snapshot.get("challenge_id") != CHALLENGE_ID:
        raise InputError("snapshot challenge_id mismatch")
    if snapshot.get("engine") != SOURCE_ENGINE:
        raise InputError("snapshot engine must be gram-tabu")
    if snapshot.get("schema_version") != 1:
        raise InputError("snapshot schema_version must be 1")
    if snapshot.get("normalization") != NORMALIZATION:
        raise InputError("snapshot normalization mismatch")
    hits = snapshot.get("hits")
    if not isinstance(hits, list) or len(hits) != 256:
        raise InputError("source snapshot must contain exactly 256 hits")
    if hasse.get("engine") != "gram-hasse" or hasse.get("schema_version") != 1:
        raise InputError("invalid Hasse report engine or schema")
    report_source = hasse.get("snapshot")
    if not isinstance(report_source, dict):
        raise InputError("Hasse report is missing snapshot binding")
    if report_source.get("sha256") != expected_snapshot_sha256:
        raise InputError("Hasse report binds a different source snapshot")
    if report_source.get("stored_hit_count") != len(hits):
        raise InputError("Hasse report source hit count mismatch")
    selection = hasse.get("selection")
    if not isinstance(selection, dict) or selection.get("hit_indices") != list(
        range(len(hits))
    ):
        raise InputError("Hasse report must select every source hit in order")
    results = hasse.get("results")
    if not isinstance(results, list) or len(results) != len(hits):
        raise InputError("Hasse report result count mismatch")
    statuses = ("rejected", "no_obstruction", "error")
    counts = {status: 0 for status in statuses}
    survivors: list[dict[str, Any]] = []
    seen: set[int] = set()
    for position, result in enumerate(results):
        if not isinstance(result, dict):
            raise InputError(f"Hasse results[{position}] is not an object")
        hit_index = result.get("hit_index")
        if type(hit_index) is not int or not 0 <= hit_index < len(hits):
            raise InputError(f"Hasse results[{position}] has invalid hit_index")
        if hit_index in seen:
            raise InputError("Hasse report contains duplicate hit indices")
        seen.add(hit_index)
        if hit_index != position:
            raise InputError("Hasse report result order is not source order")
        status = result.get("status")
        if status not in counts:
            raise InputError(f"Hasse results[{position}] has invalid status")
        counts[status] += 1
        hit = hits[hit_index]
        if not isinstance(hit, dict):
            raise InputError(f"source hits[{hit_index}] is not an object")
        for field in ("determinant", "square_root"):
            if result.get(field) != hit.get(field):
                raise InputError(
                    f"Hasse results[{position}] {field} disagrees with source hit"
                )
        if status == "no_obstruction":
            survivors.append(validate_source_hit(hit, hit_index))
    expected_summary = {"selected": len(hits), **counts}
    if hasse.get("summary") != expected_summary:
        raise InputError("Hasse report summary mismatch")
    if len(survivors) != 50:
        raise InputError(
            f"expected exactly 50 no-obstruction hits, found {len(survivors)}"
        )
    return snapshot, survivors, hasse, snapshot_bytes, hasse_bytes


def parse_enumerator_output(
    contents: str,
    order: int,
    source_edge_mask: int,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    records: list[dict[str, Any]] = []
    summary: dict[str, int] | None = None
    seen_sets: set[int] = set()
    edge_count = order * (order - 1) // 2
    expected_digits = (edge_count + 3) // 4
    for line_number, line in enumerate(contents.splitlines(), 1):
        if line.startswith("#summary\t"):
            if summary is not None:
                raise InputError("enumerator emitted multiple summary rows")
            parsed: dict[str, int] = {}
            for item in line.split("\t")[1:]:
                if "=" not in item:
                    raise InputError("malformed enumerator summary item")
                key, value = item.split("=", 1)
                if not value.isascii() or not value.isdigit():
                    raise InputError("non-decimal enumerator summary value")
                parsed[key] = int(value)
            required = {"masks_examined", "valid_sets", "labeled_mates"}
            if set(parsed) != required:
                raise InputError("enumerator summary fields mismatch")
            summary = parsed
            continue
        if summary is not None:
            raise InputError("enumerator emitted data after its summary")
        fields = line.split("\t")
        if len(fields) != 2:
            raise InputError(f"malformed enumerator row {line_number}")
        set_text, mate_text = fields
        try:
            switch_set = int(set_text, 16)
            mate_edge_mask = int(mate_text, 16)
        except ValueError as error:
            raise InputError(f"non-hex enumerator row {line_number}") from error
        if (
            not set_text
            or any(character not in "0123456789abcdef" for character in set_text)
            or len(mate_text) != expected_digits
            or any(character not in "0123456789abcdef" for character in mate_text)
        ):
            raise InputError(f"noncanonical hex in enumerator row {line_number}")
        if switch_set in seen_sets:
            raise InputError("enumerator emitted a duplicate switching set")
        seen_sets.add(switch_set)
        expected_mate = switch_edge_mask(
            order, source_edge_mask, switch_set
        )
        if mate_edge_mask != expected_mate:
            raise InputError("enumerator mate disagrees with Python reconstruction")
        records.append(
            {
                "switch_set_hex": format(switch_set, "x"),
                "switch_set_size": popcount(switch_set),
                "mate_edge_mask_hex": f"{mate_edge_mask:0{expected_digits}x}",
            }
        )
    if summary is None:
        raise InputError("enumerator output is missing its summary")
    expected_examined = (1 << (order - 1)) - 1 - (order % 2 == 0)
    if summary["masks_examined"] != expected_examined:
        raise InputError(
            "enumerator did not examine every proper nonempty even subset"
        )
    if summary["valid_sets"] != len(records):
        raise InputError("enumerator valid-set count mismatch")
    if summary["labeled_mates"] != len(records):
        raise InputError("enumerator labeled-mate count mismatch")
    return records, summary


def checkpoint_for_source(
    survivor: dict[str, Any],
    enumerator: Path,
    enumerator_sha256: str,
    source_sha256: str,
    hasse_sha256: str,
) -> dict[str, Any]:
    command = [
        str(enumerator),
        "--order",
        str(ORDER),
        "--edge-mask",
        survivor["edge_mask_hex"],
    ]
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    if completed.returncode != 0:
        raise InputError(
            f"enumerator failed for hit {survivor['hit_index']}: "
            f"{completed.stderr.strip()}"
        )
    if completed.stderr:
        raise InputError("enumerator unexpectedly wrote to standard error")
    records, summary = parse_enumerator_output(
        completed.stdout, ORDER, survivor["edge_mask"]
    )
    unique: dict[str, dict[str, Any]] = {}
    for record in records:
        edge_hex = record["mate_edge_mask_hex"]
        prior = unique.get(edge_hex)
        if prior is None or int(record["switch_set_hex"], 16) < int(
            prior["switch_set_hex"], 16
        ):
            unique[edge_hex] = record
    return {
        "schema_version": 1,
        "engine": "gram-gm-switch-checkpoint",
        "complete": True,
        "source_snapshot_sha256": source_sha256,
        "source_hasse_sha256": hasse_sha256,
        "enumerator_sha256": enumerator_sha256,
        "source_hit_index": survivor["hit_index"],
        "source_hit_identity_sha256": survivor["hit_identity_sha256"],
        "source_edge_mask_hex": survivor["edge_mask_hex"],
        "source_determinant": str(survivor["determinant"]),
        "source_square_root": str(survivor["square_root"]),
        "masks_examined": summary["masks_examined"],
        "valid_switching_sets": summary["valid_sets"],
        "unique_labeled_mates": len(unique),
        "labeled_mates": [unique[key] for key in sorted(unique)],
    }


def validate_checkpoint(
    checkpoint: dict[str, Any],
    survivor: dict[str, Any],
    enumerator_sha256: str,
    source_sha256: str,
    hasse_sha256: str,
) -> None:
    expected = {
        "schema_version": 1,
        "engine": "gram-gm-switch-checkpoint",
        "complete": True,
        "source_snapshot_sha256": source_sha256,
        "source_hasse_sha256": hasse_sha256,
        "enumerator_sha256": enumerator_sha256,
        "source_hit_index": survivor["hit_index"],
        "source_hit_identity_sha256": survivor["hit_identity_sha256"],
        "source_edge_mask_hex": survivor["edge_mask_hex"],
        "source_determinant": str(survivor["determinant"]),
        "source_square_root": str(survivor["square_root"]),
        "masks_examined": (1 << (ORDER - 1)) - 1,
    }
    for field, value in expected.items():
        if checkpoint.get(field) != value:
            raise InputError(
                f"checkpoint for hit {survivor['hit_index']} has stale {field}"
            )
    records = checkpoint.get("labeled_mates")
    if not isinstance(records, list):
        raise InputError("checkpoint labeled_mates must be a list")
    if checkpoint.get("unique_labeled_mates") != len(records):
        raise InputError("checkpoint unique-labeled-mate count mismatch")
    if (
        type(checkpoint.get("valid_switching_sets")) is not int
        or checkpoint["valid_switching_sets"] < len(records)
    ):
        raise InputError("checkpoint valid-switching-set count mismatch")
    seen_edges: set[str] = set()
    for record in records:
        if not isinstance(record, dict):
            raise InputError("checkpoint mate record must be an object")
        set_text = record.get("switch_set_hex")
        edge_text = record.get("mate_edge_mask_hex")
        if not isinstance(set_text, str) or not isinstance(edge_text, str):
            raise InputError("checkpoint mate record is missing hex masks")
        try:
            switch_set = int(set_text, 16)
            edge_mask = int(edge_text, 16)
        except ValueError as error:
            raise InputError("checkpoint mate record contains non-hex") from error
        if edge_text in seen_edges:
            raise InputError("checkpoint contains duplicate labeled mates")
        seen_edges.add(edge_text)
        expected_edge_mask = switch_edge_mask(
            ORDER, survivor["edge_mask"], switch_set
        )
        if edge_mask != expected_edge_mask:
            raise InputError("checkpoint mate failed exact reconstruction")
        if record.get("switch_set_size") != popcount(switch_set):
            raise InputError("checkpoint switch-set size mismatch")


def build_original_registry(
    snapshot: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[tuple[Any, ...], list[int]]]:
    hits = snapshot["hits"]
    representatives: list[dict[str, Any]] = []
    buckets: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for hit_index, hit in enumerate(hits):
        if not isinstance(hit, dict):
            raise InputError(f"source hits[{hit_index}] must be an object")
        edge_mask = parse_hit_edges(hit)
        graph = graph_from_mask(ORDER, edge_mask)
        key = graph_bucket_key(graph)
        matched: int | None = None
        for class_index in buckets[key]:
            if exact_isomorphic(graph, representatives[class_index]["graph"]):
                matched = class_index
                break
        if matched is None:
            matched = len(representatives)
            representatives.append(
                {
                    "graph": graph,
                    "edge_mask": edge_mask,
                    "hit_indices": [hit_index],
                }
            )
            buckets[key].append(matched)
        else:
            representatives[matched]["hit_indices"].append(hit_index)
    return representatives, buckets


def deduplicate_mates(
    checkpoints: list[dict[str, Any]],
    survivors_by_index: dict[int, dict[str, Any]],
    original_representatives: list[dict[str, Any]],
    original_buckets: dict[tuple[Any, ...], list[int]],
) -> tuple[list[dict[str, Any]], int]:
    classes: list[dict[str, Any]] = []
    buckets: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    exact_comparisons = 0
    for checkpoint in checkpoints:
        source_index = checkpoint["source_hit_index"]
        survivor = survivors_by_index[source_index]
        for record in checkpoint["labeled_mates"]:
            edge_mask = int(record["mate_edge_mask_hex"], 16)
            graph = graph_from_mask(ORDER, edge_mask)
            key = graph_bucket_key(graph)
            matched: int | None = None
            for class_index in buckets[key]:
                exact_comparisons += 1
                if exact_isomorphic(graph, classes[class_index]["graph"]):
                    matched = class_index
                    break
            witness = {
                "source_hit_index": source_index,
                "switch_set_hex": record["switch_set_hex"],
                "switch_set_size": record["switch_set_size"],
            }
            if matched is None:
                original_matches: list[int] = []
                for original_index in original_buckets.get(key, []):
                    exact_comparisons += 1
                    if exact_isomorphic(
                        graph, original_representatives[original_index]["graph"]
                    ):
                        original_matches.extend(
                            original_representatives[original_index]["hit_indices"]
                        )
                matched = len(classes)
                classes.append(
                    {
                        "graph": graph,
                        "edge_mask": edge_mask,
                        "edge_mask_hex": record["mate_edge_mask_hex"],
                        "witness": witness,
                        "source_hit_indices": {source_index},
                        "labeled_mate_multiplicity": 1,
                        "original_snapshot_hit_indices": sorted(original_matches),
                        "square_root": survivor["square_root"],
                        "determinant": survivor["determinant"],
                    }
                )
                buckets[key].append(matched)
            else:
                item = classes[matched]
                if item["square_root"] != survivor["square_root"]:
                    raise ArithmeticError(
                        "isomorphic mates have inconsistent determinants"
                    )
                item["source_hit_indices"].add(source_index)
                item["labeled_mate_multiplicity"] += 1
                if record["mate_edge_mask_hex"] < item["edge_mask_hex"]:
                    item["graph"] = graph
                    item["edge_mask"] = edge_mask
                    item["edge_mask_hex"] = record["mate_edge_mask_hex"]
                    item["witness"] = witness
    classes.sort(key=lambda item: (item["square_root"], item["edge_mask_hex"]))
    return classes, exact_comparisons


def audit_retained_class(
    item: dict[str, Any], survivors_by_index: dict[int, dict[str, Any]]
) -> dict[str, Any]:
    witness = item["witness"]
    source = survivors_by_index[witness["source_hit_index"]]
    switch_set = int(witness["switch_set_hex"], 16)
    expected_mate = switch_edge_mask(ORDER, source["edge_mask"], switch_set)
    if expected_mate != item["edge_mask"]:
        raise ArithmeticError("retained-class witness does not reconstruct mate")
    source_gram = build_gram(ORDER, source["edge_mask"])
    mate_gram = build_gram(ORDER, item["edge_mask"])
    scaled_q, size = scaled_switch_matrix(ORDER, switch_set)
    scaled_orthogonality = matrix_multiply(scaled_q, transpose(scaled_q))
    expected_identity = [
        [size * size if row == column else 0 for column in range(ORDER)]
        for row in range(ORDER)
    ]
    if scaled_orthogonality != expected_identity:
        raise ArithmeticError("scaled switching matrix is not orthogonal")
    if [sum(row) for row in scaled_q] != [size] * ORDER:
        raise ArithmeticError("scaled switching matrix does not fix all-ones")
    conjugated = matrix_multiply(
        matrix_multiply(scaled_q, source_gram), transpose(scaled_q)
    )
    expected_conjugated = [
        [size * size * value for value in row] for row in mate_gram
    ]
    if conjugated != expected_conjugated:
        raise ArithmeticError("exact scaled Gram conjugacy check failed")
    source_determinant = bareiss_determinant(source_gram)
    mate_determinant = bareiss_determinant(mate_gram)
    if source_determinant != source["determinant"]:
        raise ArithmeticError("source determinant changed during final audit")
    if mate_determinant != source_determinant:
        raise ArithmeticError("mate determinant is not preserved")
    if not exact_positive_definite(mate_gram):
        raise ArithmeticError("mate failed exact positive-definiteness")
    return {
        "scaled_conjugacy_checked": True,
        "scaled_orthogonality_checked": True,
        "all_ones_fixed_checked": True,
        "exact_determinant_checked": True,
        "exact_positive_definite_checked": True,
        "switching_set_size": size,
    }


def route_hit(item: dict[str, Any], class_index: int) -> dict[str, Any]:
    edges = [
        [first + 1, second + 1]
        for first, second in edge_pairs_from_mask(ORDER, item["edge_mask"])
    ]
    return {
        "determinant": str(item["determinant"]),
        "divisible_by_2_22": item["square_root"] % (1 << 22) == 0,
        "edge_count": len(edges),
        "edges": edges,
        "gm_class_index": class_index,
        "gm_original_snapshot_isomorphic": bool(
            item["original_snapshot_hit_indices"]
        ),
        "gm_witness": item["witness"],
        "positive_definite": True,
        "qualified": True,
        "square_root": str(item["square_root"]),
    }


def self_test(enumerator: Path) -> None:
    # A small non-isomorphic GM pair.  X={0,1,2,3} is empty and every
    # outside vertex has exactly two neighbors in X.
    source_pairs = [
        (0, 4),
        (0, 5),
        (0, 6),
        (0, 7),
        (1, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ]
    expected_pairs = [
        (1, 6),
        (1, 7),
        (2, 4),
        (2, 5),
        (2, 7),
        (3, 4),
        (3, 5),
        (3, 6),
    ]
    source_mask = edge_mask_from_pairs(8, source_pairs)
    expected_mask = edge_mask_from_pairs(8, expected_pairs)
    if switch_edge_mask(8, source_mask, 0xF) != expected_mask:
        raise AssertionError("small GM switch reconstruction failed")
    source_graph = graph_from_mask(8, source_mask)
    mate_graph = graph_from_mask(8, expected_mask)
    if exact_isomorphic(source_graph, mate_graph):
        raise AssertionError("small GM example unexpectedly isomorphic")
    source_charpoly = nx.to_numpy_array(source_graph, dtype=int)
    mate_charpoly = nx.to_numpy_array(mate_graph, dtype=int)
    # Integer traces through order 8 determine the characteristic polynomial
    # via Newton identities; checking all traces avoids floating point.
    source_power = source_charpoly.copy()
    mate_power = mate_charpoly.copy()
    for _ in range(1, 9):
        if int(source_power.trace()) != int(mate_power.trace()):
            raise AssertionError("small GM example is not exactly cospectral")
        source_power = source_power @ source_charpoly
        mate_power = mate_power @ mate_charpoly
    completed = subprocess.run(
        [
            str(enumerator),
            "--order",
            "8",
            "--edge-mask",
            f"{source_mask:07x}",
        ],
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    if completed.returncode != 0 or completed.stderr:
        raise AssertionError("small enumerator self-test invocation failed")
    records, _ = parse_enumerator_output(completed.stdout, 8, source_mask)
    selected = [
        record
        for record in records
        if int(record["switch_set_hex"], 16) == 0xF
    ]
    if len(selected) != 1 or int(selected[0]["mate_edge_mask_hex"], 16) != expected_mask:
        raise AssertionError("enumerator missed known small GM switch")

    # Every two-vertex switch is the graph obtained by transposing those
    # vertices, so it is an exact isomorphic/no-new-class calibration.
    pair_mate_mask = switch_edge_mask(8, source_mask, 0x3)
    if not exact_isomorphic(source_graph, graph_from_mask(8, pair_mate_mask)):
        raise AssertionError("two-vertex isomorphic calibration failed")

    # Malformed graph/hit/set inputs must fail closed.
    malformed_hit = {
        "edges": [[1, 2], [1, 2]],
        "edge_count": 2,
    }
    try:
        parse_hit_edges(malformed_hit)
    except InputError:
        pass
    else:
        raise AssertionError("duplicate source edge was accepted")
    try:
        validate_switch_set(8, source_mask, 0)
    except InputError:
        pass
    else:
        raise AssertionError("empty switching set was accepted")
    malformed = subprocess.run(
        [str(enumerator), "--order", "1", "--edge-mask", "0"],
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    if malformed.returncode == 0:
        raise AssertionError("enumerator accepted malformed order")

    # Exact-isomorphism dedupe calibration, including a relabeling.
    relabeled = nx.relabel_nodes(source_graph, {index: 7 - index for index in range(8)})
    if graph_bucket_key(source_graph) != graph_bucket_key(relabeled):
        raise AssertionError("isomorphic bucket keys disagree")
    if not exact_isomorphic(source_graph, relabeled):
        raise AssertionError("exact isomorphism rejected a relabeling")


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    snapshot, survivors, _, snapshot_bytes, hasse_bytes = validate_bound_sources(
        arguments.snapshot,
        arguments.hasse_report,
        arguments.expected_snapshot_sha256,
        arguments.expected_hasse_sha256,
    )
    source_sha256 = sha256_bytes(snapshot_bytes)
    hasse_sha256 = sha256_bytes(hasse_bytes)
    enumerator_sha256 = sha256_path(arguments.enumerator)
    source_code_sha256 = sha256_path(arguments.enumerator_source)
    if enumerator_sha256 != arguments.expected_enumerator_sha256:
        raise InputError("enumerator binary SHA-256 mismatch")
    if source_code_sha256 != arguments.expected_enumerator_source_sha256:
        raise InputError("enumerator source SHA-256 mismatch")
    self_test(arguments.enumerator)

    arguments.work_dir.mkdir(parents=True, exist_ok=True)
    checkpoints: list[dict[str, Any]] = []
    for survivor in survivors:
        checkpoint_path = (
            arguments.work_dir
            / f"source-hit-{survivor['hit_index']:03d}.checkpoint.json"
        )
        if checkpoint_path.exists():
            checkpoint, _ = load_json_bytes(checkpoint_path)
            validate_checkpoint(
                checkpoint,
                survivor,
                enumerator_sha256,
                source_sha256,
                hasse_sha256,
            )
        else:
            checkpoint = checkpoint_for_source(
                survivor,
                arguments.enumerator,
                enumerator_sha256,
                source_sha256,
                hasse_sha256,
            )
            atomic_write_json(checkpoint_path, checkpoint)
        checkpoints.append(checkpoint)

    survivors_by_index = {item["hit_index"]: item for item in survivors}
    original_representatives, original_buckets = build_original_registry(snapshot)
    classes, exact_comparisons = deduplicate_mates(
        checkpoints,
        survivors_by_index,
        original_representatives,
        original_buckets,
    )
    for item in classes:
        item["audit"] = audit_retained_class(item, survivors_by_index)

    route_hits = [route_hit(item, index) for index, item in enumerate(classes)]
    route_snapshot = {
        "schema_version": 1,
        "challenge_id": CHALLENGE_ID,
        "engine": ROUTE_ENGINE,
        "normalization": NORMALIZATION,
        "complete": True,
        "termination": "exhausted",
        "frontier_root": str(FRONTIER_ROOT),
        "parameters": {
            "max_stored_hits": len(route_hits),
            "source_hasse_survivor_count": len(survivors),
            "switching_scope": "all-proper-even-nontrivial-one-step-gm-sets",
        },
        "statistics": {
            "exact_squares": len(route_hits),
            "qualified_survivors": len(route_hits),
            "unrecorded_square_observations": 0,
            "source_hits_completed": len(checkpoints),
            "masks_examined": sum(item["masks_examined"] for item in checkpoints),
            "valid_switching_sets": sum(
                item["valid_switching_sets"] for item in checkpoints
            ),
        },
        "source": {
            "snapshot_path": str(arguments.snapshot),
            "snapshot_sha256": source_sha256,
            "hasse_report_path": str(arguments.hasse_report),
            "hasse_report_sha256": hasse_sha256,
            "enumerator_path": str(arguments.enumerator),
            "enumerator_sha256": enumerator_sha256,
            "enumerator_source_path": str(arguments.enumerator_source),
            "enumerator_source_sha256": source_code_sha256,
        },
        "claim_boundary": (
            "Each retained Gram class is an exact one-step GM mate and preserves "
            "the source determinant and positive-definiteness. Qualification is "
            "a necessary condition only; this snapshot constructs no sign factor."
        ),
        "hits": route_hits,
    }
    atomic_write_json(arguments.route_snapshot, route_snapshot)
    route_sha256 = sha256_path(arguments.route_snapshot)

    report_classes: list[dict[str, Any]] = []
    for class_index, item in enumerate(classes):
        report_classes.append(
            {
                "class_index": class_index,
                "representative_edge_mask_hex": item["edge_mask_hex"],
                "square_root": str(item["square_root"]),
                "determinant": str(item["determinant"]),
                "witness": item["witness"],
                "source_hit_indices": sorted(item["source_hit_indices"]),
                "labeled_mate_multiplicity": item["labeled_mate_multiplicity"],
                "original_snapshot_hit_indices": item[
                    "original_snapshot_hit_indices"
                ],
                "isomorphic_to_original_snapshot": bool(
                    item["original_snapshot_hit_indices"]
                ),
                "audit": item["audit"],
            }
        )
    per_source = [
        {
            "source_hit_index": item["source_hit_index"],
            "source_hit_identity_sha256": item["source_hit_identity_sha256"],
            "masks_examined": item["masks_examined"],
            "valid_switching_sets": item["valid_switching_sets"],
            "unique_labeled_mates": item["unique_labeled_mates"],
            "checkpoint": str(
                arguments.work_dir
                / f"source-hit-{item['source_hit_index']:03d}.checkpoint.json"
            ),
        }
        for item in checkpoints
    ]
    return {
        "schema_version": 1,
        "engine": "gram-gm-switch-audit",
        "method": "complete-one-step-godsil-mckay-switching",
        "complete": True,
        "source": route_snapshot["source"],
        "scope": {
            "source_hasse_survivors": len(survivors),
            "definition": (
                "Every proper even X with regular induced subgraph and outside "
                "neighbor counts in {0, |X|/2, |X|}; no-op sets without a "
                "|X|/2 outside vertex are excluded."
            ),
            "per_source_masks_examined": (1 << (ORDER - 1)) - 1,
        },
        "summary": {
            "source_hits_completed": len(checkpoints),
            "masks_examined": sum(item["masks_examined"] for item in checkpoints),
            "valid_switching_sets": sum(
                item["valid_switching_sets"] for item in checkpoints
            ),
            "unique_labeled_mates_across_sources_before_isomorphism": sum(
                item["unique_labeled_mates"] for item in checkpoints
            ),
            "generated_isomorphism_classes": len(classes),
            "generated_classes_isomorphic_to_original_snapshot": sum(
                bool(item["original_snapshot_hit_indices"]) for item in classes
            ),
            "novel_isomorphism_classes_relative_original_snapshot": sum(
                not item["original_snapshot_hit_indices"] for item in classes
            ),
            "original_snapshot_isomorphism_classes": len(
                original_representatives
            ),
            "exact_isomorphism_comparisons": exact_comparisons,
            "exact_scaled_conjugacy_checks": len(classes),
            "exact_determinant_checks": len(classes),
            "exact_positive_definite_checks": len(classes),
        },
        "route_snapshot": {
            "path": str(arguments.route_snapshot),
            "sha256": route_sha256,
            "hit_count": len(route_hits),
            "engine": ROUTE_ENGINE,
            "normalization": NORMALIZATION,
        },
        "per_source": per_source,
        "classes": report_classes,
        "claim_boundary": (
            "Completeness is for nontrivial one-step Godsil--McKay switching "
            "sets of the 50 source hits selected by the pinned Hasse report. "
            "Graph deduplication is exact up to ordinary unlabeled isomorphism. "
            "The experiment does not cover iterated switching, other source "
            "snapshots, or switching operations outside the GM definition. "
            "A retained Gram is not a {-1,+1} factor."
        ),
    }


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument(
        "--enumerator",
        type=Path,
        default=Path("build/research/gram_gm_switch_enum"),
    )
    parser.add_argument(
        "--enumerator-source",
        type=Path,
        default=Path("research/gram_gm_switch_enum.cpp"),
    )
    parser.add_argument("--snapshot", type=Path)
    parser.add_argument("--hasse-report", type=Path)
    parser.add_argument("--expected-snapshot-sha256")
    parser.add_argument("--expected-hasse-sha256")
    parser.add_argument("--expected-enumerator-sha256")
    parser.add_argument("--expected-enumerator-source-sha256")
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--route-snapshot", type=Path)
    parser.add_argument("--output", type=Path)
    return parser


def main() -> int:
    arguments = make_parser().parse_args()
    try:
        if arguments.self_test:
            if not arguments.enumerator.is_file():
                raise InputError(f"enumerator does not exist: {arguments.enumerator}")
            self_test(arguments.enumerator)
            print("gram_gm_switch self-test: ok")
            return 0
        required = (
            "snapshot",
            "hasse_report",
            "expected_snapshot_sha256",
            "expected_hasse_sha256",
            "expected_enumerator_sha256",
            "expected_enumerator_source_sha256",
            "work_dir",
            "route_snapshot",
            "output",
        )
        missing = [field for field in required if getattr(arguments, field) is None]
        if missing:
            raise InputError(
                "missing required run arguments: " + ", ".join(missing)
            )
        assert arguments.output is not None
        assert arguments.route_snapshot is not None
        if arguments.output.resolve() == arguments.route_snapshot.resolve():
            raise InputError("--output and --route-snapshot must differ")
        report = run(arguments)
        atomic_write_json(arguments.output, report)
        print(
            json.dumps(
                {
                    "status": "complete",
                    "report": str(arguments.output),
                    "report_sha256": sha256_path(arguments.output),
                    "route_snapshot": str(arguments.route_snapshot),
                    "route_snapshot_sha256": sha256_path(arguments.route_snapshot),
                    "summary": report["summary"],
                },
                sort_keys=True,
            )
        )
        return 0
    except (
        ArithmeticError,
        InputError,
        OSError,
        subprocess.SubprocessError,
    ) as error:
        print(f"gram_gm_switch: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
