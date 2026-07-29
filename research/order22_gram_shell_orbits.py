#!/usr/bin/env python3
"""Audit signed-Gram symmetry orbits of an exact order-22 column shell.

The signed row automorphism group of a target Gram acts on normalized sign
columns.  This tool computes the complete shell-orbit partition and exact
linear constraints on the number of selected columns in each orbit.

The count constraints come from signed unordered-pair orbit sums.  They are
necessary for every factor of the target Gram.  They need not determine every
orbit count, so observed count vectors from known factors are reported
separately and are not claimed to be forced.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time
from typing import Any

import pynauty

from order22_align_gram_factors import (
    ORDER,
    gram,
    normalized_column_masks,
    read_matrix,
    signed_gram_graph,
)
from order22_gram_factor_cpsat import parse_shell, read_factor_list, read_json


PINNED_PYNAUTY = "2.8.8.1"


class AuditError(ValueError):
    """An invalid input or failed exact consistency check."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write(path: Path, contents: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise AuditError(f"refusing to overwrite {path}")
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(contents)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def fraction_json(value: Fraction) -> int | str:
    if value.denominator == 1:
        return value.numerator
    return f"{value.numerator}/{value.denominator}"


def signed_automorphism_generators(
    target_gram: list[list[int]],
) -> tuple[list[tuple[list[int], list[int]]], int]:
    """Return generators as source-to-target row maps and source-row signs."""

    graph, _ = signed_gram_graph(target_gram)
    raw = pynauty.autgrp(graph)
    group_order = int(round(raw[1] * (10 ** raw[2])))
    generators: list[tuple[list[int], list[int]]] = []
    for generator_index, automorphism in enumerate(raw[0]):
        permutation = []
        signs = []
        for source_row in range(ORDER):
            mapped_fiber = automorphism[2 * ORDER + source_row]
            target_row = mapped_fiber - 2 * ORDER
            if not 0 <= target_row < ORDER:
                raise AuditError(
                    f"generator {generator_index} maps a fiber outside its color"
                )
            positive = automorphism[2 * source_row]
            negative = automorphism[2 * source_row + 1]
            if (
                positive // 2 != target_row
                or negative // 2 != target_row
                or {positive % 2, negative % 2} != {0, 1}
            ):
                raise AuditError(
                    f"generator {generator_index} breaks a signed row fiber"
                )
            permutation.append(target_row)
            signs.append(1 if positive % 2 == 0 else -1)
        if set(permutation) != set(range(ORDER)):
            raise AuditError(
                f"generator {generator_index} row map is not a permutation"
            )
        for first in range(ORDER):
            for second in range(ORDER):
                if (
                    target_gram[permutation[first]][permutation[second]]
                    != signs[first] * signs[second] * target_gram[first][second]
                ):
                    raise AuditError(
                        f"generator {generator_index} failed signed Gram equality"
                    )
        generators.append((permutation, signs))
    return generators, group_order


def act_mask(
    mask: int, generator: tuple[list[int], list[int]]
) -> int:
    permutation, signs = generator
    source = [1 if mask >> row & 1 else -1 for row in range(ORDER)]
    target = [0] * ORDER
    for row in range(ORDER):
        target[permutation[row]] = signs[row] * source[row]
    if target[0] < 0:
        target = [-value for value in target]
    return sum((value == 1) << row for row, value in enumerate(target))


def shell_orbits(
    masks: list[int],
    generators: list[tuple[list[int], list[int]]],
) -> list[list[int]]:
    mask_set = set(masks)
    unseen = set(masks)
    result: list[list[int]] = []
    while unseen:
        orbit = {min(unseen)}
        queue = list(orbit)
        while queue:
            mask = queue.pop()
            for generator in generators:
                image = act_mask(mask, generator)
                if image not in mask_set:
                    raise AuditError(
                        "a signed Gram automorphism left the exact column shell"
                    )
                if image not in orbit:
                    orbit.add(image)
                    queue.append(image)
        unseen.difference_update(orbit)
        result.append(sorted(orbit))
    result.sort(key=lambda orbit: (len(orbit), orbit))
    if sum(map(len, result)) != len(masks):
        raise AuditError("shell orbit partition has the wrong cardinality")
    return result


def signed_pair_orbits(
    generators: list[tuple[list[int], list[int]]],
) -> list[dict[str, Any]]:
    """Partition row pairs and track the sign character of each image.

    If a pair is reachable with both signs, its Reynolds orbit sum vanishes
    and supplies no orbit-count equation.  Such an orbit is marked
    ``sign_inconsistent`` and omitted from the invariant system.
    """

    unseen = {
        (first, second)
        for first in range(ORDER)
        for second in range(first + 1, ORDER)
    }
    result: list[dict[str, Any]] = []
    while unseen:
        seed = min(unseen)
        weights = {seed: 1}
        queue = [seed]
        inconsistent = False
        while queue:
            first, second = queue.pop()
            weight = weights[(first, second)]
            for permutation, signs in generators:
                image = tuple(
                    sorted((permutation[first], permutation[second]))
                )
                image_weight = weight * signs[first] * signs[second]
                if image in weights:
                    if weights[image] != image_weight:
                        inconsistent = True
                else:
                    weights[image] = image_weight
                    queue.append(image)
        unseen.difference_update(weights)
        result.append(
            {
                "weights": weights,
                "sign_inconsistent": inconsistent,
            }
        )
    result.sort(
        key=lambda record: (
            len(record["weights"]),
            sorted(record["weights"]),
        )
    )
    if sum(len(record["weights"]) for record in result) != ORDER * (ORDER - 1) // 2:
        raise AuditError("signed pair orbits do not partition unordered pairs")
    return result


def invariant_system(
    target_gram: list[list[int]],
    partition: list[list[int]],
    pair_partition: list[dict[str, Any]],
) -> tuple[list[list[int]], list[int], list[dict[str, Any]]]:
    equations: list[list[int]] = []
    targets: list[int] = []
    records: list[dict[str, Any]] = []
    for pair_index, record in enumerate(pair_partition):
        if record["sign_inconsistent"]:
            continue
        weights: dict[tuple[int, int], int] = record["weights"]
        coefficients = []
        for orbit in partition:
            values = {
                sum(
                    weight
                    * (
                        1
                        if ((mask >> first) & 1) == ((mask >> second) & 1)
                        else -1
                    )
                    for (first, second), weight in weights.items()
                )
                for mask in orbit
            }
            if len(values) != 1:
                raise AuditError(
                    "a signed pair sum is not constant on a shell orbit"
                )
            coefficients.append(values.pop())
        target = sum(
            weight * target_gram[first][second]
            for (first, second), weight in weights.items()
        )
        equations.append(coefficients)
        targets.append(target)
        records.append(
            {
                "kind": "signed_pair_orbit_sum",
                "pair_orbit_index": pair_index,
                "coefficients": coefficients,
                "rhs": target,
            }
        )
    equations.append([1] * len(partition))
    targets.append(ORDER)
    records.append(
        {
            "kind": "selected_column_count",
            "coefficients": equations[-1],
            "rhs": ORDER,
        }
    )
    return equations, targets, records


def reduced_system(
    equations: list[list[int]], targets: list[int]
) -> tuple[list[list[Fraction]], list[int]]:
    variable_count = len(equations[0])
    augmented = [
        [Fraction(value) for value in row] + [Fraction(target)]
        for row, target in zip(equations, targets)
    ]
    pivot_row = 0
    pivot_columns = []
    for column in range(variable_count):
        selected = next(
            (
                row
                for row in range(pivot_row, len(augmented))
                if augmented[row][column]
            ),
            None,
        )
        if selected is None:
            continue
        augmented[pivot_row], augmented[selected] = (
            augmented[selected],
            augmented[pivot_row],
        )
        pivot = augmented[pivot_row][column]
        augmented[pivot_row] = [
            value / pivot for value in augmented[pivot_row]
        ]
        for row in range(len(augmented)):
            if row == pivot_row or not augmented[row][column]:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                value - factor * base
                for value, base in zip(augmented[row], augmented[pivot_row])
            ]
        pivot_columns.append(column)
        pivot_row += 1
    for row in augmented:
        if (
            all(value == 0 for value in row[:variable_count])
            and row[-1] != 0
        ):
            raise AuditError("exact orbit-count equations are inconsistent")
    return augmented, pivot_columns


def forced_counts(
    reduced: list[list[Fraction]], variable_count: int
) -> dict[str, int]:
    result = {}
    for row in reduced:
        nonzero = [
            index for index, value in enumerate(row[:variable_count]) if value
        ]
        if len(nonzero) != 1:
            continue
        index = nonzero[0]
        value = row[-1] / row[index]
        if value.denominator == 1:
            result[str(index)] = int(value)
    return result


def resolve_known_factors(list_path: Path) -> list[Path]:
    return read_factor_list(list_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", required=True)
    parser.add_argument("--shell-report", required=True, type=Path)
    parser.add_argument("--factor", required=True, type=Path)
    parser.add_argument("--known-factor-list", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()

    if getattr(pynauty, "__version__", None) != PINNED_PYNAUTY:
        raise RuntimeError(f"requires pynauty=={PINNED_PYNAUTY}")
    started = time.monotonic()
    shell_path = arguments.shell_report.expanduser().resolve()
    factor_path = arguments.factor.expanduser().resolve()
    list_path = arguments.known_factor_list.expanduser().resolve()
    output_path = arguments.output.expanduser().resolve()

    shell = parse_shell(read_json(shell_path))
    factor = read_matrix(factor_path)
    target_gram = gram(factor)
    generators, group_order = signed_automorphism_generators(target_gram)
    partition = shell_orbits(shell, generators)
    orbit_by_mask = {
        mask: orbit_index
        for orbit_index, orbit in enumerate(partition)
        for mask in orbit
    }
    pair_partition = signed_pair_orbits(generators)
    equations, targets, equation_records = invariant_system(
        target_gram, partition, pair_partition
    )
    reduced, pivots = reduced_system(equations, targets)

    known_records = []
    vector_histogram: dict[tuple[int, ...], int] = {}
    for path in resolve_known_factors(list_path):
        known = read_matrix(path)
        if gram(known) != target_gram:
            raise AuditError(f"known factor has the wrong target Gram: {path}")
        support = normalized_column_masks(known)
        if any(mask not in orbit_by_mask for mask in support):
            raise AuditError(f"known factor leaves the shell: {path}")
        vector = [0] * len(partition)
        for mask in support:
            vector[orbit_by_mask[mask]] += 1
        for coefficients, rhs in zip(equations, targets):
            if sum(a * b for a, b in zip(coefficients, vector)) != rhs:
                raise AuditError(f"known factor violates an invariant: {path}")
        vector_key = tuple(vector)
        vector_histogram[vector_key] = vector_histogram.get(vector_key, 0) + 1
        known_records.append(
            {
                "path": str(path),
                "sha256": sha256(path),
                "orbit_count_vector": vector,
            }
        )

    pair_records = []
    for index, record in enumerate(pair_partition):
        weights: dict[tuple[int, int], int] = record["weights"]
        pair_records.append(
            {
                "index": index,
                "size": len(weights),
                "sign_inconsistent": record["sign_inconsistent"],
                "signed_pairs": [
                    [first, second, weight]
                    for (first, second), weight in sorted(weights.items())
                ],
            }
        )

    result = {
        "engine": "order22-signed-gram-shell-orbits-v1",
        "order": ORDER,
        "name": arguments.name,
        "claim": (
            "The complete exact shell is partitioned under the full "
            "pynauty-generated signed-Gram group. Reported invariant equations "
            "and forced counts are necessary for every factor of this target Gram."
        ),
        "claim_boundary": (
            "Known-factor orbit vectors are observations, not exhaustive or "
            "forced unless their entries also appear in forced_orbit_counts."
        ),
        "inputs": {
            "shell_report": {
                "path": str(shell_path),
                "sha256": sha256(shell_path),
                "shell_size": len(shell),
            },
            "factor": {
                "path": str(factor_path),
                "sha256": sha256(factor_path),
            },
            "known_factor_list": {
                "path": str(list_path),
                "sha256": sha256(list_path),
                "factor_count": len(known_records),
            },
        },
        "group": {
            "signed_double_cover_order": group_order,
            "generator_count": len(generators),
            "generators": [
                {
                    "source_to_target_row_permutation": permutation,
                    "source_row_signs": signs,
                }
                for permutation, signs in generators
            ],
        },
        "shell_partition": {
            "orbit_count": len(partition),
            "orbit_sizes": [len(orbit) for orbit in partition],
            "orbits": [
                {
                    "index": index,
                    "size": len(orbit),
                    "minimum_mask": orbit[0],
                    "masks": orbit,
                }
                for index, orbit in enumerate(partition)
            ],
        },
        "signed_pair_partition": {
            "orbit_count": len(pair_partition),
            "consistent_orbit_count": sum(
                not record["sign_inconsistent"] for record in pair_partition
            ),
            "orbits": pair_records,
        },
        "orbit_count_system": {
            "variable_count": len(partition),
            "variable_convention": (
                "k_i is the number of selected masks in shell orbit i"
            ),
            "equation_count": len(equations),
            "rank": len(pivots),
            "pivot_columns": pivots,
            "equations": equation_records,
            "rref": [
                [fraction_json(value) for value in row]
                for row in reduced
                if any(row)
            ],
            "forced_orbit_counts": forced_counts(reduced, len(partition)),
        },
        "known_factors": known_records,
        "known_orbit_count_vectors": [
            {"vector": list(vector), "factor_count": count}
            for vector, count in sorted(vector_histogram.items())
        ],
        "dependencies": {
            "pynauty": getattr(pynauty, "__version__", "unknown"),
        },
        "elapsed_seconds": time.monotonic() - started,
    }
    contents = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode(
        "utf-8"
    )
    atomic_write(output_path, contents)
    print(
        f"name={arguments.name} shell={len(shell)} "
        f"orbits={len(partition)} rank={len(pivots)} "
        f"known_vectors={len(vector_histogram)} "
        f"elapsed={result['elapsed_seconds']:.3f}s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
