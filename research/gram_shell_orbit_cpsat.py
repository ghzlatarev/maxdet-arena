#!/usr/bin/env python3
"""Solve one symmetry-fixed factor slice of an exported order-23 Gram shell.

The shell contains normalized sign columns s with s[0] = +1.  Since the
target Gram is nonsingular, no normalized column can repeat, so factorization
is the Boolean problem

    sum_s x_s s s^T = G,  x_s in {0, 1}.

For the published Gram, Aut(G) splits the 1,382 shell columns into orbits of
sizes 6, 432, 432, and 512.  This tool fixes the selected three members of the
size-six orbit, which removes the first and cheapest source of symmetry.
OR-Tools proposes a factor; exact integer reconstruction and determinant
checks remain authoritative.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tempfile
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

import ortools
from ortools.sat.python import cp_model
import pynauty

from gram_hasse import ORDER, InputError, bareiss_determinant


FRONTIER = 2_779_447_296_000_000
PUBLISHED_SMALL_ORBIT = (
    8161,
    26521,
    8356323,
    8357403,
    8363525,
    8380805,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise InputError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise InputError(f"{path} must contain a JSON object")
    return value


def read_matrix(path: Path) -> list[list[int]]:
    try:
        matrix = [
            [int(token) for token in line.split()]
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    except (OSError, ValueError) as error:
        raise InputError(f"cannot read matrix {path}: {error}") from error
    if len(matrix) != ORDER or any(
        len(row) != ORDER or any(value not in (-1, 1) for value in row)
        for row in matrix
    ):
        raise InputError(f"{path} must contain exactly 23x23 signs")
    return matrix


def gram(matrix: list[list[int]]) -> list[list[int]]:
    return [
        [
            sum(
                matrix[row][column] * matrix[other][column]
                for column in range(ORDER)
            )
            for other in range(ORDER)
        ]
        for row in range(ORDER)
    ]


def normalized_column_masks(matrix: list[list[int]]) -> list[int]:
    result = []
    for column in range(ORDER):
        sign = matrix[0][column]
        result.append(
            sum(
                (matrix[row][column] * sign == 1) << row
                for row in range(ORDER)
            )
        )
    return result


def gram_automorphism_generators(
    target_gram: list[list[int]],
) -> list[list[int]]:
    adjacency = {
        row: [
            other
            for other in range(ORDER)
            if other != row and target_gram[row][other] == 3
        ]
        for row in range(ORDER)
    }
    graph = pynauty.Graph(
        number_of_vertices=ORDER,
        directed=False,
        adjacency_dict=adjacency,
    )
    automorphisms = pynauty.autgrp(graph)
    group_order = int(
        round(automorphisms[1] * (10 ** automorphisms[2]))
    )
    if group_order != 442_368 or len(automorphisms[0]) != 15:
        raise InputError(
            "unexpected published-Gram automorphism group "
            f"(generators={len(automorphisms[0])}, order={group_order})"
        )
    return automorphisms[0]


def act_mask(mask: int, permutation: list[int]) -> int:
    old = [1 if mask >> row & 1 else -1 for row in range(ORDER)]
    new = [0] * ORDER
    for row in range(ORDER):
        new[permutation[row]] = old[row]
    if new[0] < 0:
        new = [-value for value in new]
    return sum((value == 1) << row for row, value in enumerate(new))


def shell_orbits(
    masks: list[int], generators: list[list[int]]
) -> list[list[int]]:
    """Return normalized-shell orbits under exact row-Gram automorphisms."""
    mask_set = set(masks)

    unseen = set(masks)
    result = []
    while unseen:
        orbit = {min(unseen)}
        queue = list(orbit)
        while queue:
            mask = queue.pop()
            for generator in generators:
                image = act_mask(mask, generator)
                if image not in mask_set:
                    raise InputError(
                        "Gram automorphism did not preserve the sign shell"
                    )
                if image not in orbit:
                    orbit.add(image)
                    queue.append(image)
        unseen.difference_update(orbit)
        result.append(sorted(orbit))
    result.sort(key=lambda orbit: (len(orbit), orbit))
    if [len(orbit) for orbit in result] != [6, 432, 432, 512]:
        raise InputError(
            "unexpected shell orbit sizes: "
            + ",".join(str(len(orbit)) for orbit in result)
        )
    return result


def pair_orbits(generators: list[list[int]]) -> list[list[tuple[int, int]]]:
    unseen = {
        (row, other)
        for row in range(ORDER)
        for other in range(row + 1, ORDER)
    }
    result = []
    while unseen:
        orbit = {min(unseen)}
        queue = list(orbit)
        while queue:
            row, other = queue.pop()
            for generator in generators:
                image = tuple(
                    sorted((generator[row], generator[other]))
                )
                if image not in orbit:
                    orbit.add(image)
                    queue.append(image)
        unseen.difference_update(orbit)
        result.append(sorted(orbit))
    result.sort(key=lambda orbit: (len(orbit), orbit))
    if len(result) != 16 or sum(map(len, result)) != ORDER * (ORDER - 1) // 2:
        raise InputError(
            f"unexpected unordered-pair orbit partition ({len(result)} orbits)"
        )
    return result


def exact_forced_orbit_counts(
    shell_partition: list[list[int]],
    pair_partition: list[list[tuple[int, int]]],
    target_gram: list[list[int]],
) -> tuple[list[int], dict[str, object]]:
    """Derive shell-orbit counts from exact sums of the Gram equations."""
    equations: list[list[int]] = []
    targets: list[int] = []
    for pair_orbit in pair_partition:
        row = []
        for shell_orbit in shell_partition:
            mask = shell_orbit[0]
            coefficient = 0
            for first, second in pair_orbit:
                coefficient += (
                    1
                    if ((mask >> first) & 1) == ((mask >> second) & 1)
                    else -1
                )
            # The sum is constant on a shell orbit by construction.
            for check_mask in shell_orbit:
                check = sum(
                    1
                    if ((check_mask >> first) & 1)
                    == ((check_mask >> second) & 1)
                    else -1
                    for first, second in pair_orbit
                )
                if check != coefficient:
                    raise InputError(
                        "pair-orbit equation is not constant on shell orbit"
                    )
            row.append(coefficient)
        equations.append(row)
        targets.append(
            sum(target_gram[first][second] for first, second in pair_orbit)
        )
    equations.append([1] * len(shell_partition))
    targets.append(ORDER)

    augmented = [
        [Fraction(value) for value in row] + [Fraction(target)]
        for row, target in zip(equations, targets)
    ]
    variable_count = len(shell_partition)
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
        augmented[pivot_row] = [value / pivot for value in augmented[pivot_row]]
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
        if all(value == 0 for value in row[:variable_count]) and row[-1] != 0:
            raise InputError("shell-orbit incidence equations are inconsistent")
    if pivot_columns != list(range(variable_count)):
        raise InputError(
            "shell-orbit incidence equations do not have a unique solution"
        )
    solution = [Fraction(0) for _ in range(variable_count)]
    for row, column in enumerate(pivot_columns):
        solution[column] = augmented[row][-1]
    if any(value.denominator != 1 or value < 0 for value in solution):
        raise InputError("forced shell-orbit counts are not nonnegative integers")
    counts = [int(value) for value in solution]
    if counts != [3, 6, 6, 8]:
        raise InputError(f"unexpected forced shell-orbit counts: {counts}")
    return counts, {
        "equation_count": len(equations),
        "pair_orbit_count": len(pair_partition),
        "rank": len(pivot_columns),
        "shell_orbit_sizes": [len(orbit) for orbit in shell_partition],
        "unique_counts": counts,
    }


def align_support_to_small_triple(
    support: set[int],
    selected_small: set[int],
    generators: list[list[int]],
) -> tuple[set[int], bool]:
    """Move a factor support to the requested equivalent small-orbit triple."""
    initial_triple = frozenset(support.intersection(PUBLISHED_SMALL_ORBIT))
    target = frozenset(selected_small)
    if initial_triple == target:
        return support, False
    queue = [(initial_triple, support)]
    seen = {initial_triple}
    while queue:
        triple, current_support = queue.pop(0)
        for generator in generators:
            image_support = {
                act_mask(mask, generator) for mask in current_support
            }
            image_triple = frozenset(
                image_support.intersection(PUBLISHED_SMALL_ORBIT)
            )
            if image_triple == target:
                return image_support, True
            if image_triple not in seen:
                seen.add(image_triple)
                queue.append((image_triple, image_support))
    raise InputError(
        "factor's small-orbit triple is not equivalent to --small-triple"
    )


def reconstruct_factor(masks: list[int], selected: list[int]) -> list[list[int]]:
    columns = [
        [1 if masks[index] >> row & 1 else -1 for row in range(ORDER)]
        for index in selected
    ]
    if len(columns) != ORDER:
        raise ArithmeticError(
            f"candidate uses {len(columns)} columns instead of {ORDER}"
        )
    return [
        [columns[column][row] for column in range(ORDER)]
        for row in range(ORDER)
    ]


def matrix_bytes(matrix: list[list[int]]) -> str:
    return "".join(" ".join(map(str, row)) + "\n" for row in matrix)


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


def parse_triple(text: str) -> tuple[int, int, int]:
    try:
        values = tuple(int(value) for value in text.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "--small-triple must be three comma-separated indices"
        ) from error
    if len(values) != 3 or len(set(values)) != 3 or any(
        value < 0 or value >= len(PUBLISHED_SMALL_ORBIT) for value in values
    ):
        raise argparse.ArgumentTypeError(
            "--small-triple must select three distinct indices in [0,5]"
        )
    return tuple(sorted(values))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shell-report", required=True, type=Path)
    parser.add_argument("--source-matrix", required=True, type=Path)
    parser.add_argument("--small-triple", required=True, type=parse_triple)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--hint-factor", type=Path)
    parser.add_argument(
        "--exclude-factor",
        action="append",
        default=[],
        type=Path,
        help="exclude one exact shell subset (repeatable)",
    )
    parser.add_argument(
        "--max-hint-overlap",
        type=int,
        help="require at most this many columns from --hint-factor",
    )
    parser.add_argument(
        "--randomize-search",
        action="store_true",
        help="enable seeded CP-SAT randomized search",
    )
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument("--time-limit", type=float, default=60.0)
    parser.add_argument("--workers", type=int, default=1)
    arguments = parser.parse_args()

    if (
        arguments.time_limit <= 0
        or not math.isfinite(arguments.time_limit)
        or arguments.workers <= 0
    ):
        parser.error("--time-limit and --workers must be positive")
    if arguments.output.resolve() == arguments.metadata.resolve():
        parser.error("--output and --metadata must be distinct")
    if arguments.max_hint_overlap is not None and (
        arguments.hint_factor is None
        or arguments.max_hint_overlap < 0
        or arguments.max_hint_overlap >= ORDER
    ):
        parser.error("--max-hint-overlap requires --hint-factor and must be 0..22")
    if arguments.output.exists() or arguments.metadata.exists():
        parser.error("refusing to overwrite output or metadata")

    report = read_json(arguments.shell_report)
    results = report.get("results")
    if (
        not isinstance(results, list)
        or len(results) != 1
        or not isinstance(results[0], dict)
    ):
        raise InputError("shell report must contain exactly one result")
    result = results[0]
    if result.get("rejected") is not False:
        raise InputError("selected shell result is not an accepted shell")
    masks_value = result.get("shell_sign_masks")
    if not isinstance(masks_value, list) or any(
        type(mask) is not int
        or mask < 0
        or mask >= 1 << ORDER
        or mask & 1 == 0
        for mask in masks_value
    ):
        raise InputError("shell report has invalid shell_sign_masks")
    masks = sorted(masks_value)
    if (
        len(masks) != 1382
        or result.get("shell_size") != len(masks)
        or len(set(masks)) != len(masks)
    ):
        raise InputError("expected 1,382 unique published-Gram shell masks")
    if tuple(mask for mask in masks if mask in PUBLISHED_SMALL_ORBIT) != (
        PUBLISHED_SMALL_ORBIT
    ):
        raise InputError("published size-six shell orbit is incomplete")

    source = read_matrix(arguments.source_matrix)
    if result.get("source_sha256") != sha256(arguments.source_matrix):
        raise InputError(
            "shell report source SHA-256 does not bind to --source-matrix"
        )
    target_gram = gram(source)
    gram_determinant = bareiss_determinant(target_gram)
    gram_root = math.isqrt(gram_determinant)
    if (
        gram_root * gram_root != gram_determinant
        or gram_root != FRONTIER
        or str(gram_determinant) != result.get("determinant")
    ):
        raise InputError("source and shell report do not bind to the frontier Gram")

    model = cp_model.CpModel()
    variables = [
        model.new_bool_var(f"shell_{index}") for index in range(len(masks))
    ]
    model.add(sum(variables) == ORDER)
    for row in range(ORDER):
        for other in range(row + 1, ORDER):
            same = [
                variables[index]
                for index, mask in enumerate(masks)
                if ((mask >> row) & 1) == ((mask >> other) & 1)
            ]
            model.add(
                sum(same) == (ORDER + target_gram[row][other]) // 2
            )

    generators = gram_automorphism_generators(target_gram)
    orbits = shell_orbits(masks, generators)
    forced_orbit_counts, orbit_count_certificate = exact_forced_orbit_counts(
        orbits, pair_orbits(generators), target_gram
    )
    mask_to_index = {mask: index for index, mask in enumerate(masks)}
    for orbit, count in zip(orbits, forced_orbit_counts):
        model.add(
            sum(variables[mask_to_index[mask]] for mask in orbit)
            == count
        )

    selected_small = {
        PUBLISHED_SMALL_ORBIT[index] for index in arguments.small_triple
    }
    for mask in PUBLISHED_SMALL_ORBIT:
        model.add(variables[mask_to_index[mask]] == (mask in selected_small))

    hint_support: set[int] = set()
    hint_was_aligned = False
    if arguments.hint_factor is not None:
        hint = read_matrix(arguments.hint_factor)
        if gram(hint) != target_gram:
            raise InputError("--hint-factor does not have the exact target Gram")
        hint_support = set(normalized_column_masks(hint))
        if len(hint_support) != ORDER or not hint_support.issubset(mask_to_index):
            raise InputError("--hint-factor columns do not form a shell subset")
        hint_support, hint_was_aligned = align_support_to_small_triple(
            hint_support, selected_small, generators
        )
        for index, variable in enumerate(variables):
            model.add_hint(variable, masks[index] in hint_support)
        if arguments.max_hint_overlap is not None:
            model.add(
                sum(
                    variables[mask_to_index[mask]]
                    for mask in hint_support
                )
                <= arguments.max_hint_overlap
            )

    excluded_factor_hashes = []
    excluded_factors_aligned = 0
    for path in arguments.exclude_factor:
        excluded = read_matrix(path)
        if gram(excluded) != target_gram:
            raise InputError(f"--exclude-factor does not have target Gram: {path}")
        support = set(normalized_column_masks(excluded))
        if len(support) != ORDER or not support.issubset(mask_to_index):
            raise InputError(f"--exclude-factor is not a shell subset: {path}")
        support, was_aligned = align_support_to_small_triple(
            support, selected_small, generators
        )
        excluded_factors_aligned += was_aligned
        model.add(
            sum(variables[mask_to_index[mask]] for mask in support) <= ORDER - 1
        )
        excluded_factor_hashes.append(sha256(path))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = arguments.time_limit
    solver.parameters.num_search_workers = arguments.workers
    solver.parameters.random_seed = arguments.seed
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 2
    solver.parameters.randomize_search = arguments.randomize_search

    started = time.monotonic()
    status = solver.solve(model)
    elapsed = time.monotonic() - started
    status_name = solver.status_name(status)
    factor: list[list[int]] | None = None
    selected: list[int] = []
    verified = False
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        selected = [
            index
            for index, variable in enumerate(variables)
            if solver.value(variable)
        ]
        factor = reconstruct_factor(masks, selected)
        score = abs(bareiss_determinant(factor))
        verified = gram(factor) == target_gram and score == gram_root
    output_contents = matrix_bytes(factor) if verified and factor is not None else None
    output_sha256 = (
        hashlib.sha256(output_contents.encode("ascii")).hexdigest()
        if output_contents is not None
        else None
    )

    metadata = {
        "claim_boundary": (
            "factor_found=true means exact Gram reconstruction and Bareiss "
            "determinant checks passed. UNKNOWN is a bounded search result, "
            "not an infeasibility proof."
        ),
        "engine": "gram-shell-orbit-cpsat",
        "engine_source": str(Path(__file__).resolve()),
        "engine_source_sha256": sha256(Path(__file__).resolve()),
        "dependencies": {
            "ortools": getattr(ortools, "__version__", "unknown"),
            "pynauty": getattr(pynauty, "__version__", "unknown"),
        },
        "factor_found": verified,
        "factor_score": str(gram_root) if verified else None,
        "fixed_small_orbit_masks": sorted(selected_small),
        "fixed_small_orbit_triple": list(arguments.small_triple),
        "gram_determinant": str(gram_determinant),
        "orbit_count_exact_certificate": orbit_count_certificate,
        "hint_factor": (
            str(arguments.hint_factor)
            if arguments.hint_factor is not None
            else None
        ),
        "hint_factor_sha256": (
            sha256(arguments.hint_factor)
            if arguments.hint_factor is not None
            else None
        ),
        "hint_factor_was_automorphism_aligned": hint_was_aligned,
        "max_hint_overlap": arguments.max_hint_overlap,
        "excluded_factor_count": len(arguments.exclude_factor),
        "excluded_factor_sha256": excluded_factor_hashes,
        "excluded_factors_automorphism_aligned": excluded_factors_aligned,
        "randomize_search": arguments.randomize_search,
        "runtime_seconds": elapsed,
        "requested_limits": {
            "time_limit_seconds": arguments.time_limit,
        },
        "seed": arguments.seed,
        "selected_shell_masks": (
            [masks[index] for index in selected] if verified else []
        ),
        "output": str(arguments.output),
        "output_sha256": output_sha256,
        "shell_report": str(arguments.shell_report),
        "shell_report_sha256": sha256(arguments.shell_report),
        "shell_size": len(masks),
        "solver": {
            "best_objective_bound": solver.best_objective_bound,
            "branches": solver.num_branches,
            "conflicts": solver.num_conflicts,
            "num_workers": arguments.workers,
            "status": status_name,
            "wall_time": solver.wall_time,
        },
        "source_matrix": str(arguments.source_matrix),
        "source_matrix_sha256": sha256(arguments.source_matrix),
    }
    if not verified or factor is None:
        atomic_write_text(
            arguments.metadata,
            json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        )
        print(
            f"no verified factor: triple={arguments.small_triple} "
            f"status={status_name} elapsed={elapsed:.3f}s"
        )
        return 3

    assert output_contents is not None
    atomic_write_text(arguments.output, output_contents)
    atomic_write_text(
        arguments.metadata,
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
    )
    print(
        f"verified factor triple={arguments.small_triple} "
        f"|det|={gram_root} elapsed={elapsed:.3f}s"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (InputError, ArithmeticError, OSError) as error:
        raise SystemExit(f"gram_shell_orbit_cpsat: {error}") from error
