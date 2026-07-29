#!/usr/bin/env python3
"""Sample an exact sign factor from an exported order-22 Gram shell.

For a nonsingular sign factor R with row Gram G, normalize every column so
its first sign is +1.  Every normalized column lies in the exported shell,
no two columns repeat, and factorization is the Boolean system

    sum_s x_s s s^T = G,  x_s in {0, 1}.

OR-Tools proposes one solution.  Exact Gram equality and a fraction-free
integer determinant check are authoritative before any matrix is written.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import tempfile
import time
from pathlib import Path
from typing import Any

import ortools
from ortools.sat.python import cp_model


ORDER = 22


class InputError(ValueError):
    """An invalid or inconsistent research input."""


def sha256_bytes(contents: bytes) -> str:
    return hashlib.sha256(contents).hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise InputError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise InputError(f"{path} must contain a JSON object")
    return value


def read_matrix(path: Path) -> list[list[int]]:
    try:
        rows = [
            [int(token) for token in line.split()]
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    except (OSError, UnicodeDecodeError, ValueError) as error:
        raise InputError(f"cannot read matrix {path}: {error}") from error
    if len(rows) != ORDER or any(
        len(row) != ORDER or any(value not in (-1, 1) for value in row)
        for row in rows
    ):
        raise InputError(f"{path} must contain exactly {ORDER}x{ORDER} signs")
    return rows


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


def bareiss_determinant(source: list[list[int]]) -> int:
    work = [row[:] for row in source]
    previous = 1
    sign = 1
    for column in range(len(work) - 1):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[pivot], work[column] = work[column], work[pivot]
            sign = -sign
        pivot_value = work[column][column]
        for row in range(column + 1, len(work)):
            for inner in range(column + 1, len(work)):
                numerator = (
                    work[row][inner] * pivot_value
                    - work[row][column] * work[column][inner]
                )
                if column and numerator % previous:
                    raise ArithmeticError("non-exact Bareiss division")
                work[row][inner] = numerator if not column else numerator // previous
            work[row][column] = 0
        previous = pivot_value
    return sign * work[-1][-1]


def normalized_column_masks(matrix: list[list[int]]) -> list[int]:
    masks = []
    for column in range(ORDER):
        column_sign = matrix[0][column]
        masks.append(
            sum(
                (matrix[row][column] * column_sign == 1) << row
                for row in range(ORDER)
            )
        )
    if len(set(masks)) != ORDER:
        raise InputError("nonsingular factor unexpectedly has repeated columns")
    return masks


def parse_shell(report: dict[str, Any]) -> list[int]:
    if (
        report.get("engine") != "order22-gram-shell-v1"
        or report.get("order") != ORDER
        or report.get("complete") is not True
        or report.get("assignments_completed") != 1 << (ORDER - 1)
        or report.get("factor_columns_in_shell") is not True
    ):
        raise InputError("shell report is incomplete or has an unexpected schema")
    raw_masks = report.get("shell_masks")
    shell_size = report.get("shell_size")
    if not isinstance(raw_masks, list) or shell_size != len(raw_masks):
        raise InputError("shell mask count does not match shell_size")
    masks: list[int] = []
    for index, mask in enumerate(raw_masks):
        if (
            type(mask) is not int
            or mask < 0
            or mask >= 1 << ORDER
            or mask & 1 == 0
        ):
            raise InputError(f"invalid normalized shell mask at index {index}")
        masks.append(mask)
    if len(set(masks)) != len(masks):
        raise InputError("shell masks are not unique")
    return masks


def reconstruct_factor(selected_masks: list[int]) -> list[list[int]]:
    if len(selected_masks) != ORDER:
        raise ArithmeticError(f"selected {len(selected_masks)} columns")
    columns = [
        [1 if mask >> row & 1 else -1 for row in range(ORDER)]
        for mask in selected_masks
    ]
    return [
        [columns[column][row] for column in range(ORDER)]
        for row in range(ORDER)
    ]


def matrix_bytes(matrix: list[list[int]]) -> bytes:
    return "".join(" ".join(map(str, row)) + "\n" for row in matrix).encode(
        "ascii"
    )


def atomic_write(path: Path, contents: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise InputError(f"refusing to overwrite {path}")
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


def parse_mask(value: str) -> int:
    try:
        result = int(value, 0)
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"invalid integer mask: {value}") from error
    if result < 0 or result >= 1 << ORDER or result & 1 == 0:
        raise argparse.ArgumentTypeError(
            f"mask must be in [0, 2^{ORDER}) with bit zero set"
        )
    return result


def read_factor_list(path: Path) -> list[Path]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeDecodeError) as error:
        raise InputError(f"cannot read factor list {path}: {error}") from error
    factors = []
    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        candidate = Path(line).expanduser()
        if not candidate.is_absolute():
            candidate = path.parent / candidate
        candidate = candidate.resolve()
        if not candidate.is_file():
            raise InputError(
                f"{path}:{line_number}: missing factor {candidate}"
            )
        factors.append(candidate)
    if not factors:
        raise InputError(f"factor list is empty: {path}")
    return factors


def parse_orbit_counts(value: str) -> list[int]:
    try:
        counts = [int(token) for token in value.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "--orbit-counts must be comma-separated nonnegative integers"
        ) from error
    if not counts or any(count < 0 for count in counts):
        raise argparse.ArgumentTypeError(
            "--orbit-counts must be comma-separated nonnegative integers"
        )
    return counts


def parse_orbit_report(
    report_path: Path,
    shell_path: Path,
    factor_path: Path,
    masks: list[int],
) -> tuple[list[list[int]], dict[int, int], dict[str, Any]]:
    report = read_json(report_path)
    if (
        report.get("engine") != "order22-signed-gram-shell-orbits-v1"
        or report.get("order") != ORDER
    ):
        raise InputError("orbit report has an unexpected schema")
    inputs = report.get("inputs")
    if not isinstance(inputs, dict):
        raise InputError("orbit report is missing input bindings")
    shell_input = inputs.get("shell_report")
    factor_input = inputs.get("factor")
    if (
        not isinstance(shell_input, dict)
        or shell_input.get("sha256") != sha256(shell_path)
        or not isinstance(factor_input, dict)
        or factor_input.get("sha256") != sha256(factor_path)
    ):
        raise InputError("orbit report is bound to different shell/factor bytes")

    partition_record = report.get("shell_partition")
    raw_orbits = (
        partition_record.get("orbits")
        if isinstance(partition_record, dict)
        else None
    )
    if not isinstance(raw_orbits, list) or not raw_orbits:
        raise InputError("orbit report has no shell partition")
    partition: list[list[int]] = []
    for expected_index, raw_orbit in enumerate(raw_orbits):
        if (
            not isinstance(raw_orbit, dict)
            or raw_orbit.get("index") != expected_index
            or not isinstance(raw_orbit.get("masks"), list)
        ):
            raise InputError("orbit report has malformed orbit records")
        orbit = raw_orbit["masks"]
        if (
            raw_orbit.get("size") != len(orbit)
            or any(type(mask) is not int for mask in orbit)
            or len(set(orbit)) != len(orbit)
        ):
            raise InputError("orbit report has a malformed shell orbit")
        partition.append(orbit)
    flattened = [mask for orbit in partition for mask in orbit]
    if len(flattened) != len(set(flattened)) or set(flattened) != set(masks):
        raise InputError("orbit report does not exactly partition this shell")

    count_system = report.get("orbit_count_system")
    raw_forced = (
        count_system.get("forced_orbit_counts")
        if isinstance(count_system, dict)
        else None
    )
    if not isinstance(raw_forced, dict):
        raise InputError("orbit report has no forced-count record")
    forced: dict[int, int] = {}
    for raw_index, count in raw_forced.items():
        try:
            index = int(raw_index)
        except (TypeError, ValueError) as error:
            raise InputError("invalid forced orbit index") from error
        if (
            not 0 <= index < len(partition)
            or type(count) is not int
            or not 0 <= count <= len(partition[index])
        ):
            raise InputError("invalid forced orbit count")
        forced[index] = count
    return partition, forced, report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shell-report", required=True, type=Path)
    parser.add_argument("--factor", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument(
        "--exclude-factor",
        action="append",
        default=[],
        type=Path,
        help="exclude one normalized column set (repeatable)",
    )
    parser.add_argument(
        "--exclude-factor-list",
        action="append",
        default=[],
        type=Path,
        help=(
            "UTF-8 file containing one factor path per line; relative paths "
            "are resolved from the list file (repeatable)"
        ),
    )
    parser.add_argument(
        "--hint-factor",
        type=Path,
        help="exact target-Gram factor used as a complete CP-SAT hint",
    )
    parser.add_argument(
        "--minimize-overlap-with",
        type=Path,
        help="minimize the number of columns shared with this factor",
    )
    parser.add_argument(
        "--random-objective-seed",
        type=int,
        help=(
            "minimize reproducible independent random mask weights; mutually "
            "exclusive with --minimize-overlap-with"
        ),
    )
    parser.add_argument(
        "--random-objective-maximum",
        type=int,
        default=1_000_000,
    )
    parser.add_argument(
        "--require-mask",
        action="append",
        default=[],
        type=parse_mask,
    )
    parser.add_argument(
        "--forbid-mask",
        action="append",
        default=[],
        type=parse_mask,
    )
    parser.add_argument(
        "--orbit-report",
        type=Path,
        help=(
            "exact signed-Gram shell-orbit audit; mathematically forced counts "
            "are added automatically"
        ),
    )
    parser.add_argument(
        "--orbit-counts",
        type=parse_orbit_counts,
        help=(
            "restrict to one complete comma-separated orbit-count vector; "
            "requires --orbit-report and defines a search slice"
        ),
    )
    parser.add_argument("--time-limit", default=600.0, type=float)
    parser.add_argument("--workers", default=8, type=int)
    parser.add_argument("--seed", default=1, type=int)
    parser.add_argument("--log-search-progress", action="store_true")
    arguments = parser.parse_args()

    if arguments.time_limit <= 0 or not math.isfinite(arguments.time_limit):
        parser.error("--time-limit must be finite and positive")
    if arguments.workers <= 0:
        parser.error("--workers must be positive")
    if arguments.output.resolve() == arguments.metadata.resolve():
        parser.error("--output and --metadata must differ")
    if arguments.output.exists() or arguments.metadata.exists():
        parser.error("refusing to overwrite output or metadata")
    if set(arguments.require_mask) & set(arguments.forbid_mask):
        parser.error("the same mask cannot be both required and forbidden")
    if arguments.orbit_counts is not None and arguments.orbit_report is None:
        parser.error("--orbit-counts requires --orbit-report")
    if (
        arguments.random_objective_seed is not None
        and arguments.minimize_overlap_with is not None
    ):
        parser.error(
            "--random-objective-seed and --minimize-overlap-with are exclusive"
        )
    if arguments.random_objective_maximum <= 0:
        parser.error("--random-objective-maximum must be positive")

    started = time.monotonic()
    shell_path = arguments.shell_report.expanduser().absolute()
    factor_path = arguments.factor.expanduser().absolute()
    report = read_json(shell_path)
    masks = parse_shell(report)
    mask_to_index = {mask: index for index, mask in enumerate(masks)}
    source_factor = read_matrix(factor_path)
    target_gram = gram(source_factor)
    source_determinant = bareiss_determinant(source_factor)
    if source_determinant == 0:
        raise InputError("source factor is singular")
    source_support = normalized_column_masks(source_factor)
    if any(mask not in mask_to_index for mask in source_support):
        raise InputError("source factor contains a column outside the shell")
    reported_support = report.get("known_factor_column_masks")
    if (
        not isinstance(reported_support, list)
        or sorted(reported_support) != sorted(source_support)
    ):
        raise InputError("shell report is not bound to the supplied factor")

    orbit_path: Path | None = None
    orbit_partition: list[list[int]] = []
    forced_orbit_counts: dict[int, int] = {}
    orbit_report: dict[str, Any] | None = None
    if arguments.orbit_report is not None:
        orbit_path = arguments.orbit_report.expanduser().absolute()
        (
            orbit_partition,
            forced_orbit_counts,
            orbit_report,
        ) = parse_orbit_report(
            orbit_path, shell_path, factor_path, masks
        )
    requested_orbit_counts = arguments.orbit_counts
    if requested_orbit_counts is not None:
        if len(requested_orbit_counts) != len(orbit_partition):
            raise InputError(
                "--orbit-counts length does not match the orbit partition"
            )
        if sum(requested_orbit_counts) != ORDER:
            raise InputError("--orbit-counts must sum to the matrix order")
        for index, count in enumerate(requested_orbit_counts):
            if count > len(orbit_partition[index]):
                raise InputError(
                    f"--orbit-counts[{index}] exceeds its orbit size"
                )
            forced = forced_orbit_counts.get(index)
            if forced is not None and count != forced:
                raise InputError(
                    f"--orbit-counts[{index}] conflicts with forced count {forced}"
                )
        count_system = orbit_report.get("orbit_count_system")
        equations = (
            count_system.get("equations")
            if isinstance(count_system, dict)
            else None
        )
        if not isinstance(equations, list):
            raise InputError("orbit report has no exact count equations")
        for equation in equations:
            if (
                not isinstance(equation, dict)
                or not isinstance(equation.get("coefficients"), list)
                or type(equation.get("rhs")) is not int
                or len(equation["coefficients"]) != len(orbit_partition)
                or any(type(value) is not int for value in equation["coefficients"])
            ):
                raise InputError("orbit report has a malformed count equation")
            if (
                sum(
                    coefficient * count
                    for coefficient, count in zip(
                        equation["coefficients"], requested_orbit_counts
                    )
                )
                != equation["rhs"]
            ):
                raise InputError("--orbit-counts violates an exact invariant")

    factor_list_records = []
    exclusion_paths = [
        path.expanduser().absolute() for path in arguments.exclude_factor
    ]
    for raw_list_path in arguments.exclude_factor_list:
        list_path = raw_list_path.expanduser().absolute()
        listed = read_factor_list(list_path)
        factor_list_records.append(
            {
                "path": str(list_path),
                "sha256": sha256(list_path),
                "factor_count": len(listed),
            }
        )
        exclusion_paths.extend(listed)
    unique_exclusion_paths: dict[Path, None] = {}
    for path in exclusion_paths:
        unique_exclusion_paths[path.resolve()] = None

    exclusion_records: list[dict[str, Any]] = []
    exclusion_supports: list[list[int]] = []
    for path in unique_exclusion_paths:
        matrix = read_matrix(path)
        if gram(matrix) != target_gram:
            raise InputError(f"excluded factor has a different row Gram: {path}")
        support = normalized_column_masks(matrix)
        if any(mask not in mask_to_index for mask in support):
            raise InputError(f"excluded factor leaves the shell: {path}")
        exclusion_supports.append(support)
        exclusion_records.append({"path": str(path), "sha256": sha256(path)})

    hint_path = (
        arguments.hint_factor.expanduser().absolute()
        if arguments.hint_factor is not None
        else factor_path
    )
    hint_matrix = read_matrix(hint_path)
    if gram(hint_matrix) != target_gram:
        raise InputError(f"hint factor has a different row Gram: {hint_path}")
    hint_support = normalized_column_masks(hint_matrix)
    if any(mask not in mask_to_index for mask in hint_support):
        raise InputError(f"hint factor leaves the shell: {hint_path}")

    overlap_support: list[int] | None = None
    overlap_record: dict[str, Any] | None = None
    if arguments.minimize_overlap_with is not None:
        overlap_path = arguments.minimize_overlap_with.expanduser().absolute()
        overlap_matrix = read_matrix(overlap_path)
        if gram(overlap_matrix) != target_gram:
            raise InputError("overlap factor has a different row Gram")
        overlap_support = normalized_column_masks(overlap_matrix)
        if any(mask not in mask_to_index for mask in overlap_support):
            raise InputError("overlap factor leaves the shell")
        overlap_record = {
            "path": str(overlap_path),
            "sha256": sha256(overlap_path),
        }

    for mask in [*arguments.require_mask, *arguments.forbid_mask]:
        if mask not in mask_to_index:
            raise InputError(f"requested mask is outside the shell: {mask}")

    model_started = time.monotonic()
    model = cp_model.CpModel()
    variables = [
        model.NewBoolVar(f"x_{index}") for index in range(len(masks))
    ]
    model.Add(sum(variables) == ORDER)
    constraint_term_count = len(variables)
    for row in range(ORDER):
        for column in range(row + 1, ORDER):
            target = target_gram[row][column]
            if (ORDER + target) % 2:
                raise InputError("Gram parity is incompatible with sign columns")
            same_indices = [
                index
                for index, mask in enumerate(masks)
                if ((mask >> row) & 1) == ((mask >> column) & 1)
            ]
            different_count = len(masks) - len(same_indices)
            if len(same_indices) <= different_count:
                model.Add(
                    sum(variables[index] for index in same_indices)
                    == (ORDER + target) // 2
                )
                constraint_term_count += len(same_indices)
            else:
                same = set(same_indices)
                different_indices = [
                    index for index in range(len(masks)) if index not in same
                ]
                model.Add(
                    sum(variables[index] for index in different_indices)
                    == (ORDER - target) // 2
                )
                constraint_term_count += len(different_indices)

    for support in exclusion_supports:
        model.Add(
            sum(variables[mask_to_index[mask]] for mask in support) <= ORDER - 1
        )
    for mask in arguments.require_mask:
        model.Add(variables[mask_to_index[mask]] == 1)
    for mask in arguments.forbid_mask:
        model.Add(variables[mask_to_index[mask]] == 0)
    for orbit_index, count in sorted(forced_orbit_counts.items()):
        model.Add(
            sum(
                variables[mask_to_index[mask]]
                for mask in orbit_partition[orbit_index]
            )
            == count
        )
        constraint_term_count += len(orbit_partition[orbit_index])
    if requested_orbit_counts is not None:
        for orbit, count in zip(orbit_partition, requested_orbit_counts):
            model.Add(
                sum(variables[mask_to_index[mask]] for mask in orbit) == count
            )
            constraint_term_count += len(orbit)

    hint_support_set = set(hint_support)
    for mask, variable in zip(masks, variables):
        model.AddHint(variable, int(mask in hint_support_set))
    random_objective_weights: list[int] | None = None
    if overlap_support is not None:
        model.Minimize(
            sum(variables[mask_to_index[mask]] for mask in overlap_support)
        )
    elif arguments.random_objective_seed is not None:
        objective_rng = random.Random(arguments.random_objective_seed)
        random_objective_weights = [
            objective_rng.randint(1, arguments.random_objective_maximum)
            for _ in masks
        ]
        model.Minimize(
            sum(
                weight * variable
                for weight, variable in zip(
                    random_objective_weights, variables
                )
            )
        )

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = arguments.time_limit
    solver.parameters.num_search_workers = arguments.workers
    solver.parameters.random_seed = arguments.seed
    solver.parameters.randomize_search = True
    solver.parameters.symmetry_level = 2
    solver.parameters.log_search_progress = arguments.log_search_progress
    model_seconds = time.monotonic() - model_started

    status = solver.Solve(model)
    status_name = solver.StatusName(status)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    selected_masks: list[int] = []
    output_record: dict[str, Any] | None = None
    exact_overlap: int | None = None
    exact_random_objective: int | None = None
    if feasible:
        selected_masks = [
            mask
            for mask, variable in zip(masks, variables)
            if solver.BooleanValue(variable)
        ]
        factor = reconstruct_factor(selected_masks)
        if gram(factor) != target_gram:
            raise ArithmeticError("solver factor failed exact Gram reconstruction")
        determinant = bareiss_determinant(factor)
        if abs(determinant) != abs(source_determinant):
            raise ArithmeticError("solver factor determinant changed")
        contents = matrix_bytes(factor)
        atomic_write(arguments.output.expanduser().absolute(), contents)
        output_record = {
            "path": str(arguments.output.expanduser().absolute()),
            "sha256": sha256_bytes(contents),
            "determinant": determinant,
            "selected_masks": selected_masks,
        }
        if overlap_support is not None:
            exact_overlap = len(set(selected_masks) & set(overlap_support))
        if random_objective_weights is not None:
            exact_random_objective = sum(
                random_objective_weights[mask_to_index[mask]]
                for mask in selected_masks
            )

    metadata = {
        "engine": "order22-gram-factor-cpsat-v3",
        "order": ORDER,
        "claim": (
            "solver proposal independently checked by exact Gram equality and "
            "Bareiss determinant; a timeout or solver status is not an "
            "arena optimality proof"
        ),
        "claim_boundary": (
            "A requested orbit-count vector restricts the search to that "
            "family and cannot prove absence outside it. Solver UNKNOWN is not "
            "infeasibility, and exact-support exclusions do not quotient "
            "signed-Gram automorphism or H-equivalence orbits."
        ),
        "environment": {
            "ortools": getattr(ortools, "__version__", "unknown"),
        },
        "inputs": {
            "shell_report": {
                "path": str(shell_path),
                "sha256": sha256(shell_path),
                "shell_size": len(masks),
            },
            "factor": {
                "path": str(factor_path),
                "sha256": sha256(factor_path),
                "determinant": source_determinant,
            },
            "excluded_factors": exclusion_records,
            "excluded_factor_lists": factor_list_records,
            "hint_factor": {
                "path": str(hint_path),
                "sha256": sha256(hint_path),
                "used_source_factor_default": arguments.hint_factor is None,
            },
            "minimize_overlap_with": overlap_record,
            "random_objective": (
                {
                    "seed": arguments.random_objective_seed,
                    "maximum_weight": arguments.random_objective_maximum,
                    "mask_order": "shell_report.shell_masks",
                }
                if arguments.random_objective_seed is not None
                else None
            ),
            "required_masks": arguments.require_mask,
            "forbidden_masks": arguments.forbid_mask,
            "orbit_report": (
                {
                    "path": str(orbit_path),
                    "sha256": sha256(orbit_path),
                    "forced_orbit_counts": {
                        str(index): count
                        for index, count in sorted(forced_orbit_counts.items())
                    },
                }
                if orbit_path is not None
                else None
            ),
            "requested_orbit_count_slice": requested_orbit_counts,
        },
        "model": {
            "boolean_variables": len(variables),
            "gram_constraints": ORDER * (ORDER - 1) // 2,
            "constraint_literal_occurrences": constraint_term_count,
            "build_seconds": model_seconds,
            "orbit_count_slice": requested_orbit_counts is not None,
        },
        "solver": {
            "status": status_name,
            "feasible": feasible,
            "seed": arguments.seed,
            "workers": arguments.workers,
            "time_limit_seconds": arguments.time_limit,
            "wall_time_seconds": solver.WallTime(),
            "branches": solver.NumBranches(),
            "conflicts": solver.NumConflicts(),
            "objective_value": solver.ObjectiveValue() if feasible else None,
            "best_objective_bound": (
                solver.BestObjectiveBound() if feasible else None
            ),
            "exact_overlap": exact_overlap,
            "exact_random_objective": exact_random_objective,
        },
        "output": output_record,
        "elapsed_seconds": time.monotonic() - started,
    }
    metadata_bytes = (
        json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    atomic_write(arguments.metadata.expanduser().absolute(), metadata_bytes)
    print(
        f"status={status_name} shell={len(masks)} "
        f"selected={len(selected_masks)} elapsed={metadata['elapsed_seconds']:.3f}s"
    )
    return 0 if feasible else 1


if __name__ == "__main__":
    raise SystemExit(main())
