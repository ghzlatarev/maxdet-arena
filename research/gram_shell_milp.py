#!/usr/bin/env python3
"""Search an exact sign factor from an exported Gram sign-column shell.

For a positive-definite Gram G, every normalized sign column s of a factor
R lies on s^T G^-1 s = 1.  If the shell vectors are s_j, factorization is
equivalent to finding nonnegative integer multiplicities x_j satisfying

    sum_j x_j s_j s_j^T = G.

SciPy/HiGHS supplies candidate integer multiplicities.  A result is promoted
only after exact reconstruction, exact Gram equality, and Bareiss determinant
checks.  Solver infeasibility or timeout is exploratory, not a proof.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tempfile
import time
from pathlib import Path
from typing import Any

import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, milp

from gram_hasse import (
    ORDER,
    InputError,
    bareiss_determinant,
    build_gram,
)


FRONTIER = 2_779_447_296_000_000


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


def read_sign_matrix(path: Path) -> list[list[int]]:
    try:
        rows = [
            [int(value) for value in line.split()]
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    except (OSError, ValueError) as error:
        raise InputError(f"cannot read sign matrix {path}: {error}") from error
    if len(rows) != ORDER or any(
        len(row) != ORDER or any(value not in (-1, 1) for value in row)
        for row in rows
    ):
        raise InputError(f"{path} must contain exactly 23x23 signs")
    original_gram = gram(rows)
    supported_levels = {-5, -1, 3, 7}
    for row in range(1, ORDER):
        value = original_gram[0][row]
        if value in supported_levels:
            continue
        if -value not in supported_levels:
            raise InputError(
                f"{path} Gram cannot be row-switched into supported levels"
            )
        rows[row] = [-entry for entry in rows[row]]
    return rows


def gram(matrix: list[list[int]]) -> list[list[int]]:
    return [
        [
            sum(matrix[row][index] * matrix[column][index] for index in range(ORDER))
            for column in range(ORDER)
        ]
        for row in range(ORDER)
    ]


def resolve_source_path(report_path: Path, source: str) -> Path:
    path = Path(source)
    if path.is_absolute():
        return path
    cwd_path = Path.cwd() / path
    if cwd_path.exists():
        return cwd_path
    return report_path.parent / path


def selected_gram(
    report_path: Path,
    report: dict[str, Any],
    result: dict[str, Any],
) -> tuple[list[list[int]], dict[str, Any]]:
    if "hit_index" in result:
        snapshot_value = report.get("snapshot")
        if not isinstance(snapshot_value, str):
            raise InputError("snapshot result is missing report.snapshot")
        snapshot_path = resolve_source_path(report_path, snapshot_value)
        expected_sha = report.get("snapshot_sha256")
        actual_sha = sha256(snapshot_path)
        if expected_sha != actual_sha:
            raise InputError("shell report snapshot SHA-256 no longer matches")
        snapshot = read_json(snapshot_path)
        hits = snapshot.get("hits")
        hit_index = result.get("hit_index")
        if (
            not isinstance(hits, list)
            or type(hit_index) is not int
            or not (0 <= hit_index < len(hits))
            or not isinstance(hits[hit_index], dict)
        ):
            raise InputError("invalid snapshot hit selection")
        normalization = snapshot.get("normalization")
        if not isinstance(normalization, str):
            raise InputError("snapshot normalization is missing")
        return build_gram(hits[hit_index], normalization), {
            "kind": "snapshot-hit",
            "path": str(snapshot_path),
            "sha256": actual_sha,
            "hit_index": hit_index,
            "normalization": normalization,
        }

    source = result.get("source")
    if not isinstance(source, str) or source == "snapshot-hit":
        raise InputError("matrix result is missing a usable source path")
    matrix_path = resolve_source_path(report_path, source)
    expected_sha = result.get("source_sha256")
    actual_sha = sha256(matrix_path)
    if expected_sha != actual_sha:
        raise InputError("shell report matrix SHA-256 no longer matches")
    matrix = read_sign_matrix(matrix_path)
    return gram(matrix), {
        "kind": "matrix",
        "path": str(matrix_path),
        "sha256": actual_sha,
        "normalization": result.get("source_normalization"),
    }


def parse_shell_masks(result: dict[str, Any]) -> list[int]:
    if result.get("rejected") is True:
        reason = result.get("reason", "shell obstruction")
        raise InputError(
            f"selected Gram is already exactly rejected: {reason}"
        )
    masks = result.get("shell_sign_masks")
    shell_size = result.get("shell_size")
    if not isinstance(masks, list) or type(shell_size) is not int:
        raise InputError(
            "shell vectors missing; rerun gram_shell_filter "
            "with --include-shell-vectors"
        )
    if len(masks) != shell_size:
        raise InputError("exported shell count does not match shell_size")
    if shell_size == 0:
        raise InputError("empty sign-column shell cannot contain a factor")
    parsed: list[int] = []
    for index, mask in enumerate(masks):
        if (
            type(mask) is not int
            or mask < 0
            or mask >= 1 << ORDER
            or mask & 1 == 0
        ):
            raise InputError(f"invalid normalized shell mask at index {index}")
        parsed.append(mask)
    if len(set(parsed)) != len(parsed):
        raise InputError("exported shell masks are not unique")
    return parsed


def feature_vector(mask: int) -> list[int]:
    signs = [1 if mask >> index & 1 else -1 for index in range(ORDER)]
    features = [1]
    for row in range(ORDER):
        for column in range(row + 1, ORDER):
            features.append(signs[row] * signs[column])
    return features


def target_vector(gram_matrix: list[list[int]]) -> list[int]:
    target = [ORDER]
    for row in range(ORDER):
        if len(gram_matrix[row]) != ORDER or gram_matrix[row][row] != ORDER:
            raise InputError("Gram matrix has invalid shape or diagonal")
        for column in range(row + 1, ORDER):
            if gram_matrix[row][column] != gram_matrix[column][row]:
                raise InputError("Gram matrix must be symmetric")
            target.append(gram_matrix[row][column])
    return target


def reconstruct_factor(
    masks: list[int], multiplicities: list[int]
) -> list[list[int]]:
    columns: list[list[int]] = []
    for mask, count in zip(masks, multiplicities):
        signs = [1 if mask >> index & 1 else -1 for index in range(ORDER)]
        columns.extend([signs] * count)
    if len(columns) != ORDER:
        raise ArithmeticError(
            f"candidate uses {len(columns)} columns instead of {ORDER}"
        )
    return [
        [columns[column][row] for column in range(ORDER)]
        for row in range(ORDER)
    ]


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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shell-report", required=True, type=Path)
    parser.add_argument("--result-index", type=int, default=0)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--node-limit", type=int)
    parser.add_argument(
        "--seed",
        type=int,
        help="use a reproducible random linear objective to sample another factor",
    )
    parser.add_argument(
        "--positive-objective",
        action="store_true",
        help=(
            "draw strictly positive seeded objective coefficients; useful with "
            "a nonzero MIP gap when only the first exact feasible factor matters"
        ),
    )
    parser.add_argument(
        "--binary",
        action="store_true",
        help=(
            "bound multiplicities by one (exact for a nonsingular normalized "
            "sign factor, whose columns cannot repeat up to sign)"
        ),
    )
    parser.add_argument("--mip-rel-gap", type=float, default=0.0)
    parser.add_argument(
        "--exclude-factor",
        action="append",
        default=[],
        type=Path,
        help="exclude one exact normalized-column subset (repeatable)",
    )
    arguments = parser.parse_args()

    if arguments.time_limit <= 0 or not math.isfinite(arguments.time_limit):
        parser.error("--time-limit must be finite and positive")
    if arguments.node_limit is not None and arguments.node_limit <= 0:
        parser.error("--node-limit must be positive")
    if not 0.0 <= arguments.mip_rel_gap <= 1.0:
        parser.error("--mip-rel-gap must be in [0,1]")
    if arguments.positive_objective and arguments.seed is None:
        parser.error("--positive-objective requires --seed")
    if arguments.output.resolve() == arguments.metadata.resolve():
        parser.error("--output and --metadata must be distinct")
    if arguments.output.exists() or arguments.metadata.exists():
        parser.error("refusing to overwrite output or metadata")

    report = read_json(arguments.shell_report)
    results = report.get("results")
    if (
        not isinstance(results, list)
        or not (0 <= arguments.result_index < len(results))
        or not isinstance(results[arguments.result_index], dict)
    ):
        raise InputError("--result-index is outside shell report results")
    result = results[arguments.result_index]
    masks = parse_shell_masks(result)
    gram_matrix, source = selected_gram(
        arguments.shell_report, report, result
    )
    gram_determinant = bareiss_determinant(gram_matrix)
    if str(gram_determinant) != result.get("determinant"):
        raise InputError("shell report determinant does not match source Gram")
    gram_root = math.isqrt(gram_determinant)
    if gram_root * gram_root != gram_determinant:
        raise InputError("source Gram determinant is not a square")

    constraint_matrix = np.asarray(
        [feature_vector(mask) for mask in masks], dtype=np.float64
    ).T
    target = np.asarray(target_vector(gram_matrix), dtype=np.float64)
    constraints: list[LinearConstraint] = [
        LinearConstraint(constraint_matrix, target, target)
    ]
    mask_to_index = {mask: index for index, mask in enumerate(masks)}
    excluded_factor_hashes = []
    for path in arguments.exclude_factor:
        excluded = read_sign_matrix(path)
        if gram(excluded) != gram_matrix:
            raise InputError(f"--exclude-factor does not have target Gram: {path}")
        support = set()
        for column in range(ORDER):
            sign = excluded[0][column]
            mask = sum(
                (excluded[row][column] * sign == 1) << row
                for row in range(ORDER)
            )
            support.add(mask)
        if len(support) != ORDER or not support.issubset(mask_to_index):
            raise InputError(f"--exclude-factor is not a shell subset: {path}")
        no_good = np.zeros(len(masks), dtype=np.float64)
        for mask in support:
            no_good[mask_to_index[mask]] = 1.0
        constraints.append(
            LinearConstraint(no_good, -np.inf, float(ORDER - 1))
        )
        excluded_factor_hashes.append(sha256(path))
    options: dict[str, float | int | bool] = {
        "time_limit": arguments.time_limit,
        "mip_rel_gap": arguments.mip_rel_gap,
        "presolve": True,
    }
    if arguments.node_limit is not None:
        options["node_limit"] = arguments.node_limit
    objective = np.zeros(len(masks), dtype=np.float64)
    if arguments.seed is not None:
        randomizer = np.random.default_rng(arguments.seed)
        if arguments.positive_objective:
            objective = randomizer.integers(
                1,
                1_000_001,
                size=len(masks),
            ).astype(np.float64)
        else:
            objective = randomizer.integers(
                -1_000_000,
                1_000_001,
                size=len(masks),
            ).astype(np.float64)
    solve_started = time.monotonic()
    solution = milp(
        c=objective,
        integrality=np.ones(len(masks), dtype=np.uint8),
        bounds=Bounds(
            np.zeros(len(masks), dtype=np.float64),
            np.full(
                len(masks),
                1 if arguments.binary else ORDER,
                dtype=np.float64,
            ),
        ),
        constraints=constraints,
        options=options,
    )
    runtime_seconds = time.monotonic() - solve_started

    verified = False
    multiplicities: list[int] = []
    factor: list[list[int]] | None = None
    factor_score: int | None = None
    if solution.x is not None:
        multiplicities = [int(round(value)) for value in solution.x]
        residual = np.max(
            np.abs(
                constraint_matrix
                @ np.asarray(multiplicities, dtype=np.float64)
                - target
            )
        )
        if residual == 0:
            factor = reconstruct_factor(masks, multiplicities)
            factor_gram = gram(factor)
            factor_score = abs(bareiss_determinant(factor))
            verified = (
                factor_gram == gram_matrix
                and factor_score == gram_root
                and factor_score * factor_score == gram_determinant
            )

    metadata = {
        "claim_boundary": (
            "A factor_found=true result is independently reconstructed and "
            "checked in exact integers. Other solver statuses are exploratory "
            "and do not prove non-factorability. Solver status and a nonzero "
            "MIP gap are not used as an objective-optimality claim."
        ),
        "engine": "gram-shell-milp",
        "engine_source": str(Path(__file__).resolve()),
        "engine_source_sha256": sha256(Path(__file__).resolve()),
        "dependencies": {
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "factor_found": verified,
        "factor_score": str(factor_score) if verified else None,
        "frontier": str(FRONTIER),
        "gram_determinant": str(gram_determinant),
        "gram_square_root": str(gram_root),
        "multiplicities": [
            {"count": count, "shell_mask": mask}
            for mask, count in zip(masks, multiplicities)
            if count
        ],
        "result_index": arguments.result_index,
        "seed": arguments.seed,
        "binary": arguments.binary,
        "positive_objective": arguments.positive_objective,
        "mip_rel_gap": arguments.mip_rel_gap,
        "excluded_factor_count": len(arguments.exclude_factor),
        "excluded_factor_sha256": excluded_factor_hashes,
        "runtime_seconds": runtime_seconds,
        "requested_limits": {
            "time_limit_seconds": arguments.time_limit,
            "node_limit": arguments.node_limit,
        },
        "shell_report": str(arguments.shell_report),
        "shell_report_sha256": sha256(arguments.shell_report),
        "shell_size": len(masks),
        "solver": {
            "message": solution.message,
            "mip_gap": getattr(solution, "mip_gap", None),
            "mip_node_count": getattr(solution, "mip_node_count", None),
            "status": int(solution.status),
            "success": bool(solution.success),
        },
        "source": source,
    }
    output_contents = (
        "".join(" ".join(map(str, row)) + "\n" for row in factor)
        if verified and factor is not None
        else None
    )
    metadata["output"] = str(arguments.output)
    metadata["output_sha256"] = (
        hashlib.sha256(output_contents.encode("ascii")).hexdigest()
        if output_contents is not None
        else None
    )
    if not verified or factor is None:
        atomic_write_text(
            arguments.metadata,
            json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        )
        print(
            f"no verified factor: status={solution.status} "
            f"message={solution.message}"
        )
        return 3

    assert output_contents is not None
    atomic_write_text(arguments.output, output_contents)
    atomic_write_text(
        arguments.metadata,
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
    )
    print(
        f"verified factor |det|={factor_score} shell={len(masks)} "
        f"nonzero_multiplicities={sum(count != 0 for count in multiplicities)}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (InputError, ArithmeticError, OSError) as error:
        raise SystemExit(f"gram_shell_milp: {error}") from error
