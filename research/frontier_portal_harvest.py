#!/usr/bin/env python3
"""Harvest exact frontier factors from novelty-biased Gram-shell slices.

The published order-23 frontier Gram has a complete normalized sign-column
shell of size 1,382.  Exact orbit equations force every factor to choose three
of the six columns in the smallest shell orbit.  Those triples have four Gram
automorphism orbits; every currently known HT class aligns to orbit (0,2,4).

This tool deterministically shards randomized positive-objective HiGHS runs
over all four triple orbits, prioritizing the three unobserved orbits.  It does
not add raw-support no-goods.  Every emitted factor receives exact Gram and
Bareiss checks, an arena receipt, and pinned pynauty H/HT/Gram certificates.
Novelty is always reported relative to an explicit, hashed baseline audit.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import math
import os
import platform
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
ORDER = 23
FRONTIER = 2_779_447_296_000_000
PINNED_PYNAUTY = "2.8.8.1"
SMALL_ORBIT = (8161, 26521, 8356323, 8357403, 8363525, 8380805)
TRIPLE_ORBITS = (
    (0, 1, 2),
    (0, 1, 3),
    (0, 2, 5),
    (0, 2, 4),
)
DEFAULT_SHELL = (
    ROOT / "runs/direct-search/gram-shell-reference-columns-29760.json"
)
DEFAULT_BASELINE = (
    ROOT
    / (
        "runs/direct-search/frontier-factor-class-expansion-20260728/"
        "final-h-equivalence-audit.json"
    )
)
DEFAULT_ALIGNED = (
    ROOT
    / (
        "runs/direct-search/frontier-factor-class-expansion-20260728/"
        "known-ht-aligned"
    )
)
DEFAULT_EB138A = (
    ROOT
    / (
        "runs/direct-search/frontier-factor-class-expansion-20260728/"
        "seeds/"
        "h-eb138a06ec638735c34bdacf77bd1cdd869c5d2fbc3450be25d63af0cde1a134"
        ".matrix.txt"
    )
)
DEFAULT_OUTPUT = (
    ROOT / "runs/direct-search/frontier-portal-harvest-20260729"
)

Matrix = list[list[int]]


def sha256_bytes(contents: bytes) -> str:
    return hashlib.sha256(contents).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def display_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def atomic_write(path: Path, contents: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(contents)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_json(path: Path, value: Any) -> None:
    atomic_write(
        path,
        (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )


def read_matrix(path: Path) -> Matrix:
    matrix = [
        [int(token) for token in line.split()]
        for line in path.read_text(encoding="ascii").splitlines()
        if line.strip()
    ]
    if len(matrix) != ORDER or any(
        len(row) != ORDER or any(value not in (-1, 1) for value in row)
        for row in matrix
    ):
        raise ValueError(f"{path}: expected exactly {ORDER}x{ORDER} signs")
    return matrix


def matrix_bytes(matrix: Matrix) -> bytes:
    return "".join(
        " ".join(map(str, row)) + "\n" for row in matrix
    ).encode("ascii")


def gram(matrix: Matrix) -> Matrix:
    return [
        [
            sum(
                matrix[row][index] * matrix[column][index]
                for index in range(ORDER)
            )
            for column in range(ORDER)
        ]
        for row in range(ORDER)
    ]


def determinant(matrix: Matrix) -> int:
    work = [row[:] for row in matrix]
    sign = 1
    denominator = 1
    for pivot_index in range(len(work) - 1):
        pivot_row = next(
            (
                row
                for row in range(pivot_index, len(work))
                if work[row][pivot_index]
            ),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = (
                work[pivot_row],
                work[pivot_index],
            )
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, len(work)):
            for column in range(pivot_index + 1, len(work)):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                if numerator % denominator:
                    raise ArithmeticError("Bareiss division was not exact")
                work[row][column] = numerator // denominator
            work[row][pivot_index] = 0
        denominator = pivot
    return sign * work[-1][-1]


def normalized_column_masks(matrix: Matrix) -> list[int]:
    masks = []
    for column in range(ORDER):
        switch = matrix[0][column]
        masks.append(
            sum(
                (
                    matrix[row][column] * switch == 1
                )
                << row
                for row in range(ORDER)
            )
        )
    return masks


def feature_vector(mask: int) -> list[int]:
    signs = [1 if mask >> row & 1 else -1 for row in range(ORDER)]
    features = [1]
    for row in range(ORDER):
        for column in range(row + 1, ORDER):
            features.append(signs[row] * signs[column])
    return features


def target_vector(target_gram: Matrix) -> list[int]:
    target = [ORDER]
    for row in range(ORDER):
        for column in range(row + 1, ORDER):
            target.append(target_gram[row][column])
    return target


def reconstruct_factor(
    masks: list[int], multiplicities: list[int]
) -> Matrix:
    selected = [
        [1 if mask >> row & 1 else -1 for row in range(ORDER)]
        for mask, count in zip(masks, multiplicities)
        if count == 1
    ]
    if len(selected) != ORDER:
        raise ArithmeticError("factor does not select exactly 23 columns")
    return [
        [selected[column][row] for column in range(ORDER)]
        for row in range(ORDER)
    ]


def shell_inputs(path: Path) -> tuple[list[int], Path]:
    report = json.loads(path.read_text(encoding="utf-8"))
    results = report.get("results")
    if (
        not isinstance(results, list)
        or len(results) != 1
        or not isinstance(results[0], dict)
    ):
        raise ValueError("shell report must contain exactly one result")
    result = results[0]
    masks = result.get("shell_sign_masks")
    if (
        not isinstance(masks, list)
        or len(masks) != 1382
        or len(set(masks)) != 1382
        or not set(SMALL_ORBIT).issubset(masks)
    ):
        raise ValueError("shell report is not the complete 1,382-column shell")
    source_value = result.get("source")
    if not isinstance(source_value, str):
        raise ValueError("shell report lacks a matrix source")
    source = Path(source_value)
    if not source.is_absolute():
        source = ROOT / source
    source = source.resolve()
    if sha256_file(source) != result.get("source_sha256"):
        raise ValueError("shell source hash mismatch")
    if result.get("determinant") != str(FRONTIER * FRONTIER):
        raise ValueError("shell Gram determinant is not the frontier square")
    return [int(mask) for mask in masks], source


def objective_spec(global_index: int) -> tuple[str, tuple[int, int, int]]:
    """Return a deterministic mix favoring the three unobserved slices."""
    schedule = (
        ("rarity", TRIPLE_ORBITS[0]),
        ("rarity", TRIPLE_ORBITS[1]),
        ("rarity", TRIPLE_ORBITS[2]),
        ("anchor", TRIPLE_ORBITS[0]),
        ("anchor", TRIPLE_ORBITS[1]),
        ("anchor", TRIPLE_ORBITS[2]),
        ("rarity", TRIPLE_ORBITS[3]),
        ("control", TRIPLE_ORBITS[3]),
    )
    return schedule[global_index % len(schedule)]


def solve_hidden(spec_path: Path) -> int:
    import numpy as np
    import scipy
    from scipy.optimize import Bounds, LinearConstraint, milp

    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    trial_dir = spec_path.parent
    shell_path = Path(spec["shell_report"])
    masks, source = shell_inputs(shell_path)
    target_matrix = read_matrix(source)
    target_gram = gram(target_matrix)
    target = np.asarray(target_vector(target_gram), dtype=np.float64)
    constraint_matrix = np.asarray(
        [feature_vector(mask) for mask in masks], dtype=np.float64
    ).T
    known_supports = [
        {int(mask) for mask in support}
        for support in spec["known_supports"]
    ]
    frequency = {
        mask: sum(mask in support for support in known_supports)
        for mask in masks
    }
    seed = int(spec["seed"])
    randomizer = np.random.default_rng(seed)
    noise = randomizer.integers(
        1, 1_000_001, size=len(masks), dtype=np.int64
    )
    mode = str(spec["objective_mode"])
    global_index_mod = int(spec["global_index"]) % len(known_supports)
    anchor = known_supports[global_index_mod]
    if mode == "rarity":
        coefficients = (
            1
            + noise // 4
            + np.asarray(
                [450_000 * frequency[mask] for mask in masks],
                dtype=np.int64,
            )
        )
    elif mode == "anchor":
        coefficients = (
            1
            + noise // 2
            + np.asarray(
                [
                    125_000 * frequency[mask]
                    + 900_000 * (mask in anchor)
                    for mask in masks
                ],
                dtype=np.int64,
            )
        )
    elif mode == "control":
        coefficients = noise
    else:
        raise ValueError(f"unknown objective mode {mode}")
    objective_hash = sha256_bytes(
        np.asarray(coefficients, dtype="<i8").tobytes()
    )

    mask_to_index = {mask: index for index, mask in enumerate(masks)}
    selected = {
        SMALL_ORBIT[index] for index in spec["small_orbit_triple"]
    }
    lower = np.zeros(len(masks), dtype=np.float64)
    upper = np.ones(len(masks), dtype=np.float64)
    for mask in SMALL_ORBIT:
        index = mask_to_index[mask]
        if mask in selected:
            lower[index] = 1.0
        else:
            upper[index] = 0.0

    started = time.monotonic()
    solution = milp(
        c=np.asarray(coefficients, dtype=np.float64),
        integrality=np.ones(len(masks), dtype=np.uint8),
        bounds=Bounds(lower, upper),
        constraints=[
            LinearConstraint(constraint_matrix, target, target)
        ],
        options={
            "time_limit": float(spec["time_limit_seconds"]),
            "mip_rel_gap": 1.0,
            "presolve": True,
        },
    )
    elapsed = time.monotonic() - started

    factor: Matrix | None = None
    exact_determinant: int | None = None
    multiplicities: list[int] = []
    if solution.x is not None:
        multiplicities = [int(round(value)) for value in solution.x]
        if all(value in (0, 1) for value in multiplicities):
            residual = (
                constraint_matrix
                @ np.asarray(multiplicities, dtype=np.float64)
                - target
            )
            if np.max(np.abs(residual)) == 0:
                candidate = reconstruct_factor(masks, multiplicities)
                exact_determinant = abs(determinant(candidate))
                if (
                    gram(candidate) == target_gram
                    and exact_determinant == FRONTIER
                ):
                    factor = candidate

    metadata = {
        "schema_version": 1,
        "engine": "frontier-portal-gram-shell-milp-v1",
        "claim_boundary": (
            "factor_found=true is backed by exact Gram and Bareiss checks; "
            "a timeout or other non-factor result is not an infeasibility proof"
        ),
        "factor_found": factor is not None,
        "exact_absolute_determinant": (
            str(exact_determinant) if factor is not None else None
        ),
        "global_index": spec["global_index"],
        "shard_index": spec["shard_index"],
        "shard_count": spec["shard_count"],
        "seed": seed,
        "small_orbit_triple": spec["small_orbit_triple"],
        "small_orbit_masks": sorted(selected),
        "objective_mode": mode,
        "objective_sha256": objective_hash,
        "objective_min": int(np.min(coefficients)),
        "objective_max": int(np.max(coefficients)),
        "anchor_known_support_index": global_index_mod,
        "raw_support_no_goods": 0,
        "known_support_count": len(known_supports),
        "requested_time_limit_seconds": spec["time_limit_seconds"],
        "elapsed_seconds": elapsed,
        "solver": {
            "status": int(solution.status),
            "success": bool(solution.success),
            "message": solution.message,
            "mip_gap": getattr(solution, "mip_gap", None),
            "mip_node_count": getattr(solution, "mip_node_count", None),
        },
        "dependencies": {
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "shell_report": display_path(shell_path),
        "shell_report_sha256": sha256_file(shell_path),
        "source_matrix": display_path(source),
        "source_matrix_sha256": sha256_file(source),
        "selected_shell_masks": (
            [mask for mask, value in zip(masks, multiplicities) if value]
            if factor is not None
            else []
        ),
    }
    if factor is not None:
        output = trial_dir / "factor.matrix.txt"
        atomic_write(output, matrix_bytes(factor))
        metadata["output"] = display_path(output)
        metadata["output_sha256"] = sha256_file(output)
    write_json(trial_dir / "solver.json", metadata)
    return 0 if factor is not None else 3


def classify_hidden(matrix_path: Path) -> int:
    import pynauty
    from h_equivalence_audit import (
        h_certificate,
        normalized_gram_graph,
        sha256_hex,
        transpose,
    )

    if getattr(pynauty, "__version__", None) != PINNED_PYNAUTY:
        raise RuntimeError(f"requires pynauty=={PINNED_PYNAUTY}")
    matrix = read_matrix(matrix_path)
    direct = h_certificate(matrix)
    transposed = h_certificate(transpose(matrix))
    gram_certificate, group_order, edges, degrees = normalized_gram_graph(
        matrix
    )
    print(
        json.dumps(
            {
                "pynauty_version": pynauty.__version__,
                "exact_absolute_determinant": str(abs(determinant(matrix))),
                "h_certificate_sha256": sha256_hex(direct),
                "transpose_h_certificate_sha256": sha256_hex(transposed),
                "ht_certificate_sha256": sha256_hex(
                    min(direct, transposed)
                ),
                "self_dual": direct == transposed,
                "gram_certificate_sha256": sha256_hex(gram_certificate),
                "gram_automorphism_group_order": group_order,
                "gram_edge_count": edges,
                "gram_degree_multiset": list(degrees),
            },
            sort_keys=True,
        )
    )
    return 0


def classifier_python(explicit: Path | None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit)
    environment = os.environ.get("MAXDET_CLASSIFIER_PYTHON")
    if environment:
        candidates.append(Path(environment))
    candidates.extend(
        [
            Path(sys.executable),
            Path("/tmp/maxdet-h-audit-20260729/bin/python"),
        ]
    )
    seen: set[Path] = set()
    for candidate in candidates:
        executable = Path(os.path.abspath(candidate.expanduser()))
        if executable in seen or not executable.is_file():
            continue
        seen.add(executable)
        probe = subprocess.run(
            [
                str(executable),
                "-c",
                (
                    "import pynauty; "
                    "print(getattr(pynauty, '__version__', 'unknown'))"
                ),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if probe.returncode == 0 and probe.stdout.strip() == PINNED_PYNAUTY:
            return executable
    raise RuntimeError(
        "no pynauty==2.8.8.1 interpreter; pass --classifier-python"
    )


def execute_trial(
    python: Path, script: Path, spec_path: Path
) -> tuple[int, float]:
    started = time.monotonic()
    completed = subprocess.run(
        [str(python), str(script), "_solve", str(spec_path)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    atomic_write(
        spec_path.parent / "stdout.txt", completed.stdout.encode("utf-8")
    )
    atomic_write(
        spec_path.parent / "stderr.txt", completed.stderr.encode("utf-8")
    )
    return completed.returncode, time.monotonic() - started


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--shell-report", type=Path, default=DEFAULT_SHELL)
    parser.add_argument("--baseline-audit", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--aligned-dir", type=Path, default=DEFAULT_ALIGNED)
    parser.add_argument("--classifier-python", type=Path)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--trials-per-shard", type=int, default=8)
    parser.add_argument("--seed-base", type=int, default=730_000)
    parser.add_argument("--time-limit", type=float, default=90.0)
    parser.add_argument("--jobs", type=int, default=2)
    arguments = parser.parse_args()
    if (
        arguments.shard_count <= 0
        or not 0 <= arguments.shard_index < arguments.shard_count
        or arguments.trials_per_shard <= 0
        or arguments.jobs <= 0
        or arguments.time_limit <= 0
        or not math.isfinite(arguments.time_limit)
    ):
        parser.error("invalid shard, trial, job, or time-limit argument")

    output_dir = arguments.output_dir.expanduser().resolve()
    if output_dir.exists():
        raise FileExistsError(output_dir)
    output_dir.mkdir(parents=True)
    shell_path = arguments.shell_report.expanduser().resolve()
    baseline_path = arguments.baseline_audit.expanduser().resolve()
    aligned_dir = arguments.aligned_dir.expanduser().resolve()
    classifier = classifier_python(arguments.classifier_python)
    masks, source = shell_inputs(shell_path)
    target = read_matrix(source)
    target_gram = gram(target)
    if abs(determinant(target)) != FRONTIER:
        raise ArithmeticError("shell source is not an exact frontier matrix")

    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    if (
        baseline.get("expected_absolute_determinant") != FRONTIER
        or baseline.get("ht_class_count") != 8
        or baseline.get("gram_class_count") != 1
        or baseline.get("pynauty_version") != PINNED_PYNAUTY
    ):
        raise ValueError("baseline is not the pinned 8-HT / 1-Gram audit")
    baseline_h = {
        item["certificate_sha256"] for item in baseline["h_classes"]
    }
    baseline_ht = {
        item["certificate_sha256"] for item in baseline["ht_classes"]
    }
    baseline_gram = {
        item["certificate_sha256"] for item in baseline["gram_classes"]
    }

    known_paths = sorted(aligned_dir.glob("*.matrix.txt")) + [DEFAULT_EB138A]
    known_supports = []
    known_records = []
    shell_set = set(masks)
    expected_observed_small_masks = {
        SMALL_ORBIT[index] for index in (0, 2, 4)
    }
    for path in known_paths:
        matrix = read_matrix(path)
        support = sorted(set(normalized_column_masks(matrix)))
        if (
            gram(matrix) != target_gram
            or abs(determinant(matrix)) != FRONTIER
            or len(support) != ORDER
            or not set(support).issubset(shell_set)
        ):
            raise ArithmeticError(f"invalid aligned known factor: {path}")
        selected_small_masks = set(support) & set(SMALL_ORBIT)
        if selected_small_masks != expected_observed_small_masks:
            raise ArithmeticError(
                f"{path}: aligned known support does not select exactly "
                "the observed small-orbit triple (0,2,4)"
            )
        known_supports.append(support)
        known_records.append(
            {
                "path": display_path(path),
                "sha256": sha256_file(path),
                "support_sha256": sha256_bytes(
                    ",".join(map(str, support)).encode("ascii")
                ),
                "small_orbit_triple": [0, 2, 4],
                "small_orbit_masks": sorted(selected_small_masks),
            }
        )
    if len(known_supports) != 8:
        raise ArithmeticError("expected exactly eight aligned HT supports")
    print(
        "startup assertion passed: 8/8 aligned known HT supports select "
        "exactly small-orbit triple (0,2,4)",
        flush=True,
    )

    script = Path(__file__).resolve()
    specs = []
    for local_index in range(arguments.trials_per_shard):
        global_index = (
            arguments.shard_index + local_index * arguments.shard_count
        )
        mode, triple = objective_spec(global_index)
        trial_dir = output_dir / f"trial-{global_index:05d}"
        trial_dir.mkdir()
        spec = {
            "schema_version": 1,
            "global_index": global_index,
            "local_index": local_index,
            "shard_index": arguments.shard_index,
            "shard_count": arguments.shard_count,
            "seed": arguments.seed_base + global_index,
            "time_limit_seconds": arguments.time_limit,
            "objective_mode": mode,
            "small_orbit_triple": list(triple),
            "shell_report": str(shell_path),
            "shell_report_sha256": sha256_file(shell_path),
            "known_supports": known_supports,
            "known_support_count": len(known_supports),
            "raw_support_no_goods": 0,
        }
        spec_path = trial_dir / "spec.json"
        write_json(spec_path, spec)
        specs.append((global_index, spec_path))

    started = time.monotonic()
    executions: dict[int, tuple[int, float]] = {}
    with ThreadPoolExecutor(max_workers=arguments.jobs) as executor:
        futures = {
            executor.submit(
                execute_trial, Path(sys.executable), script, spec_path
            ): global_index
            for global_index, spec_path in specs
        }
        for future in as_completed(futures):
            global_index = futures[future]
            executions[global_index] = future.result()
            print(
                f"trial={global_index} returncode="
                f"{executions[global_index][0]} elapsed="
                f"{executions[global_index][1]:.1f}s",
                flush=True,
            )

    records = []
    observed_h = set(baseline_h)
    observed_ht = set(baseline_ht)
    observed_gram = set(baseline_gram)
    for global_index, spec_path in sorted(specs):
        trial_dir = spec_path.parent
        returncode, process_elapsed = executions[global_index]
        if returncode not in (0, 3):
            error = (trial_dir / "stderr.txt").read_text(encoding="utf-8")
            raise RuntimeError(f"trial {global_index} failed: {error.strip()}")
        solver_path = trial_dir / "solver.json"
        if not solver_path.is_file():
            raise RuntimeError(f"trial {global_index} omitted solver.json")
        solver = json.loads(solver_path.read_text(encoding="utf-8"))
        factor_found = bool(solver["factor_found"])
        classification = None
        receipt_record = None
        if factor_found:
            matrix_path = trial_dir / "factor.matrix.txt"
            classify = subprocess.run(
                [
                    str(classifier),
                    str(script),
                    "_classify",
                    str(matrix_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            atomic_write(
                trial_dir / "classifier.stderr.txt",
                classify.stderr.encode("utf-8"),
            )
            if classify.returncode != 0:
                raise RuntimeError(
                    f"classifier failed: {classify.stderr.strip()}"
                )
            classification = json.loads(classify.stdout)
            support = set(normalized_column_masks(read_matrix(matrix_path)))
            overlaps = [
                len(support & set(known)) for known in known_supports
            ]
            h = classification["h_certificate_sha256"]
            ht = classification["ht_certificate_sha256"]
            gram_sha = classification["gram_certificate_sha256"]
            classification.update(
                {
                    "baseline_known_h": h in baseline_h,
                    "baseline_known_ht": ht in baseline_ht,
                    "baseline_known_gram": gram_sha in baseline_gram,
                    "stream_new_h": h not in observed_h,
                    "stream_new_ht": ht not in observed_ht,
                    "stream_new_gram": gram_sha not in observed_gram,
                    "maximum_aligned_support_overlap": max(overlaps),
                    "minimum_aligned_support_symmetric_difference": (
                        2 * (ORDER - max(overlaps))
                    ),
                    "aligned_support_overlaps": overlaps,
                }
            )
            observed_h.add(h)
            observed_ht.add(ht)
            observed_gram.add(gram_sha)
            write_json(trial_dir / "classification.json", classification)

            receipt_path = trial_dir / "arena.receipt.json"
            arena = subprocess.run(
                [
                    str(ROOT / "arena"),
                    "verify",
                    str(matrix_path),
                    "--json",
                    str(receipt_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            atomic_write(
                trial_dir / "arena.stdout.txt", arena.stdout.encode("utf-8")
            )
            atomic_write(
                trial_dir / "arena.stderr.txt", arena.stderr.encode("utf-8")
            )
            if arena.returncode != 0:
                raise RuntimeError(f"arena verify failed: {arena.stderr}")
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            if (
                receipt["score"]["absolute_determinant"] != str(FRONTIER)
                or receipt["checks"]["exact_bareiss"] != "passed"
            ):
                raise ArithmeticError("arena receipt is not a frontier tie")
            receipt_record = {
                "path": display_path(receipt_path),
                "file_sha256": sha256_file(receipt_path),
                "internal_receipt_sha256": receipt["receipt_sha256"],
            }

        records.append(
            {
                "global_index": global_index,
                "returncode": returncode,
                "process_elapsed_seconds": process_elapsed,
                "factor_found": factor_found,
                "solver": display_path(solver_path),
                "solver_sha256": sha256_file(solver_path),
                "solver_dependencies": solver["dependencies"],
                "classification": classification,
                "arena_receipt": receipt_record,
            }
        )

    feasible = [record for record in records if record["factor_found"]]
    solver_dependencies = {
        (
            record["solver_dependencies"]["numpy"],
            record["solver_dependencies"]["scipy"],
        )
        for record in records
    }
    if len(solver_dependencies) != 1:
        raise ArithmeticError("solver subprocess dependency versions disagree")
    solver_numpy, solver_scipy = next(iter(solver_dependencies))
    manifest = {
        "schema_version": 1,
        "campaign": "frontier-portal-harvest-20260729",
        "claim_boundary": (
            "Novelty is relative to the pinned local baseline only. "
            "Non-feasible bounded runs do not prove empty slices. The fixed "
            "known Gram makes a new row-Gram class impossible in this model; "
            "transpose Gram is still classified explicitly."
        ),
        "method": (
            "exact 1,382-column frontier-Gram shell MILP; theorem-forced "
            "size-six-orbit triples; soft rarity objectives; no raw-support "
            "no-goods; exact Gram/Bareiss/arena and H/HT/Gram classification"
        ),
        "frontier": str(FRONTIER),
        "shard": {
            "index": arguments.shard_index,
            "count": arguments.shard_count,
            "trials_per_shard": arguments.trials_per_shard,
            "seed_base": arguments.seed_base,
        },
        "requested_time_limit_seconds_per_trial": arguments.time_limit,
        "jobs": arguments.jobs,
        "trials_completed": len(records),
        "verified_frontier_ties": len(feasible),
        "new_h_classes": sum(
            record["classification"]["stream_new_h"] for record in feasible
        ),
        "new_ht_classes": sum(
            record["classification"]["stream_new_ht"] for record in feasible
        ),
        "new_gram_classes": sum(
            record["classification"]["stream_new_gram"]
            for record in feasible
        ),
        "baseline": {
            "path": display_path(baseline_path),
            "sha256": sha256_file(baseline_path),
            "h_classes": len(baseline_h),
            "ht_classes": len(baseline_ht),
            "gram_classes": len(baseline_gram),
        },
        "shell": {
            "path": display_path(shell_path),
            "sha256": sha256_file(shell_path),
            "size": len(masks),
            "source": display_path(source),
            "source_sha256": sha256_file(source),
        },
        "known_aligned_factors": known_records,
        "known_support_small_orbit_assertion": {
            "checked_supports": len(known_supports),
            "all_select_exactly_observed_triple": True,
            "observed_triple": [0, 2, 4],
            "observed_masks": sorted(expected_observed_small_masks),
        },
        "solver_runtime": {
            "python_executable": sys.executable,
            "python_version": platform.python_version(),
            "python_full_version": sys.version,
            "numpy": solver_numpy,
            "scipy": solver_scipy,
        },
        "classifier": {
            "python": str(classifier),
            "pynauty": PINNED_PYNAUTY,
        },
        "source": {
            "path": display_path(script),
            "sha256": sha256_file(script),
        },
        "records": records,
        "elapsed_seconds": time.monotonic() - started,
    }
    canonical = json.dumps(manifest, sort_keys=True, separators=(",", ":"))
    manifest["campaign_receipt_sha256"] = sha256_bytes(
        canonical.encode("utf-8")
    )
    write_json(output_dir / "manifest.json", manifest)
    print(
        f"complete ties={len(feasible)} new_HT="
        f"{manifest['new_ht_classes']} new_Gram="
        f"{manifest['new_gram_classes']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "_solve":
        raise SystemExit(solve_hidden(Path(sys.argv[2]).resolve()))
    if len(sys.argv) >= 2 and sys.argv[1] == "_classify":
        raise SystemExit(classify_hidden(Path(sys.argv[2]).resolve()))
    raise SystemExit(main())
