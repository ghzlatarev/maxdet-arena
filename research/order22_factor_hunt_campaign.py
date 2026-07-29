#!/usr/bin/env python3
"""Run bounded alternate-factor hunts against both exact order-22 Grams.

Each run excludes every currently aligned exact support, uses a distinct
known factor as both a complete CP-SAT hint and overlap-minimization anchor,
and receives a distinct randomized-search seed.  Any feasible factor is
immediately checked by exact Gram equality and Bareiss determinant and
classified under H and HT equivalence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

import pynauty

from h_equivalence_audit import determinant, h_certificate, transpose


ROOT = Path(__file__).resolve().parents[1]
ORDER = 22
TARGET_DETERMINANT = 409_600_000_000_000
PINNED_PYNAUTY = "2.8.8.1"
DEFAULT_ALIGNMENT = (
    ROOT
    / "runs/direct-search/order22-gram-aligned-factors-20260729/manifest.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / (
        "runs/direct-search/"
        "order22-gram-factor-cpsat-wave3-20260729"
    )
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


def resolve_path(value: str, base: Path = ROOT) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


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


def gram(matrix: Matrix) -> list[list[int]]:
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


def certificate_sha256(matrix: Matrix) -> str:
    return sha256_bytes(h_certificate(matrix))


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


def select_diverse_hints(
    records: list[dict[str, Any]],
    count: int,
) -> list[dict[str, Any]]:
    if count >= len(records):
        return sorted(records, key=lambda record: record["h_certificate_sha256"])
    supports = {
        str(record["h_certificate_sha256"]): set(
            int(mask)
            for mask in record["normalized_column_masks_decimal"]
        )
        for record in records
    }
    by_h = {
        str(record["h_certificate_sha256"]): record
        for record in records
    }
    certificates = sorted(by_h)
    if count == 1:
        return [by_h[certificates[0]]]
    first, second = min(
        (
            (
                left,
                right,
                len(supports[left] & supports[right]),
            )
            for index, left in enumerate(certificates)
            for right in certificates[index + 1 :]
        ),
        key=lambda item: (item[2], item[0], item[1]),
    )[:2]
    selected = [first, second]
    while len(selected) < count:
        candidate = min(
            (
                certificate
                for certificate in certificates
                if certificate not in selected
            ),
            key=lambda certificate: (
                max(
                    len(supports[certificate] & supports[other])
                    for other in selected
                ),
                sum(
                    len(supports[certificate] & supports[other])
                    for other in selected
                ),
                certificate,
            ),
        )
        selected.append(candidate)
    return [by_h[certificate] for certificate in selected]


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def classify_factor(
    path: Path,
    target_gram: list[list[int]],
    known_factors: list[dict[str, Any]],
) -> dict[str, Any]:
    matrix = read_matrix(path)
    checked = determinant(matrix)
    if abs(checked) != TARGET_DETERMINANT:
        raise ArithmeticError(f"{path}: unexpected determinant {checked}")
    if gram(matrix) != target_gram:
        raise ArithmeticError(f"{path}: exact target Gram mismatch")
    direct = certificate_sha256(matrix)
    transposed = certificate_sha256(transpose(matrix))
    ht = min(direct, transposed)
    support = set(normalized_column_masks(matrix))
    known_by_h = {
        str(record["h_certificate_sha256"]): record
        for record in known_factors
    }
    exact_support_matches = sorted(
        str(record["h_certificate_sha256"])
        for record in known_factors
        if support
        == set(
            int(mask)
            for mask in record["normalized_column_masks_decimal"]
        )
    )
    return {
        "path": display_path(path),
        "raw_sha256": sha256_file(path),
        "determinant": checked,
        "h_certificate_sha256": direct,
        "transpose_h_certificate_sha256": transposed,
        "ht_certificate_sha256": ht,
        "known_h_class": direct in known_by_h,
        "known_ht_class": ht
        in {
            str(record["ht_certificate_sha256"])
            for record in known_factors
        },
        "exact_known_support_matches_h": exact_support_matches,
        "normalized_column_masks_decimal": sorted(support),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--alignment-manifest",
        type=Path,
        default=DEFAULT_ALIGNMENT,
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--runs-per-target", type=int, default=3)
    parser.add_argument("--time-limit", type=float, default=60.0)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--seed-base", type=int, default=592_000)
    arguments = parser.parse_args()
    if getattr(pynauty, "__version__", None) != PINNED_PYNAUTY:
        raise RuntimeError(f"requires pynauty=={PINNED_PYNAUTY}")
    if (
        arguments.runs_per_target <= 0
        or arguments.workers <= 0
        or arguments.time_limit <= 0
        or not math.isfinite(arguments.time_limit)
    ):
        parser.error("run count, workers, and time limit must be positive")

    alignment_path = arguments.alignment_manifest.expanduser().resolve()
    output_dir = arguments.output_dir.expanduser().resolve()
    if output_dir.exists():
        raise FileExistsError(output_dir)
    output_dir.mkdir(parents=True)
    alignment = json.loads(alignment_path.read_text(encoding="utf-8"))
    if (
        alignment.get("h_class_count") != 30
        or alignment.get("ht_class_count") != 26
        or alignment.get("target_gram_class_count") != 2
    ):
        raise ValueError("alignment manifest is not the 30-H / 26-HT package")
    factors_by_target: dict[str, list[dict[str, Any]]] = {}
    for record in alignment["h_factors"]:
        factors_by_target.setdefault(str(record["target"]), []).append(record)
    targets = {
        str(record["name"]): record
        for record in alignment["targets"]
    }
    if set(factors_by_target) != set(targets):
        raise ArithmeticError("target/factor allocation mismatch")

    solver_source = ROOT / "research/order22_gram_factor_cpsat.py"
    run_specs = []
    for target_name in sorted(targets):
        known_factors = factors_by_target[target_name]
        hints = select_diverse_hints(
            known_factors,
            min(arguments.runs_per_target, len(known_factors)),
        )
        target = targets[target_name]
        shell_report = resolve_path(target["complete_shell_report"])
        target_factor = resolve_path(target["source_factor_path"])
        exclusion_list = resolve_path(
            str(
                Path(alignment_path).parent
                / target["cpsat_known_factor_list_path"]
            ),
            Path("/"),
        )
        target_gram = gram(read_matrix(target_factor))
        for target_run_index, hint in enumerate(hints):
            run_specs.append(
                {
                    "target": target_name,
                    "target_run_index": target_run_index,
                    "hint": hint,
                    "shell_report": shell_report,
                    "target_factor": target_factor,
                    "exclusion_list": exclusion_list,
                    "target_gram": target_gram,
                    "known_factors": known_factors,
                }
            )

    started = time.monotonic()
    run_records = []
    feasible_records = []
    for run_index, spec in enumerate(run_specs):
        target_name = str(spec["target"])
        hint = spec["hint"]
        h = str(hint["h_certificate_sha256"])
        hint_path = resolve_path(
            str(alignment_path.parent / hint["aligned_factor_path"]),
            Path("/"),
        )
        run_dir = (
            output_dir
            / f"run-{run_index:02d}-{target_name}-hint-{h[:12]}"
        )
        run_dir.mkdir()
        output_matrix = run_dir / "factor.matrix.txt"
        metadata_path = run_dir / "solver.json"
        stdout_path = run_dir / "stdout.txt"
        stderr_path = run_dir / "stderr.txt"
        seed = arguments.seed_base + run_index
        command = [
            sys.executable,
            str(solver_source),
            "--shell-report",
            str(spec["shell_report"]),
            "--factor",
            str(spec["target_factor"]),
            "--hint-factor",
            str(hint_path),
            "--minimize-overlap-with",
            str(hint_path),
            "--exclude-factor-list",
            str(spec["exclusion_list"]),
            "--output",
            str(output_matrix),
            "--metadata",
            str(metadata_path),
            "--time-limit",
            str(arguments.time_limit),
            "--workers",
            str(arguments.workers),
            "--seed",
            str(seed),
        ]
        print(
            f"starting {run_index + 1}/{len(run_specs)} "
            f"target={target_name} hint={h[:12]} seed={seed}",
            flush=True,
        )
        run_started = time.monotonic()
        completed = subprocess.run(
            command,
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        stdout_path.write_text(completed.stdout, encoding="utf-8")
        stderr_path.write_text(completed.stderr, encoding="utf-8")
        if completed.returncode not in (0, 1):
            raise RuntimeError(
                f"solver failed with {completed.returncode}: "
                f"{completed.stderr.strip()}"
            )
        if not metadata_path.is_file():
            raise RuntimeError("solver did not write metadata")
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        feasible = bool(metadata["solver"]["feasible"])
        if feasible != output_matrix.is_file():
            raise ArithmeticError("solver feasible/output binding failed")
        classification = None
        if feasible:
            classification = classify_factor(
                output_matrix,
                spec["target_gram"],
                spec["known_factors"],
            )
            write_json(run_dir / "classification.json", classification)
            feasible_records.append(classification)
            print(
                f"FEASIBLE target={target_name} "
                f"H={classification['h_certificate_sha256']} "
                f"known_H={classification['known_h_class']} "
                f"known_HT={classification['known_ht_class']}",
                flush=True,
            )
        run_records.append(
            {
                "run_index": run_index,
                "target": target_name,
                "seed": seed,
                "time_limit_seconds": arguments.time_limit,
                "workers": arguments.workers,
                "hint_h_certificate_sha256": h,
                "hint_factor": display_path(hint_path),
                "exclude_factor_list": display_path(
                    spec["exclusion_list"]
                ),
                "known_supports_excluded": len(
                    spec["known_factors"]
                ),
                "command": command,
                "returncode": completed.returncode,
                "elapsed_seconds": time.monotonic() - run_started,
                "solver_status": metadata["solver"]["status"],
                "solver_feasible": feasible,
                "solver_wall_time_seconds": metadata["solver"][
                    "wall_time_seconds"
                ],
                "solver_branches": metadata["solver"]["branches"],
                "solver_conflicts": metadata["solver"]["conflicts"],
                "solver_exact_overlap": metadata["solver"][
                    "exact_overlap"
                ],
                "metadata_path": display_path(metadata_path),
                "metadata_sha256": sha256_file(metadata_path),
                "classification": classification,
            }
        )
        print(
            f"finished target={target_name} status="
            f"{metadata['solver']['status']} "
            f"elapsed={run_records[-1]['elapsed_seconds']:.1f}s",
            flush=True,
        )

    manifest = {
        "schema_version": 1,
        "engine": "order22-gram-factor-cpsat-wave2-campaign-v1",
        "claim": (
            "Each run excluded every known exact aligned support for its "
            "target Gram and used a distinct diversity-selected known hint. "
            "Any feasible output was independently checked and H/HT "
            "classified exactly."
        ),
        "claim_boundary": (
            "UNKNOWN is a bounded CP-SAT result, not an infeasibility proof. "
            "Exact support exclusions do not eliminate full signed-Gram "
            "automorphism orbits."
        ),
        "alignment_manifest": {
            "path": display_path(alignment_path),
            "sha256": sha256_file(alignment_path),
            "h_class_count": alignment["h_class_count"],
            "ht_class_count": alignment["ht_class_count"],
        },
        "solver": {
            "path": display_path(solver_source),
            "sha256": sha256_file(solver_source),
            "python": sys.executable,
        },
        "runs_per_target": arguments.runs_per_target,
        "requested_time_limit_seconds_per_run": arguments.time_limit,
        "workers_per_run": arguments.workers,
        "runs_completed": len(run_records),
        "feasible_factor_count": len(feasible_records),
        "novel_h_factor_count": sum(
            record["known_h_class"] is False
            for record in feasible_records
        ),
        "novel_ht_factor_count": sum(
            record["known_ht_class"] is False
            for record in feasible_records
        ),
        "runs": run_records,
        "elapsed_seconds": time.monotonic() - started,
        "dependencies": {
            "pynauty": getattr(pynauty, "__version__", "unknown"),
        },
        "source": {
            "path": display_path(Path(__file__)),
            "sha256": sha256_file(Path(__file__)),
        },
    }
    write_json(output_dir / "manifest.json", manifest)
    print(
        f"campaign complete runs={len(run_records)} "
        f"feasible={len(feasible_records)} "
        f"novel_H={manifest['novel_h_factor_count']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
