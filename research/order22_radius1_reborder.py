#!/usr/bin/env python3
"""Exhaust radius-one core perturbations and globally reoptimize each border.

For every entry of each supplied 22x22 core, flip that entry and invoke the
exact order22_border engine, which exhausts all 2^21 normalized border
columns and chooses the optimal border row analytically.  Every retained
per-core score is independently recomputed with exact Bareiss elimination.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from order22_border_structure import bareiss_determinant, read_sign_matrix


ORDER = 22
ASSIGNMENTS_PER_CORE = 1 << (ORDER - 1)


def parse_labeled_path(text: str) -> tuple[str, Path]:
    if "=" not in text:
        raise argparse.ArgumentTypeError("base must be LABEL=PATH")
    label, raw_path = text.split("=", 1)
    if not label or not raw_path:
        raise argparse.ArgumentTypeError("base must be LABEL=PATH")
    return label, Path(raw_path)


def matrix_payload(matrix: list[list[int]]) -> bytes:
    return "".join(
        " ".join(map(str, row)) + "\n" for row in matrix
    ).encode("ascii")


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_bytes(payload)
    temporary.replace(path)


def read_events(path: Path) -> list[dict[str, object]]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def increment(histogram: dict[str, int], value: int) -> None:
    key = str(value)
    histogram[key] = histogram.get(key, 0) + 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base",
        action="append",
        required=True,
        type=parse_labeled_path,
        metavar="LABEL=PATH",
    )
    parser.add_argument("--engine", required=True, type=Path)
    parser.add_argument("--arena", type=Path, default=Path("./arena"))
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--progress-every", type=int, default=25)
    arguments = parser.parse_args()

    labels = [label for label, _ in arguments.base]
    if len(labels) != len(set(labels)):
        parser.error("base labels must be unique")
    if arguments.progress_every < 0:
        parser.error("--progress-every must be non-negative")

    engine = arguments.engine.resolve()
    arena = arguments.arena.resolve()
    output_dir = arguments.output_dir.resolve()
    if not engine.is_file():
        raise FileNotFoundError(engine)
    if not arena.is_file():
        raise FileNotFoundError(arena)
    output_dir.mkdir(parents=True, exist_ok=True)
    results_path = output_dir / "results.jsonl"
    report_path = output_dir / "report.json"
    if results_path.exists() or report_path.exists():
        raise FileExistsError("refusing to overwrite an existing campaign")

    bases: list[tuple[str, Path, bytes, list[list[int]]]] = []
    for label, path in arguments.base:
        resolved = path.resolve()
        payload, matrix = read_sign_matrix(resolved, ORDER)
        determinant = bareiss_determinant(matrix)
        if determinant == 0:
            raise ArithmeticError(f"{resolved}: singular base core")
        bases.append((label, resolved, payload, matrix))

    scratch = output_dir / "scratch"
    scratch.mkdir()
    scratch_core = scratch / "core.matrix.txt"
    scratch_best = scratch / "best.matrix.txt"
    scratch_log = scratch / "border.jsonl"

    started = time.monotonic()
    completed = 0
    expected = len(bases) * ORDER * ORDER
    global_best = -1
    global_records: list[dict[str, object]] = []
    score_histogram: dict[str, int] = {}
    core_determinant_histogram: dict[str, int] = {}

    try:
        with results_path.open("x", encoding="utf-8") as results:
            for label, base_path, base_payload, base in bases:
                base_sha256 = hashlib.sha256(base_payload).hexdigest()
                for row in range(ORDER):
                    for column in range(ORDER):
                        core = [source_row[:] for source_row in base]
                        core[row][column] = -core[row][column]
                        core_payload = matrix_payload(core)
                        atomic_write(scratch_core, core_payload)
                        process = subprocess.run(
                            (
                                str(engine),
                                "--start",
                                str(scratch_core),
                                "--output",
                                str(scratch_best),
                                "--log",
                                str(scratch_log),
                            ),
                            check=False,
                            capture_output=True,
                            text=True,
                        )
                        if process.returncode:
                            raise RuntimeError(
                                f"{label} ({row + 1},{column + 1}) failed: "
                                f"{process.stderr.strip()}"
                            )
                        events = read_events(scratch_log)
                        if (
                            not events
                            or events[0].get("event") != "start"
                            or events[-1].get("event") != "finished"
                        ):
                            raise RuntimeError("malformed border-engine log")
                        first = events[0]
                        finished = events[-1]
                        if (
                            finished.get("complete") is not True
                            or finished.get("assignments_completed")
                            != ASSIGNMENTS_PER_CORE
                        ):
                            raise RuntimeError("incomplete border enumeration")

                        best_payload, best_matrix = read_sign_matrix(
                            scratch_best, ORDER + 1
                        )
                        exact_score = abs(bareiss_determinant(best_matrix))
                        logged_score = int(
                            str(finished["absolute_determinant"])
                        )
                        if exact_score != logged_score:
                            raise ArithmeticError(
                                "border score failed independent Bareiss check"
                            )
                        exact_core_determinant = abs(
                            bareiss_determinant(core)
                        )
                        logged_core_determinant = int(
                            str(first["absolute_core_determinant"])
                        )
                        if exact_core_determinant != logged_core_determinant:
                            raise ArithmeticError(
                                "core determinant failed independent check"
                            )

                        record: dict[str, object] = {
                            "assignments_completed": ASSIGNMENTS_PER_CORE,
                            "base": label,
                            "base_path": str(base_path),
                            "base_sha256": base_sha256,
                            "best_absolute_determinant": str(exact_score),
                            "best_border_columns_up_to_global_sign": finished[
                                "best_border_columns_up_to_global_sign"
                            ],
                            "best_matrix_sha256": hashlib.sha256(
                                best_payload
                            ).hexdigest(),
                            "column_1_based": column + 1,
                            "core_absolute_determinant": str(
                                exact_core_determinant
                            ),
                            "core_sha256": hashlib.sha256(
                                core_payload
                            ).hexdigest(),
                            "engine_elapsed_seconds": finished[
                                "elapsed_seconds"
                            ],
                            "exact_bareiss_checked": True,
                            "row_1_based": row + 1,
                        }
                        results.write(
                            json.dumps(
                                record,
                                sort_keys=True,
                                separators=(",", ":"),
                            )
                            + "\n"
                        )
                        results.flush()
                        completed += 1
                        increment(score_histogram, exact_score)
                        increment(
                            core_determinant_histogram,
                            exact_core_determinant,
                        )
                        if exact_score > global_best:
                            global_best = exact_score
                            global_records = [record]
                            atomic_write(
                                output_dir / "global-best-core.matrix.txt",
                                core_payload,
                            )
                            atomic_write(
                                output_dir / "global-best.matrix.txt",
                                best_payload,
                            )
                            shutil.copyfile(
                                scratch_log,
                                output_dir / "global-best-border.jsonl",
                            )
                        elif exact_score == global_best:
                            global_records.append(record)

                        if arguments.progress_every and (
                            completed % arguments.progress_every == 0
                            or completed == expected
                        ):
                            print(
                                f"completed={completed}/{expected} "
                                f"global_best={global_best} "
                                f"elapsed_seconds="
                                f"{time.monotonic() - started:.3f}",
                                flush=True,
                            )

        receipt_path = output_dir / "global-best.receipt.json"
        verification = subprocess.run(
            (
                str(arena),
                "verify",
                "--quiet",
                "--json",
                str(receipt_path),
                str(output_dir / "global-best.matrix.txt"),
            ),
            check=False,
            capture_output=True,
            text=True,
        )
        if verification.returncode:
            raise RuntimeError(
                f"arena verification failed: {verification.stderr.strip()}"
            )
        verified_score = int(verification.stdout.strip())
        if verified_score != global_best:
            raise ArithmeticError("arena score disagrees with campaign best")
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))

        elapsed = time.monotonic() - started
        report = {
            "arithmetic": "exact-integer-only",
            "assignments_completed": completed * ASSIGNMENTS_PER_CORE,
            "assignments_per_core": ASSIGNMENTS_PER_CORE,
            "bases": [
                {
                    "label": label,
                    "path": str(path),
                    "raw_sha256": hashlib.sha256(payload).hexdigest(),
                    "absolute_determinant": str(
                        abs(bareiss_determinant(matrix))
                    ),
                }
                for label, path, payload, matrix in bases
            ],
            "complete": completed == expected,
            "core_determinant_histogram": core_determinant_histogram,
            "cores_completed": completed,
            "cores_expected": expected,
            "elapsed_seconds": elapsed,
            "engine": "order22-radius1-exact-reborder-v1",
            "border_engine": {
                "path": str(engine),
                "sha256": hashlib.sha256(engine.read_bytes()).hexdigest(),
            },
            "global_best_absolute_determinant": str(global_best),
            "global_best_records": global_records,
            "global_best_receipt": {
                "path": str(receipt_path),
                "receipt_sha256": receipt["receipt_sha256"],
                "matrix_raw_sha256": receipt["matrix"]["raw_sha256"],
            },
            "radius": 1,
            "score_histogram": score_histogram,
        }
        atomic_write(
            report_path,
            (json.dumps(report, indent=2, sort_keys=True) + "\n").encode(
                "utf-8"
            ),
        )
        print(
            f"complete cores={completed} "
            f"assignments={completed * ASSIGNMENTS_PER_CORE} "
            f"global_best={global_best} elapsed_seconds={elapsed:.3f}",
            flush=True,
        )
    finally:
        if scratch.exists():
            shutil.rmtree(scratch)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
