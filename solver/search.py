#!/usr/bin/env python3
"""Small dependency-free starter search.

This is deliberately simple and agent-editable. The trusted verifier never
imports or executes it. Better solvers may use any language or dependency; they
only need to emit candidate/matrix.txt.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
import time
from pathlib import Path
from typing import TextIO

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from maxdet.contract import load_contract, load_matrix
from maxdet.exact import bareiss_determinant, matrix_text
from maxdet.frontier import effective_frontier

ORDER = 23


def random_matrix(randomizer: random.Random) -> list[list[int]]:
    return [
        [randomizer.choice((-1, 1)) for _ in range(ORDER)]
        for _ in range(ORDER)
    ]


def read_matrix(path: Path) -> list[list[int]]:
    contract = load_contract(REPOSITORY_ROOT / "challenge.json")
    _, matrix = load_matrix(path, contract)
    return matrix


def score(matrix: list[list[int]]) -> int:
    return abs(bareiss_determinant(matrix))


def first_improvement(
    matrix: list[list[int]],
    current_score: int,
    randomizer: random.Random,
    deadline: float,
) -> tuple[int, bool]:
    coordinates = [(row, column) for row in range(ORDER) for column in range(ORDER)]
    randomizer.shuffle(coordinates)
    for row, column in coordinates:
        if time.monotonic() >= deadline:
            return current_score, False
        matrix[row][column] *= -1
        changed_score = score(matrix)
        if changed_score > current_score:
            return changed_score, True
        matrix[row][column] *= -1
    return current_score, False


def atomic_write(path: Path, matrix: list[list[int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(matrix_text(matrix), encoding="ascii")
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=60.0)
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("candidate/matrix.txt"),
    )
    parser.add_argument(
        "--start",
        type=Path,
        default=None,
        help="start matrix; defaults to the better of candidate and trusted frontier",
    )
    parser.add_argument(
        "--log",
        type=Path,
        default=Path("runs/starter-search.jsonl"),
    )
    parser.add_argument(
        "--report-interval",
        type=float,
        default=1.0,
        help="minimum seconds between console improvement reports; 0 reports all",
    )
    parser.add_argument(
        "--heartbeat-seconds",
        type=float,
        default=30.0,
        help="seconds between JSONL progress events; 0 disables heartbeats",
    )
    return parser.parse_args()


def write_event(log: TextIO, event: dict[str, object]) -> None:
    log.write(json.dumps(event, sort_keys=True) + "\n")
    log.flush()


def main() -> int:
    args = parse_args()
    if not math.isfinite(args.seconds) or args.seconds <= 0:
        raise SystemExit("--seconds must be finite and positive")
    if not math.isfinite(args.report_interval) or args.report_interval < 0:
        raise SystemExit("--report-interval must be finite and non-negative")
    if not math.isfinite(args.heartbeat_seconds) or args.heartbeat_seconds < 0:
        raise SystemExit("--heartbeat-seconds must be finite and non-negative")
    randomizer = random.Random(args.seed)
    started = time.monotonic()
    deadline = started + args.seconds
    run_started_unix_ns = time.time_ns()
    args.log.parent.mkdir(parents=True, exist_ok=True)

    start_path = args.start
    if start_path is None:
        candidate_path = REPOSITORY_ROOT / "candidate" / "matrix.txt"
        candidate_matrix = read_matrix(candidate_path)
        candidate_score = score(candidate_matrix)
        contract = load_contract(REPOSITORY_ROOT / "challenge.json")
        frontier = effective_frontier(REPOSITORY_ROOT, contract)
        if candidate_score >= frontier.absolute_determinant:
            start_path = candidate_path
        else:
            start_path = REPOSITORY_ROOT / frontier.source / "matrix.txt"

    incumbent = read_matrix(start_path)
    incumbent_score = score(incumbent)
    best = [row[:] for row in incumbent]
    best_score = incumbent_score
    try:
        start_label = start_path.resolve().relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError:
        start_label = start_path.name
    restart = 0
    total_improvements = 0
    last_report = -math.inf
    next_heartbeat = (
        started + args.heartbeat_seconds
        if args.heartbeat_seconds
        else math.inf
    )

    with args.log.open("a", encoding="utf-8") as log:
        write_event(
            log,
            {
                "event": "start",
                "seed": args.seed,
                "run_started_unix_ns": run_started_unix_ns,
                "start_matrix": start_label,
                "absolute_determinant": str(best_score),
            },
        )
        while time.monotonic() < deadline:
            matrix = [row[:] for row in incumbent] if restart == 0 else random_matrix(randomizer)
            current_score = score(matrix)
            improvements = 0
            while time.monotonic() < deadline:
                current_score, changed = first_improvement(
                    matrix,
                    current_score,
                    randomizer,
                    deadline,
                )
                if not changed:
                    break
                improvements += 1
                total_improvements += 1
                if current_score > best_score:
                    best_score = current_score
                    best = [row[:] for row in matrix]
                    atomic_write(args.output, best)
                    event = {
                        "event": "new_best",
                        "seed": args.seed,
                        "run_started_unix_ns": run_started_unix_ns,
                        "restart": restart,
                        "improvements": improvements,
                        "total_improvements": total_improvements,
                        "absolute_determinant": str(best_score),
                    }
                    write_event(log, event)
                    now = time.monotonic()
                    if now - last_report >= args.report_interval:
                        print(f"new best |det|={best_score}", flush=True)
                        last_report = now
                now = time.monotonic()
                if now >= next_heartbeat:
                    write_event(
                        log,
                        {
                            "event": "heartbeat",
                            "seed": args.seed,
                            "run_started_unix_ns": run_started_unix_ns,
                            "restart": restart,
                            "total_improvements": total_improvements,
                            "absolute_determinant": str(best_score),
                            "elapsed_seconds": round(now - started, 3),
                        },
                    )
                    next_heartbeat = now + args.heartbeat_seconds
            restart += 1
            now = time.monotonic()
            if now >= next_heartbeat:
                write_event(
                    log,
                    {
                        "event": "heartbeat",
                        "seed": args.seed,
                        "run_started_unix_ns": run_started_unix_ns,
                        "restart": restart,
                        "total_improvements": total_improvements,
                        "absolute_determinant": str(best_score),
                        "elapsed_seconds": round(now - started, 3),
                    },
                )
                next_heartbeat = now + args.heartbeat_seconds

    atomic_write(args.output, best)
    elapsed = time.monotonic() - started
    with args.log.open("a", encoding="utf-8") as log:
        write_event(
            log,
            {
                "event": "finished",
                "seed": args.seed,
                "run_started_unix_ns": run_started_unix_ns,
                "restarts": restart,
                "total_improvements": total_improvements,
                "absolute_determinant": str(best_score),
                "elapsed_seconds": round(elapsed, 3),
            },
        )
    print(f"finished |det|={best_score} restarts={restart}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
