#!/usr/bin/env python3
"""Publish compact live telemetry for the local search-space visualization."""

from __future__ import annotations

import argparse
import json
import math
import os
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DEFAULT_FRONTIER = 2_779_447_296_000_000


def read_events(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []

    events: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                continue
            if isinstance(event, dict):
                events.append(event)
    return events


def first_number(events: list[dict[str, Any]], *keys: str) -> float:
    for event in reversed(events):
        for key in keys:
            value = event.get(key)
            if isinstance(value, (int, float)) and math.isfinite(float(value)):
                return float(value)
    return 0.0


def first_integer(events: list[dict[str, Any]], *keys: str) -> int:
    for event in reversed(events):
        for key in keys:
            value = event.get(key)
            if isinstance(value, int):
                return value
            if isinstance(value, str) and value.isdigit():
                return int(value)
    return 0


def arm_payload(
    arm_dir: Path, index: int, frontier: int, default_budget: float
) -> dict[str, Any]:
    events = read_events(arm_dir / "run.jsonl")
    elapsed = first_number(events, "elapsed_seconds")
    budget = first_number(events, "seconds") or default_budget
    best = first_integer(
        events, "best_absolute_determinant", "absolute_determinant"
    )
    engine = next(
        (
            str(event["engine"])
            for event in reversed(events)
            if isinstance(event.get("engine"), str)
        ),
        "reactive-tabu",
    )
    final_events = {"summary", "complete", "result", "finished"}
    final = any(str(event.get("event", "")).lower() in final_events for event in events)
    progress = min(100.0, 100.0 * elapsed / budget) if budget > 0 else 0.0
    status = "complete" if final or progress >= 100.0 else "active"

    return {
        "id": f"arm{index}",
        "label": f"probe {index + 1}",
        "engine": engine,
        "status": status,
        "elapsed_seconds": round(elapsed, 3),
        "budget_seconds": round(budget, 3),
        "progress_percent": round(progress, 2),
        "best_absolute_determinant": str(best),
        "best_ratio_percent": round(100.0 * best / frontier, 4)
        if frontier
        else 0.0,
    }


def write_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent, text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")
        os.replace(temporary_name, path)
    finally:
        if os.path.exists(temporary_name):
            os.unlink(temporary_name)


def publish(args: argparse.Namespace) -> dict[str, Any]:
    arms = [
        arm_payload(
            args.run_root / f"arm{index}",
            index,
            args.frontier,
            args.budget_seconds,
        )
        for index in range(args.arms)
    ]
    payload = {
        "schema_version": 1,
        "updated_at": datetime.now(timezone.utc).isoformat(),
        "frontier_absolute_determinant": str(args.frontier),
        "arms": arms,
    }
    write_atomic(args.output, payload)
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument(
        "--output", type=Path, default=Path("public/search-progress.json")
    )
    parser.add_argument("--frontier", type=int, default=DEFAULT_FRONTIER)
    parser.add_argument("--arms", type=int, default=4)
    parser.add_argument("--budget-seconds", type=float, default=900.0)
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--interval", type=float, default=2.0)
    parser.add_argument("--settle-seconds", type=float, default=6.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    while True:
        payload = publish(args)
        if not args.watch:
            return 0
        if payload["arms"] and all(
            arm["status"] == "complete" for arm in payload["arms"]
        ):
            time.sleep(max(0.0, args.settle_seconds))
            publish(args)
            return 0
        time.sleep(max(0.25, args.interval))


if __name__ == "__main__":
    raise SystemExit(main())
