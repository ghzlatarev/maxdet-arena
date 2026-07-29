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
from decimal import Decimal, InvalidOperation, localcontext
from math import gcd
from pathlib import Path
from typing import Any


DEFAULT_FRONTIER = 2_779_447_296_000_000
DEFAULT_STALE_SECONDS = 30.0
DEFAULT_HISTORY_POINTS = 72
CORE_SCALE = 1 << 22
COMPLETE_EVENTS = {"summary", "complete", "result", "finished"}
STOPPED_EVENTS = {"stopped"}


def read_events(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []

    events: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                # The search process may be in the middle of appending its
                # final JSONL record while the publisher reads the file.
                continue
            if isinstance(event, dict):
                events.append(event)
    return events


def read_object(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def number_value(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    parsed = float(value)
    return parsed if math.isfinite(parsed) else None


def integer_value(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.isdigit():
        return int(value)
    return None


def first_number(events: list[dict[str, Any]], *keys: str) -> float:
    for event in reversed(events):
        for key in keys:
            value = number_value(event.get(key))
            if value is not None:
                return value
    return 0.0


def first_integer(events: list[dict[str, Any]], *keys: str) -> int:
    for event in reversed(events):
        for key in keys:
            value = integer_value(event.get(key))
            if value is not None:
                return value
    return 0


def object_integer(source: dict[str, Any], key: str, fallback: int = 0) -> int:
    value = integer_value(source.get(key))
    return fallback if value is None else value


def signed_integer_value(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if (
        isinstance(value, str)
        and value
        and (value.isdigit() or (value[0] == "-" and value[1:].isdigit()))
    ):
        return int(value)
    return None


def safe_ratio_decimal(numerator: int, denominator: int) -> float | None:
    if denominator <= 0:
        return None
    try:
        with localcontext() as context:
            context.prec = 30
            value = float(Decimal(numerator) / Decimal(denominator))
    except (InvalidOperation, OverflowError, ValueError):
        return None
    return value if math.isfinite(value) else None


def source_modified_at(paths: list[Path]) -> float | None:
    timestamps: list[float] = []
    for path in paths:
        try:
            timestamps.append(path.stat().st_mtime)
        except OSError:
            continue
    return max(timestamps) if timestamps else None


def iso_timestamp(timestamp: float) -> str:
    return datetime.fromtimestamp(timestamp, timezone.utc).isoformat()


def compact_history(
    events: list[dict[str, Any]], frontier: int, maximum_points: int
) -> list[dict[str, Any]]:
    points: list[dict[str, Any]] = []
    best = 0
    archive_size = 0
    discoveries = 0
    epochs = 0
    for event in events:
        elapsed = number_value(event.get("elapsed_seconds"))
        if elapsed is None:
            continue
        candidate_best = integer_value(
            event.get("best_absolute_determinant")
        )
        if candidate_best is not None:
            best = candidate_best
        candidate_archive = integer_value(event.get("archive_size"))
        if candidate_archive is not None:
            archive_size = candidate_archive
        candidate_discoveries = integer_value(
            event.get("search_added_sketch_discoveries")
        )
        if candidate_discoveries is None:
            candidate_discoveries = integer_value(
                event.get("sketch_discoveries")
            )
        if candidate_discoveries is not None:
            discoveries = candidate_discoveries
        candidate_epochs = integer_value(event.get("epochs_completed"))
        if candidate_epochs is not None:
            epochs = candidate_epochs
        points.append(
            {
                "elapsed_seconds": round(max(0.0, elapsed), 3),
                "best_absolute_determinant": str(best),
                "best_ratio_percent": (
                    round(100.0 * best / frontier, 4) if frontier else 0.0
                ),
                "archive_size": archive_size,
                "sketch_discoveries": discoveries,
                "epochs_completed": epochs,
            }
        )

    if len(points) <= maximum_points:
        return points
    # Uniform index sampling is deterministic and preserves both endpoints.
    # Hopper epoch records are already dense and nearly uniform in time.
    indices = {
        round(index * (len(points) - 1) / (maximum_points - 1))
        for index in range(maximum_points)
    }
    return [points[index] for index in sorted(indices)]


def legacy_arm_payload(
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
    final = any(
        str(event.get("event", "")).lower() in COMPLETE_EVENTS
        for event in events
    )
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
        "best_ratio_percent": (
            round(100.0 * best / frontier, 4) if frontier else 0.0
        ),
    }


def arm_status(
    events: list[dict[str, Any]],
    summary: dict[str, Any],
    modified_at: float | None,
    stale_seconds: float,
    now: float,
) -> str:
    event_names = {
        str(event.get("event", "")).lower() for event in events
    }
    if event_names & COMPLETE_EVENTS or summary.get("complete") is True:
        return "complete"
    if (
        event_names & STOPPED_EVENTS
        or str(summary.get("reason", "")).lower() == "signal"
    ):
        return "stopped"
    if not events:
        return "pending"
    if modified_at is None or now - modified_at > stale_seconds:
        return "stale"
    return "active"


def campaign_arm_payload(
    run_root: Path,
    arm_config: dict[str, Any],
    seed_ids: list[str],
    frontier: int,
    default_budget: float,
    stale_seconds: float,
    maximum_history_points: int,
    now: float,
) -> dict[str, Any]:
    arm_id = str(arm_config["id"])
    arm_dir = run_root / arm_id
    log_path = arm_dir / "run.jsonl"
    summary_path = arm_dir / "summary.json"
    events = read_events(log_path)
    summary = read_object(summary_path)
    modified_at = source_modified_at([log_path, summary_path])
    elapsed = first_number(events, "elapsed_seconds")
    if not elapsed:
        elapsed_number = number_value(summary.get("elapsed_seconds"))
        elapsed = 0.0 if elapsed_number is None else elapsed_number
    budget_number = number_value(arm_config.get("budget_seconds"))
    budget = default_budget if budget_number is None else budget_number
    best = first_integer(
        events, "best_absolute_determinant", "absolute_determinant"
    )
    if not best:
        best = object_integer(summary, "best_absolute_determinant")
    engine = next(
        (
            str(event["engine"])
            for event in reversed(events)
            if isinstance(event.get("engine"), str)
        ),
        str(summary.get("engine", "gram-sketch-basin-hopper")),
    )
    status = arm_status(
        events, summary, modified_at, stale_seconds, now
    )
    progress = min(100.0, 100.0 * elapsed / budget) if budget > 0 else 0.0
    kick_flips = object_integer(arm_config, "kick_flips")
    label = arm_config.get("label")
    if not isinstance(label, str) or not label.strip():
        label = f"k={kick_flips}" if kick_flips else arm_id
    configured_seed_ids = arm_config.get("seed_basin_ids")
    if (
        isinstance(configured_seed_ids, list)
        and all(isinstance(item, str) for item in configured_seed_ids)
    ):
        arm_seed_ids = configured_seed_ids
    else:
        arm_seed_ids = seed_ids
    search_discoveries = integer_value(
        summary.get("search_added_sketch_discoveries")
    )
    if search_discoveries is None:
        search_discoveries = integer_value(summary.get("sketch_discoveries"))
    if search_discoveries is None:
        search_discoveries = first_integer(
            events,
            "search_added_sketch_discoveries",
            "sketch_discoveries",
        )

    return {
        "id": arm_id,
        "label": label,
        "engine": engine,
        "status": status,
        "source_updated_at": (
            iso_timestamp(modified_at) if modified_at is not None else None
        ),
        "elapsed_seconds": round(elapsed, 3),
        "budget_seconds": round(budget, 3),
        "progress_percent": round(progress, 2),
        "best_absolute_determinant": str(best),
        "best_ratio_percent": (
            round(100.0 * best / frontier, 4) if frontier else 0.0
        ),
        "best_core_quotient": object_integer(
            summary,
            "best_core_quotient",
            first_integer(events, "best_core_quotient"),
        ),
        "seed_basin_ids": arm_seed_ids,
        "random_seed": object_integer(arm_config, "random_seed"),
        "kick_flips": kick_flips,
        "epoch_moves": object_integer(arm_config, "epoch_moves"),
        "archive_size": object_integer(
            summary, "archive_size", first_integer(events, "archive_size")
        ),
        "archive_capacity": object_integer(
            summary, "archive_capacity"
        ),
        "sketch_discoveries": search_discoveries,
        "epochs_completed": object_integer(
            summary,
            "epochs_completed",
            first_integer(events, "epochs_completed"),
        ),
        "strict_target_states": object_integer(
            summary,
            "strict_target_states",
            first_integer(events, "strict_target_states"),
        ),
        "history": compact_history(
            events, frontier, maximum_history_points
        ),
    }


def campaign_status(arms: list[dict[str, Any]]) -> str:
    statuses = {str(arm["status"]) for arm in arms}
    if statuses == {"complete"}:
        return "complete"
    if "active" in statuses:
        return "active"
    if "stale" in statuses:
        return "stale"
    if "stopped" in statuses:
        return "stopped"
    return "pending"


def sanitized_seed_basins(campaign: dict[str, Any]) -> list[dict[str, str]]:
    raw_seeds = campaign.get("seed_basins")
    if not isinstance(raw_seeds, list):
        raise ValueError("campaign seed_basins must be a list")
    seeds: list[dict[str, str]] = []
    for index, seed in enumerate(raw_seeds):
        if not isinstance(seed, dict):
            raise ValueError(f"campaign seed_basins[{index}] must be an object")
        seed_id = seed.get("id")
        label = seed.get("label")
        score = integer_value(seed.get("absolute_determinant"))
        basin_key = seed.get("ht_gram_basin_key_sha256")
        if (
            not isinstance(seed_id, str)
            or not isinstance(label, str)
            or score is None
            or not isinstance(basin_key, str)
        ):
            raise ValueError(
                f"campaign seed_basins[{index}] is missing required fields"
            )
        seeds.append(
            {
                "id": seed_id,
                "label": label,
                "absolute_determinant": str(score),
                "gram_basin_key_sha256": basin_key,
            }
        )
    return seeds


def coronal_identity(
    elite: dict[str, Any],
    seeds: list[dict[str, str]],
    frontier: int,
    score: int,
) -> dict[str, str]:
    matching_seeds = [
        seed for seed in seeds if int(seed["absolute_determinant"]) == score
    ]
    identity: dict[str, str] = {}
    if score == frontier:
        identity["identity"] = "frontier"
    elif elite.get("origin") == "seed" and len(matching_seeds) == 1:
        identity["identity"] = "seed"
    if matching_seeds and (
        identity.get("identity") == "frontier"
        or elite.get("origin") == "seed"
    ):
        # A determinant alone identifies the seed only when it is unique in
        # this campaign. Ambiguous ties remain intentionally unlabeled.
        if len(matching_seeds) == 1:
            identity["seed_id"] = matching_seeds[0]["id"]
            identity["seed_label"] = matching_seeds[0]["label"]
    return identity


def coronal_dominates(
    candidate: dict[str, Any], other: dict[str, Any]
) -> bool:
    candidate_det_m = int(candidate["det_m"])
    other_det_m = int(other["det_m"])
    candidate_numerator = int(candidate["kappa_numerator"])
    candidate_denominator = int(candidate["kappa_denominator"])
    other_numerator = int(other["kappa_numerator"])
    other_denominator = int(other["kappa_denominator"])
    kappa_comparison = (
        candidate_numerator * other_denominator
        - other_numerator * candidate_denominator
    )
    return (
        candidate_det_m >= other_det_m
        and kappa_comparison <= 0
        and (candidate_det_m > other_det_m or kappa_comparison < 0)
    )


def sanitized_coronal_projection(
    run_root: Path,
    campaign: dict[str, Any],
    arms: list[dict[str, Any]],
    seeds: list[dict[str, str]],
    frontier: int,
) -> tuple[dict[str, Any] | None, list[Path]]:
    if campaign.get("parent_policy") != "exact_coronal_pareto":
        return None, []

    points: list[dict[str, Any]] = []
    manifest_paths: list[Path] = []
    for arm in arms:
        arm_id = str(arm["id"])
        manifest_path = run_root / arm_id / "archive" / "manifest.json"
        manifest_paths.append(manifest_path)
        manifest = read_object(manifest_path)
        raw_elites = manifest.get("elites")
        if not isinstance(raw_elites, list):
            continue

        for elite in raw_elites:
            if not isinstance(elite, dict):
                continue
            rank = integer_value(elite.get("rank"))
            core_quotient = integer_value(elite.get("core_quotient"))
            absolute_determinant = integer_value(
                elite.get("absolute_determinant")
            )
            raw_coronal_points = elite.get("coronal_points")
            if (
                rank is None
                or rank < 0
                or core_quotient is None
                or core_quotient <= 0
                or absolute_determinant is None
                or absolute_determinant <= 0
                or not isinstance(raw_coronal_points, list)
            ):
                continue

            # A symmetric row/column descriptor commonly emits the same exact
            # point twice. Collapse only duplicates within this retained elite;
            # identical coordinates from different arms remain separate so the
            # public projection retains arm provenance and color.
            seen_orientations: set[tuple[int, int, int]] = set()
            for orientation, raw_point in enumerate(raw_coronal_points):
                if not isinstance(raw_point, dict):
                    continue
                det_m = integer_value(raw_point.get("det_m"))
                numerator = signed_integer_value(
                    raw_point.get("kappa_numerator")
                )
                denominator = integer_value(
                    raw_point.get("kappa_denominator")
                )
                if (
                    det_m is None
                    or det_m <= 0
                    or numerator is None
                    or denominator is None
                    or denominator <= 0
                ):
                    continue
                common = gcd(abs(numerator), denominator)
                normalized = (
                    det_m,
                    numerator // common,
                    denominator // common,
                )
                if normalized in seen_orientations:
                    continue
                seen_orientations.add(normalized)
                decimal = safe_ratio_decimal(numerator, denominator)
                if decimal is None:
                    continue
                point: dict[str, Any] = {
                    "arm_id": arm_id,
                    "rank": rank,
                    "orientation": orientation,
                    "core_quotient": str(core_quotient),
                    "absolute_determinant": str(absolute_determinant),
                    "det_m": str(det_m),
                    "kappa_numerator": str(numerator),
                    "kappa_denominator": str(denominator),
                    "kappa_decimal": decimal,
                }
                point.update(
                    coronal_identity(
                        elite,
                        seeds,
                        frontier,
                        absolute_determinant,
                    )
                )
                points.append(point)

    if not points:
        return None, manifest_paths

    for point in points:
        point["pareto_front"] = not any(
            coronal_dominates(other, point)
            for other in points
            if other is not point
        )

    points.sort(
        key=lambda point: (
            str(point["arm_id"]),
            int(point["rank"]),
            int(point["orientation"]),
        )
    )
    return (
        {
            "kind": "retained_exact_coronal",
            "x_axis": "det_m",
            "x_preference": "larger",
            "y_axis": "kappa",
            "y_preference": "smaller",
            "points": points,
        },
        manifest_paths,
    )


def campaign_payload(
    args: argparse.Namespace, campaign_path: Path
) -> dict[str, Any]:
    campaign = read_object(campaign_path)
    if not campaign:
        raise ValueError(f"cannot read campaign config: {campaign_path}")
    campaign_id = campaign.get("campaign_id")
    if not isinstance(campaign_id, str) or not campaign_id:
        raise ValueError("campaign config has no campaign_id")

    frontier = object_integer(
        campaign, "frontier_absolute_determinant", args.frontier
    )
    frontier_core = object_integer(
        campaign, "known_core_frontier", frontier // CORE_SCALE
    )
    strict_core = object_integer(
        campaign, "strict_core_target", frontier_core + 1
    )
    strict = object_integer(
        campaign,
        "first_strict_absolute_determinant",
        strict_core * CORE_SCALE,
    )
    gate_core = object_integer(campaign, "quotient_gate")
    gate = object_integer(
        campaign, "gate_absolute_determinant", gate_core * CORE_SCALE
    )
    near_core = (frontier_core * 99 + 99) // 100
    near = near_core * CORE_SCALE
    budget_number = number_value(campaign.get("budget_seconds_per_arm"))
    budget = (
        args.budget_seconds if budget_number is None else budget_number
    )
    seeds = sanitized_seed_basins(campaign)
    seed_ids = [seed["id"] for seed in seeds]
    raw_arms = campaign.get("arms")
    if not isinstance(raw_arms, list) or not raw_arms:
        raise ValueError("campaign arms must be a nonempty list")
    if not all(
        isinstance(arm, dict) and isinstance(arm.get("id"), str)
        for arm in raw_arms
    ):
        raise ValueError("every campaign arm must have a string id")

    now = time.time()
    arms = [
        campaign_arm_payload(
            args.run_root,
            arm,
            seed_ids,
            frontier,
            budget,
            args.stale_seconds,
            args.max_history_points,
            now,
        )
        for arm in raw_arms
    ]
    coronal_projection, manifest_paths = sanitized_coronal_projection(
        args.run_root,
        campaign,
        raw_arms,
        seeds,
        frontier,
    )
    source_timestamps = [
        str(arm["source_updated_at"])
        for arm in arms
        if isinstance(arm.get("source_updated_at"), str)
    ]
    manifest_modified_at = source_modified_at(manifest_paths)
    if manifest_modified_at is not None:
        source_timestamps.append(iso_timestamp(manifest_modified_at))
    best = max(
        (int(str(arm["best_absolute_determinant"])) for arm in arms),
        default=0,
    )
    label = campaign.get("label")
    if not isinstance(label, str) or not label.strip():
        label = "Gram-basin hopper pilot"
    claim_boundary = campaign.get("claim_boundary")
    if not isinstance(claim_boundary, str):
        claim_boundary = (
            "Live Gram sketches are invariant search descriptors, not "
            "canonical classes."
        )

    payload: dict[str, Any] = {
        "schema_version": 2,
        "updated_at": datetime.now(timezone.utc).isoformat(),
        "source_updated_at": (
            max(source_timestamps) if source_timestamps else None
        ),
        "campaign": {
            "id": campaign_id,
            "label": label,
            "status": campaign_status(arms),
        },
        "frontier_absolute_determinant": str(frontier),
        "strict_target_absolute_determinant": str(strict),
        "best_absolute_determinant": str(best),
        "best_ratio_percent": (
            round(100.0 * best / frontier, 4) if frontier else 0.0
        ),
        "thresholds": [
            {
                "id": "archive-gate",
                "label": "Archive gate",
                "core_quotient": str(gate_core),
                "absolute_determinant": str(gate),
            },
            {
                "id": "near-frontier",
                "label": "99% corridor",
                "core_quotient": str(near_core),
                "absolute_determinant": str(near),
            },
            {
                "id": "frontier",
                "label": "Frontier tie",
                "core_quotient": str(frontier_core),
                "absolute_determinant": str(frontier),
            },
            {
                "id": "strict",
                "label": "First strict score",
                "core_quotient": str(strict_core),
                "absolute_determinant": str(strict),
            },
        ],
        "seed_basins": seeds,
        "totals": {
            "archive_cells": sum(int(arm["archive_size"]) for arm in arms),
            "sketch_discoveries": sum(
                int(arm["sketch_discoveries"]) for arm in arms
            ),
            "epochs_completed": sum(
                int(arm["epochs_completed"]) for arm in arms
            ),
            "strict_target_states": sum(
                int(arm["strict_target_states"]) for arm in arms
            ),
        },
        "arms": arms,
        "claim_boundary": claim_boundary,
    }
    if coronal_projection is not None:
        payload["coronal_projection"] = coronal_projection
    return payload


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
    campaign_path = args.campaign_config
    if campaign_path is None:
        default_campaign_path = args.run_root / "campaign.json"
        if default_campaign_path.exists():
            campaign_path = default_campaign_path
    if campaign_path is not None:
        payload = campaign_payload(args, campaign_path)
    else:
        arms = [
            legacy_arm_payload(
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
        "--campaign-config",
        type=Path,
        help=(
            "schema-v2 campaign metadata; defaults to RUN_ROOT/campaign.json "
            "when that file exists"
        ),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("public/search-progress.json")
    )
    parser.add_argument("--frontier", type=int, default=DEFAULT_FRONTIER)
    parser.add_argument("--arms", type=int, default=4)
    parser.add_argument("--budget-seconds", type=float, default=900.0)
    parser.add_argument(
        "--stale-seconds", type=float, default=DEFAULT_STALE_SECONDS
    )
    parser.add_argument(
        "--max-history-points",
        type=int,
        default=DEFAULT_HISTORY_POINTS,
    )
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--interval", type=float, default=2.0)
    parser.add_argument("--settle-seconds", type=float, default=6.0)
    args = parser.parse_args()
    if args.arms <= 0:
        parser.error("--arms must be positive")
    if args.budget_seconds <= 0:
        parser.error("--budget-seconds must be positive")
    if args.stale_seconds <= 0:
        parser.error("--stale-seconds must be positive")
    if args.max_history_points < 2:
        parser.error("--max-history-points must be at least 2")
    return args


def main() -> int:
    args = parse_args()
    while True:
        payload = publish(args)
        if not args.watch:
            return 0
        terminal_statuses = {"complete", "stopped"}
        if payload["arms"] and all(
            str(arm["status"]) in terminal_statuses
            for arm in payload["arms"]
        ):
            time.sleep(max(0.0, args.settle_seconds))
            publish(args)
            return 0
        time.sleep(max(0.25, args.interval))


if __name__ == "__main__":
    raise SystemExit(main())
