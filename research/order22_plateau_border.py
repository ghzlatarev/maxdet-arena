#!/usr/bin/env python3
"""Exhaust fixed borders for one representative per order-22 HT class."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any


ORDER = 22
ASSIGNMENTS_PER_CORE = 1 << (ORDER - 1)


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
        raise ValueError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def atomic_write(path: Path, contents: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
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


def resolve_bound_path(report_path: Path, value: object) -> Path:
    if not isinstance(value, str):
        raise ValueError("class path is missing")
    path = Path(value)
    if not path.is_absolute():
        path = report_path.parent / path
    return path.expanduser().absolute()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--classes-report", required=True, type=Path)
    parser.add_argument("--engine", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    arguments = parser.parse_args()

    report_path = arguments.classes_report.expanduser().absolute()
    engine = arguments.engine.expanduser().absolute()
    output_dir = arguments.output_dir.expanduser().absolute()
    if not engine.is_file():
        parser.error("--engine must be a regular file")
    if output_dir.exists() and any(output_dir.iterdir()):
        parser.error("--output-dir must be absent or empty")
    output_dir.mkdir(parents=True, exist_ok=True)

    classes_report = read_json(report_path)
    if (
        classes_report.get("engine")
        != "order22-maxdet-plateau-h-class-bfs-v1"
        or classes_report.get("order") != ORDER
        or classes_report.get("complete_seeded_component_closure") is not True
    ):
        raise ValueError("plateau class report is incomplete or incompatible")
    raw_classes = classes_report.get("classes")
    if not isinstance(raw_classes, list):
        raise ValueError("plateau report is missing classes")

    representatives: dict[str, dict[str, Any]] = {}
    for record in raw_classes:
        if not isinstance(record, dict):
            raise ValueError("invalid plateau class record")
        certificate = record.get("ht_certificate_sha256")
        index = record.get("index")
        if not isinstance(certificate, str) or type(index) is not int:
            raise ValueError("invalid HT certificate or class index")
        existing = representatives.get(certificate)
        if existing is None or index < existing["index"]:
            representatives[certificate] = record
    selected = sorted(representatives.values(), key=lambda item: item["index"])
    expected_count = classes_report.get("ht_class_count")
    if len(selected) != expected_count:
        raise ValueError("HT representative count does not match report")

    started = time.monotonic()
    results: list[dict[str, Any]] = []
    global_best = -1
    global_best_result: dict[str, Any] | None = None
    for position, class_record in enumerate(selected):
        index = int(class_record["index"])
        h_certificate = str(class_record["h_certificate_sha256"])
        ht_certificate = str(class_record["ht_certificate_sha256"])
        core_path = resolve_bound_path(report_path, class_record.get("path"))
        expected_core_sha = class_record.get("raw_sha256")
        actual_core_sha = sha256(core_path)
        if actual_core_sha != expected_core_sha:
            raise ValueError(f"plateau core SHA-256 changed: {core_path}")

        destination = (
            output_dir / f"class-{index:03d}-{h_certificate[:12]}"
        )
        destination.mkdir()
        best_path = destination / "best.matrix.txt"
        log_path = destination / "run.jsonl"
        completed = subprocess.run(
            (
                str(engine),
                "--start",
                str(core_path),
                "--output",
                str(best_path),
                "--log",
                str(log_path),
            ),
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode:
            raise RuntimeError(
                f"border engine failed for class {index}: "
                f"{completed.stderr.strip()}"
            )
        events = [
            json.loads(line)
            for line in log_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        if (
            len(events) < 2
            or events[0].get("event") != "start"
            or events[-1].get("event") != "finished"
            or events[-1].get("complete") is not True
            or events[-1].get("assignments_completed")
            != ASSIGNMENTS_PER_CORE
        ):
            raise RuntimeError(f"incomplete border log for class {index}")
        score = int(events[-1]["absolute_determinant"])
        result = {
            "class_index": index,
            "h_certificate_sha256": h_certificate,
            "ht_certificate_sha256": ht_certificate,
            "core": {
                "path": str(core_path),
                "sha256": actual_core_sha,
                "absolute_determinant": str(
                    abs(int(class_record["determinant"]))
                ),
            },
            "best": {
                "path": str(best_path),
                "sha256": sha256(best_path),
                "absolute_determinant": str(score),
            },
            "log": {
                "path": str(log_path),
                "sha256": sha256(log_path),
            },
            "assignments_completed": ASSIGNMENTS_PER_CORE,
            "best_border_columns_up_to_global_sign": events[-1][
                "best_border_columns_up_to_global_sign"
            ],
            "engine_elapsed_seconds": events[-1]["elapsed_seconds"],
        }
        results.append(result)
        if score > global_best:
            global_best = score
            global_best_result = result
            shutil.copyfile(best_path, output_dir / ".global-best.matrix.tmp")
            os.replace(
                output_dir / ".global-best.matrix.tmp",
                output_dir / "global-best.matrix.txt",
            )
        print(
            f"completed={position + 1}/{len(selected)} "
            f"class={index} score={score} global_best={global_best}",
            flush=True,
        )

    if global_best_result is None:
        raise ArithmeticError("no HT class was screened")
    manifest = {
        "schema_version": 1,
        "engine": "order22-plateau-border-campaign-v1",
        "claim": (
            "All normalized border columns were scored exactly for one "
            "representative of every HT class in the bound seeded plateau "
            "component report."
        ),
        "claim_boundary": (
            "The class report closes only components reached from its seeds; "
            "this is not a global order-22 or order-23 optimality claim."
        ),
        "order22_classes_report": {
            "path": str(report_path),
            "sha256": sha256(report_path),
        },
        "border_engine": {
            "path": str(engine),
            "sha256": sha256(engine),
        },
        "ht_classes_completed": len(results),
        "ht_classes_expected": expected_count,
        "assignments_per_core": ASSIGNMENTS_PER_CORE,
        "assignments_completed": len(results) * ASSIGNMENTS_PER_CORE,
        "complete": len(results) == expected_count,
        "global_best_absolute_determinant": str(global_best),
        "global_best": global_best_result,
        "results": results,
        "elapsed_seconds": time.monotonic() - started,
        "source": {
            "path": str(Path(__file__).absolute()),
            "sha256": sha256(Path(__file__).absolute()),
        },
    }
    atomic_write(
        output_dir / "manifest.json",
        (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    print(
        f"complete={manifest['complete']} "
        f"ht_classes={len(results)} best={global_best} "
        f"elapsed={manifest['elapsed_seconds']:.3f}s",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
