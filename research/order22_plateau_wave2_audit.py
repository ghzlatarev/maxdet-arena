#!/usr/bin/env python3
"""Audit and durably package the order-22 plateau wave-2 campaign.

The large raw harvest corpus is independently checked in place and bound by
an ordered SHA-256 corpus digest.  Compact summaries, exact classifications,
one raw representative per reported H class, the novel seed, the complete
seeded-component closure, and every exhaustive border result are retained.

This audit independently replays all 24x484 one-entry determinant tests with
exact SymPy arithmetic, reclassifies every neutral neighbor, checks every
retained H/HT certificate, recomputes all harvest and border determinants by
Bareiss elimination, and sends all 20 border maxima through ``arena verify``.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Any

import pynauty
import sympy

from h_equivalence_audit import determinant, h_certificate, transpose


ROOT = Path(__file__).resolve().parents[1]
ORDER = 22
TARGET_DETERMINANT = 409_600_000_000_000
ASSIGNMENTS_PER_BORDER = 1 << (ORDER - 1)
PINNED_PYNAUTY = "2.8.8.1"
EXPECTED_CLOSURE_RAW_SHA256 = (
    "7382f37b2df9ba68404d30d4617a0a819652446fcb5c7911bc81a074d30b61c8"
)
EXPECTED_BORDER_RAW_SHA256 = (
    "5ddf151d6257caac156a839cc8c171659abd35a051fe0e94b7ceeb4d4e5de0a8"
)
NOVEL_H = "704a5cfc16dab41953be7b0efbcc04dae093492f7122b6213d724c1fe8e40064"
NOVEL_HT = "3fd66050ccd7f0483110b0d20658b508364d7056d8bd338f969f4e9c24f4cb07"
NOVEL_RAW_SHA256 = (
    "9129d158cb0d07de90543048b6cb1aaf557e7444125db647c17bf2af0a38ecfa"
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


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return value


def read_matrix(path: Path, order: int) -> Matrix:
    try:
        matrix = [
            [int(token) for token in line.split()]
            for line in path.read_text(encoding="ascii").splitlines()
            if line.strip()
        ]
    except (OSError, UnicodeDecodeError, ValueError) as error:
        raise ValueError(f"cannot parse {path}: {error}") from error
    if len(matrix) != order or any(
        len(row) != order or any(value not in (-1, 1) for value in row)
        for row in matrix
    ):
        raise ValueError(f"{path}: expected exactly {order}x{order} signs")
    return matrix


def certificate_sha256(matrix: Matrix) -> str:
    return sha256_bytes(h_certificate(matrix))


def trace_gram_squared(matrix: Matrix) -> int:
    size = len(matrix)
    gram = [
        [
            sum(
                matrix[inner][row] * matrix[inner][column]
                for inner in range(size)
            )
            for column in range(size)
        ]
        for row in range(size)
    ]
    return sum(value * value for row in gram for value in row)


def copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def summary_integer(contents: str, key: str) -> int:
    match = re.search(
        rf"(?:^|\s){re.escape(key)}\s+(-?[0-9]+)(?:\s|$)",
        contents,
    )
    if match is None:
        raise ValueError(f"summary is missing integer field {key}")
    return int(match.group(1))


def rewrite_strings(value: Any, replacements: dict[str, str]) -> Any:
    if isinstance(value, str):
        result = value
        for source, destination in replacements.items():
            if result == source:
                return destination
            if result.startswith(source + os.sep):
                return destination + result[len(source) :]
        return result
    if isinstance(value, list):
        return [rewrite_strings(item, replacements) for item in value]
    if isinstance(value, dict):
        return {
            key: rewrite_strings(item, replacements)
            for key, item in value.items()
        }
    return value


def ordered_corpus_digest(records: list[tuple[str, str]]) -> str:
    digest = hashlib.sha256()
    for relative, raw_sha256 in records:
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(raw_sha256.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def audit_harvest(
    source: Path,
    destination: Path,
    final_destination: Path,
    known_report: Path,
) -> dict[str, Any]:
    known = read_json(known_report)
    known_h = {
        str(record["h_certificate_sha256"])
        for record in known["classes"]
    }
    known_sha256 = sha256_file(known_report)
    arm_records = []
    corpus_records: list[tuple[str, str]] = []
    total_matrices = 0
    total_exact_determinants = 0
    novel_records: list[dict[str, Any]] = []

    for arm in range(8):
        arm_name = f"arm{arm}"
        arm_source = source / arm_name
        arm_destination = destination / arm_name
        summary_path = arm_source / "summary.txt"
        stdout_path = source / f"{arm_name}.stdout.txt"
        if not summary_path.is_file() or not stdout_path.is_file():
            raise FileNotFoundError(f"incomplete wave-2 arm: {arm_name}")
        summary_contents = summary_path.read_text(encoding="utf-8")
        copy_file(summary_path, arm_destination / "summary.txt")
        copy_file(stdout_path, arm_destination / "stdout.txt")
        matrix_paths = sorted(arm_source.glob("target-*.matrix.txt"))
        saved_raw = summary_integer(summary_contents, "saved_raw")
        target_hits = summary_integer(summary_contents, "target_hits")
        attempts = summary_integer(summary_contents, "attempts")
        completed = summary_integer(summary_contents, "completed")
        best = summary_integer(summary_contents, "best")
        if saved_raw != len(matrix_paths):
            raise ArithmeticError(
                f"{arm_name}: saved_raw={saved_raw}, files={len(matrix_paths)}"
            )

        matrix_hashes: dict[Path, str] = {}
        for path in matrix_paths:
            matrix = read_matrix(path, ORDER)
            checked = determinant(matrix)
            if abs(checked) != TARGET_DETERMINANT:
                raise ArithmeticError(
                    f"{path}: exact determinant is {checked}"
                )
            raw_sha256 = sha256_file(path)
            matrix_hashes[path.resolve()] = raw_sha256
            corpus_records.append(
                (f"{arm_name}/{path.name}", raw_sha256)
            )
        total_matrices += len(matrix_paths)
        total_exact_determinants += len(matrix_paths)

        classification_path = arm_source / "classification.json"
        classification_record = None
        if matrix_paths:
            if not classification_path.is_file():
                raise FileNotFoundError(
                    f"{arm_name}: matrices exist without classification"
                )
            classification = read_json(classification_path)
            if (
                classification.get("engine")
                != "order22-gradient-harvest-classification-v1"
                or classification.get("order") != ORDER
                or classification.get("input_matrix_count")
                != len(matrix_paths)
                or classification.get("exact_determinants_checked")
                != len(matrix_paths)
                or classification.get("all_exact_determinants_on_target")
                is not True
                or classification.get("known_report", {}).get("sha256")
                != known_sha256
            ):
                raise ArithmeticError(
                    f"{arm_name}: classification binding failed"
                )
            classes = classification.get("classes")
            if not isinstance(classes, list):
                raise ValueError(f"{arm_name}: missing classification classes")
            if (
                sum(
                    int(record["raw_representative_count"])
                    for record in classes
                )
                != len(matrix_paths)
                or len(classes) != classification["h_class_count"]
                or len(
                    {
                        str(record["ht_certificate_sha256"])
                        for record in classes
                    }
                )
                != classification["ht_class_count"]
            ):
                raise ArithmeticError(
                    f"{arm_name}: classification counts failed"
                )

            normalized = json.loads(json.dumps(classification))
            normalized_classes = normalized["classes"]
            for index, (record, normalized_record) in enumerate(
                zip(classes, normalized_classes)
            ):
                representative = record["representative"]
                representative_path = Path(
                    str(representative["path"])
                ).resolve()
                if representative_path not in matrix_hashes:
                    raise ArithmeticError(
                        f"{arm_name}: representative leaves source corpus"
                    )
                if (
                    matrix_hashes[representative_path]
                    != representative["raw_sha256"]
                ):
                    raise ArithmeticError(
                        f"{arm_name}: representative SHA-256 mismatch"
                    )
                matrix = read_matrix(representative_path, ORDER)
                checked = determinant(matrix)
                direct = certificate_sha256(matrix)
                transposed = certificate_sha256(transpose(matrix))
                ht = min(direct, transposed)
                if (
                    checked != representative["determinant"]
                    or direct != record["h_certificate_sha256"]
                    or transposed
                    != record["transpose_h_certificate_sha256"]
                    or ht != record["ht_certificate_sha256"]
                    or ht != representative["ht_certificate_sha256"]
                    or trace_gram_squared(matrix)
                    != record["trace_gram_squared"]
                    or (direct in known_h)
                    is not bool(record["known_h_class"])
                ):
                    raise ArithmeticError(
                        f"{arm_name}: representative class audit failed"
                    )
                representative_destination = (
                    arm_destination
                    / "representatives"
                    / f"class-{index:03d}-{direct[:12]}.matrix.txt"
                )
                copy_file(
                    representative_path, representative_destination
                )
                normalized_record["representative"]["path"] = display_path(
                    final_destination
                    / arm_name
                    / "representatives"
                    / representative_destination.name
                )
                if not record["known_h_class"]:
                    novel_records.append(
                        {
                            "arm": arm,
                            "h_certificate_sha256": direct,
                            "transpose_h_certificate_sha256": transposed,
                            "ht_certificate_sha256": ht,
                            "path": representative_path,
                            "raw_sha256": representative[
                                "raw_sha256"
                            ],
                        }
                    )

            copy_file(
                classification_path,
                arm_destination / "classification.raw.json",
            )
            normalized["input_dir"] = str(
                final_destination / arm_name / "raw-corpus-not-retained"
            )
            normalized["source_input_dir"] = classification["input_dir"]
            write_json(
                arm_destination / "classification.json", normalized
            )
            classify_stdout = source / f"{arm_name}.classify.stdout.txt"
            if not classify_stdout.is_file():
                raise FileNotFoundError(
                    f"{arm_name}: classification stdout is missing"
                )
            copy_file(
                classify_stdout,
                arm_destination / "classification.stdout.txt",
            )
            classification_record = {
                "raw_path": display_path(
                    final_destination
                    / arm_name
                    / "classification.raw.json"
                ),
                "raw_sha256": sha256_file(classification_path),
                "h_class_count": classification["h_class_count"],
                "ht_class_count": classification["ht_class_count"],
                "known_h_class_count": classification[
                    "known_h_class_count"
                ],
                "novel_h_class_count": classification[
                    "novel_h_class_count"
                ],
                "representatives_independently_rechecked": len(classes),
            }
        elif classification_path.exists():
            raise ArithmeticError(
                f"{arm_name}: empty arm unexpectedly has a classification"
            )
        if not matrix_paths and (target_hits != 0 or saved_raw != 0):
            raise ArithmeticError(f"{arm_name}: empty-arm counters disagree")

        arm_records.append(
            {
                "arm": arm,
                "attempts": attempts,
                "completed": completed,
                "target_hits": target_hits,
                "saved_raw": saved_raw,
                "best_absolute_determinant": best,
                "summary_path": display_path(
                    final_destination / arm_name / "summary.txt"
                ),
                "summary_sha256": sha256_file(summary_path),
                "raw_matrix_count": len(matrix_paths),
                "exact_determinants_independently_rechecked": len(
                    matrix_paths
                ),
                "classification": classification_record,
            }
        )

    if (
        len(novel_records) != 1
        or novel_records[0]["h_certificate_sha256"] != NOVEL_H
        or novel_records[0]["ht_certificate_sha256"] != NOVEL_HT
        or novel_records[0]["raw_sha256"] != NOVEL_RAW_SHA256
    ):
        raise ArithmeticError("wave-2 novel-class claim failed")
    novel_destination = destination / f"novel-h-{NOVEL_H[:12]}.matrix.txt"
    copy_file(novel_records[0]["path"], novel_destination)

    return {
        "engine": "order22-gradient-harvest-v1",
        "engine_source": {
            "path": "research/order22_gradient_harvest.cpp",
            "sha256": sha256_file(
                ROOT / "research/order22_gradient_harvest.cpp"
            ),
        },
        "engine_binary": {
            "path": "build/research/order22_gradient_harvest",
            "sha256": sha256_file(
                ROOT / "build/research/order22_gradient_harvest"
            ),
        },
        "known_report": {
            "path": display_path(known_report),
            "sha256": known_sha256,
            "h_class_count": len(known_h),
        },
        "arms": arm_records,
        "raw_matrix_count": total_matrices,
        "exact_determinants_independently_rechecked": (
            total_exact_determinants
        ),
        "raw_corpus_persisted": False,
        "raw_corpus_ordered_sha256": ordered_corpus_digest(
            sorted(corpus_records)
        ),
        "retention_boundary": (
            "The 6,882 raw harvest matrices are bound by ordered name/hash "
            "digest but are not copied. All summaries, classifications, "
            "classification representatives, and the novel raw seed are "
            "retained."
        ),
        "novel_class": {
            "h_certificate_sha256": NOVEL_H,
            "transpose_h_certificate_sha256": novel_records[0][
                "transpose_h_certificate_sha256"
            ],
            "ht_certificate_sha256": NOVEL_HT,
            "raw_sha256": NOVEL_RAW_SHA256,
            "path": display_path(
                final_destination / novel_destination.name
            ),
        },
    }


def exact_neighbor_multiplicities(matrix: Matrix) -> Counter[str]:
    symbolic = sympy.Matrix(matrix)
    signed_determinant = int(symbolic.det(method="domain-ge"))
    if abs(signed_determinant) != TARGET_DETERMINANT:
        raise ArithmeticError("closure representative is off target")
    inverse = symbolic.inv(method="DM")
    if symbolic * inverse != sympy.eye(ORDER):
        raise ArithmeticError("exact inverse identity failed")
    result: Counter[str] = Counter()
    for row in range(ORDER):
        for column in range(ORDER):
            proposed = sympy.Integer(signed_determinant) * (
                1
                - 2
                * matrix[row][column]
                * inverse[column, row]
            )
            if proposed.q != 1:
                raise ArithmeticError(
                    "entry-flip determinant was not integral"
                )
            if abs(int(proposed)) != TARGET_DETERMINANT:
                continue
            neighbor = [source_row[:] for source_row in matrix]
            neighbor[row][column] = -neighbor[row][column]
            if determinant(neighbor) != int(proposed):
                raise ArithmeticError(
                    "independent flip determinant check failed"
                )
            result[certificate_sha256(neighbor)] += 1
    return result


def component_summary(
    classes: list[dict[str, Any]],
    adjacency: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    class_by_h = {
        str(record["h_certificate_sha256"]): record
        for record in classes
    }
    graph: dict[str, set[str]] = {
        certificate: set() for certificate in class_by_h
    }
    for record in adjacency:
        source = str(record["source_h_certificate_sha256"])
        for target in record["neighbor_h_class_multiplicities"]:
            graph[source].add(target)
            graph[target].add(source)
    unseen = set(graph)
    components = []
    while unseen:
        first = min(unseen)
        members = {first}
        queue = deque([first])
        while queue:
            source = queue.popleft()
            for target in graph[source]:
                if target not in members:
                    members.add(target)
                    queue.append(target)
        unseen.difference_update(members)
        components.append(
            {
                "h_class_count": len(members),
                "ht_class_count": len(
                    {
                        str(
                            class_by_h[certificate][
                                "ht_certificate_sha256"
                            ]
                        )
                        for certificate in members
                    }
                ),
                "h_certificate_sha256": sorted(members),
            }
        )
    components.sort(
        key=lambda record: (
            record["h_class_count"],
            record["h_certificate_sha256"],
        )
    )
    return components


def audit_closure(
    source: Path,
    destination: Path,
    final_destination: Path,
    novel_path: Path,
) -> tuple[dict[str, Any], dict[str, Path]]:
    raw_report_path = source / "report.json"
    raw_report_sha256 = sha256_file(raw_report_path)
    if raw_report_sha256 != EXPECTED_CLOSURE_RAW_SHA256:
        raise ArithmeticError("unexpected v3 closure report SHA-256")
    report = read_json(raw_report_path)
    if (
        report.get("engine")
        != "order22-maxdet-plateau-h-class-bfs-v1"
        or report.get("order") != ORDER
        or report.get("h_class_count") != 24
        or report.get("ht_class_count") != 20
        or report.get("classes_explored") != 24
        or report.get("entry_flip_determinants_checked") != 24 * ORDER * ORDER
        or report.get("complete_seeded_component_closure") is not True
    ):
        raise ArithmeticError("v3 closure report header failed")
    classes = report.get("classes")
    adjacency = report.get("adjacency")
    if not isinstance(classes, list) or not isinstance(adjacency, list):
        raise ValueError("v3 closure classes or adjacency are missing")
    if len(classes) != 24 or len(adjacency) != 24:
        raise ArithmeticError("v3 closure record counts failed")

    class_paths: dict[str, Path] = {}
    matrices: dict[str, Matrix] = {}
    seen_h: set[str] = set()
    seen_ht: set[str] = set()
    for record in classes:
        h = str(record["h_certificate_sha256"])
        source_path = Path(str(record["path"]))
        if not source_path.is_file():
            source_path = source / "classes" / source_path.name
        matrix = read_matrix(source_path, ORDER)
        checked = determinant(matrix)
        direct = certificate_sha256(matrix)
        transposed = certificate_sha256(transpose(matrix))
        ht = min(direct, transposed)
        if (
            sha256_file(source_path) != record["raw_sha256"]
            or checked != record["determinant"]
            or abs(checked) != TARGET_DETERMINANT
            or direct != h
            or transposed
            != record["transpose_h_certificate_sha256"]
            or ht != record["ht_certificate_sha256"]
            or record.get("explored") is not True
        ):
            raise ArithmeticError(
                f"closure class audit failed for {source_path}"
            )
        destination_path = destination / "classes" / source_path.name
        copy_file(source_path, destination_path)
        class_paths[h] = destination_path
        matrices[h] = matrix
        seen_h.add(h)
        seen_ht.add(ht)
    if len(seen_h) != 24 or len(seen_ht) != 20:
        raise ArithmeticError("closure H/HT distinct counts failed")

    adjacency_by_h = {
        str(record["source_h_certificate_sha256"]): record
        for record in adjacency
    }
    neutral_total = 0
    for h, matrix in matrices.items():
        observed = exact_neighbor_multiplicities(matrix)
        expected_record = adjacency_by_h[h]
        expected = Counter(
            {
                str(key): int(value)
                for key, value in expected_record[
                    "neighbor_h_class_multiplicities"
                ].items()
            }
        )
        if observed != expected:
            raise ArithmeticError(
                f"exact neighbor replay disagrees for H class {h}"
            )
        if sum(observed.values()) != expected_record[
            "neutral_entry_flip_count"
        ]:
            raise ArithmeticError(
                f"neutral count disagrees for H class {h}"
            )
        neutral_total += sum(observed.values())
    if neutral_total != report["neutral_labeled_neighbors"]:
        raise ArithmeticError("closure neutral total failed")

    components = component_summary(classes, adjacency)
    if sorted(
        (record["h_class_count"], record["ht_class_count"])
        for record in components
    ) != [(4, 4), (4, 4), (6, 5), (10, 7)]:
        raise ArithmeticError("unexpected closure component structure")

    copy_file(raw_report_path, destination / "report.raw.json")
    replacements = {
        str(source): str(final_destination),
        str(
            Path(
                "/tmp/order22-harvest-wave2-20260729/arm1/"
                "target-000828.matrix.txt"
            )
        ): str(novel_path),
    }
    normalized = rewrite_strings(report, replacements)
    normalized["raw_source_report"] = {
        "path": display_path(final_destination / "report.raw.json"),
        "sha256": raw_report_sha256,
    }
    write_json(destination / "report.json", normalized)
    return (
        {
            "raw_report_path": display_path(
                final_destination / "report.raw.json"
            ),
            "raw_report_sha256": raw_report_sha256,
            "normalized_report_path": display_path(
                final_destination / "report.json"
            ),
            "normalized_report_sha256": sha256_file(
                destination / "report.json"
            ),
            "h_class_count": 24,
            "ht_class_count": 20,
            "entry_flip_determinants_independently_rechecked": (
                24 * ORDER * ORDER
            ),
            "neutral_neighbors_independently_reclassified": neutral_total,
            "components": components,
        },
        class_paths,
    )


def audit_border(
    source: Path,
    destination: Path,
    final_destination: Path,
    closure_source: Path,
    final_closure: Path,
    class_paths: dict[str, Path],
    arena: Path,
) -> dict[str, Any]:
    raw_manifest_path = source / "manifest.json"
    raw_manifest_sha256 = sha256_file(raw_manifest_path)
    if raw_manifest_sha256 != EXPECTED_BORDER_RAW_SHA256:
        raise ArithmeticError("unexpected v3 border manifest SHA-256")
    manifest = read_json(raw_manifest_path)
    results = manifest.get("results")
    if (
        manifest.get("engine")
        != "order22-plateau-border-campaign-v1"
        or manifest.get("complete") is not True
        or manifest.get("ht_classes_completed") != 20
        or manifest.get("ht_classes_expected") != 20
        or manifest.get("assignments_per_core")
        != ASSIGNMENTS_PER_BORDER
        or manifest.get("assignments_completed")
        != 20 * ASSIGNMENTS_PER_BORDER
        or not isinstance(results, list)
        or len(results) != 20
    ):
        raise ArithmeticError("v3 border manifest header failed")

    normalized = json.loads(json.dumps(manifest))
    normalized_results = normalized["results"]
    seen_ht: set[str] = set()
    exact_scores = []
    receipt_records = []
    for result, normalized_result in zip(results, normalized_results):
        h = str(result["h_certificate_sha256"])
        ht = str(result["ht_certificate_sha256"])
        if ht in seen_ht:
            raise ArithmeticError(f"duplicate border HT class {ht}")
        seen_ht.add(ht)
        core_source = Path(str(result["core"]["path"]))
        if not core_source.is_file():
            core_source = closure_source / "classes" / core_source.name
        core = read_matrix(core_source, ORDER)
        if (
            sha256_file(core_source) != result["core"]["sha256"]
            or abs(determinant(core))
            != int(result["core"]["absolute_determinant"])
            or certificate_sha256(core) != h
            or min(
                h, certificate_sha256(transpose(core))
            ) != ht
            or class_paths[h].name != core_source.name
        ):
            raise ArithmeticError(f"border core audit failed for {h}")

        best_source = Path(str(result["best"]["path"]))
        log_source = Path(str(result["log"]["path"]))
        if (
            not best_source.is_file()
            or not log_source.is_file()
            or sha256_file(best_source) != result["best"]["sha256"]
            or sha256_file(log_source) != result["log"]["sha256"]
        ):
            raise ArithmeticError(f"border output hash failed for {h}")
        best = read_matrix(best_source, ORDER + 1)
        if [row[:ORDER] for row in best[:ORDER]] != core:
            raise ArithmeticError(f"border output does not contain core {h}")
        exact_score = abs(determinant(best))
        if exact_score != int(result["best"]["absolute_determinant"]):
            raise ArithmeticError(f"border exact score failed for {h}")

        events = [
            json.loads(line)
            for line in log_source.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        if (
            len(events) < 2
            or events[0].get("event") != "start"
            or events[-1].get("event") != "finished"
            or events[-1].get("complete") is not True
            or events[-1].get("assignments_completed")
            != ASSIGNMENTS_PER_BORDER
            or int(events[-1]["absolute_determinant"]) != exact_score
            or result["assignments_completed"]
            != ASSIGNMENTS_PER_BORDER
        ):
            raise ArithmeticError(f"border log completion failed for {h}")
        improvements = [
            int(event["absolute_determinant"])
            for event in events
            if event.get("event") == "new_best"
        ]
        if improvements != sorted(set(improvements)):
            raise ArithmeticError(f"border improvements are not strict for {h}")
        if not improvements or improvements[-1] != exact_score:
            raise ArithmeticError(f"border final improvement failed for {h}")

        relative_directory = Path(best_source.parent.name)
        result_destination = destination / relative_directory
        best_destination = result_destination / "best.matrix.txt"
        log_destination = result_destination / "run.jsonl"
        receipt_destination = result_destination / "best.receipt.json"
        copy_file(best_source, best_destination)
        copy_file(log_source, log_destination)
        completed = subprocess.run(
            [
                str(arena),
                "verify",
                str(best_destination),
                "--json",
                str(receipt_destination),
                "--quiet",
            ],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        if int(completed.stdout.strip()) != exact_score:
            raise ArithmeticError(f"arena verify score failed for {h}")

        final_result_directory = final_destination / relative_directory
        normalized_result["core"]["path"] = display_path(
            final_closure / "classes" / core_source.name
        )
        normalized_result["best"]["path"] = display_path(
            final_result_directory / "best.matrix.txt"
        )
        normalized_result["log"]["path"] = display_path(
            final_result_directory / "run.jsonl"
        )
        normalized_result["arena_receipt"] = {
            "path": display_path(
                final_result_directory / "best.receipt.json"
            ),
            "sha256": sha256_file(receipt_destination),
        }
        exact_scores.append(exact_score)
        receipt_records.append(
            {
                "h_certificate_sha256": h,
                "ht_certificate_sha256": ht,
                "absolute_determinant": exact_score,
                "receipt_path": display_path(
                    final_result_directory / "best.receipt.json"
                ),
                "receipt_file_sha256": sha256_file(
                    receipt_destination
                ),
            }
        )

    if len(seen_ht) != 20:
        raise ArithmeticError("border HT coverage failed")
    global_best_source = source / "global-best.matrix.txt"
    global_best_destination = destination / "global-best.matrix.txt"
    copy_file(global_best_source, global_best_destination)
    global_score = abs(
        determinant(read_matrix(global_best_destination, ORDER + 1))
    )
    if (
        global_score != max(exact_scores)
        or global_score
        != int(manifest["global_best_absolute_determinant"])
        or global_score != 2_465_792_000_000_000
    ):
        raise ArithmeticError("border global best audit failed")
    global_receipt = destination / "global-best.receipt.json"
    completed = subprocess.run(
        [
            str(arena),
            "verify",
            str(global_best_destination),
            "--json",
            str(global_receipt),
            "--quiet",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    if int(completed.stdout.strip()) != global_score:
        raise ArithmeticError("global arena verification failed")

    copy_file(raw_manifest_path, destination / "manifest.raw.json")
    normalized["order22_classes_report"] = {
        "path": display_path(final_closure / "report.raw.json"),
        "sha256": EXPECTED_CLOSURE_RAW_SHA256,
    }
    normalized["global_best"] = next(
        record
        for record in normalized_results
        if record["class_index"]
        == manifest["global_best"]["class_index"]
    )
    normalized["global_best_matrix"] = {
        "path": display_path(
            final_destination / "global-best.matrix.txt"
        ),
        "sha256": sha256_file(global_best_destination),
        "arena_receipt_path": display_path(
            final_destination / "global-best.receipt.json"
        ),
        "arena_receipt_sha256": sha256_file(global_receipt),
    }
    write_json(destination / "manifest.json", normalized)
    return {
        "raw_manifest_path": display_path(
            final_destination / "manifest.raw.json"
        ),
        "raw_manifest_sha256": raw_manifest_sha256,
        "normalized_manifest_path": display_path(
            final_destination / "manifest.json"
        ),
        "normalized_manifest_sha256": sha256_file(
            destination / "manifest.json"
        ),
        "ht_classes_exhausted": 20,
        "assignments_per_class": ASSIGNMENTS_PER_BORDER,
        "assignments_completed": 20 * ASSIGNMENTS_PER_BORDER,
        "best_matrices_exactly_rechecked": 20,
        "best_matrices_arena_verified": 20,
        "global_best_absolute_determinant": global_score,
        "frontier_absolute_determinant": 2_779_447_296_000_000,
        "frontier_improved": False,
        "receipts": receipt_records,
    }


def file_inventory(root: Path) -> list[dict[str, Any]]:
    records = []
    for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file()):
        if path.name == "manifest.json" and path.parent == root:
            continue
        records.append(
            {
                "path": str(path.relative_to(root)),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
        )
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--harvest-source",
        type=Path,
        default=Path("/tmp/order22-harvest-wave2-20260729"),
    )
    parser.add_argument(
        "--closure-source",
        type=Path,
        default=Path("/tmp/order22-maxdet-plateau-bfs-v3"),
    )
    parser.add_argument(
        "--border-source",
        type=Path,
        default=Path("/tmp/order22-plateau-border-v3b"),
    )
    parser.add_argument(
        "--known-report",
        type=Path,
        default=(
            ROOT
            / (
                "runs/direct-search/order22-plateau-harvest-20260729/"
                "closure/report.json"
            )
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=(
            ROOT
            / (
                "runs/direct-search/"
                "order22-plateau-harvest-wave2-20260729"
            )
        ),
    )
    parser.add_argument("--arena", type=Path, default=ROOT / "arena")
    arguments = parser.parse_args()

    if getattr(pynauty, "__version__", None) != PINNED_PYNAUTY:
        raise RuntimeError(f"requires pynauty=={PINNED_PYNAUTY}")
    harvest_source = arguments.harvest_source.expanduser().resolve()
    closure_source = arguments.closure_source.expanduser().resolve()
    border_source = arguments.border_source.expanduser().resolve()
    known_report = arguments.known_report.expanduser().resolve()
    output_dir = arguments.output_dir.expanduser().resolve()
    arena = arguments.arena.expanduser().resolve()
    for path in (
        harvest_source,
        closure_source,
        border_source,
    ):
        if not path.is_dir():
            raise FileNotFoundError(path)
    if not known_report.is_file() or not arena.is_file():
        raise FileNotFoundError("known report or arena executable is missing")
    if output_dir.exists():
        raise FileExistsError(output_dir)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{output_dir.name}.",
            dir=output_dir.parent,
        )
    )
    try:
        final_harvest = output_dir / "harvest"
        final_closure = output_dir / "closure"
        final_border = output_dir / "border"
        harvest_audit = audit_harvest(
            harvest_source,
            temporary / "harvest",
            final_harvest,
            known_report,
        )
        novel_path = Path(harvest_audit["novel_class"]["path"])
        closure_audit, class_paths = audit_closure(
            closure_source,
            temporary / "closure",
            final_closure,
            novel_path,
        )
        border_audit = audit_border(
            border_source,
            temporary / "border",
            final_border,
            closure_source,
            final_closure,
            class_paths,
            arena,
        )
        audit = {
            "schema_version": 1,
            "engine": "order22-plateau-wave2-independent-audit-v1",
            "claim": (
                "All retained closure representatives, all 11,616 "
                "one-entry determinant tests, all 96 neutral-neighbor "
                "classifications, all 6,882 raw harvest determinants, and "
                "all 20 exhaustive border maxima passed independent exact "
                "checks."
            ),
            "claim_boundary": (
                "The closure is complete only for components reached from "
                "the four bound seeds. The wave-2 raw corpus is hash-bound "
                "but only its class representatives are retained."
            ),
            "harvest": harvest_audit,
            "closure": closure_audit,
            "border": border_audit,
            "dependencies": {
                "pynauty": getattr(pynauty, "__version__", "unknown"),
                "sympy": sympy.__version__,
            },
            "source": {
                "path": display_path(Path(__file__)),
                "sha256": sha256_file(Path(__file__)),
            },
        }
        write_json(temporary / "independent-audit.json", audit)
        readme = """# Order-22 plateau harvest wave 2

This package retains the exact wave-2 summaries and classifications, one raw
representative per reported H class, the novel `704a5cfc…` seed, the complete
24-H / 20-HT seeded-component closure, and exhaustive fixed-border results for
all 20 HT representatives.

`independent-audit.json` rechecks all 6,882 raw harvest determinants in the
temporary source corpus, all 11,616 one-entry closure determinants, every
neutral-neighbor H certificate, every retained class certificate, every
border log and determinant, and all 20 border maxima through `arena verify`.

The raw harvest corpus is not copied; its ordered filename/SHA-256 digest and
all compact evidence are retained. Closure is only a seeded-component claim,
not a proof that no disconnected order-22 maximal-determinant components
exist. No border exceeds the order-23 frontier.
"""
        (temporary / "README.md").write_text(readme, encoding="utf-8")
        manifest = {
            "schema_version": 1,
            "engine": "order22-plateau-wave2-package-v1",
            "claim": audit["claim"],
            "claim_boundary": audit["claim_boundary"],
            "source_artifacts": {
                "harvest_directory": str(harvest_source),
                "harvest_raw_corpus_ordered_sha256": harvest_audit[
                    "raw_corpus_ordered_sha256"
                ],
                "closure_report": {
                    "path": str(closure_source / "report.json"),
                    "sha256": EXPECTED_CLOSURE_RAW_SHA256,
                },
                "border_manifest": {
                    "path": str(border_source / "manifest.json"),
                    "sha256": EXPECTED_BORDER_RAW_SHA256,
                },
            },
            "summary": {
                "wave2_raw_matrices": harvest_audit[
                    "raw_matrix_count"
                ],
                "novel_h_certificate_sha256": NOVEL_H,
                "closure_h_classes": 24,
                "closure_ht_classes": 20,
                "closure_components_h": [4, 4, 6, 10],
                "closure_components_ht": [4, 4, 5, 7],
                "border_ht_classes_completed": 20,
                "border_assignments_completed": (
                    20 * ASSIGNMENTS_PER_BORDER
                ),
                "border_global_best_absolute_determinant": (
                    2_465_792_000_000_000
                ),
            },
            "independent_audit": {
                "path": "independent-audit.json",
                "sha256": sha256_file(
                    temporary / "independent-audit.json"
                ),
            },
            "generator": {
                "path": display_path(Path(__file__)),
                "sha256": sha256_file(Path(__file__)),
            },
            "files": file_inventory(temporary),
        }
        write_json(temporary / "manifest.json", manifest)
        os.replace(temporary, output_dir)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise

    print(
        "persisted wave2=6882 closure=24H/20HT "
        "border=20x2097152 best=2465792000000000"
    )
    print(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
