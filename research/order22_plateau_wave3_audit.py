#!/usr/bin/env python3
"""Independently audit and durably package order-22 plateau harvest wave 3.

The raw 15,246-matrix harvest corpus is checked in place and bound by two
ordered SHA-256 digests, but is not copied. The package retains campaign
summaries, exact classification evidence and representatives, the two novel
seeds, the complete 30-H/26-HT seeded closure, and exhaustive fixed-border
results for all 26 HT representatives.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any

import pynauty
import sympy

from order22_plateau_wave2_audit import (
    ASSIGNMENTS_PER_BORDER,
    ORDER,
    PINNED_PYNAUTY,
    ROOT,
    TARGET_DETERMINANT,
    certificate_sha256,
    component_summary,
    copy_file,
    determinant,
    display_path,
    exact_neighbor_multiplicities,
    file_inventory,
    ordered_corpus_digest,
    read_json,
    read_matrix,
    rewrite_strings,
    sha256_file,
    summary_integer,
    transpose,
    write_json,
)


EXPECTED_KNOWN_REPORT_SHA256 = (
    "7382f37b2df9ba68404d30d4617a0a819652446fcb5c7911bc81a074d30b61c8"
)
EXPECTED_CLASSIFICATION_SHA256 = (
    "666ea211f9aa575886f7adcad8e7ce61b361f9e8c7ca7c4a1f517d83333fe798"
)
EXPECTED_CLOSURE_RAW_SHA256 = (
    "ee4d07fa3b64f66321ff4b1e1e2a71b21c452720319d39fa6789e990ddf98ee7"
)
EXPECTED_BORDER_RAW_SHA256 = (
    "baa4dd0b8b517a7a1a9e18cf561edf9f6911c92d7bc68c95a142fcc2ac731f20"
)
FRONTIER_ABSOLUTE_DETERMINANT = 2_779_447_296_000_000
EXPECTED_RAW_MATRIX_COUNT = 15_246
EXPECTED_RAW_ARM_COUNTS = {
    "arm0": 1225,
    "arm1": 1197,
    "arm2": 3000,
    "arm3": 3000,
    "arm4": 3000,
    "arm5": 3000,
    "arm6": 753,
    "arm7": 71,
}
NOVEL_CLASSES = {
    "c85238e708440b4af9d2e1109172276238b0ff0e84db8b24c2712b7a5a8b7b17": {
        "raw_sha256": (
            "cc70a44e3e1b8c9f00a2ea00425b75ed464b6ceb6e104efdfce28005899eb560"
        ),
        "source_relative": "arm3/target-001689.matrix.txt",
    },
    "6bef9d4795a1b3501057997effe8281a62f41a869776a87e27eef58a188e153e": {
        "raw_sha256": (
            "02ff8a262b13ae1806afffdbce4a5abb1374f055a7d4bcd00d5d771928ab583b"
        ),
        "source_relative": "arm3/target-002993.matrix.txt",
    },
}
NEW_COMPONENT_H = {
    "c85238e708440b4af9d2e1109172276238b0ff0e84db8b24c2712b7a5a8b7b17",
    "6bef9d4795a1b3501057997effe8281a62f41a869776a87e27eef58a188e153e",
    "11b9e8a8cffa0f3b1eb062a4307d25911d552de07411b2e6cd0488b3299a5f43",
    "d75efa364892cb76bff7bb23d07989f93c885828ee9c3542f811266de3782897",
    "af2703da0e96d59dd6ca5066322a8369d2feef1fd580db244d7f66ecbe9219ff",
    "45cd81ce407988632171653f01221d4055871ff8e48d3377c1cd6d659e71dbc7",
}
EXPECTED_NEW_BORDER = {
    "c85238e708440b4af9d2e1109172276238b0ff0e84db8b24c2712b7a5a8b7b17": (
        2_465_792_000_000_000,
        10,
    ),
    "6bef9d4795a1b3501057997effe8281a62f41a869776a87e27eef58a188e153e": (
        2_392_064_000_000_000,
        20,
    ),
    "11b9e8a8cffa0f3b1eb062a4307d25911d552de07411b2e6cd0488b3299a5f43": (
        2_392_064_000_000_000,
        20,
    ),
    "d75efa364892cb76bff7bb23d07989f93c885828ee9c3542f811266de3782897": (
        2_392_064_000_000_000,
        20,
    ),
    "af2703da0e96d59dd6ca5066322a8369d2feef1fd580db244d7f66ecbe9219ff": (
        2_465_792_000_000_000,
        10,
    ),
    "45cd81ce407988632171653f01221d4055871ff8e48d3377c1cd6d659e71dbc7": (
        2_424_832_000_000_000,
        8,
    ),
}


Matrix = list[list[int]]


def audit_harvest(
    source: Path,
    destination: Path,
    final_destination: Path,
    known_report_path: Path,
) -> dict[str, Any]:
    known_report_sha256 = sha256_file(known_report_path)
    if known_report_sha256 != EXPECTED_KNOWN_REPORT_SHA256:
        raise ArithmeticError("unexpected v3 known-report SHA-256")
    known_report = read_json(known_report_path)
    known_h = {
        str(record["h_certificate_sha256"])
        for record in known_report["classes"]
    }
    if (
        len(known_h) != 24
        or known_report.get("h_class_count") != 24
        or known_report.get("ht_class_count") != 20
    ):
        raise ArithmeticError("unexpected v3 known-report class counts")

    arm_records = []
    corpus_records: list[tuple[str, str]] = []
    raw_paths: dict[Path, tuple[str, str, Matrix, int]] = {}
    named_payload = hashlib.sha256()
    sign_counts: Counter[str] = Counter()
    for arm in range(8):
        arm_name = f"arm{arm}"
        arm_source = source / arm_name
        summary_path = arm_source / "summary.txt"
        stdout_path = source / f"{arm_name}.stdout.txt"
        if not summary_path.is_file() or not stdout_path.is_file():
            raise FileNotFoundError(f"incomplete wave-3 arm: {arm_name}")
        summary_contents = summary_path.read_text(encoding="utf-8")
        matrix_paths = sorted(arm_source.glob("target-*.matrix.txt"))
        if (
            len(matrix_paths) != EXPECTED_RAW_ARM_COUNTS[arm_name]
            or summary_integer(summary_contents, "saved_raw")
            != len(matrix_paths)
        ):
            raise ArithmeticError(f"{arm_name}: raw-count mismatch")

        for path in matrix_paths:
            contents = path.read_bytes()
            matrix = read_matrix(path, ORDER)
            checked = determinant(matrix)
            if abs(checked) != TARGET_DETERMINANT:
                raise ArithmeticError(f"{path}: exact determinant {checked}")
            raw_sha256 = hashlib.sha256(contents).hexdigest()
            relative = f"{arm_name}/{path.name}"
            corpus_records.append((relative, raw_sha256))
            raw_paths[path.resolve()] = (
                arm_name,
                raw_sha256,
                matrix,
                checked,
            )
            named_payload.update(arm_name.encode("utf-8"))
            named_payload.update(b"/")
            named_payload.update(path.name.encode("ascii"))
            named_payload.update(b"\0")
            named_payload.update(contents)
            sign_counts["positive" if checked > 0 else "negative"] += 1

        arm_destination = destination / arm_name
        copy_file(summary_path, arm_destination / "summary.txt")
        copy_file(stdout_path, arm_destination / "stdout.txt")
        arm_records.append(
            {
                "arm": arm,
                "random_seed": summary_integer(
                    summary_contents, "random_seed"
                ),
                "minimum_kick": summary_integer(
                    summary_contents, "minimum_kick"
                ),
                "maximum_kick": summary_integer(
                    summary_contents, "maximum_kick"
                ),
                "attempts": summary_integer(summary_contents, "attempts"),
                "completed": summary_integer(
                    summary_contents, "completed"
                ),
                "target_hits": summary_integer(
                    summary_contents, "target_hits"
                ),
                "saved_raw": len(matrix_paths),
                "best_absolute_determinant": summary_integer(
                    summary_contents, "best"
                ),
                "summary_path": display_path(
                    final_destination / arm_name / "summary.txt"
                ),
                "summary_sha256": sha256_file(summary_path),
                "stdout_path": display_path(
                    final_destination / arm_name / "stdout.txt"
                ),
                "stdout_sha256": sha256_file(stdout_path),
                "exact_determinants_independently_rechecked": len(
                    matrix_paths
                ),
            }
        )

    if len(raw_paths) != EXPECTED_RAW_MATRIX_COUNT:
        raise ArithmeticError("wave-3 total raw-count mismatch")

    classification_path = source / "exact-classification.json"
    if sha256_file(classification_path) != EXPECTED_CLASSIFICATION_SHA256:
        raise ArithmeticError("unexpected wave-3 classification SHA-256")
    classification = read_json(classification_path)
    classes = classification.get("classes")
    if (
        classification.get("engine")
        != "order22-known-class-pivot-audit-v1"
        or classification.get("order") != ORDER
        or classification.get("input_matrix_count")
        != EXPECTED_RAW_MATRIX_COUNT
        or classification.get("exact_determinants_checked")
        != EXPECTED_RAW_MATRIX_COUNT
        or classification.get("all_exact_determinants_on_target") is not True
        or classification.get("known_report", {}).get("sha256")
        != EXPECTED_KNOWN_REPORT_SHA256
        or classification.get("input_arm_counts")
        != EXPECTED_RAW_ARM_COUNTS
        or classification.get("determinant_sign_counts")
        != dict(sorted(sign_counts.items()))
        or classification.get("named_payload_sha256")
        != named_payload.hexdigest()
        or classification.get("observed_h_class_count") != 26
        or classification.get("observed_ht_class_count") != 22
        or classification.get("known_h_classes_observed") != 24
        or classification.get("novel_h_class_count") != 2
        or not isinstance(classes, list)
        or len(classes) != 26
        or sum(int(record["matrix_count"]) for record in classes)
        != EXPECTED_RAW_MATRIX_COUNT
    ):
        raise ArithmeticError("wave-3 classification header failed")

    normalized = json.loads(json.dumps(classification))
    normalized_classes = normalized["classes"]
    observed_novel: dict[str, dict[str, Any]] = {}
    representative_records = []
    for index, (record, normalized_record) in enumerate(
        zip(classes, normalized_classes)
    ):
        representative = record["representative"]
        source_path = Path(str(representative["path"])).resolve()
        if source_path not in raw_paths:
            raise ArithmeticError("classification representative left corpus")
        arm_name, raw_sha256, matrix, checked = raw_paths[source_path]
        direct = certificate_sha256(matrix)
        transposed = certificate_sha256(transpose(matrix))
        ht = min(direct, transposed)
        if (
            raw_sha256 != representative["raw_sha256"]
            or checked != int(representative["determinant"])
            or direct != record["h_certificate_sha256"]
            or transposed != record["transpose_h_certificate_sha256"]
            or ht != record["ht_certificate_sha256"]
            or (direct in known_h) is not bool(record["known_h_class"])
            or sum(int(value) for value in record["arm_counts"].values())
            != int(record["matrix_count"])
            or int(record["arm_counts"].get(arm_name, 0)) < 1
        ):
            raise ArithmeticError(
                f"classification representative audit failed: {direct}"
            )
        representative_destination = (
            destination
            / "representatives"
            / f"class-{index:03d}-{direct[:12]}.matrix.txt"
        )
        copy_file(source_path, representative_destination)
        final_representative = (
            final_destination
            / "representatives"
            / representative_destination.name
        )
        normalized_record["representative"]["path"] = display_path(
            final_representative
        )
        representative_records.append(
            {
                "h_certificate_sha256": direct,
                "ht_certificate_sha256": ht,
                "path": display_path(final_representative),
                "raw_sha256": raw_sha256,
            }
        )
        if record["known_h_class"] is False:
            observed_novel[direct] = {
                "source_path": source_path,
                "raw_sha256": raw_sha256,
                "determinant": checked,
                "transpose_h_certificate_sha256": transposed,
                "ht_certificate_sha256": ht,
            }

    if set(observed_novel) != set(NOVEL_CLASSES):
        raise ArithmeticError("unexpected wave-3 novel H classes")

    novel_destinations: dict[str, Path] = {}
    for direct, expected in NOVEL_CLASSES.items():
        record = observed_novel[direct]
        if (
            record["raw_sha256"] != expected["raw_sha256"]
            or str(
                record["source_path"].relative_to(source)
            )
            != expected["source_relative"]
            or record["transpose_h_certificate_sha256"] != direct
            or record["ht_certificate_sha256"] != direct
        ):
            raise ArithmeticError(f"novel-class binding failed: {direct}")
        novel_destination = (
            destination / f"novel-h-{direct[:12]}.matrix.txt"
        )
        copy_file(record["source_path"], novel_destination)
        novel_destinations[direct] = novel_destination

    normalized["input_dirs"] = [
        display_path(final_destination / arm / "raw-corpus-not-retained")
        for arm in sorted(EXPECTED_RAW_ARM_COUNTS)
    ]
    normalized["known_report"]["path"] = display_path(known_report_path)
    normalized["source"]["path"] = "research/order22_known_class_audit.py"
    for record in normalized["novel_matrices"]:
        direct = str(record["h_certificate_sha256"])
        record["path"] = display_path(
            final_destination / novel_destinations[direct].name
        )
    normalized["raw_source_classification"] = {
        "path": display_path(
            final_destination / "exact-classification.raw.json"
        ),
        "sha256": EXPECTED_CLASSIFICATION_SHA256,
    }
    copy_file(
        classification_path,
        destination / "exact-classification.raw.json",
    )
    write_json(destination / "exact-classification.json", normalized)

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
        "classification_source": {
            "path": "research/order22_known_class_audit.py",
            "sha256": sha256_file(
                ROOT / "research/order22_known_class_audit.py"
            ),
        },
        "classification_raw_sha256": EXPECTED_CLASSIFICATION_SHA256,
        "classification_normalized_path": display_path(
            final_destination / "exact-classification.json"
        ),
        "classification_normalized_sha256": sha256_file(
            destination / "exact-classification.json"
        ),
        "known_report": {
            "path": display_path(known_report_path),
            "sha256": known_report_sha256,
            "h_class_count": 24,
            "ht_class_count": 20,
        },
        "arms": arm_records,
        "raw_matrix_count": EXPECTED_RAW_MATRIX_COUNT,
        "exact_determinants_independently_rechecked": (
            EXPECTED_RAW_MATRIX_COUNT
        ),
        "raw_corpus_persisted": False,
        "raw_corpus_ordered_name_hash_sha256": ordered_corpus_digest(
            sorted(corpus_records)
        ),
        "raw_corpus_ordered_named_payload_sha256": named_payload.hexdigest(),
        "classification_representatives_retained": representative_records,
        "novel_classes": [
            {
                "h_certificate_sha256": direct,
                "transpose_h_certificate_sha256": direct,
                "ht_certificate_sha256": direct,
                "determinant": observed_novel[direct]["determinant"],
                "raw_sha256": observed_novel[direct]["raw_sha256"],
                "source_relative": NOVEL_CLASSES[direct][
                    "source_relative"
                ],
                "path": display_path(
                    final_destination / novel_destinations[direct].name
                ),
            }
            for direct in sorted(NOVEL_CLASSES)
        ],
        "retention_boundary": (
            "The 15,246 raw matrices were independently determinant-checked "
            "and are bound by ordered digests, but only the 26 classification "
            "representatives and two novel raw seeds are copied."
        ),
    }


def audit_closure(
    source: Path,
    destination: Path,
    final_destination: Path,
    harvest_source: Path,
    final_harvest: Path,
    previous_closure_source: Path,
) -> tuple[dict[str, Any], dict[str, Path]]:
    raw_report_path = source / "report.json"
    if sha256_file(raw_report_path) != EXPECTED_CLOSURE_RAW_SHA256:
        raise ArithmeticError("unexpected v4 closure report SHA-256")
    report = read_json(raw_report_path)
    classes = report.get("classes")
    adjacency = report.get("adjacency")
    if (
        report.get("engine")
        != "order22-maxdet-plateau-h-class-bfs-v1"
        or report.get("order") != ORDER
        or report.get("h_class_count") != 30
        or report.get("ht_class_count") != 26
        or report.get("classes_explored") != 30
        or report.get("entry_flip_determinants_checked")
        != 30 * ORDER * ORDER
        or report.get("neutral_labeled_neighbors") != 120
        or report.get("complete_seeded_component_closure") is not True
        or not isinstance(classes, list)
        or not isinstance(adjacency, list)
        or len(classes) != 30
        or len(adjacency) != 30
    ):
        raise ArithmeticError("v4 closure header failed")

    class_paths: dict[str, Path] = {}
    matrices: dict[str, Matrix] = {}
    seen_ht: set[str] = set()
    for record in classes:
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
            or checked != int(record["determinant"])
            or abs(checked) != TARGET_DETERMINANT
            or direct != record["h_certificate_sha256"]
            or transposed != record["transpose_h_certificate_sha256"]
            or ht != record["ht_certificate_sha256"]
            or record.get("explored") is not True
        ):
            raise ArithmeticError(f"closure class failed: {source_path}")
        destination_path = destination / "classes" / source_path.name
        copy_file(source_path, destination_path)
        class_paths[direct] = destination_path
        matrices[direct] = matrix
        seen_ht.add(ht)
    if len(class_paths) != 30 or len(seen_ht) != 26:
        raise ArithmeticError("v4 closure distinct-count audit failed")

    adjacency_by_h = {
        str(record["source_h_certificate_sha256"]): record
        for record in adjacency
    }
    neutral_total = 0
    for direct, matrix in matrices.items():
        observed = exact_neighbor_multiplicities(matrix)
        expected_record = adjacency_by_h[direct]
        expected = Counter(
            {
                str(key): int(value)
                for key, value in expected_record[
                    "neighbor_h_class_multiplicities"
                ].items()
            }
        )
        if (
            observed != expected
            or sum(observed.values())
            != int(expected_record["neutral_entry_flip_count"])
        ):
            raise ArithmeticError(f"closure adjacency failed: {direct}")
        neutral_total += sum(observed.values())
    if neutral_total != 120:
        raise ArithmeticError("closure neutral-neighbor total failed")

    components = component_summary(classes, adjacency)
    if sorted(
        (record["h_class_count"], record["ht_class_count"])
        for record in components
    ) != [(4, 4), (4, 4), (6, 5), (6, 6), (10, 7)]:
        raise ArithmeticError("unexpected v4 component topology")
    new_components = [
        record
        for record in components
        if set(record["h_certificate_sha256"]) & NEW_COMPONENT_H
    ]
    if (
        len(new_components) != 1
        or set(new_components[0]["h_certificate_sha256"])
        != NEW_COMPONENT_H
        or new_components[0]["h_class_count"] != 6
        or new_components[0]["ht_class_count"] != 6
    ):
        raise ArithmeticError("new classes are not one 6-H/6-HT component")

    copy_file(raw_report_path, destination / "report.raw.json")
    replacements: dict[str, str] = {}
    for direct, expected in NOVEL_CLASSES.items():
        replacements[str(harvest_source / expected["source_relative"])] = (
            display_path(final_harvest / f"novel-h-{direct[:12]}.matrix.txt")
        )
    replacements[str(previous_closure_source)] = display_path(
        final_destination
    )
    replacements[str(source)] = display_path(final_destination)
    normalized = rewrite_strings(report, replacements)
    normalized["raw_source_report"] = {
        "path": display_path(final_destination / "report.raw.json"),
        "sha256": EXPECTED_CLOSURE_RAW_SHA256,
    }
    normalized["component_summary"] = components
    write_json(destination / "report.json", normalized)
    return (
        {
            "raw_report_path": display_path(
                final_destination / "report.raw.json"
            ),
            "raw_report_sha256": EXPECTED_CLOSURE_RAW_SHA256,
            "normalized_report_path": display_path(
                final_destination / "report.json"
            ),
            "normalized_report_sha256": sha256_file(
                destination / "report.json"
            ),
            "h_class_count": 30,
            "ht_class_count": 26,
            "entry_flip_determinants_independently_rechecked": (
                30 * ORDER * ORDER
            ),
            "neutral_neighbors_independently_reclassified": neutral_total,
            "components": components,
            "new_component": new_components[0],
            "orrick_count_comparison": {
                "our_seeded_closure_h_class_count": 30,
                "orrick_2006_reported_h_class_count": 30,
                "counts_equal": True,
                "claim_boundary": (
                    "Count equality does not identify our classes with "
                    "Orrick's unpublished 30-class corpus and does not prove "
                    "global exhaustion beyond the seeded components."
                ),
            },
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
    if sha256_file(raw_manifest_path) != EXPECTED_BORDER_RAW_SHA256:
        raise ArithmeticError("unexpected v4 border manifest SHA-256")
    manifest = read_json(raw_manifest_path)
    results = manifest.get("results")
    if (
        manifest.get("engine")
        != "order22-plateau-border-campaign-v1"
        or manifest.get("complete") is not True
        or manifest.get("ht_classes_completed") != 26
        or manifest.get("ht_classes_expected") != 26
        or manifest.get("assignments_per_core")
        != ASSIGNMENTS_PER_BORDER
        or manifest.get("assignments_completed")
        != 26 * ASSIGNMENTS_PER_BORDER
        or not isinstance(results, list)
        or len(results) != 26
    ):
        raise ArithmeticError("v4 border header failed")

    normalized = json.loads(json.dumps(manifest))
    normalized_results = normalized["results"]
    seen_ht: set[str] = set()
    exact_scores = []
    receipt_records = []
    new_results = []
    for result, normalized_result in zip(results, normalized_results):
        direct = str(result["h_certificate_sha256"])
        ht = str(result["ht_certificate_sha256"])
        if ht in seen_ht:
            raise ArithmeticError(f"duplicate border HT class: {ht}")
        seen_ht.add(ht)
        core_source = Path(str(result["core"]["path"]))
        if not core_source.is_file():
            core_source = closure_source / "classes" / core_source.name
        core = read_matrix(core_source, ORDER)
        if (
            sha256_file(core_source) != result["core"]["sha256"]
            or abs(determinant(core))
            != int(result["core"]["absolute_determinant"])
            or certificate_sha256(core) != direct
            or min(direct, certificate_sha256(transpose(core))) != ht
            or class_paths[direct].name != core_source.name
        ):
            raise ArithmeticError(f"border core failed: {direct}")

        best_source = Path(str(result["best"]["path"]))
        log_source = Path(str(result["log"]["path"]))
        if (
            not best_source.is_file()
            or not log_source.is_file()
            or sha256_file(best_source) != result["best"]["sha256"]
            or sha256_file(log_source) != result["log"]["sha256"]
        ):
            raise ArithmeticError(f"border file hash failed: {direct}")
        best = read_matrix(best_source, ORDER + 1)
        if [row[:ORDER] for row in best[:ORDER]] != core:
            raise ArithmeticError(f"border core embedding failed: {direct}")
        exact_score = abs(determinant(best))
        if exact_score != int(result["best"]["absolute_determinant"]):
            raise ArithmeticError(f"border determinant failed: {direct}")

        events = [
            json.loads(line)
            for line in log_source.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        improvements = [
            int(event["absolute_determinant"])
            for event in events
            if event.get("event") == "new_best"
        ]
        if (
            len(events) < 2
            or events[0].get("event") != "start"
            or events[-1].get("event") != "finished"
            or events[-1].get("complete") is not True
            or events[-1].get("assignments_completed")
            != ASSIGNMENTS_PER_BORDER
            or int(events[-1]["absolute_determinant"]) != exact_score
            or int(result["assignments_completed"])
            != ASSIGNMENTS_PER_BORDER
            or improvements != sorted(set(improvements))
            or not improvements
            or improvements[-1] != exact_score
        ):
            raise ArithmeticError(f"border log failed: {direct}")

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
            raise ArithmeticError(f"arena verification failed: {direct}")

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
        record = {
            "class_index": int(result["class_index"]),
            "h_certificate_sha256": direct,
            "ht_certificate_sha256": ht,
            "absolute_determinant": exact_score,
            "best_border_columns_up_to_global_sign": int(
                result["best_border_columns_up_to_global_sign"]
            ),
            "matrix_path": display_path(
                final_result_directory / "best.matrix.txt"
            ),
            "matrix_sha256": sha256_file(best_destination),
            "receipt_path": display_path(
                final_result_directory / "best.receipt.json"
            ),
            "receipt_sha256": sha256_file(receipt_destination),
        }
        receipt_records.append(record)
        if direct in NEW_COMPONENT_H:
            expected_score, expected_multiplicity = EXPECTED_NEW_BORDER[
                direct
            ]
            if (
                exact_score != expected_score
                or record["best_border_columns_up_to_global_sign"]
                != expected_multiplicity
            ):
                raise ArithmeticError(
                    f"new-component border changed: {direct}"
                )
            new_results.append(record)
        exact_scores.append(exact_score)

    if len(seen_ht) != 26 or len(new_results) != 6:
        raise ArithmeticError("v4 border coverage failed")
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
        raise ArithmeticError("v4 border global-best audit failed")
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
        if record["class_index"] == manifest["global_best"]["class_index"]
    )
    normalized["global_best_matrix"] = {
        "path": display_path(final_destination / "global-best.matrix.txt"),
        "sha256": sha256_file(global_best_destination),
        "arena_receipt_path": display_path(
            final_destination / "global-best.receipt.json"
        ),
        "arena_receipt_sha256": sha256_file(global_receipt),
    }
    normalized["new_component_results"] = sorted(
        new_results, key=lambda record: record["class_index"]
    )
    write_json(destination / "manifest.json", normalized)
    return {
        "raw_manifest_path": display_path(
            final_destination / "manifest.raw.json"
        ),
        "raw_manifest_sha256": EXPECTED_BORDER_RAW_SHA256,
        "normalized_manifest_path": display_path(
            final_destination / "manifest.json"
        ),
        "normalized_manifest_sha256": sha256_file(
            destination / "manifest.json"
        ),
        "ht_classes_exhausted": 26,
        "assignments_per_class": ASSIGNMENTS_PER_BORDER,
        "assignments_completed": 26 * ASSIGNMENTS_PER_BORDER,
        "best_matrices_exactly_rechecked": 26,
        "best_matrices_arena_verified": 26,
        "global_best_absolute_determinant": global_score,
        "frontier_absolute_determinant": FRONTIER_ABSOLUTE_DETERMINANT,
        "frontier_improved": global_score > FRONTIER_ABSOLUTE_DETERMINANT,
        "new_component_best_absolute_determinant": max(
            record["absolute_determinant"] for record in new_results
        ),
        "new_component_results": sorted(
            new_results, key=lambda record: record["class_index"]
        ),
        "receipts": receipt_records,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--harvest-source",
        type=Path,
        default=Path("/tmp/order22-harvest-wave3-20260729"),
    )
    parser.add_argument(
        "--closure-source",
        type=Path,
        default=Path("/tmp/order22-maxdet-plateau-bfs-v4"),
    )
    parser.add_argument(
        "--border-source",
        type=Path,
        default=Path("/tmp/order22-maxdet-plateau-border-v4"),
    )
    parser.add_argument(
        "--previous-closure-source",
        type=Path,
        default=Path("/tmp/order22-maxdet-plateau-bfs-v3"),
    )
    parser.add_argument(
        "--known-report",
        type=Path,
        default=(
            ROOT
            / (
                "runs/direct-search/order22-plateau-harvest-wave2-20260729/"
                "closure/report.raw.json"
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
                "order22-plateau-harvest-wave3-20260729"
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
    previous_closure_source = (
        arguments.previous_closure_source.expanduser().resolve()
    )
    known_report = arguments.known_report.expanduser().resolve()
    output_dir = arguments.output_dir.expanduser().resolve()
    arena = arguments.arena.expanduser().resolve()
    for path in (
        harvest_source,
        closure_source,
        border_source,
        previous_closure_source,
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
        closure_audit, class_paths = audit_closure(
            closure_source,
            temporary / "closure",
            final_closure,
            harvest_source,
            final_harvest,
            previous_closure_source,
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
            "engine": "order22-plateau-wave3-independent-audit-v1",
            "claim": (
                "All 15,246 raw harvest determinants, all 14,520 "
                "one-entry closure determinants, all 120 neutral-neighbor "
                "classifications, all retained H/HT certificates, and all "
                "26 exhaustive border maxima passed independent exact checks."
            ),
            "claim_boundary": (
                "The closure is complete only for components reached from "
                "the 26 bound seed orientations. Equality between our 30 "
                "H-class count and Orrick's historical count neither proves "
                "the class sets are identical nor proves global exhaustion."
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
            "support_source": {
                "path": "research/order22_plateau_wave2_audit.py",
                "sha256": sha256_file(
                    ROOT / "research/order22_plateau_wave2_audit.py"
                ),
            },
        }
        write_json(temporary / "independent-audit.json", audit)
        readme = """# Order-22 plateau harvest wave 3

Wave 3 exact-audited 15,246 raw order-22 maximal-determinant matrices and
found two H-classes absent from the frozen 24-H closure. Their one-entry
plateau closure adds six H/HT classes, producing a 30-H / 26-HT seeded
closure with five connected components of H sizes 4, 4, 6, 6, and 10.

The six new classes form one connected 6-H / 6-HT component. Exhaustive
fixed-border enumeration checked all 2^21 sign columns for every one of the
26 HT representatives. The new component's best order-23 start is
2,465,792,000,000,000, tied by classes `c85238e70844…` and `af2703da0e96…`;
it does not beat the 2,779,447,296,000,000 frontier.

`independent-audit.json` rechecks every raw determinant, every closure entry
flip and class certificate, every border log and determinant, and all 26
border maxima with `arena verify`. The raw harvest corpus is not copied; it
is bound by ordered filename/SHA-256 and named-payload digests.

Our 30-H count equals the count reported by Orrick (2006). This does not show
that the two class sets are identical because Orrick's 30 representatives
were not published, and seeded closure is not a proof that no disconnected
components exist.
"""
        (temporary / "README.md").write_text(readme, encoding="utf-8")
        manifest = {
            "schema_version": 1,
            "engine": "order22-plateau-wave3-package-v1",
            "claim": audit["claim"],
            "claim_boundary": audit["claim_boundary"],
            "source_artifacts": {
                "harvest_directory": str(harvest_source),
                "harvest_raw_corpus_ordered_name_hash_sha256": (
                    harvest_audit[
                        "raw_corpus_ordered_name_hash_sha256"
                    ]
                ),
                "harvest_raw_corpus_ordered_named_payload_sha256": (
                    harvest_audit[
                        "raw_corpus_ordered_named_payload_sha256"
                    ]
                ),
                "classification": {
                    "path": str(
                        harvest_source / "exact-classification.json"
                    ),
                    "sha256": EXPECTED_CLASSIFICATION_SHA256,
                },
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
                "wave3_raw_matrices": EXPECTED_RAW_MATRIX_COUNT,
                "novel_h_certificates_sha256": sorted(NOVEL_CLASSES),
                "closure_h_classes": 30,
                "closure_ht_classes": 26,
                "closure_components_h": [4, 4, 6, 6, 10],
                "closure_components_ht": [4, 4, 5, 6, 7],
                "new_component_h_classes": 6,
                "new_component_ht_classes": 6,
                "border_ht_classes_completed": 26,
                "border_assignments_completed": (
                    26 * ASSIGNMENTS_PER_BORDER
                ),
                "border_global_best_absolute_determinant": (
                    2_465_792_000_000_000
                ),
                "new_component_border_best_absolute_determinant": (
                    2_465_792_000_000_000
                ),
                "frontier_improved": False,
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
        "persisted wave3=15246 closure=30H/26HT "
        "border=26x2097152 best=2465792000000000"
    )
    print(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
