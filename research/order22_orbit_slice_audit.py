#!/usr/bin/env python3
"""Freeze and independently audit the order-22 Gram-orbit slice pilot."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Any

import pynauty

from order22_factor_hunt_campaign import (
    ORDER,
    TARGET_DETERMINANT,
    certificate_sha256,
    gram,
    normalized_column_masks,
    read_matrix,
)
from h_equivalence_audit import determinant, transpose


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ALIGNMENT = (
    ROOT
    / "runs/direct-search/order22-gram-aligned-factors-20260729/manifest.json"
)
DEFAULT_CAMPAIGN = (
    ROOT
    / "runs/direct-search/order22-gram-factor-orbit-slices-20260729"
)
DEFAULT_ORBITS = (
    ROOT / "runs/direct-search/order22-gram-shell-orbits-20260729"
)
PINNED_PYNAUTY = "2.8.8.1"

RUN_SPECS = (
    ("gsds-known-vector", "gsds", "UNKNOWN", False),
    ("mendeley-vector-b", "mendeley", "UNKNOWN", False),
    ("gsds-anchored", "gsds", "OPTIMAL", True),
    ("mendeley-b-anchored", "mendeley", "FEASIBLE", True),
    ("random-gsds-595001", "gsds", "FEASIBLE", True),
    ("random-mendeley-a-595002", "mendeley", "FEASIBLE", True),
    ("random-mendeley-b-595003", "mendeley", "FEASIBLE", True),
)
HINT_SPECS = (
    ("gsds", "gsds"),
    ("mendeley-a", "mendeley"),
    ("mendeley-b", "mendeley"),
)


class AuditError(ValueError):
    """An inconsistent frozen result."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def atomic_write(path: Path, contents: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise AuditError(f"refusing to overwrite {path}")
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


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise AuditError(f"{path}: expected a JSON object")
    return value


def exact_classification(
    path: Path,
    target_gram: list[list[int]],
    known: list[dict[str, Any]],
) -> dict[str, Any]:
    matrix = read_matrix(path)
    checked = determinant(matrix)
    if abs(checked) != TARGET_DETERMINANT or gram(matrix) != target_gram:
        raise AuditError(f"{path}: exact determinant/Gram check failed")
    direct = certificate_sha256(matrix)
    transposed = certificate_sha256(transpose(matrix))
    ht = min(direct, transposed)
    support = sorted(set(normalized_column_masks(matrix)))
    known_h = {str(record["h_certificate_sha256"]) for record in known}
    known_ht = {str(record["ht_certificate_sha256"]) for record in known}
    exact_matches = sorted(
        str(record["h_certificate_sha256"])
        for record in known
        if support
        == sorted(
            set(int(mask) for mask in record["normalized_column_masks_decimal"])
        )
    )
    return {
        "path": display_path(path),
        "sha256": sha256(path),
        "determinant": checked,
        "h_certificate_sha256": direct,
        "transpose_h_certificate_sha256": transposed,
        "ht_certificate_sha256": ht,
        "known_h_class": direct in known_h,
        "known_ht_class": ht in known_ht,
        "exact_known_support_matches_h": exact_matches,
        "normalized_column_masks_decimal": support,
    }


def projective_factor_graph(matrix: list[list[int]]) -> pynauty.Graph:
    """Encode signed row/column monomial automorphisms with colors fixed."""

    column_clone_offset = 3 * ORDER
    column_fiber_offset = 5 * ORDER
    vertex_count = 6 * ORDER
    adjacency = {vertex: set() for vertex in range(vertex_count)}

    def add_edge(left: int, right: int) -> None:
        adjacency[left].add(right)
        adjacency[right].add(left)

    for row in range(ORDER):
        add_edge(2 * ORDER + row, 2 * row)
        add_edge(2 * ORDER + row, 2 * row + 1)
    for column in range(ORDER):
        add_edge(column_fiber_offset + column, column_clone_offset + 2 * column)
        add_edge(
            column_fiber_offset + column,
            column_clone_offset + 2 * column + 1,
        )
    for row in range(ORDER):
        for column in range(ORDER):
            cross = 1 if matrix[row][column] < 0 else 0
            add_edge(
                2 * row,
                column_clone_offset + 2 * column + cross,
            )
            add_edge(
                2 * row + 1,
                column_clone_offset + 2 * column + 1 - cross,
            )
    return pynauty.Graph(
        number_of_vertices=vertex_count,
        directed=False,
        adjacency_dict={
            vertex: sorted(neighbors)
            for vertex, neighbors in adjacency.items()
        },
        vertex_coloring=[
            set(range(0, 2 * ORDER)),
            set(range(2 * ORDER, 3 * ORDER)),
            set(range(3 * ORDER, 5 * ORDER)),
            set(range(5 * ORDER, 6 * ORDER)),
        ],
    )


def group_order(graph: pynauty.Graph) -> int:
    automorphisms = pynauty.autgrp(graph)
    return int(round(automorphisms[1] * (10 ** automorphisms[2])))


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", type=Path, default=DEFAULT_ALIGNMENT)
    parser.add_argument("--campaign-dir", type=Path, default=DEFAULT_CAMPAIGN)
    parser.add_argument("--orbit-dir", type=Path, default=DEFAULT_ORBITS)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    if getattr(pynauty, "__version__", None) != PINNED_PYNAUTY:
        raise RuntimeError(f"requires pynauty=={PINNED_PYNAUTY}")

    alignment_path = arguments.alignment.expanduser().resolve()
    campaign_dir = arguments.campaign_dir.expanduser().resolve()
    orbit_dir = arguments.orbit_dir.expanduser().resolve()
    output_path = (
        arguments.output.expanduser().resolve()
        if arguments.output is not None
        else campaign_dir / "manifest.json"
    )
    alignment = read_json(alignment_path)
    if (
        alignment.get("h_class_count") != 30
        or alignment.get("ht_class_count") != 26
        or alignment.get("target_gram_class_count") != 2
    ):
        raise AuditError("alignment is not the audited 30-H / 26-HT package")
    known_by_target: dict[str, list[dict[str, Any]]] = {}
    for record in alignment["h_factors"]:
        known_by_target.setdefault(str(record["target"]), []).append(record)
    targets = {
        str(record["name"]): record for record in alignment["targets"]
    }

    orbit_records = {}
    target_grams = {}
    target_group_orders = {}
    for target_name in ("gsds", "mendeley"):
        path = orbit_dir / f"{target_name}.json"
        report = read_json(path)
        if (
            report.get("engine") != "order22-signed-gram-shell-orbits-v1"
            or report.get("name") != target_name
            or len(report.get("known_factors", []))
            != len(known_by_target[target_name])
        ):
            raise AuditError(f"{path}: unexpected orbit audit")
        factor_path = Path(report["inputs"]["factor"]["path"])
        if (
            report["inputs"]["factor"]["sha256"] != sha256(factor_path)
            or report["inputs"]["shell_report"]["sha256"]
            != sha256(Path(report["inputs"]["shell_report"]["path"]))
        ):
            raise AuditError(f"{path}: broken input hash binding")
        target_grams[target_name] = gram(read_matrix(factor_path))
        target_group_orders[target_name] = int(
            report["group"]["signed_double_cover_order"]
        )
        orbit_records[target_name] = {
            "path": display_path(path),
            "sha256": sha256(path),
            "signed_double_cover_group_order": target_group_orders[target_name],
            "generator_count": report["group"]["generator_count"],
            "shell_size": report["inputs"]["shell_report"]["shell_size"],
            "shell_orbit_sizes": report["shell_partition"]["orbit_sizes"],
            "count_system_rank": report["orbit_count_system"]["rank"],
            "forced_orbit_counts": report["orbit_count_system"][
                "forced_orbit_counts"
            ],
            "known_orbit_count_vectors": report[
                "known_orbit_count_vectors"
            ],
        }

    hint_records = []
    for hint_name, target_name in HINT_SPECS:
        metadata_path = campaign_dir / "hints" / f"{hint_name}.json"
        metadata = read_json(metadata_path)
        output_path_hint = campaign_dir / "hints" / f"{hint_name}.matrix.txt"
        if (
            metadata.get("engine") != "order22-gram-support-orbit-hint-v1"
            or metadata["output"]["sha256"] != sha256(output_path_hint)
            or metadata["output"]["exact_excluded_support_match"] is not False
        ):
            raise AuditError(f"{metadata_path}: invalid orbit hint")
        classification = exact_classification(
            output_path_hint,
            target_grams[target_name],
            known_by_target[target_name],
        )
        if (
            not classification["known_h_class"]
            or classification["exact_known_support_matches_h"]
        ):
            raise AuditError(f"{output_path_hint}: unexpected hint classification")
        hint_records.append(
            {
                "name": hint_name,
                "target": target_name,
                "metadata_path": display_path(metadata_path),
                "metadata_sha256": sha256(metadata_path),
                "found_step": metadata["walk"]["found_step"],
                "classification": classification,
            }
        )

    run_records = []
    feasible_classifications = []
    for run_name, target_name, expected_status, expected_feasible in RUN_SPECS:
        run_dir = campaign_dir / run_name
        metadata_path = run_dir / "solver.json"
        metadata = read_json(metadata_path)
        solver = metadata.get("solver", {})
        if (
            metadata.get("engine") != "order22-gram-factor-cpsat-v3"
            or solver.get("status") != expected_status
            or solver.get("feasible") is not expected_feasible
            or (metadata.get("output") is None)
            != (not expected_feasible)
        ):
            raise AuditError(f"{metadata_path}: unexpected solver result")
        record = {
            "name": run_name,
            "target": target_name,
            "metadata_path": display_path(metadata_path),
            "metadata_sha256": sha256(metadata_path),
            "status": solver["status"],
            "feasible": solver["feasible"],
            "branches": solver["branches"],
            "conflicts": solver["conflicts"],
            "wall_time_seconds": solver["wall_time_seconds"],
            "requested_orbit_count_slice": metadata["inputs"][
                "requested_orbit_count_slice"
            ],
            "required_masks": metadata["inputs"]["required_masks"],
            "objective_value": solver["objective_value"],
            "best_objective_bound": solver["best_objective_bound"],
            "exact_overlap": solver["exact_overlap"],
            "exact_random_objective": solver.get("exact_random_objective"),
        }
        if expected_feasible:
            matrix_path = run_dir / "factor.matrix.txt"
            if metadata["output"]["sha256"] != sha256(matrix_path):
                raise AuditError(f"{metadata_path}: output hash mismatch")
            classification = exact_classification(
                matrix_path,
                target_grams[target_name],
                known_by_target[target_name],
            )
            if (
                not classification["known_h_class"]
                or not classification["known_ht_class"]
                or classification["exact_known_support_matches_h"]
            ):
                raise AuditError(
                    f"{matrix_path}: expected a transformed known-class support"
                )
            record["classification"] = classification
            feasible_classifications.append(classification)
        run_records.append(record)

    stabilizer_records = []
    histograms: dict[str, Counter[tuple[int, int]]] = {
        "gsds": Counter(),
        "mendeley": Counter(),
    }
    for target_name, records in sorted(known_by_target.items()):
        gram_group_order = target_group_orders[target_name]
        for record in records:
            factor_path = (
                alignment_path.parent / record["aligned_factor_path"]
            ).resolve()
            stabilizer_order = group_order(
                projective_factor_graph(read_matrix(factor_path))
            )
            if gram_group_order % stabilizer_order:
                raise AuditError("factor stabilizer does not divide Gram group")
            support_orbit_size = gram_group_order // stabilizer_order
            histograms[target_name][
                (stabilizer_order, support_orbit_size)
            ] += 1
            stabilizer_records.append(
                {
                    "target": target_name,
                    "h_certificate_sha256": record[
                        "h_certificate_sha256"
                    ],
                    "factor_path": display_path(factor_path),
                    "factor_sha256": sha256(factor_path),
                    "factor_stabilizer_order": stabilizer_order,
                    "exact_support_orbit_size": support_orbit_size,
                }
            )

    status_histogram = Counter(record["status"] for record in run_records)
    result = {
        "engine": "order22-gram-orbit-slice-audit-v1",
        "order": ORDER,
        "claim": (
            "Two complete signed-Gram shell orbit partitions, three exact "
            "automorphic hints, and seven bounded CP-SAT slice runs were hash-"
            "bound and independently checked. All five feasible solver outputs "
            "have exact target Gram/determinant and belong to known H/HT classes."
        ),
        "claim_boundary": (
            "Observed known-factor orbit vectors define restricted search "
            "families unless forced by the invariant system. UNKNOWN is not "
            "infeasibility. Exact-support cuts do not quotient H-equivalence; "
            "no absence, global classification, or frontier claim follows."
        ),
        "inputs": {
            "alignment_manifest": {
                "path": display_path(alignment_path),
                "sha256": sha256(alignment_path),
            },
            "orbit_reports": orbit_records,
        },
        "hints": hint_records,
        "runs": run_records,
        "support_orbit_audit": {
            "method": (
                "full colored signed row/column factor graph; exact orbit-"
                "stabilizer ratio against the signed-double-cover Gram group"
            ),
            "records": stabilizer_records,
            "histograms": {
                target: [
                    {
                        "factor_stabilizer_order": stabilizer,
                        "exact_support_orbit_size": orbit_size,
                        "h_class_count": count,
                    }
                    for (stabilizer, orbit_size), count in sorted(histogram.items())
                ]
                for target, histogram in histograms.items()
            },
        },
        "summary": {
            "run_count": len(run_records),
            "solver_status_histogram": dict(sorted(status_histogram.items())),
            "feasible_output_count": len(feasible_classifications),
            "novel_h_class_count": sum(
                not record["known_h_class"]
                for record in feasible_classifications
            ),
            "novel_ht_class_count": sum(
                not record["known_ht_class"]
                for record in feasible_classifications
            ),
            "exact_known_support_match_count": sum(
                bool(record["exact_known_support_matches_h"])
                for record in feasible_classifications
            ),
        },
        "dependencies": {
            "pynauty": getattr(pynauty, "__version__", "unknown"),
        },
        "sources": [
            {
                "path": display_path(path),
                "sha256": sha256(path),
            }
            for path in (
                ROOT / "research/order22_gram_shell_orbits.py",
                ROOT / "research/order22_gram_factor_cpsat.py",
                ROOT / "research/order22_gram_support_orbit_hint.py",
                Path(__file__).resolve(),
            )
        ],
    }
    contents = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode(
        "utf-8"
    )
    atomic_write(output_path, contents)
    print(
        f"runs={len(run_records)} feasible={len(feasible_classifications)} "
        f"novel_h={result['summary']['novel_h_class_count']} "
        f"novel_ht={result['summary']['novel_ht_class_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
