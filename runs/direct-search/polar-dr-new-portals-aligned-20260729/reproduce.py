#!/usr/bin/env python3
"""Align the two new polar/DR frontier classes to the published Gram frame.

This deliberately reuses the pinned exact alignment primitives from the prior
``known-ht-aligned/reproduce.py`` audit.  In addition to that audit's Gram
isomorphism witness, this driver records and replays the complete second-stage
row-automorphism, per-column sign, and column-permutation transport into the
sorted canonical shell-factor representation.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any


ORDER = 23
FRONTIER = 2_779_447_296_000_000
REPRESENTATIVES = (
    {
        "label": "polar-dr-de764",
        "path": (
            "runs/direct-search/polar-dr-portals-20260729/frontier-ties/"
            "8241af155359ed39fc39ea8c53533ba63f974c4b179aa890628fc5e27b211208"
            ".matrix.txt"
        ),
        "expected_h": (
            "de7642266b6997198a6892fe26c1e080d0c673f095f9710b5003bbb8b5de2345"
        ),
        "expected_transpose_h": (
            "de7642266b6997198a6892fe26c1e080d0c673f095f9710b5003bbb8b5de2345"
        ),
        "expected_ht": (
            "de7642266b6997198a6892fe26c1e080d0c673f095f9710b5003bbb8b5de2345"
        ),
    },
    {
        "label": "polar-dr-1e4b",
        "path": (
            "runs/direct-search/polar-dr-portals-20260729/frontier-ties/"
            "85e5efc9c4b5bfa8df951ee344f32b212e83384f5f568560c9f0cbedcb235e03"
            ".matrix.txt"
        ),
        "expected_h": (
            "1e4b14334f15ec751275a8383b5b7e8aa8270c951e3ff5a3d342c8fe203d9efe"
        ),
        "expected_transpose_h": (
            "411ffbb1f7635c339ca3de9d970062026e50b371571426515c98e867d29e4f3a"
        ),
        "expected_ht": (
            "1e4b14334f15ec751275a8383b5b7e8aa8270c951e3ff5a3d342c8fe203d9efe"
        ),
    },
)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_base(path: Path) -> Any:
    specification = importlib.util.spec_from_file_location(
        "known_ht_alignment_base", path
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load alignment base: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def replay_initial_alignment(
    source: list[list[int]],
    row_signs: list[int],
    row_mapping: list[int],
    column_signs: list[int],
) -> list[list[int]]:
    aligned = [[0] * ORDER for _ in range(ORDER)]
    for source_row, target_row in enumerate(row_mapping):
        aligned[target_row] = [
            row_signs[source_row] * value for value in source[source_row]
        ]
    for column, sign in enumerate(column_signs):
        for row in range(ORDER):
            aligned[row][column] *= sign
    return aligned


def full_aut_transport(
    aligned: list[list[int]],
    row_permutation: list[int],
) -> tuple[list[list[int]], list[int], list[int], list[int]]:
    """Replay Aut(G), column normalization, and sorted-mask reindexing."""
    permuted = [[0] * ORDER for _ in range(ORDER)]
    for source_row, target_row in enumerate(row_permutation):
        permuted[target_row] = aligned[source_row][:]
    column_signs = permuted[0][:]
    masks: list[int] = []
    normalized_columns: list[list[int]] = []
    for column, sign in enumerate(column_signs):
        values = [sign * permuted[row][column] for row in range(ORDER)]
        normalized_columns.append(values)
        masks.append(
            sum((value == 1) << row for row, value in enumerate(values))
        )
    sorted_masks = sorted(masks)
    if len(set(sorted_masks)) != ORDER:
        raise ArithmeticError("Aut(G) replay produced duplicate shell columns")
    mask_to_column = {mask: index for index, mask in enumerate(sorted_masks)}
    column_mapping = [mask_to_column[mask] for mask in masks]
    factor = [[0] * ORDER for _ in range(ORDER)]
    for source_column, target_column in enumerate(column_mapping):
        for row in range(ORDER):
            factor[row][target_column] = normalized_columns[source_column][row]
    return factor, column_signs, column_mapping, sorted_masks


def certificate_triplet(
    matrix: list[list[int]], h_certificate: Any, transpose: Any, sha256_hex: Any
) -> tuple[str, str, str]:
    direct = h_certificate(matrix)
    transposed = h_certificate(transpose(matrix))
    return (
        sha256_hex(direct),
        sha256_hex(transposed),
        sha256_hex(min(direct, transposed)),
    )


def run(arguments: argparse.Namespace) -> int:
    script_path = Path(__file__).resolve()
    repo = script_path.parents[3]
    os.chdir(repo)
    output_dir = arguments.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "alignment-report.json"
    if report_path.exists():
        raise FileExistsError(f"refusing to overwrite {report_path}")

    base_path = (
        repo
        / "runs/direct-search/frontier-factor-class-expansion-20260728"
        / "known-ht-aligned/reproduce.py"
    )
    base = load_base(base_path)
    if getattr(base.pynauty, "__version__", "unknown") != "2.8.8.1":
        raise RuntimeError("alignment requires pynauty==2.8.8.1")
    if base.nx.__version__ != "3.2.1":
        raise RuntimeError("alignment requires networkx==3.2.1")

    sys.path.insert(0, str(repo / "research"))
    from h_equivalence_audit import (  # pylint: disable=import-outside-toplevel
        determinant,
        h_certificate,
        read_matrix,
        sha256_hex,
        transpose,
    )

    target_path = arguments.target.resolve()
    target_raw, target = read_matrix(target_path)
    target_gram = base.gram(target)
    target_graph = base.defect_graph(target_gram)
    target_determinant = abs(determinant(target))
    target_gram_determinant = determinant(target_gram)
    if (
        target_determinant != FRONTIER
        or target_gram_determinant != FRONTIER * FRONTIER
    ):
        raise ArithmeticError("target determinant/Gram identity failed")

    shell_path = arguments.shell_report.resolve()
    shell_payload = json.loads(shell_path.read_text(encoding="utf-8"))
    shell_results = shell_payload.get("results")
    if (
        not isinstance(shell_results, list)
        or len(shell_results) != 1
        or not isinstance(shell_results[0], dict)
    ):
        raise ValueError("shell report must contain exactly one result")
    shell_masks_value = shell_results[0].get("shell_sign_masks")
    if (
        not isinstance(shell_masks_value, list)
        or len(shell_masks_value) != 1382
        or len(set(shell_masks_value)) != 1382
    ):
        raise ValueError("shell report lacks 1,382 unique masks")
    shell_masks = set(shell_masks_value)
    if shell_results[0].get("source_sha256") != sha256_bytes(target_raw):
        raise ValueError("shell report is not bound to target")
    if shell_results[0].get("determinant") != str(FRONTIER * FRONTIER):
        raise ValueError("shell report has wrong Gram determinant")

    generators, group_order = base.gram_automorphism_generators(target_gram)
    entries: list[dict[str, Any]] = []
    for representative in REPRESENTATIVES:
        source_path = (repo / representative["path"]).resolve()
        source_raw, source = read_matrix(source_path)
        source_score = abs(determinant(source))
        if source_score != FRONTIER:
            raise ArithmeticError(f"{source_path}: source is not a frontier tie")
        direct_h, transpose_h, ht = certificate_triplet(
            source, h_certificate, transpose, sha256_hex
        )
        if (
            direct_h != representative["expected_h"]
            or transpose_h != representative["expected_transpose_h"]
            or ht != representative["expected_ht"]
        ):
            raise ArithmeticError(
                f"{source_path}: source H/HT certificates changed"
            )

        aligned, row_signs, row_mapping, initial_column_signs = (
            base.align_to_target_gram(source, target_gram, target_graph)
        )
        if (
            replay_initial_alignment(
                source, row_signs, row_mapping, initial_column_signs
            )
            != aligned
        ):
            raise ArithmeticError("initial alignment witness failed replay")
        aligned_support = set(base.normalized_column_masks(aligned))
        if (
            len(aligned_support) != ORDER
            or not aligned_support.issubset(shell_masks)
        ):
            raise ArithmeticError("aligned factor is outside exact shell")
        aligned_triple = base.triple_indices(aligned_support)

        (
            canonical_support,
            canonical_triple,
            generator_word,
            automorphism_permutation,
        ) = base.canonicalize_support(aligned_support, generators)
        factor = base.reconstruct_factor(canonical_support)
        (
            replayed_factor,
            aut_column_signs,
            aut_column_mapping,
            sorted_masks,
        ) = full_aut_transport(aligned, automorphism_permutation)
        if replayed_factor != factor:
            raise ArithmeticError("full Aut(G)/column witness failed replay")
        if sorted_masks != sorted(canonical_support):
            raise ArithmeticError("canonical support witness disagrees")
        if base.gram(factor) != target_gram:
            raise ArithmeticError("canonical factor misses exact target Gram")
        factor_score = abs(determinant(factor))
        if factor_score != FRONTIER:
            raise ArithmeticError("canonical factor changed determinant")
        factor_h, factor_transpose_h, factor_ht = certificate_triplet(
            factor, h_certificate, transpose, sha256_hex
        )
        if (factor_h, factor_transpose_h, factor_ht) != (
            direct_h,
            transpose_h,
            ht,
        ):
            raise ArithmeticError("canonicalization changed H/HT certificates")

        full_row_mapping = [
            automorphism_permutation[row_mapping[source_row]]
            for source_row in range(ORDER)
        ]
        stem = f"ht-{ht[:12]}"
        matrix_path = output_dir / f"{stem}.matrix.txt"
        receipt_path = output_dir / f"{stem}.receipt.json"
        if matrix_path.exists() or receipt_path.exists():
            raise FileExistsError(f"refusing to overwrite {stem}")
        base.atomic_write(matrix_path, base.matrix_bytes(factor))
        completed = subprocess.run(
            [
                str(arguments.arena.resolve()),
                "verify",
                str(matrix_path),
                "--json",
                str(receipt_path),
                "--quiet",
            ],
            cwd=repo,
            check=True,
            capture_output=True,
            text=True,
        )
        if completed.stdout.strip() != str(FRONTIER):
            raise ArithmeticError("arena verifier returned wrong score")
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if (
            receipt["score"]["absolute_determinant"] != str(FRONTIER)
            or receipt["matrix"]["raw_sha256"] != sha256_file(matrix_path)
        ):
            raise ArithmeticError("arena receipt failed matrix/score binding")

        entries.append(
            {
                "label": representative["label"],
                "source_path": base.display_path(repo, source_path),
                "source_matrix_sha256": sha256_bytes(source_raw),
                "source_absolute_determinant": str(source_score),
                "source_h_certificate_sha256": direct_h,
                "source_transpose_h_certificate_sha256": transpose_h,
                "ht_certificate_sha256": ht,
                "normalized_gram_pynauty_certificate_sha256": (
                    sha256_bytes(
                        base.pynauty_graph_certificate(
                            base.defect_graph(base.normalize_gram(source)[1])
                        )
                    )
                ),
                "initial_exact_graph_isomorphism_source_to_target": (
                    row_mapping
                ),
                "initial_source_row_signs": row_signs,
                "initial_column_signs": initial_column_signs,
                "initial_alignment_witness_replayed_exactly": True,
                "aligned_small_orbit_triple": list(aligned_triple),
                "aligned_small_orbit_masks": [
                    base.PUBLISHED_SMALL_ORBIT[index]
                    for index in aligned_triple
                ],
                "canonical_small_orbit_triple": list(canonical_triple),
                "canonical_small_orbit_masks": [
                    base.PUBLISHED_SMALL_ORBIT[index]
                    for index in canonical_triple
                ],
                "aut_g_generator_word": generator_word,
                "aut_g_row_permutation": automorphism_permutation,
                "full_source_row_mapping_to_output": full_row_mapping,
                "aut_stage_column_signs": aut_column_signs,
                "full_source_column_mapping_to_output": aut_column_mapping,
                "selected_shell_masks": sorted(canonical_support),
                "full_transport_witness_replayed_exactly": True,
                "matrix_path": base.display_path(repo, matrix_path),
                "matrix_sha256": sha256_file(matrix_path),
                "matrix_h_certificate_sha256": factor_h,
                "matrix_transpose_h_certificate_sha256": (
                    factor_transpose_h
                ),
                "matrix_ht_certificate_sha256": factor_ht,
                "exact_row_gram_matches_target": True,
                "exact_absolute_determinant": str(factor_score),
                "arena_verified_score": completed.stdout.strip(),
                "arena_receipt_path": base.display_path(repo, receipt_path),
                "arena_receipt_file_sha256": sha256_file(receipt_path),
                "arena_receipt_sha256": receipt["receipt_sha256"],
            }
        )

    report = {
        "schema_version": 1,
        "complete": True,
        "claim_boundary": (
            "These are two new classes relative to the frozen local frontier "
            "audit, not a literature-novelty, optimality, or world-record claim."
        ),
        "method": (
            "reuse of pinned known-ht exact defect-graph isomorphism and "
            "Aut(G) shell canonicalization, plus explicit full row/column "
            "monomial witness replay"
        ),
        "command": [sys.executable, *sys.argv],
        "driver_path": base.display_path(repo, script_path),
        "driver_sha256": sha256_file(script_path),
        "alignment_base_path": base.display_path(repo, base_path),
        "alignment_base_sha256": sha256_file(base_path),
        "python": {
            "executable": sys.executable,
            "version": sys.version,
            "platform": platform.platform(),
        },
        "pynauty_version": getattr(base.pynauty, "__version__", "unknown"),
        "networkx_version": base.nx.__version__,
        "target_matrix_path": base.display_path(repo, target_path),
        "target_matrix_sha256": sha256_bytes(target_raw),
        "target_absolute_determinant": str(target_determinant),
        "target_gram_determinant": str(target_gram_determinant),
        "target_gram_sha256": sha256_bytes(
            json.dumps(target_gram, separators=(",", ":")).encode("ascii")
        ),
        "target_gram_pynauty_certificate_sha256": sha256_bytes(
            base.pynauty_graph_certificate(target_graph)
        ),
        "target_gram_automorphism_group_order": group_order,
        "target_gram_automorphism_generator_count": len(generators),
        "shell_report_path": base.display_path(repo, shell_path),
        "shell_report_sha256": sha256_file(shell_path),
        "shell_size": len(shell_masks),
        "expected_ht_class_count": len(REPRESENTATIVES),
        "aligned_ht_class_count": len(entries),
        "all_initial_witnesses_replayed_exactly": all(
            entry["initial_alignment_witness_replayed_exactly"]
            for entry in entries
        ),
        "all_full_transport_witnesses_replayed_exactly": all(
            entry["full_transport_witness_replayed_exactly"]
            for entry in entries
        ),
        "all_aligned_factors_have_exact_target_gram": all(
            entry["exact_row_gram_matches_target"] for entry in entries
        ),
        "all_aligned_factors_arena_verified": all(
            entry["arena_verified_score"] == str(FRONTIER)
            for entry in entries
        ),
        "observed_canonical_triples": sorted(
            {
                tuple(entry["canonical_small_orbit_triple"])
                for entry in entries
            }
        ),
        "representatives": entries,
    }
    base.atomic_write(
        report_path,
        (json.dumps(report, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    report_sidecar = report_path.with_suffix(".json.sha256")
    base.atomic_write(
        report_sidecar,
        f"{sha256_file(report_path)}  {report_path.name}\n".encode("ascii"),
    )

    artifacts = [report_path, report_sidecar, script_path]
    for entry in entries:
        artifacts.extend(
            [
                repo / entry["matrix_path"],
                repo / entry["arena_receipt_path"],
            ]
        )
    provenance = {
        "schema_version": 1,
        "complete": True,
        "command": [sys.executable, *sys.argv],
        "artifacts": [
            {
                "path": base.display_path(repo, path),
                "sha256": sha256_file(path),
            }
            for path in artifacts
        ],
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip(),
    }
    provenance_path = output_dir / "provenance.json"
    base.atomic_write(
        provenance_path,
        (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode(
            "utf-8"
        ),
    )
    base.atomic_write(
        provenance_path.with_suffix(".json.sha256"),
        f"{sha256_file(provenance_path)}  {provenance_path.name}\n".encode(
            "ascii"
        ),
    )
    print(
        json.dumps(
            {
                "aligned_ht_classes": len(entries),
                "canonical_triples": report["observed_canonical_triples"],
                "report": base.display_path(repo, report_path),
                "report_sha256": sha256_file(report_path),
            },
            sort_keys=True,
        )
    )
    return 0


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "runs/direct-search/polar-dr-new-portals-aligned-20260729"
        ),
    )
    parser.add_argument(
        "--target",
        type=Path,
        default=Path("references/orrick-et-al-2003/matrix.txt"),
    )
    parser.add_argument(
        "--shell-report",
        type=Path,
        default=Path(
            "runs/direct-search/gram-shell-reference-columns-29760.json"
        ),
    )
    parser.add_argument("--arena", type=Path, default=Path("./arena"))
    return parser.parse_args()


if __name__ == "__main__":
    try:
        raise SystemExit(run(parse_arguments()))
    except (
        ArithmeticError,
        FileExistsError,
        OSError,
        RuntimeError,
        subprocess.SubprocessError,
        ValueError,
    ) as error:
        raise SystemExit(f"polar_dr_new_portal_align: {error}") from error
