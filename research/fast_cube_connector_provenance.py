#!/usr/bin/env python3
"""Create or verify a durable sidecar for a completed exact connector cube."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Sequence

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from maxdet.exact import bareiss_determinant
from research.fast_cube_batch import (
    atomic_json,
    atomic_write,
    canonical_matrix_bytes,
    read_matrix,
    sha256_bytes,
)
from research.fast_cube_lnps import (
    affine_fingerprint,
    apply_mask,
    parse_coordinate_file,
)

PROVENANCE_NAME = "provenance.json"
SIDECAR_HASH_NAME = "provenance.sha256"
INVENTORY_TAG = b"maxdet-connector-artifact-inventory-v1\0"

Coordinate = tuple[int, int]
Matrix = tuple[tuple[int, ...], ...]


def resolve_repository_path(path: Path) -> Path:
    resolved = (
        path.resolve()
        if path.is_absolute()
        else (REPOSITORY_ROOT / path).resolve()
    )
    resolved.relative_to(REPOSITORY_ROOT)
    return resolved


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPOSITORY_ROOT))


def read_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return value


def file_record(path: Path) -> dict[str, str | int]:
    return {
        "path": relative(path),
        "raw_sha256": sha256_bytes(path.read_bytes()),
        "size_bytes": path.stat().st_size,
    }


def inventory(directory: Path) -> tuple[list[dict[str, str | int]], str]:
    excluded = {PROVENANCE_NAME, SIDECAR_HASH_NAME}
    records = [
        file_record(path)
        for path in sorted(directory.rglob("*"))
        if path.is_file() and path.name not in excluded
    ]
    encoded = json.dumps(
        records, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    digest = sha256_bytes(INVENTORY_TAG + encoded)
    return records, digest


def differing_coordinates(first: Matrix, second: Matrix) -> tuple[Coordinate, ...]:
    if len(first) != len(second):
        raise ValueError("connector endpoints have different orders")
    return tuple(
        (row, column)
        for row in range(len(first))
        for column in range(len(first))
        if first[row][column] != second[row][column]
    )


def validate_completed_connector(
    run_directory: Path,
    input_directory: Path,
    target_path: Path,
    source_path: Path,
    alignment_path: Path,
) -> dict[str, Any]:
    aggregate_path = run_directory / "aggregate-report.json"
    manifest_path = run_directory / "manifest.json"
    aggregate = read_object(aggregate_path)
    manifest = read_object(manifest_path)
    if (
        aggregate.get("complete") is not True
        or aggregate.get("reason") != "complete"
        or int(aggregate.get("completed_leaves", 0)) != 32
        or int(aggregate.get("total_assignments", 0)) != 1 << 32
        or int(aggregate.get("tie_mask_truncation_leaf_count", -1)) != 0
    ):
        raise ValueError("connector aggregate is not a complete exact 32-cube")
    if aggregate.get("full_affine_fingerprint_sha256") != manifest.get(
        "full_affine_fingerprint_sha256"
    ):
        raise ValueError("manifest and aggregate fingerprints disagree")

    alignment = read_object(alignment_path)
    if alignment.get("complete") is not True:
        raise ValueError("alignment metadata is incomplete")
    aligned_path = resolve_repository_path(Path(str(alignment["output"])))
    if resolve_repository_path(
        Path(str(alignment["target"]))
    ) != target_path:
        raise ValueError("alignment target does not match requested target")
    if resolve_repository_path(
        Path(str(alignment["source"]))
    ) != source_path:
        raise ValueError("alignment source does not match requested source")
    for path, hash_key in (
        (target_path, "target_raw_sha256"),
        (source_path, "source_raw_sha256"),
        (aligned_path, "output_raw_sha256"),
    ):
        if sha256_bytes(path.read_bytes()) != alignment[hash_key]:
            raise ValueError(f"alignment hash mismatch: {path}")

    target = read_matrix(target_path)
    source = read_matrix(source_path)
    aligned = read_matrix(aligned_path)
    target_determinant = bareiss_determinant(target)
    source_determinant = bareiss_determinant(source)
    aligned_determinant = bareiss_determinant(aligned)
    if abs(source_determinant) != abs(aligned_determinant):
        raise ValueError("alignment did not preserve the source determinant")
    if (
        abs(source_determinant)
        != int(alignment["source_absolute_determinant"])
        or abs(aligned_determinant)
        != int(alignment["aligned_absolute_determinant"])
    ):
        raise ValueError("alignment determinant metadata disagrees")

    inner_path = input_directory / "inner27.coords.txt"
    outer_path = input_directory / "outer5.coords.txt"
    support_path = input_directory / "support32.coords.txt"
    inner = parse_coordinate_file(inner_path)
    outer = parse_coordinate_file(outer_path)
    support = parse_coordinate_file(support_path)
    if (
        len(inner) != 27
        or len(outer) != 5
        or support != inner + outer
        or len(set(support)) != 32
    ):
        raise ValueError("input supports do not form one unique 27+5 support")
    if support != differing_coordinates(target, aligned):
        raise ValueError("support is not the ordered endpoint difference set")
    if apply_mask(target, support, (1 << 32) - 1) != aligned:
        raise ValueError("all-ones connector mask missed the aligned endpoint")
    manifest_support = tuple(
        (int(row) - 1, int(column) - 1)
        for row, column in manifest["full_support"]
    )
    if manifest_support != support:
        raise ValueError("manifest support disagrees with connector inputs")
    run_support = parse_coordinate_file(
        run_directory / "support32.coords.txt"
    )
    if run_support != support:
        raise ValueError("run support disagrees with connector inputs")
    fingerprint = affine_fingerprint(target, support)
    if fingerprint != manifest["full_affine_fingerprint_sha256"]:
        raise ValueError("recomputed connector fingerprint disagrees")

    leaves = sorted((run_directory / "leaves").glob("leaf-*"))
    if len(leaves) != 32:
        raise ValueError("connector does not contain exactly 32 leaf directories")
    outer_masks: set[int] = set()
    leaf_fingerprints: set[str] = set()
    assignment_visits = 0
    for leaf in leaves:
        report = read_object(leaf / "report.json")
        summary = read_object(leaf / "partition-summary.json")
        if (
            report.get("complete") is not True
            or report.get("all_assignments_bound_checked") is not True
            or int(report.get("dimension", 0)) != 27
            or int(report.get("assignments", 0)) != 1 << 27
        ):
            raise ValueError(f"incomplete evaluator leaf: {leaf}")
        leaf_start = read_matrix(leaf / "start.matrix.txt")
        leaf_support = parse_coordinate_file(leaf / "support.coords.txt")
        if leaf_support != inner:
            raise ValueError(f"leaf support changed: {leaf}")
        if report["start_parsed_matrix_sha256"] != sha256_bytes(
            canonical_matrix_bytes(leaf_start)
        ):
            raise ValueError(f"leaf start hash mismatch: {leaf}")
        if report["coordinate_file_raw_sha256"] != sha256_bytes(
            (leaf / "support.coords.txt").read_bytes()
        ):
            raise ValueError(f"leaf support hash mismatch: {leaf}")
        outer_mask = int(summary["outer_mask_decimal"])
        if outer_mask in outer_masks:
            raise ValueError("duplicate outer-mask leaf")
        outer_masks.add(outer_mask)
        expected_start = apply_mask(target, outer, outer_mask)
        if int(summary["reroot_xor_mask_decimal"]) != 0:
            raise ValueError("this pinned connector unexpectedly rerooted")
        if leaf_start != expected_start:
            raise ValueError(f"leaf does not match fixed outer mask: {leaf}")
        leaf_fingerprint = affine_fingerprint(leaf_start, inner)
        if leaf_fingerprint != summary["leaf_affine_fingerprint_sha256"]:
            raise ValueError(f"leaf fingerprint mismatch: {leaf}")
        leaf_fingerprints.add(leaf_fingerprint)
        best = read_matrix(leaf / "best.matrix.txt")
        if abs(bareiss_determinant(best)) != int(
            report["best_absolute_determinant"]
        ):
            raise ValueError(f"leaf best determinant mismatch: {leaf}")
        assignment_visits += int(report["assignments"])
    if outer_masks != set(range(32)) or len(leaf_fingerprints) != 32:
        raise ValueError("leaf partition is not a unique full outer-mask cover")
    if assignment_visits != 1 << 32:
        raise ValueError("leaf assignment total is not 2^32")

    for record in aggregate["global_top_k"]:
        matrix = read_matrix(resolve_repository_path(Path(record["artifact"])))
        if bareiss_determinant(matrix) != int(record["signed_determinant"]):
            raise ValueError("global top-K determinant mismatch")
    for record in aggregate["global_ties"]:
        artifact = resolve_repository_path(Path(record["artifact"]))
        matrix = read_matrix(artifact)
        if abs(bareiss_determinant(matrix)) != int(
            aggregate["frontier"]
        ):
            raise ValueError("global frontier tie determinant mismatch")
        verification = artifact.with_suffix(".arena-verify.txt")
        if not verification.read_text().startswith("VERIFIED\n"):
            raise ValueError("global frontier tie lacks arena verification")

    return {
        "aligned_endpoint": file_record(aligned_path),
        "alignment_metadata": file_record(alignment_path),
        "assignment_visits": assignment_visits,
        "endpoint_difference": len(support),
        "endpoint_exact_reconstruction": True,
        "full_affine_fingerprint_sha256": fingerprint,
        "global_frontier_ties": len(aggregate["global_ties"]),
        "global_top_k_replayed": len(aggregate["global_top_k"]),
        "leaf_fingerprints_unique": len(leaf_fingerprints),
        "source": {
            **file_record(source_path),
            "signed_determinant": str(source_determinant),
        },
        "target": {
            **file_record(target_path),
            "signed_determinant": str(target_determinant),
        },
    }


def build_sidecar(arguments: argparse.Namespace) -> dict[str, Any]:
    run_directory = resolve_repository_path(arguments.run_dir)
    input_directory = resolve_repository_path(arguments.input_dir)
    target_path = resolve_repository_path(arguments.target)
    source_path = resolve_repository_path(arguments.source)
    alignment_path = resolve_repository_path(arguments.alignment)
    sidecar_path = run_directory / PROVENANCE_NAME
    sidecar_hash_path = run_directory / SIDECAR_HASH_NAME
    if sidecar_path.exists() or sidecar_hash_path.exists():
        if not arguments.refresh:
            raise FileExistsError(
                "connector provenance sidecar already exists"
            )
        if not sidecar_path.is_file() or not sidecar_hash_path.is_file():
            raise ValueError("connector provenance sidecar is incomplete")
        expected = sidecar_hash_path.read_text().split()[0]
        if sha256_bytes(sidecar_path.read_bytes()) != expected:
            raise ValueError(
                "refusing to refresh a sidecar with a stale hash"
            )

    connector = validate_completed_connector(
        run_directory,
        input_directory,
        target_path,
        source_path,
        alignment_path,
    )
    run_records, run_digest = inventory(run_directory)
    input_records, input_digest = inventory(input_directory)
    manifest = read_object(run_directory / "manifest.json")
    sidecar = {
        "artifact_inventory": {
            "algorithm":
                "SHA256(tag || canonical compact JSON of path/hash/size)",
            "count": len(run_records),
            "sha256": run_digest,
        },
        "connector_validation": connector,
        "input_inventory": {
            "algorithm":
                "SHA256(tag || canonical compact JSON of path/hash/size)",
            "count": len(input_records),
            "sha256": input_digest,
        },
        "manifest": file_record(run_directory / "manifest.json"),
        "aggregate_report": file_record(
            run_directory / "aggregate-report.json"
        ),
        "pinned_run_provenance": manifest["provenance"],
        "result_changed": False,
        "run_directory": relative(run_directory),
        "schema_version": 1,
        "tool": file_record(Path(__file__).resolve()),
    }
    atomic_json(sidecar_path, sidecar)
    sidecar_hash = sha256_bytes(sidecar_path.read_bytes())
    atomic_write(
        sidecar_hash_path,
        f"{sidecar_hash}  {PROVENANCE_NAME}\n".encode("ascii"),
    )
    return sidecar


def verify_sidecar(arguments: argparse.Namespace) -> dict[str, Any]:
    run_directory = resolve_repository_path(arguments.run_dir)
    input_directory = resolve_repository_path(arguments.input_dir)
    sidecar_path = run_directory / PROVENANCE_NAME
    sidecar_hash_path = run_directory / SIDECAR_HASH_NAME
    sidecar = read_object(sidecar_path)
    expected_hash = sidecar_hash_path.read_text().split()[0]
    actual_hash = sha256_bytes(sidecar_path.read_bytes())
    if expected_hash != actual_hash:
        raise ValueError("connector provenance sidecar hash mismatch")
    run_records, run_digest = inventory(run_directory)
    input_records, input_digest = inventory(input_directory)
    if (
        len(run_records) != int(sidecar["artifact_inventory"]["count"])
        or run_digest != sidecar["artifact_inventory"]["sha256"]
        or len(input_records) != int(sidecar["input_inventory"]["count"])
        or input_digest != sidecar["input_inventory"]["sha256"]
    ):
        raise ValueError("connector artifact/input inventory changed")
    validate_completed_connector(
        run_directory,
        input_directory,
        resolve_repository_path(arguments.target),
        resolve_repository_path(arguments.source),
        resolve_repository_path(arguments.alignment),
    )
    return {
        "artifact_inventory_sha256": run_digest,
        "input_inventory_sha256": input_digest,
        "sidecar_raw_sha256": actual_hash,
        "verified": True,
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--create", action="store_true")
    mode.add_argument("--verify", action="store_true")
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="atomically replace an existing valid sidecar in create mode",
    )
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--alignment", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.refresh and not arguments.create:
        parser.error("--refresh requires --create")
    return arguments


def main() -> int:
    arguments = parse_arguments()
    result = (
        build_sidecar(arguments)
        if arguments.create
        else verify_sidecar(arguments)
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
