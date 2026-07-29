#!/usr/bin/env python3
"""Close a frontier-factor corpus under transpose and retain new H-classes.

This is deliberately separate from the trusted arena verifier.  H and H-plus-
transpose classification uses the complete 23^2-pivot pynauty certificate from
``h_equivalence_audit.py``.  Every retained transpose is then independently
accepted by ``arena verify`` before its temporary file is promoted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path

from h_equivalence_audit import (
    FRONTIER,
    determinant,
    h_certificate,
    normalized_gram_graph,
    read_matrix,
    sha256_hex,
    transpose,
)


def matrix_bytes(matrix: list[list[int]]) -> bytes:
    return "".join(
        " ".join(map(str, row)) + "\n" for row in matrix
    ).encode("ascii")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def expand_globs(patterns: list[str]) -> list[Path]:
    paths = []
    for pattern in patterns:
        matches = sorted(Path(".").glob(pattern))
        if not matches:
            raise ValueError(f"glob matched no files: {pattern}")
        paths.extend(matches)
    return paths


def unique_paths(paths: list[Path]) -> list[Path]:
    result: dict[Path, None] = {}
    for path in paths:
        resolved = path.resolve()
        if not resolved.is_file():
            raise ValueError(f"missing matrix: {path}")
        result[resolved] = None
    return list(result)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--baseline", action="append", default=[], type=Path
    )
    parser.add_argument("--baseline-glob", action="append", default=[])
    parser.add_argument("--source", action="append", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--arena", default=Path("./arena"), type=Path)
    arguments = parser.parse_args()

    if arguments.report.exists():
        parser.error("refusing to overwrite --report")
    baseline_paths = unique_paths(
        arguments.baseline + expand_globs(arguments.baseline_glob)
    )
    source_paths = unique_paths(arguments.source)
    if not baseline_paths:
        parser.error("provide --baseline or --baseline-glob")
    arguments.output_dir.mkdir(parents=True, exist_ok=True)

    started = time.monotonic()
    baseline_h: set[bytes] = set()
    baseline_ht: set[bytes] = set()
    baseline_raw: set[str] = set()
    for path in baseline_paths:
        raw, matrix = read_matrix(path)
        score = abs(determinant(matrix))
        if score != FRONTIER:
            raise ValueError(f"{path}: |det|={score}, expected {FRONTIER}")
        direct = h_certificate(matrix)
        transposed = h_certificate(transpose(matrix))
        baseline_h.add(direct)
        baseline_ht.add(min(direct, transposed))
        baseline_raw.add(sha256_hex(raw))

    retained_h = set(baseline_h)
    retained = []
    observations = []
    for source_path in source_paths:
        source_raw, source = read_matrix(source_path)
        source_score = abs(determinant(source))
        if source_score != FRONTIER:
            raise ValueError(
                f"{source_path}: |det|={source_score}, expected {FRONTIER}"
            )
        candidate = transpose(source)
        candidate_raw = matrix_bytes(candidate)
        candidate_raw_sha256 = sha256_hex(candidate_raw)
        candidate_h = h_certificate(candidate)
        candidate_back = h_certificate(transpose(candidate))
        candidate_ht = min(candidate_h, candidate_back)
        gram_certificate, group_order, edges, degrees = (
            normalized_gram_graph(candidate)
        )
        novel_h_at_start = candidate_h not in baseline_h
        novel_ht_at_start = candidate_ht not in baseline_ht
        novel_h_online = candidate_h not in retained_h

        observation: dict[str, object] = {
            "source": str(source_path),
            "source_matrix_sha256": sha256_hex(source_raw),
            "transpose_matrix_sha256": candidate_raw_sha256,
            "h_certificate_sha256": sha256_hex(candidate_h),
            "ht_certificate_sha256": sha256_hex(candidate_ht),
            "novel_h_vs_baseline": novel_h_at_start,
            "novel_ht_vs_baseline": novel_ht_at_start,
            "retained_online": novel_h_online,
            "gram_certificate_sha256": sha256_hex(gram_certificate),
            "gram_automorphism_group_order": group_order,
            "gram_edge_count": edges,
            "gram_degree_multiset": list(degrees),
        }
        if novel_h_online:
            stem = f"h-{sha256_hex(candidate_h)}"
            output_path = arguments.output_dir / f"{stem}.matrix.txt"
            receipt_path = arguments.output_dir / f"{stem}.receipt.json"
            if output_path.exists() or receipt_path.exists():
                raise ValueError(f"refusing to overwrite retained seed {stem}")

            matrix_fd, matrix_name = tempfile.mkstemp(
                prefix=f".{stem}.", suffix=".matrix.tmp",
                dir=arguments.output_dir,
            )
            os.close(matrix_fd)
            receipt_fd, receipt_name = tempfile.mkstemp(
                prefix=f".{stem}.", suffix=".receipt.tmp",
                dir=arguments.output_dir,
            )
            os.close(receipt_fd)
            temporary_matrix = Path(matrix_name)
            temporary_receipt = Path(receipt_name)
            try:
                temporary_matrix.write_bytes(candidate_raw)
                verified = subprocess.run(
                    [
                        str(arguments.arena.resolve()),
                        "verify",
                        str(temporary_matrix),
                        "--json",
                        str(temporary_receipt),
                        "--quiet",
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                )
                if verified.stdout.strip() != str(FRONTIER):
                    raise RuntimeError(
                        "arena verify returned unexpected score: "
                        + verified.stdout.strip()
                    )
                os.replace(temporary_matrix, output_path)
                os.replace(temporary_receipt, receipt_path)
            finally:
                temporary_matrix.unlink(missing_ok=True)
                temporary_receipt.unlink(missing_ok=True)

            retained_h.add(candidate_h)
            retained.append(str(output_path))
            observation.update(
                {
                    "retained_matrix": str(output_path),
                    "retained_matrix_sha256": file_sha256(output_path),
                    "arena_receipt": str(receipt_path),
                    "arena_receipt_sha256": file_sha256(receipt_path),
                    "arena_verified_score": str(FRONTIER),
                }
            )
        observations.append(observation)

    report = {
        "schema_version": 1,
        "claim_boundary": (
            "H novelty is relative to the explicit baseline paths. "
            "It is not a literature novelty claim. Every retained matrix "
            "passed the trusted arena verifier."
        ),
        "method": (
            "deterministic transpose closure; minimum color-preserving "
            "pynauty certificate over all 23^2 dephased pivots"
        ),
        "pynauty_version": "2.8.8.1",
        "baseline_input_count": len(baseline_paths),
        "baseline_unique_raw_matrix_count": len(baseline_raw),
        "baseline_h_class_count": len(baseline_h),
        "baseline_ht_class_count": len(baseline_ht),
        "source_count": len(source_paths),
        "transpose_observation_count": len(observations),
        "new_h_class_count": len(retained_h) - len(baseline_h),
        "new_ht_class_count": 0,
        "final_h_class_count": len(retained_h),
        "final_ht_class_count": len(baseline_ht),
        "retained_seed_count": len(retained),
        "retained_seeds": retained,
        "elapsed_seconds": time.monotonic() - started,
        "observations": observations,
    }
    arguments.report.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{arguments.report.name}.",
        suffix=".tmp",
        dir=arguments.report.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(report, output, indent=2, sort_keys=True)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, arguments.report)
    finally:
        temporary.unlink(missing_ok=True)

    print(
        f"baseline H/HT={len(baseline_h)}/{len(baseline_ht)} "
        f"new H/HT={report['new_h_class_count']}/0 "
        f"retained={len(retained)} elapsed={report['elapsed_seconds']:.3f}s"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        raise SystemExit(f"frontier_factor_transpose_closure: {error}") from error
