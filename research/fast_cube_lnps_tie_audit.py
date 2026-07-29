#!/usr/bin/env python3
"""Audit frontier ties preserved by the iterative fast-cube LNPS pilot."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from maxdet.exact import bareiss_determinant, normalize_signs
from research.fast_cube_batch import (
    FRONTIER,
    atomic_json,
    atomic_write,
    canonical_matrix_bytes,
    read_matrix,
    sha256_bytes,
)

DEFAULT_LNPS = Path(
    "runs/direct-search/fast-principal-cube/"
    "lnps-beam-h012-20260728"
)
DEFAULT_KNOWN24 = Path(
    "runs/direct-search/neutral-cycle/two-cycle-union-29952"
)
DEFAULT_QUBO = Path(
    "runs/qubo-trust-pilot-20260728-seed31003/"
    "best-proposal.matrix.txt"
)


def normalized_hash(path: Path) -> str:
    matrix = read_matrix(path)
    normalized = tuple(
        tuple(value for value in row) for row in normalize_signs(matrix)
    )
    return sha256_bytes(canonical_matrix_bytes(normalized))


def parse_arena(output: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for key, label in (
        ("determinant", "determinant"),
        ("matrix_sha256", "matrix sha256"),
        ("normalized_sha256", "normalized sha256"),
        ("receipt_sha256", "receipt sha256"),
    ):
        match = re.search(
            rf"^{re.escape(label)}:\s*(\S+)\s*$",
            output,
            flags=re.MULTILINE,
        )
        if not match:
            raise RuntimeError(f"arena output omitted {label}")
        result[key] = match.group(1)
    return result


def relative(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root))
    except ValueError:
        return str(path.resolve())


def run(arguments: argparse.Namespace) -> int:
    root = REPOSITORY_ROOT
    lnps = (root / arguments.lnps_dir).resolve()
    known24_directory = (root / arguments.known24_dir).resolve()
    qubo_path = (root / arguments.qubo).resolve()
    tie_directory = lnps / "frontier-ties"
    tie_paths = sorted(tie_directory.glob("*.matrix.txt"))
    if not tie_paths:
        raise ValueError(f"no preserved ties under {tie_directory}")
    known24_paths = sorted(known24_directory.glob("*.matrix.txt"))
    if len(known24_paths) != 24:
        raise ValueError(
            f"expected 24 neutral-network matrices, found {len(known24_paths)}"
        )

    known_by_normalized: dict[str, list[str]] = {}
    for path in known24_paths:
        if abs(bareiss_determinant(read_matrix(path))) != FRONTIER:
            raise ValueError(f"known24 input is not frontier: {path}")
        known_by_normalized.setdefault(normalized_hash(path), []).append(
            relative(path, root)
        )
    if abs(bareiss_determinant(read_matrix(qubo_path))) != FRONTIER:
        raise ValueError("QUBO comparison matrix is not frontier")
    qubo_normalized = normalized_hash(qubo_path)

    records: list[dict[str, Any]] = []
    for path in tie_paths:
        matrix = read_matrix(path)
        determinant = bareiss_determinant(matrix)
        if abs(determinant) != FRONTIER:
            raise RuntimeError(f"preserved tie is not frontier: {path}")
        verification = subprocess.run(
            ["./arena", "verify", str(path)],
            cwd=root,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        verification_path = path.with_suffix("").with_suffix(
            ".arena-verify.txt"
        )
        atomic_write(verification_path, verification.stdout.encode("utf-8"))
        if verification.returncode != 0:
            raise RuntimeError(f"arena rejected {path}")
        fields = parse_arena(verification.stdout)
        raw_hash = sha256_bytes(canonical_matrix_bytes(matrix))
        local_normalized = normalized_hash(path)
        if (
            fields["matrix_sha256"] != raw_hash
            or fields["normalized_sha256"] != local_normalized
        ):
            raise RuntimeError("local hashes disagree with arena verifier")
        matches_known24 = known_by_normalized.get(local_normalized, [])
        matches_qubo = local_normalized == qubo_normalized
        records.append(
            {
                "arena_verified": True,
                "artifact": relative(path, root),
                "matches_known24_by_sign_normalized_identity":
                    matches_known24,
                "matches_qubo_by_sign_normalized_identity": matches_qubo,
                "normalized_new_vs_known24_plus_qubo":
                    not matches_known24 and not matches_qubo,
                "normalized_sha256": local_normalized,
                "raw_sha256": raw_hash,
                "receipt_sha256": fields["receipt_sha256"],
                "signed_determinant": str(determinant),
                "verification_artifact": relative(
                    verification_path, root
                ),
            }
        )

    output = (
        (root / arguments.output).resolve()
        if arguments.output is not None
        else lnps / "tie-audit" / "normalized-comparison.json"
    )
    report: dict[str, Any] = {
        "all_arena_verified": all(
            record["arena_verified"] for record in records
        ),
        "complete": True,
        "frontier": str(FRONTIER),
        "known24_directory": relative(known24_directory, root),
        "known24_file_count": len(known24_paths),
        "method":
            "exact-Bareiss-plus-arena-and-sign-normalized-hash-comparison",
        "normalized_new_count_vs_known24_plus_qubo": sum(
            bool(record["normalized_new_vs_known24_plus_qubo"])
            for record in records
        ),
        "preserved_tie_count": len(records),
        "qubo": {
            "normalized_sha256": qubo_normalized,
            "path": relative(qubo_path, root),
        },
        "schema_version": 1,
        "ties": records,
    }
    if arguments.h_audit_python is not None:
        audit_python = (
            arguments.h_audit_python
            if arguments.h_audit_python.is_absolute()
            else (Path.cwd() / arguments.h_audit_python).absolute()
        )
        if not audit_python.is_file():
            raise ValueError(f"missing H-audit Python: {audit_python}")
        command = [
            str(audit_python),
            "research/h_equivalence_audit.py",
            "--local-corpus",
            str(qubo_path),
            *(str(path) for path in tie_paths),
        ]
        audit = subprocess.run(
            command,
            cwd=root,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if audit.returncode != 0:
            raise RuntimeError(
                "pinned H/HT audit failed: " + audit.stderr.strip()
            )
        h_report = json.loads(audit.stdout)
        h_output = lnps / "tie-audit" / "h-equivalence-audit.json"
        atomic_json(h_output, h_report)
        atomic_json(
            lnps / "tie-audit" / "h-equivalence-command.json",
            {
                "command": command,
                "pynauty_version": h_report["pynauty_version"],
                "returncode": audit.returncode,
            },
        )
        tie_raw_hashes = {record["raw_sha256"] for record in records}

        def certificate_map(key: str) -> dict[str, str]:
            result: dict[str, str] = {}
            for group in h_report[key]:
                certificate = group["certificate_sha256"]
                for member in group["members"]:
                    result[member["matrix_sha256"]] = certificate
            return result

        h_by_raw = certificate_map("h_classes")
        ht_by_raw = certificate_map("ht_classes")
        baseline_h = {
            group["certificate_sha256"]
            for group in h_report["h_classes"]
            if any(
                member["matrix_sha256"] not in tie_raw_hashes
                for member in group["members"]
            )
        }
        baseline_ht = {
            group["certificate_sha256"]
            for group in h_report["ht_classes"]
            if any(
                member["matrix_sha256"] not in tie_raw_hashes
                for member in group["members"]
            )
        }
        for record in records:
            raw_hash = record["raw_sha256"]
            record["h_certificate_sha256"] = h_by_raw[raw_hash]
            record["ht_certificate_sha256"] = ht_by_raw[raw_hash]
            record["adds_local_h_class"] = (
                h_by_raw[raw_hash] not in baseline_h
            )
            record["adds_local_ht_class"] = (
                ht_by_raw[raw_hash] not in baseline_ht
            )
        report["h_ht_audit"] = {
            "artifact": relative(h_output, root),
            "baseline_h_class_count": len(baseline_h),
            "baseline_ht_class_count": len(baseline_ht),
            "candidate_adds_h_classes": len(
                {
                    record["h_certificate_sha256"]
                    for record in records
                    if record["adds_local_h_class"]
                }
            ),
            "candidate_adds_ht_classes": len(
                {
                    record["ht_certificate_sha256"]
                    for record in records
                    if record["adds_local_ht_class"]
                }
            ),
            "combined_gram_class_count": h_report["gram_class_count"],
            "combined_h_class_count": h_report["h_class_count"],
            "combined_ht_class_count": h_report["ht_class_count"],
            "combined_unique_matrix_count": h_report[
                "unique_matrix_count"
            ],
            "pynauty_version": h_report["pynauty_version"],
            "scope":
                "novel only relative to frozen local corpus plus QUBO tie",
        }
    atomic_json(output, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lnps-dir", type=Path, default=DEFAULT_LNPS)
    parser.add_argument(
        "--known24-dir", type=Path, default=DEFAULT_KNOWN24
    )
    parser.add_argument("--qubo", type=Path, default=DEFAULT_QUBO)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--h-audit-python",
        type=Path,
        help=(
            "run the full local-corpus H/HT audit with this Python "
            "(must contain pinned pynauty==2.8.8.1)"
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    try:
        raise SystemExit(run(parse_arguments()))
    except (OSError, ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
