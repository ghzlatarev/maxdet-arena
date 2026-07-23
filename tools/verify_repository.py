#!/usr/bin/env python3
"""Reproduce every trusted score and cross-check frontier metadata."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from maxdet.contract import load_contract  # noqa: E402
from maxdet.frontier import (  # noqa: E402
    effective_frontier,
    load_frontier_data,
    trusted_artifacts,
)
from maxdet.json_tools import loads_strict_json  # noqa: E402
from maxdet.receipt import canonical_json_bytes, verify_matrix  # noqa: E402
from maxdet.submission import (  # noqa: E402
    discover_submission_directories,
    verify_submission,
)


def check_receipt(matrix_path: Path, receipt_path: Path) -> dict:
    contract = load_contract(ROOT / "challenge.json")
    verified = verify_matrix(matrix_path, contract)
    expected = canonical_json_bytes(verified.receipt)
    observed = receipt_path.read_bytes()
    if observed != expected:
        raise RuntimeError(f"stale trusted receipt: {receipt_path.relative_to(ROOT)}")
    return verified.receipt


def check_reference_source(reference: dict) -> None:
    source_path = ROOT / "references/orrick-et-al-2003/source.json"
    source = loads_strict_json(source_path.read_bytes())
    if not isinstance(source, dict):
        raise RuntimeError("reference source metadata must be an object")
    if (
        source.get("reported_absolute_determinant")
        != reference["score"]["absolute_determinant"]
    ):
        raise RuntimeError("reference source determinant does not match receipt")

    archive = source.get("source_archive")
    if not isinstance(archive, dict):
        raise RuntimeError("reference source archive metadata is missing")
    if archive.get("byte_match_after_conversion") is not True:
        raise RuntimeError("reference source does not record an exact matrix match")
    if archive.get("matrix_raw_sha256") != reference["matrix"]["raw_sha256"]:
        raise RuntimeError("reference source matrix hash does not match receipt")
    archive_hash = archive.get("archive_sha256")
    if (
        not isinstance(archive_hash, str)
        or len(archive_hash) != 64
        or any(character not in "0123456789abcdef" for character in archive_hash)
    ):
        raise RuntimeError("reference source archive hash is not lowercase SHA-256")

    upper_bound = source.get("reported_order_23_upper_bound")
    if not isinstance(upper_bound, dict):
        raise RuntimeError("reference source upper-bound metadata is missing")
    if (
        upper_bound.get("determinant_squared")
        != reference["score"]["ehlich_bound_squared"]
    ):
        raise RuntimeError("reference source upper bound does not match contract")

    survey = source.get("survey_cross_check")
    if not isinstance(survey, dict) or survey.get("doi") != "10.37236/10367":
        raise RuntimeError("reference survey cross-check is missing or altered")


def main() -> int:
    genesis = check_receipt(
        ROOT / "records/genesis/matrix.txt",
        ROOT / "records/genesis/receipt.json",
    )
    reference = check_receipt(
        ROOT / "references/orrick-et-al-2003/matrix.txt",
        ROOT / "references/orrick-et-al-2003/receipt.json",
    )
    check_reference_source(reference)

    contract = load_contract(ROOT / "challenge.json")
    frontier, floor = load_frontier_data(ROOT, contract)
    target = frontier["target_to_beat"]
    if target["absolute_determinant"] != reference["score"]["absolute_determinant"]:
        raise RuntimeError("frontier target does not match published reference receipt")
    if target["receipt_sha256"] != reference["receipt_sha256"]:
        raise RuntimeError("frontier target receipt hash is stale")
    if target["source"] != "references/orrick-et-al-2003":
        raise RuntimeError("frontier target source is not the published reference")

    arena_best = frontier["arena_best"]
    arena_source = ROOT / arena_best["source"]
    arena_receipt = check_receipt(
        arena_source / "matrix.txt",
        arena_source / "receipt.json",
    )
    if (
        arena_best["absolute_determinant"]
        != arena_receipt["score"]["absolute_determinant"]
    ):
        raise RuntimeError("arena_best score does not match its source matrix")
    if arena_best["receipt_sha256"] != arena_receipt["receipt_sha256"]:
        raise RuntimeError("arena_best receipt hash is stale")
    verify_matrix(ROOT / "candidate/matrix.txt", contract)

    trusted_receipts: set[str] = set()
    normalized_hashes: set[str] = set()
    for artifact in trusted_artifacts(ROOT, contract):
        if artifact.source.startswith("submissions/"):
            continue
        if artifact.receipt_sha256 in trusted_receipts:
            raise RuntimeError(
                f"duplicate receipt in trusted artifacts: {artifact.source}"
            )
        if artifact.sign_normalized_sha256 in normalized_hashes:
            raise RuntimeError(
                f"sign-normalized duplicate in trusted artifacts: {artifact.source}"
            )
        trusted_receipts.add(artifact.receipt_sha256)
        normalized_hashes.add(artifact.sign_normalized_sha256)
    submission_results: list[dict] = []
    for directory in discover_submission_directories(ROOT / "submissions"):
        result = verify_submission(directory, contract)
        relative = directory.relative_to(ROOT)
        if relative.parts[1] != result["handle"]:
            raise RuntimeError(f"submission handle path mismatch: {relative}")
        if relative.parts[2] != result["submission_id"]:
            raise RuntimeError(f"submission id path mismatch: {relative}")
        if int(result["absolute_determinant"]) <= floor.absolute_determinant:
            raise RuntimeError(f"accepted submission does not beat floor: {relative}")
        if result["receipt_sha256"] in trusted_receipts:
            raise RuntimeError(f"duplicate receipt in trusted repository: {relative}")
        if result["sign_normalized_sha256"] in normalized_hashes:
            raise RuntimeError(f"sign-normalized duplicate in repository: {relative}")
        trusted_receipts.add(result["receipt_sha256"])
        normalized_hashes.add(result["sign_normalized_sha256"])
        submission_results.append(result)

    for result in submission_results:
        parent = result["parent_receipt_sha256"]
        if parent is not None and parent not in trusted_receipts:
            raise RuntimeError(
                f"submission {result['handle']}/{result['submission_id']} "
                "references an unknown parent"
            )
        if parent == result["receipt_sha256"]:
            raise RuntimeError(
                f"submission {result['handle']}/{result['submission_id']} "
                "references itself as parent"
            )

    parent_by_receipt = {
        result["receipt_sha256"]: result["parent_receipt_sha256"]
        for result in submission_results
    }
    visit_state: dict[str, int] = {}

    def visit(receipt_sha256: str) -> None:
        state = visit_state.get(receipt_sha256, 0)
        if state == 1:
            raise RuntimeError("submission parent graph contains a cycle")
        if state == 2:
            return
        visit_state[receipt_sha256] = 1
        parent = parent_by_receipt[receipt_sha256]
        if parent in parent_by_receipt:
            visit(parent)
        visit_state[receipt_sha256] = 2

    for receipt_sha256 in parent_by_receipt:
        visit(receipt_sha256)

    effective = effective_frontier(ROOT, contract)

    print("TRUSTED REPOSITORY VERIFIED")
    print(f"contract sha256: {contract.sha256}")
    print(f"genesis |det|: {genesis['score']['absolute_determinant']}")
    print(f"published floor |det|: {target['absolute_determinant']}")
    print(f"effective frontier |det|: {effective.absolute_determinant}")
    print(f"accepted submissions: {len(submission_results)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
