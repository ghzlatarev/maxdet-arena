#!/usr/bin/env python3
"""Verify one untrusted pull-request submission with trusted base code."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path, PurePosixPath

TRUSTED_ROOT = Path(__file__).resolve().parent.parent
if str(TRUSTED_ROOT) not in sys.path:
    sys.path.insert(0, str(TRUSTED_ROOT))

from maxdet.contract import load_contract  # noqa: E402
from maxdet.errors import SubmissionError  # noqa: E402
from maxdet.frontier import effective_frontier, trusted_artifacts  # noqa: E402
from maxdet.submission import ALLOWED_FILES, verify_submission  # noqa: E402


def tracked_entries(root: Path) -> dict[str, tuple[str, str]]:
    process = subprocess.run(
        ["git", "-C", str(root), "ls-files", "-s", "-z"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    entries: dict[str, tuple[str, str]] = {}
    for record in process.stdout.split(b"\0"):
        if not record:
            continue
        header, raw_path = record.split(b"\t", 1)
        mode, object_id, stage = header.decode("ascii").split()
        if stage != "0":
            raise SubmissionError("unmerged index entry in pull request")
        path = raw_path.decode("utf-8")
        entries[path] = (mode, object_id)
    return entries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--untrusted-root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    untrusted_root = args.untrusted_root.resolve()
    base = tracked_entries(TRUSTED_ROOT)
    head = tracked_entries(untrusted_root)
    changed = {
        path
        for path in base.keys() | head.keys()
        if base.get(path) != head.get(path)
    }
    if not changed:
        raise SubmissionError("pull request contains no tracked changes")

    bundle_paths: set[PurePosixPath] = set()
    for raw_path in changed:
        path = PurePosixPath(raw_path)
        if (
            len(path.parts) != 4
            or path.parts[0] != "submissions"
            or path.name not in ALLOWED_FILES
        ):
            raise SubmissionError(
                f"pull request may change only one submission bundle: {raw_path}"
            )
        if head.get(raw_path, (None,))[0] != "100644":
            raise SubmissionError(f"submission file must be a regular 0644 blob: {raw_path}")
        bundle_paths.add(path.parent)

    if len(bundle_paths) != 1:
        raise SubmissionError("pull request must contain exactly one submission bundle")
    bundle = next(iter(bundle_paths))
    prefix = bundle.as_posix() + "/"
    if any(path.startswith(prefix) for path in base):
        raise SubmissionError("existing submissions are immutable")

    directory = untrusted_root.joinpath(*bundle.parts)
    contract = load_contract(TRUSTED_ROOT / "challenge.json")
    result = verify_submission(directory, contract)
    if bundle.parts[1] != result["handle"]:
        raise SubmissionError("submission path handle does not match metadata")
    if bundle.parts[2] != result["submission_id"]:
        raise SubmissionError("submission path id does not match metadata")

    target = effective_frontier(TRUSTED_ROOT, contract).absolute_determinant
    score = int(result["absolute_determinant"])
    if score <= target:
        raise SubmissionError(
            f"score {score} does not beat current frontier {target}"
        )
    artifacts = trusted_artifacts(TRUSTED_ROOT, contract)
    normalized_hashes = {
        artifact.sign_normalized_sha256 for artifact in artifacts
    }
    if result["sign_normalized_sha256"] in normalized_hashes:
        raise SubmissionError("sign-normalized matrix duplicates a trusted artifact")
    parent = result["parent_receipt_sha256"]
    receipt_hashes = {artifact.receipt_sha256 for artifact in artifacts}
    if parent is not None and parent not in receipt_hashes:
        raise SubmissionError("parent receipt is not present on the trusted base branch")

    print("PULL REQUEST SUBMISSION VERIFIED")
    print(f"path: {bundle}")
    print(f"handle: {result['handle']}")
    print(f"score |det|: {score}")
    print(f"improvement: {score - target}")
    print(f"receipt sha256: {result['receipt_sha256']}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SubmissionError as exc:
        print(f"REJECTED: {exc}", file=sys.stderr)
        raise SystemExit(2)
