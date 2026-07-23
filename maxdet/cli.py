"""Command-line interface for the trusted local verifier."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import stat
import sys
import tempfile
import unittest
from pathlib import Path

from .contract import load_contract
from .errors import ArenaError
from .frontier import effective_frontier, trusted_artifacts
from .receipt import canonical_json_bytes, verify_matrix, verify_matrix_bytes
from .submission import (
    HANDLE_PATTERN,
    MAX_FILE_BYTES,
    SUBMISSION_ID_PATTERN,
    parse_metadata_bytes,
    verify_submission,
)


def repository_root() -> Path:
    return Path(__file__).resolve().parent.parent


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="arena",
        description="Exact verifier and research harness for MaxDet Arena.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    verify = subparsers.add_parser("verify", help="verify a candidate matrix")
    verify.add_argument(
        "matrix",
        nargs="?",
        type=Path,
        default=Path("candidate/matrix.txt"),
    )
    verify.add_argument(
        "--contract",
        type=Path,
        default=None,
        help="challenge contract (defaults to repository challenge.json)",
    )
    verify.add_argument("--json", type=Path, default=None, dest="json_path")
    verify.add_argument(
        "--quiet",
        action="store_true",
        help="print only the exact absolute determinant",
    )

    tests = subparsers.add_parser("test", help="run the dependency-free test suite")
    tests.add_argument("--verbose", action="store_true")

    check = subparsers.add_parser(
        "check-submission",
        help="strictly verify a complete submission bundle",
    )
    check.add_argument("directory", type=Path)

    prepare = subparsers.add_parser(
        "prepare",
        help="prepare candidate files as a submission bundle",
    )
    prepare.add_argument("submission_id")
    prepare.add_argument("--handle", required=True)
    prepare.add_argument("--method", required=True)
    prepare.add_argument("--parent", default=None)

    return parser


def command_verify(args: argparse.Namespace) -> int:
    root = repository_root()
    contract_path = args.contract or root / "challenge.json"
    contract = load_contract(contract_path)
    verified = verify_matrix(args.matrix, contract)
    receipt_bytes = canonical_json_bytes(verified.receipt)

    if args.json_path:
        args.json_path.parent.mkdir(parents=True, exist_ok=True)
        args.json_path.write_bytes(receipt_bytes)

    if args.quiet:
        print(verified.abs_determinant)
        return 0

    score = verified.receipt["score"]
    matrix = verified.receipt["matrix"]
    print("VERIFIED")
    print(f"challenge: {verified.receipt['challenge_id']}")
    print(f"order: {matrix['order']}")
    print(f"determinant: {score['determinant']}")
    print(f"score |det|: {score['absolute_determinant']}")
    print(
        "Hadamard ratio: "
        f"{score['hadamard_ratio_parts_per_billion'] / 10_000_000:.7f}%"
    )
    barba_ratio_percent = (
        math.sqrt(
            score["barba_ratio_squared_parts_per_billion"] / 1_000_000_000
        )
        * 100
    )
    print(f"Barba ratio: {barba_ratio_percent:.7f}%")
    print(f"matrix sha256: {matrix['raw_sha256']}")
    print(f"normalized sha256: {matrix['sign_normalized_sha256']}")
    print(f"receipt sha256: {verified.receipt['receipt_sha256']}")
    frontier = effective_frontier(root, contract)
    if verified.abs_determinant > frontier.absolute_determinant:
        print(
            "frontier: strict improvement by "
            f"{verified.abs_determinant - frontier.absolute_determinant}"
        )
    else:
        print(
            "frontier gap: "
            f"{frontier.absolute_determinant - verified.abs_determinant}"
        )
    if args.json_path:
        print(f"receipt: {args.json_path}")
    print("claim: verified arena score; not a world-record or optimality claim")
    return 0


def command_test(args: argparse.Namespace) -> int:
    root = repository_root()
    suite = unittest.defaultTestLoader.discover(str(root / "tests"))
    result = unittest.TextTestRunner(verbosity=2 if args.verbose else 1).run(suite)
    return 0 if result.wasSuccessful() else 1


def command_check_submission(args: argparse.Namespace) -> int:
    root = repository_root()
    contract = load_contract(root / "challenge.json")
    result = verify_submission(args.directory, contract)
    submissions_root = (root / "submissions").absolute()
    try:
        relative = args.directory.absolute().relative_to(submissions_root)
    except ValueError:
        relative = None
    if relative is not None:
        if len(relative.parts) != 2:
            raise ValueError("submission path must be submissions/HANDLE/RESULT_ID")
        if relative.parts[0] != result["handle"]:
            raise ValueError("submission path handle does not match metadata")
        if relative.parts[1] != result["submission_id"]:
            raise ValueError("submission path id does not match metadata")
    frontier = effective_frontier(
        root,
        contract,
        exclude_directory=args.directory,
    )
    score = int(result["absolute_determinant"])
    if score <= frontier.absolute_determinant:
        raise ValueError(
            f"score {score} does not beat current frontier "
            f"{frontier.absolute_determinant}"
        )
    parent = result["parent_receipt_sha256"]
    trusted_receipts = {
        artifact.receipt_sha256 for artifact in trusted_artifacts(root, contract)
    }
    if parent == result["receipt_sha256"]:
        raise ValueError("submission cannot reference itself as parent")
    if parent is not None and parent not in trusted_receipts:
        raise ValueError("parent receipt is not present in trusted artifacts")
    print("SUBMISSION VERIFIED")
    print(f"id: {result['submission_id']}")
    print(f"handle: {result['handle']}")
    print(f"score |det|: {result['absolute_determinant']}")
    print(f"receipt sha256: {result['receipt_sha256']}")
    return 0


def _lstat_exists(path: Path) -> bool:
    try:
        path.lstat()
    except FileNotFoundError:
        return False
    return True


def _ensure_real_directory(path: Path) -> None:
    try:
        metadata = path.lstat()
    except FileNotFoundError:
        path.mkdir()
        metadata = path.lstat()
    if not stat.S_ISDIR(metadata.st_mode):
        raise ValueError(f"path must be a real directory: {path}")


def command_prepare(args: argparse.Namespace) -> int:
    root = repository_root()
    if not HANDLE_PATTERN.fullmatch(args.handle):
        raise ValueError("--handle must be a safe lowercase handle")
    if not SUBMISSION_ID_PATTERN.fullmatch(args.submission_id):
        raise ValueError("submission_id must be a safe lowercase slug")
    if args.parent is not None and not re.fullmatch(r"[0-9a-f]{64}", args.parent):
        raise ValueError("--parent must be a lowercase SHA-256")
    contract = load_contract(root / "challenge.json")
    candidate = root / "candidate" / "matrix.txt"
    try:
        candidate_bytes = candidate.read_bytes()
    except OSError as exc:
        raise ValueError(f"cannot read candidate matrix: {exc}") from exc
    verified = verify_matrix_bytes(candidate_bytes, contract)

    frontier = effective_frontier(root, contract)
    if verified.abs_determinant <= frontier.absolute_determinant:
        raise ValueError(
            f"candidate score {verified.abs_determinant} does not beat current "
            f"frontier {frontier.absolute_determinant}"
        )
    if args.parent is not None:
        trusted_receipts = {
            artifact.receipt_sha256
            for artifact in trusted_artifacts(root, contract)
        }
        if args.parent not in trusted_receipts:
            raise ValueError("--parent is not present in trusted artifacts")

    metadata = {
        "schema_version": 1,
        "submission_id": args.submission_id,
        "handle": args.handle,
        "method": args.method,
        "parent_receipt_sha256": args.parent,
        "artifact_license": "CC0-1.0",
    }
    metadata_bytes = (
        json.dumps(metadata, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    )
    if len(metadata_bytes) > MAX_FILE_BYTES["metadata.json"]:
        raise ValueError("generated metadata exceeds size limit")
    parse_metadata_bytes(metadata_bytes)

    submissions_root = root / "submissions"
    directory = submissions_root / args.handle / args.submission_id
    if _lstat_exists(directory):
        raise ValueError(f"submission path already exists: {directory}")

    staging_root = root / "runs"
    _ensure_real_directory(staging_root)
    staging = Path(tempfile.mkdtemp(prefix="prepare-", dir=staging_root))
    created_handle = False
    try:
        (staging / "matrix.txt").write_bytes(candidate_bytes)
        (staging / "metadata.json").write_bytes(metadata_bytes)
        (staging / "receipt.json").write_bytes(canonical_json_bytes(verified.receipt))
        result = verify_submission(staging, contract)

        _ensure_real_directory(submissions_root)
        handle_directory = submissions_root / args.handle
        if not _lstat_exists(handle_directory):
            handle_directory.mkdir()
            created_handle = True
        _ensure_real_directory(handle_directory)
        if _lstat_exists(directory):
            raise ValueError(f"submission path already exists: {directory}")
        staging.rename(directory)
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        if created_handle:
            try:
                (submissions_root / args.handle).rmdir()
            except OSError:
                pass
        raise

    print(f"prepared: {directory.relative_to(root)}")
    print(f"score |det|: {result['absolute_determinant']}")
    print("next: commit this directory and open a pull request")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "verify":
            return command_verify(args)
        if args.command == "test":
            return command_test(args)
        if args.command == "check-submission":
            return command_check_submission(args)
        if args.command == "prepare":
            return command_prepare(args)
        parser.error(f"unknown command: {args.command}")
    except ArenaError as exc:
        print(f"REJECTED: {exc}", file=sys.stderr)
        return 2
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
