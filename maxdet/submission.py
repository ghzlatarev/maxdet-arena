"""Strict inspection of untrusted submission artifacts.

The verifier reads data files only. It never imports or executes contributor
solver code.
"""

from __future__ import annotations

import os
import re
import stat
from pathlib import Path
from typing import Any

from .contract import Contract
from .errors import SubmissionError
from .json_tools import StrictJsonError, loads_strict_json
from .receipt import canonical_json_bytes, verify_matrix_bytes

ALLOWED_FILES = {"matrix.txt", "metadata.json", "receipt.json", "notes.md"}
REQUIRED_FILES = {"matrix.txt", "metadata.json", "receipt.json"}
MAX_FILE_BYTES = {
    "matrix.txt": 8192,
    "metadata.json": 8192,
    "receipt.json": 32768,
    "notes.md": 65536,
}
HANDLE_PATTERN = re.compile(r"^[a-z0-9](?:[a-z0-9_-]{0,38}[a-z0-9])?$")
SUBMISSION_ID_PATTERN = re.compile(
    r"^[a-z0-9](?:[a-z0-9_-]{0,62}[a-z0-9])?$"
)


def _read_regular_file(path: Path, maximum: int) -> bytes:
    try:
        metadata = path.lstat()
    except OSError as exc:
        raise SubmissionError(f"cannot inspect {path.name}: {exc}") from exc
    if not stat.S_ISREG(metadata.st_mode):
        raise SubmissionError(f"{path.name} must be a regular file, not a link")
    if metadata.st_size > maximum:
        raise SubmissionError(f"{path.name} exceeds {maximum}-byte limit")
    try:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    except OSError as exc:
        raise SubmissionError(f"cannot safely open {path.name}: {exc}") from exc
    try:
        with os.fdopen(descriptor, "rb") as stream:
            data = stream.read(maximum + 1)
    except OSError as exc:
        raise SubmissionError(f"cannot read {path.name}: {exc}") from exc
    if len(data) > maximum:
        raise SubmissionError(f"{path.name} exceeds {maximum}-byte limit")
    return data


def parse_metadata_bytes(raw: bytes) -> dict[str, Any]:
    try:
        data = loads_strict_json(raw)
    except StrictJsonError as exc:
        raise SubmissionError(
            f"metadata.json is not strict UTF-8 JSON: {exc}"
        ) from exc
    if not isinstance(data, dict):
        raise SubmissionError("metadata.json root must be an object")

    required = {
        "schema_version",
        "submission_id",
        "handle",
        "method",
        "parent_receipt_sha256",
        "artifact_license",
    }
    optional = {"agent", "notes", "runtime_seconds", "seed"}
    missing = required - data.keys()
    unknown = data.keys() - required - optional
    if missing:
        raise SubmissionError(
            "metadata.json missing keys: " + ", ".join(sorted(missing))
        )
    if unknown:
        raise SubmissionError(
            "metadata.json has unknown keys: " + ", ".join(sorted(unknown))
        )
    if data["schema_version"] != 1:
        raise SubmissionError("metadata schema_version must be 1")

    submission_id = data["submission_id"]
    if not isinstance(submission_id, str) or not SUBMISSION_ID_PATTERN.fullmatch(
        submission_id
    ):
        raise SubmissionError("submission_id must be a safe lowercase slug")
    handle = data["handle"]
    if not isinstance(handle, str) or not HANDLE_PATTERN.fullmatch(handle):
        raise SubmissionError("handle must be a safe lowercase handle")
    method = data["method"]
    if (
        not isinstance(method, str)
        or not (1 <= len(method) <= 500)
        or not method.strip()
    ):
        raise SubmissionError(
            "method must contain 1 to 500 characters and not be blank"
        )
    parent = data["parent_receipt_sha256"]
    if parent is not None and (
        not isinstance(parent, str)
        or not re.fullmatch(r"[0-9a-f]{64}", parent)
    ):
        raise SubmissionError("parent_receipt_sha256 must be null or lowercase SHA-256")
    if data["artifact_license"] != "CC0-1.0":
        raise SubmissionError("artifact_license must be CC0-1.0")

    for key in ("agent", "notes"):
        if key in data and (
            not isinstance(data[key], str) or len(data[key]) > 2000
        ):
            raise SubmissionError(f"{key} must be a string of at most 2000 characters")
    if "runtime_seconds" in data and (
        not isinstance(data["runtime_seconds"], int)
        or isinstance(data["runtime_seconds"], bool)
        or data["runtime_seconds"] < 0
    ):
        raise SubmissionError("runtime_seconds must be a non-negative integer")
    if "seed" in data and (
        not isinstance(data["seed"], int) or isinstance(data["seed"], bool)
    ):
        raise SubmissionError("seed must be an integer")
    return data


def discover_submission_directories(submissions_root: Path) -> list[Path]:
    """Return every bundle in a strict handle/id directory tree."""

    try:
        root_metadata = submissions_root.lstat()
    except FileNotFoundError:
        return []
    except OSError as exc:
        raise SubmissionError(f"cannot inspect submissions directory: {exc}") from exc
    if not stat.S_ISDIR(root_metadata.st_mode):
        raise SubmissionError("submissions path must be a real directory")

    directories: list[Path] = []
    try:
        handle_paths = sorted(submissions_root.iterdir(), key=lambda path: path.name)
    except OSError as exc:
        raise SubmissionError(f"cannot list submissions directory: {exc}") from exc
    for handle_path in handle_paths:
        if not HANDLE_PATTERN.fullmatch(handle_path.name):
            raise SubmissionError(
                f"invalid entry under submissions/: {handle_path.name}"
            )
        try:
            handle_metadata = handle_path.lstat()
        except OSError as exc:
            raise SubmissionError(
                f"cannot inspect submission handle {handle_path.name}: {exc}"
            ) from exc
        if not stat.S_ISDIR(handle_metadata.st_mode):
            raise SubmissionError(
                f"submission handle must be a real directory: {handle_path.name}"
            )
        try:
            submission_paths = sorted(
                handle_path.iterdir(),
                key=lambda path: path.name,
            )
        except OSError as exc:
            raise SubmissionError(
                f"cannot list submission handle {handle_path.name}: {exc}"
            ) from exc
        if not submission_paths:
            raise SubmissionError(
                f"submission handle directory is empty: {handle_path.name}"
            )
        for submission_path in submission_paths:
            if not SUBMISSION_ID_PATTERN.fullmatch(submission_path.name):
                raise SubmissionError(
                    "invalid submission id under "
                    f"{handle_path.name}/: {submission_path.name}"
                )
            try:
                submission_metadata = submission_path.lstat()
            except OSError as exc:
                raise SubmissionError(
                    f"cannot inspect submission {submission_path}: {exc}"
                ) from exc
            if not stat.S_ISDIR(submission_metadata.st_mode):
                raise SubmissionError(
                    f"submission must be a real directory: {submission_path}"
                )
            directories.append(submission_path)
    return directories


def verify_submission(directory: Path, contract: Contract) -> dict[str, Any]:
    try:
        directory_metadata = directory.lstat()
    except OSError as exc:
        raise SubmissionError(f"cannot inspect submission directory: {exc}") from exc
    if not stat.S_ISDIR(directory_metadata.st_mode):
        raise SubmissionError("submission path must be a real directory")

    try:
        entries = list(directory.iterdir())
    except OSError as exc:
        raise SubmissionError(f"cannot list submission directory: {exc}") from exc
    names = {entry.name for entry in entries}
    unknown = names - ALLOWED_FILES
    missing = REQUIRED_FILES - names
    if unknown:
        raise SubmissionError("unexpected submission files: " + ", ".join(sorted(unknown)))
    if missing:
        raise SubmissionError("missing submission files: " + ", ".join(sorted(missing)))

    raw_files = {
        name: _read_regular_file(directory / name, MAX_FILE_BYTES[name])
        for name in names
    }
    metadata = parse_metadata_bytes(raw_files["metadata.json"])

    verified = verify_matrix_bytes(raw_files["matrix.txt"], contract)
    expected_receipt = canonical_json_bytes(verified.receipt)
    if raw_files["receipt.json"] != expected_receipt:
        raise SubmissionError(
            "receipt.json is stale or altered; regenerate it with ./arena verify"
        )

    if "notes.md" in raw_files:
        try:
            raw_files["notes.md"].decode("utf-8")
        except UnicodeDecodeError as exc:
            raise SubmissionError("notes.md must be valid UTF-8") from exc

    return {
        "submission_id": metadata["submission_id"],
        "handle": metadata["handle"],
        "method": metadata["method"],
        "parent_receipt_sha256": metadata["parent_receipt_sha256"],
        "absolute_determinant": str(verified.abs_determinant),
        "receipt_sha256": verified.receipt["receipt_sha256"],
        "sign_normalized_sha256": verified.receipt["matrix"][
            "sign_normalized_sha256"
        ],
    }
