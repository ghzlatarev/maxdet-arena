"""Load the declared frontier floor and derive the effective trusted frontier."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

from .contract import Contract
from .errors import SubmissionError
from .json_tools import StrictJsonError, loads_strict_json
from .submission import discover_submission_directories, verify_submission

SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
DECIMAL_PATTERN = re.compile(r"^(?:0|[1-9][0-9]*)$")
DATE_PATTERN = re.compile(r"^[0-9]{4}-[0-9]{2}-[0-9]{2}$")


@dataclass(frozen=True)
class Frontier:
    absolute_determinant: int
    receipt_sha256: str
    source: str
    label: str


def _validate_artifact(data: Any, label: str) -> Frontier:
    if not isinstance(data, dict):
        raise SubmissionError(f"{label} must be an object")
    required = {"absolute_determinant", "label", "receipt_sha256", "source"}
    if set(data) != required:
        raise SubmissionError(f"{label} fields do not match schema")
    score = data["absolute_determinant"]
    if not isinstance(score, str) or not DECIMAL_PATTERN.fullmatch(score):
        raise SubmissionError(f"{label}.absolute_determinant must be a decimal string")
    receipt = data["receipt_sha256"]
    if not isinstance(receipt, str) or not SHA256_PATTERN.fullmatch(receipt):
        raise SubmissionError(f"{label}.receipt_sha256 must be a lowercase SHA-256")
    source = data["source"]
    if not isinstance(source, str) or not source:
        raise SubmissionError(f"{label}.source must be a non-empty relative path")
    source_path = PurePosixPath(source)
    if source_path.is_absolute() or ".." in source_path.parts:
        raise SubmissionError(f"{label}.source must be a safe relative path")
    artifact_label = data["label"]
    if not isinstance(artifact_label, str) or not (1 <= len(artifact_label) <= 200):
        raise SubmissionError(f"{label}.label must contain 1 to 200 characters")
    return Frontier(
        absolute_determinant=int(score),
        receipt_sha256=receipt,
        source=source,
        label=artifact_label,
    )


def load_frontier_data(root: Path, contract: Contract) -> tuple[dict[str, Any], Frontier]:
    path = root / "data" / "frontier.json"
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise SubmissionError(f"cannot read frontier data: {exc}") from exc
    try:
        data = loads_strict_json(raw)
    except StrictJsonError as exc:
        raise SubmissionError(f"frontier data is not strict UTF-8 JSON: {exc}") from exc
    if not isinstance(data, dict):
        raise SubmissionError("frontier data root must be an object")
    required = {
        "schema_version",
        "challenge_id",
        "status",
        "target_to_beat",
        "arena_best",
        "claim_boundary",
        "updated",
    }
    if set(data) != required:
        raise SubmissionError("frontier data fields do not match schema")
    if data["schema_version"] != 1:
        raise SubmissionError("frontier schema_version must be 1")
    if data["challenge_id"] != contract.challenge_id:
        raise SubmissionError("frontier challenge_id does not match contract")
    if data["status"] not in {"private-dogfooding", "open"}:
        raise SubmissionError("frontier status must be private-dogfooding or open")
    if (
        not isinstance(data["claim_boundary"], str)
        or not (1 <= len(data["claim_boundary"]) <= 500)
    ):
        raise SubmissionError("frontier claim_boundary must contain 1 to 500 characters")
    if not isinstance(data["updated"], str) or not DATE_PATTERN.fullmatch(
        data["updated"]
    ):
        raise SubmissionError("frontier updated must be YYYY-MM-DD")
    floor = _validate_artifact(data["target_to_beat"], "target_to_beat")
    _validate_artifact(data["arena_best"], "arena_best")
    return data, floor


def effective_frontier(
    root: Path,
    contract: Contract,
    *,
    exclude_directory: Path | None = None,
) -> Frontier:
    """Return max(declared floor, every exactly verified accepted submission)."""

    _, best = load_frontier_data(root, contract)
    for directory in discover_submission_directories(root / "submissions"):
        if exclude_directory is not None and directory.absolute() == exclude_directory.absolute():
            continue
        result = verify_submission(directory, contract)
        score = int(result["absolute_determinant"])
        if score > best.absolute_determinant:
            best = Frontier(
                absolute_determinant=score,
                receipt_sha256=result["receipt_sha256"],
                source=directory.relative_to(root).as_posix(),
                label=f"{result['handle']}/{result['submission_id']}",
            )
    return best
