"""Load and validate the immutable JSON challenge contract."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .errors import ContractError, MatrixFormatError
from .json_tools import StrictJsonError, loads_strict_json

SUPPORTED_SCHEMA_VERSION = 1
SUPPORTED_ALGORITHM = "fraction-free-bareiss-v1"
EXPECTED_CONTRACT_DATA = {
    "schema_version": 1,
    "challenge_id": "maxdet-23-v1",
    "title": "MaxDet Arena — Order 23",
    "order": 23,
    "allowed_entries": [-1, 1],
    "objective": {
        "direction": "maximize",
        "quantity": "absolute_determinant",
        "ranking_type": "exact_integer",
    },
    "matrix_format": {
        "encoding": "utf-8",
        "rows": 23,
        "columns": 23,
        "tokens": ["-1", "1"],
        "separator": "ascii_whitespace",
        "max_bytes": 8192,
    },
    "verification": {
        "algorithm": "fraction-free-bareiss-v1",
        "independent_checks": [
            "gram-determinant-identity",
            "modular-determinants",
            "hadamard-bound",
            "power-of-two-divisibility",
        ],
        "modular_primes": [998244353, 1000000007, 1000000009],
        "sign_normalization": "first-row-and-column-positive-v1",
    },
    "claims": {
        "verified_means": (
            "The submitted matrix is valid and has the reported exact determinant."
        ),
        "verified_does_not_mean": (
            "The matrix is globally optimal or a world record."
        ),
        "world_record_policy": (
            "Any world-record claim requires a separate literature review and "
            "independent expert confirmation."
        ),
    },
}


@dataclass(frozen=True)
class Contract:
    path: Path
    raw_bytes: bytes
    data: dict[str, Any]

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.raw_bytes).hexdigest()

    @property
    def challenge_id(self) -> str:
        return str(self.data["challenge_id"])

    @property
    def order(self) -> int:
        return int(self.data["order"])

    @property
    def max_matrix_bytes(self) -> int:
        return int(self.data["matrix_format"]["max_bytes"])

    @property
    def modular_primes(self) -> tuple[int, ...]:
        return tuple(int(value) for value in self.data["verification"]["modular_primes"])


def load_contract(path: Path) -> Contract:
    try:
        raw_bytes = path.read_bytes()
    except OSError as exc:
        raise ContractError(f"cannot read challenge contract: {exc}") from exc

    try:
        data = loads_strict_json(raw_bytes)
    except StrictJsonError as exc:
        raise ContractError(f"challenge contract is not strict UTF-8 JSON: {exc}") from exc

    if not isinstance(data, dict):
        raise ContractError("challenge contract root must be an object")
    if data != EXPECTED_CONTRACT_DATA:
        raise ContractError("v1 challenge contract fields are immutable")

    return Contract(path=path, raw_bytes=raw_bytes, data=data)


def parse_matrix_bytes(raw_bytes: bytes, contract: Contract) -> list[list[int]]:
    if len(raw_bytes) > contract.max_matrix_bytes:
        raise MatrixFormatError(
            f"matrix exceeds {contract.max_matrix_bytes}-byte limit"
        )
    if b"\x00" in raw_bytes:
        raise MatrixFormatError("matrix contains a NUL byte")
    allowed_bytes = set(b"-1 \t\r\n\v\f")
    if any(value not in allowed_bytes for value in raw_bytes):
        raise MatrixFormatError(
            "matrix may contain only -1, 1, and ASCII whitespace"
        )
    text = raw_bytes.decode("ascii")
    lines = text.strip(" \t\r\n\v\f").splitlines()
    if len(lines) != contract.order:
        raise MatrixFormatError(
            f"expected {contract.order} rows, found {len(lines)}"
        )

    matrix: list[list[int]] = []
    for row_number, line in enumerate(lines, start=1):
        tokens = line.split()
        if len(tokens) != contract.order:
            raise MatrixFormatError(
                f"row {row_number}: expected {contract.order} entries, "
                f"found {len(tokens)}"
            )
        if any(token not in {"-1", "1"} for token in tokens):
            raise MatrixFormatError(
                f"row {row_number}: entries must be literal -1 or 1"
            )
        matrix.append([int(token) for token in tokens])

    return matrix


def load_matrix(path: Path, contract: Contract) -> tuple[bytes, list[list[int]]]:
    try:
        raw_bytes = path.read_bytes()
    except OSError as exc:
        raise MatrixFormatError(f"cannot read matrix: {exc}") from exc
    return raw_bytes, parse_matrix_bytes(raw_bytes, contract)
