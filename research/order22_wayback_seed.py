#!/usr/bin/env python3
"""Recover the two archived Indiana MaxDet order-22 factors exactly."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path


ORDER = 22
EXPECTED_ABS_DETERMINANT = 409_600_000_000_000
MAX_SOURCE_BYTES = 1_000_000


@dataclass(frozen=True)
class Capture:
    source_sha256: str
    matrix_sha256: str


CAPTURES = {
    "20030308032721": Capture(
        source_sha256=(
            "457fd2a4b76c1fabbb64b238a2fcd39e3824fb2305e3ee177bfa9d5e5e131fbd"
        ),
        matrix_sha256=(
            "d919dac4566093f0676ddc0a304d8c9526a85f8ab3333af1a0e26bc41eefaec6"
        ),
    ),
    "20030713162905": Capture(
        source_sha256=(
            "6ef5551084a64b460ed9a4fcd049b7fa2a2dc17dd2b616f8718fbf4d62cd3357"
        ),
        matrix_sha256=(
            "c7684c22aa4dc37fa18617e186c292268d2cb5f70b2a04268dd1c4b008c04f93"
        ),
    ),
}

FACTOR_PATTERN = re.compile(
    rb"(?:^|\n)R:\s*<pre>\s*(.*?)\s*</pre>",
    flags=re.IGNORECASE | re.DOTALL,
)


def read_source(path: Path) -> bytes:
    if not path.is_file():
        raise ValueError("source path is not a regular file")
    size = path.stat().st_size
    if size > MAX_SOURCE_BYTES:
        raise ValueError(f"source exceeds {MAX_SOURCE_BYTES} bytes")
    return path.read_bytes()


def parse_factor(source: bytes) -> list[list[int]]:
    matches = FACTOR_PATTERN.findall(source)
    if len(matches) != 1:
        raise ValueError(f"expected one R factor block, found {len(matches)}")

    encoded_rows = [row.strip() for row in matches[0].splitlines() if row.strip()]
    if len(encoded_rows) != ORDER:
        raise ValueError(f"expected {ORDER} factor rows, found {len(encoded_rows)}")

    matrix: list[list[int]] = []
    for index, encoded in enumerate(encoded_rows, start=1):
        if len(encoded) != ORDER or any(value not in b"+-" for value in encoded):
            raise ValueError(f"invalid factor row {index}")
        matrix.append([1 if value == ord("+") else -1 for value in encoded])
    return matrix


def bareiss_determinant(source: list[list[int]]) -> int:
    work = [row[:] for row in source]
    previous = 1
    sign = 1
    for column in range(len(work) - 1):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[pivot], work[column] = work[column], work[pivot]
            sign = -sign
        pivot_value = work[column][column]
        for row in range(column + 1, len(work)):
            for inner in range(column + 1, len(work)):
                numerator = (
                    work[row][inner] * pivot_value
                    - work[row][column] * work[column][inner]
                )
                if column and numerator % previous:
                    raise ArithmeticError("non-exact Bareiss division")
                work[row][inner] = numerator if not column else numerator // previous
            work[row][column] = 0
        previous = pivot_value
    return sign * work[-1][-1]


def payload(matrix: list[list[int]]) -> bytes:
    return "".join(" ".join(map(str, row)) + "\n" for row in matrix).encode("ascii")


def atomic_write(path: Path, contents: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not path.is_file():
        raise ValueError("output path is not a regular file")
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(contents)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capture", choices=sorted(CAPTURES), required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    expected = CAPTURES[arguments.capture]
    source = read_source(arguments.source.expanduser().absolute())
    source_digest = hashlib.sha256(source).hexdigest()
    if source_digest != expected.source_sha256:
        raise ValueError(f"unexpected source SHA-256: {source_digest}")

    matrix = parse_factor(source)
    determinant = bareiss_determinant(matrix)
    if abs(determinant) != EXPECTED_ABS_DETERMINANT:
        raise ArithmeticError(f"unexpected determinant: {determinant}")

    contents = payload(matrix)
    matrix_digest = hashlib.sha256(contents).hexdigest()
    if matrix_digest != expected.matrix_sha256:
        raise ValueError(f"unexpected matrix SHA-256: {matrix_digest}")

    atomic_write(arguments.output.expanduser().absolute(), contents)
    print(
        f"recovered capture={arguments.capture} order={ORDER} "
        f"det={determinant} sha256={matrix_digest}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
