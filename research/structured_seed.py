#!/usr/bin/env python3
"""Create a deterministic whole-line perturbation of an order-23 matrix."""

from __future__ import annotations

import argparse
import os
import random
import tempfile
from pathlib import Path

ORDER = 23
MAX_BYTES = 8192
ASCII_WHITESPACE = " \t\r\n\v\f"
ALLOWED_BYTES = set(b"-1 \t\r\n\v\f")


def parse_matrix(path: Path) -> list[list[int]]:
    try:
        if not path.is_file():
            raise ValueError("start path is not a regular file")
        raw = path.read_bytes()
    except OSError as exc:
        raise ValueError(f"cannot read start matrix: {exc}") from exc
    if len(raw) > MAX_BYTES:
        raise ValueError(f"start matrix exceeds {MAX_BYTES} bytes")
    if any(byte not in ALLOWED_BYTES for byte in raw):
        raise ValueError("start matrix contains invalid bytes")

    lines = raw.decode("ascii").strip(ASCII_WHITESPACE).splitlines()
    if len(lines) != ORDER:
        raise ValueError(f"expected {ORDER} rows, found {len(lines)}")
    matrix: list[list[int]] = []
    for row_number, line in enumerate(lines, start=1):
        tokens = line.split()
        if len(tokens) != ORDER or any(token not in {"-1", "1"} for token in tokens):
            raise ValueError(f"row {row_number} is not a valid {ORDER}-entry sign row")
        matrix.append([int(token) for token in tokens])
    return matrix


def fresh_pattern(rng: random.Random, old: list[int]) -> list[int]:
    while True:
        pattern = [rng.choice((-1, 1)) for _ in range(ORDER)]
        if pattern != old:
            return pattern


def perturb(
    matrix: list[list[int]], seed: int, lines: int, orientation: str
) -> list[list[int]]:
    if orientation == "mixed":
        if not 2 <= lines <= 2 * ORDER:
            raise ValueError(f"mixed --lines must be between 2 and {2 * ORDER}")
        row_count = (lines + 1) // 2
        column_count = lines // 2
    else:
        if not 1 <= lines <= ORDER:
            raise ValueError(f"{orientation} --lines must be between 1 and {ORDER}")
        row_count = lines if orientation == "rows" else 0
        column_count = lines if orientation == "columns" else 0

    rng = random.Random(seed)
    result = [row[:] for row in matrix]
    for row in rng.sample(range(ORDER), row_count):
        result[row] = fresh_pattern(rng, result[row])
    for column in rng.sample(range(ORDER), column_count):
        old = [result[row][column] for row in range(ORDER)]
        pattern = fresh_pattern(rng, old)
        for row, value in enumerate(pattern):
            result[row][column] = value
    return result


def write_matrix(path: Path, matrix: list[list[int]]) -> None:
    parent = path.parent.resolve()
    if not parent.is_dir():
        raise ValueError("output parent is not a directory")
    if path.exists() and not path.is_file():
        raise ValueError("output path is not a regular file")

    payload = "".join(" ".join(map(str, row)) + "\n" for row in matrix)
    temporary: Path | None = None
    try:
        descriptor, name = tempfile.mkstemp(
            dir=parent, prefix=f".{path.name}.", suffix=".tmp"
        )
        temporary = Path(name)
        with os.fdopen(descriptor, "w", encoding="ascii", newline="\n") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--lines", required=True, type=int)
    parser.add_argument(
        "--orientation", required=True, choices=("rows", "columns", "mixed")
    )
    args = parser.parse_args()

    try:
        start = args.start.resolve(strict=True)
        output = args.output.resolve(strict=False)
        if start == output:
            raise ValueError("start and output paths must differ")
        matrix = parse_matrix(start)
        perturbed = perturb(matrix, args.seed, args.lines, args.orientation)
        write_matrix(output, perturbed)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
