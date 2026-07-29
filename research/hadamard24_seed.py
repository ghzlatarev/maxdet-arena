#!/usr/bin/env python3
"""Extract an order-23 seed from the Mendeley order-24 Hadamard corpus."""

from __future__ import annotations

import argparse
import os
import re
import tempfile
from pathlib import Path

HADAMARD_ORDER = 24
SEED_ORDER = 23
CLASS_COUNT = 60
MAX_INPUT_BYTES = 4 * 1024 * 1024
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    REPOSITORY_ROOT
    / "runs"
    / "direct-search"
    / "reference-data"
    / "mendeley-hzf94h43c5-v1"
    / "Text"
    / "Hadamard_24.txt"
)
ASSIGNMENT_PATTERN = re.compile(
    r"^[ \t]*Mat\[24(?:,([0-9]+))?\][ \t]*=",
    re.MULTILINE,
)


def read_text(path: Path) -> str:
    """Read a bounded ASCII dataset file."""

    try:
        if not path.is_file():
            raise ValueError("input path is not a regular file")
        with path.open("rb") as stream:
            raw = stream.read(MAX_INPUT_BYTES + 1)
    except OSError as exc:
        raise ValueError(f"cannot read input: {exc}") from exc

    if len(raw) > MAX_INPUT_BYTES:
        raise ValueError(f"input exceeds {MAX_INPUT_BYTES} bytes")
    try:
        return raw.decode("ascii")
    except UnicodeDecodeError as exc:
        raise ValueError("input is not ASCII") from exc


def matrix_assignments(text: str) -> dict[int, int]:
    """Map each corpus class to the offset after its assignment operator."""

    assignments: dict[int, int] = {}
    for match in ASSIGNMENT_PATTERN.finditer(text):
        matrix_class = int(match.group(1) or "1")
        if matrix_class in assignments:
            raise ValueError(f"duplicate Mat[24] class {matrix_class}")
        assignments[matrix_class] = match.end()
    return assignments


def matrix_literal(text: str, start: int, matrix_class: int) -> str:
    """Extract one balanced Mathematica list literal."""

    cursor = start
    while cursor < len(text) and text[cursor].isspace():
        cursor += 1
    if cursor >= len(text) or text[cursor] != "{":
        raise ValueError(f"class {matrix_class} has no matrix literal")

    literal_start = cursor
    depth = 0
    while cursor < len(text):
        character = text[cursor]
        if character == "{":
            depth += 1
        elif character == "}":
            depth -= 1
            if depth < 0:
                raise ValueError(f"class {matrix_class} has unbalanced braces")
            if depth == 0:
                cursor += 1
                break
        cursor += 1
    else:
        raise ValueError(f"class {matrix_class} has an unterminated matrix literal")

    literal_end = cursor
    while cursor < len(text) and text[cursor].isspace():
        cursor += 1
    if cursor >= len(text) or text[cursor] != ";":
        raise ValueError(f"class {matrix_class} matrix is not terminated by ';'")
    return text[literal_start:literal_end]


class LiteralParser:
    """Small parser for nested Mathematica lists containing only -1 and 1."""

    def __init__(self, literal: str) -> None:
        self.literal = literal
        self.cursor = 0

    def skip_whitespace(self) -> None:
        while (
            self.cursor < len(self.literal)
            and self.literal[self.cursor].isspace()
        ):
            self.cursor += 1

    def accept(self, character: str) -> bool:
        self.skip_whitespace()
        if (
            self.cursor < len(self.literal)
            and self.literal[self.cursor] == character
        ):
            self.cursor += 1
            return True
        return False

    def expect(self, character: str) -> None:
        if not self.accept(character):
            raise ValueError(
                f"expected {character!r} at matrix-literal offset {self.cursor}"
            )

    def value(self) -> int:
        self.skip_whitespace()
        if self.literal.startswith("-1", self.cursor):
            self.cursor += 2
            return -1
        if self.literal.startswith("1", self.cursor):
            self.cursor += 1
            return 1
        raise ValueError(
            f"expected -1 or 1 at matrix-literal offset {self.cursor}"
        )

    def row(self) -> list[int]:
        self.expect("{")
        row = [self.value()]
        while self.accept(","):
            if len(row) >= HADAMARD_ORDER:
                raise ValueError(f"matrix row exceeds {HADAMARD_ORDER} entries")
            row.append(self.value())
        self.expect("}")
        return row

    def matrix(self) -> list[list[int]]:
        self.expect("{")
        matrix = [self.row()]
        while self.accept(","):
            if len(matrix) >= HADAMARD_ORDER:
                raise ValueError(f"matrix exceeds {HADAMARD_ORDER} rows")
            matrix.append(self.row())
        self.expect("}")
        self.skip_whitespace()
        if self.cursor != len(self.literal):
            raise ValueError(
                f"unexpected data at matrix-literal offset {self.cursor}"
            )
        return matrix


def validate_hadamard(matrix: list[list[int]], matrix_class: int) -> None:
    """Validate shape, entries, and H H^T = 24 I exactly."""

    if len(matrix) != HADAMARD_ORDER:
        raise ValueError(
            f"class {matrix_class} has {len(matrix)} rows, "
            f"expected {HADAMARD_ORDER}"
        )
    for row_number, row in enumerate(matrix, start=1):
        if len(row) != HADAMARD_ORDER:
            raise ValueError(
                f"class {matrix_class} row {row_number} has {len(row)} entries, "
                f"expected {HADAMARD_ORDER}"
            )
        if any(value not in (-1, 1) for value in row):
            raise ValueError(
                f"class {matrix_class} row {row_number} has a non-sign entry"
            )

    for left in range(HADAMARD_ORDER):
        for right in range(left, HADAMARD_ORDER):
            inner_product = sum(
                matrix[left][column] * matrix[right][column]
                for column in range(HADAMARD_ORDER)
            )
            expected = HADAMARD_ORDER if left == right else 0
            if inner_product != expected:
                raise ValueError(
                    f"class {matrix_class} is not Hadamard: "
                    f"(H H^T)[{left + 1},{right + 1}]={inner_product}, "
                    f"expected {expected}"
                )


def load_hadamard(path: Path, matrix_class: int) -> list[list[int]]:
    text = read_text(path)
    assignments = matrix_assignments(text)
    if matrix_class not in assignments:
        available = ", ".join(map(str, sorted(assignments))) or "none"
        raise ValueError(
            f"class {matrix_class} is absent; available classes: {available}"
        )
    literal = matrix_literal(text, assignments[matrix_class], matrix_class)
    matrix = LiteralParser(literal).matrix()
    validate_hadamard(matrix, matrix_class)
    return matrix


def delete_row_and_column(
    matrix: list[list[int]], delete_row: int, delete_column: int
) -> list[list[int]]:
    """Delete one one-based row and column."""

    row_index = delete_row - 1
    column_index = delete_column - 1
    seed = [
        [
            value
            for current_column, value in enumerate(row)
            if current_column != column_index
        ]
        for current_row, row in enumerate(matrix)
        if current_row != row_index
    ]
    if len(seed) != SEED_ORDER or any(len(row) != SEED_ORDER for row in seed):
        raise ValueError("internal error: deletion did not produce a 23x23 matrix")
    return seed


def write_matrix(path: Path, matrix: list[list[int]]) -> None:
    """Atomically replace path with a strict ASCII sign matrix."""

    parent = path.parent
    if not parent.is_dir():
        raise ValueError("output parent is not a directory")
    if path.exists() and not path.is_file():
        raise ValueError("output path is not a regular file")

    payload = "".join(" ".join(map(str, row)) + "\n" for row in matrix)
    temporary: Path | None = None
    try:
        descriptor, name = tempfile.mkstemp(
            dir=parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
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


def bounded_integer(name: str, minimum: int, maximum: int):
    def parse(value: str) -> int:
        try:
            number = int(value)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"{name} must be an integer") from exc
        if not minimum <= number <= maximum:
            raise argparse.ArgumentTypeError(
                f"{name} must be between {minimum} and {maximum}"
            )
        return number

    return parse


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Mendeley Hadamard_24.txt (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--class",
        dest="matrix_class",
        type=bounded_integer("class", 1, CLASS_COUNT),
        default=1,
        help="Hadamard equivalence-class representative, 1..60 (default: 1)",
    )
    parser.add_argument(
        "--delete-row",
        type=bounded_integer("delete row", 1, HADAMARD_ORDER),
        default=1,
        help="one-based row to delete, 1..24 (default: 1)",
    )
    parser.add_argument(
        "--delete-column",
        type=bounded_integer("delete column", 1, HADAMARD_ORDER),
        default=1,
        help="one-based column to delete, 1..24 (default: 1)",
    )
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    try:
        input_path = args.input.expanduser().resolve(strict=True)
        output_path = args.output.expanduser().absolute()
        if input_path == output_path.resolve(strict=False):
            raise ValueError("input and output paths must differ")
        matrix = load_hadamard(input_path, args.matrix_class)
        seed = delete_row_and_column(
            matrix,
            args.delete_row,
            args.delete_column,
        )
        write_matrix(output_path, seed)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    print(
        f"wrote class {args.matrix_class}, deleted row {args.delete_row} "
        f"and column {args.delete_column}: {output_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
