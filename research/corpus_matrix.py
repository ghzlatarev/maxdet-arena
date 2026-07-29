#!/usr/bin/env python3
"""Extract one strict sign matrix from the Mendeley maximal-determinant corpus."""

from __future__ import annotations

import argparse
import os
import re
import tempfile
from pathlib import Path

MAX_INPUT_BYTES = 4 * 1024 * 1024
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    REPOSITORY_ROOT
    / "runs"
    / "direct-search"
    / "reference-data"
    / "mendeley-hzf94h43c5-v1"
    / "Text"
    / "Matrices_Plus_Minus_One_MaxDet.txt"
)


def read_text(path: Path) -> str:
    """Read the bounded legacy text file without interpreting its prose."""

    if not path.is_file():
        raise ValueError("input path is not a regular file")
    try:
        with path.open("rb") as stream:
            raw = stream.read(MAX_INPUT_BYTES + 1)
    except OSError as exc:
        raise ValueError(f"cannot read input: {exc}") from exc
    if len(raw) > MAX_INPUT_BYTES:
        raise ValueError(f"input exceeds {MAX_INPUT_BYTES} bytes")
    return raw.decode("latin-1")


def assignment_offsets(text: str, order: int) -> dict[int, int]:
    """Map corpus matrix indices to offsets after their assignment operators."""

    pattern = re.compile(
        rf"^[ \t]*Mat\[{order}(?:,([0-9]+))?\][ \t]*=",
        re.MULTILINE,
    )
    offsets: dict[int, int] = {}
    for match in pattern.finditer(text):
        matrix_index = int(match.group(1) or "1")
        if matrix_index in offsets:
            raise ValueError(
                f"duplicate Mat[{order}] matrix index {matrix_index}"
            )
        offsets[matrix_index] = match.end()
    return offsets


def matrix_literal(text: str, start: int, label: str) -> str:
    """Extract a balanced Mathematica list literal."""

    cursor = start
    while cursor < len(text) and text[cursor].isspace():
        cursor += 1
    if cursor >= len(text) or text[cursor] != "{":
        raise ValueError(f"{label} has no matrix literal")

    literal_start = cursor
    depth = 0
    while cursor < len(text):
        character = text[cursor]
        if character == "{":
            depth += 1
        elif character == "}":
            depth -= 1
            if depth < 0:
                raise ValueError(f"{label} has unbalanced braces")
            if depth == 0:
                cursor += 1
                break
        cursor += 1
    else:
        raise ValueError(f"{label} has an unterminated matrix literal")

    literal_end = cursor
    while cursor < len(text) and text[cursor].isspace():
        cursor += 1
    if cursor >= len(text) or text[cursor] != ";":
        raise ValueError(f"{label} matrix is not terminated by ';'")
    return text[literal_start:literal_end]


class LiteralParser:
    """Parser for nested Mathematica lists containing only -1 and 1."""

    def __init__(self, literal: str, order: int) -> None:
        self.literal = literal
        self.order = order
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
            if len(row) >= self.order:
                raise ValueError(f"matrix row exceeds {self.order} entries")
            row.append(self.value())
        self.expect("}")
        if len(row) != self.order:
            raise ValueError(
                f"matrix row has {len(row)} entries, expected {self.order}"
            )
        return row

    def matrix(self) -> list[list[int]]:
        self.expect("{")
        matrix = [self.row()]
        while self.accept(","):
            if len(matrix) >= self.order:
                raise ValueError(f"matrix exceeds {self.order} rows")
            matrix.append(self.row())
        self.expect("}")
        self.skip_whitespace()
        if self.cursor != len(self.literal):
            raise ValueError(
                f"unexpected data at matrix-literal offset {self.cursor}"
            )
        if len(matrix) != self.order:
            raise ValueError(
                f"matrix has {len(matrix)} rows, expected {self.order}"
            )
        return matrix


def write_matrix(path: Path, matrix: list[list[int]]) -> None:
    """Atomically replace path with a strict ASCII sign matrix."""

    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not path.is_file():
        raise ValueError("output path is not a regular file")
    payload = "".join(" ".join(map(str, row)) + "\n" for row in matrix)
    temporary: Path | None = None
    try:
        descriptor, name = tempfile.mkstemp(
            dir=path.parent,
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


def positive_integer(name: str):
    def parse(value: str) -> int:
        try:
            number = int(value)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"{name} must be an integer") from exc
        if number < 1:
            raise argparse.ArgumentTypeError(f"{name} must be positive")
        return number

    return parse


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "--order",
        type=positive_integer("order"),
        required=True,
        help="matrix order in Mat[order] (for example, 22)",
    )
    parser.add_argument(
        "--matrix-index",
        type=positive_integer("matrix index"),
        default=1,
        help="one-based corpus index; an unsuffixed Mat[order] is index 1",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    try:
        input_path = args.input.expanduser().resolve(strict=True)
        output_path = args.output.expanduser().absolute()
        if input_path == output_path.resolve(strict=False):
            raise ValueError("input and output paths must differ")
        text = read_text(input_path)
        offsets = assignment_offsets(text, args.order)
        if args.matrix_index not in offsets:
            available = ", ".join(map(str, sorted(offsets))) or "none"
            raise ValueError(
                f"Mat[{args.order}] index {args.matrix_index} is absent; "
                f"available indices: {available}"
            )
        label = f"Mat[{args.order},{args.matrix_index}]"
        literal = matrix_literal(
            text,
            offsets[args.matrix_index],
            label,
        )
        matrix = LiteralParser(literal, args.order).matrix()
        write_matrix(output_path, matrix)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    print(
        f"wrote Mat[{args.order}] index {args.matrix_index}: {output_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
