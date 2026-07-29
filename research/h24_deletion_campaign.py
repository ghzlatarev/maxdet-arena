#!/usr/bin/env python3
"""Screen many row/column deletions of order-24 Hadamard representatives."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

from hadamard24_seed import (
    CLASS_COUNT,
    DEFAULT_INPUT,
    HADAMARD_ORDER,
    delete_row_and_column,
    load_hadamard,
    write_matrix,
)

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BINARY = REPOSITORY_ROOT / "build" / "research" / "reactive_tabu"
DEFAULT_SOURCE = REPOSITORY_ROOT / "research" / "reactive_tabu.cpp"
DEFAULT_DESIGN_CATALOG = (
    REPOSITORY_ROOT
    / "runs"
    / "direct-search"
    / "reference-data"
    / "spence-hadamard.24"
)
CHECKPOINT_INTERVAL = 25


@dataclass
class Result:
    matrix_class: int
    delete_row: int
    delete_column: int
    seed: int
    score: int
    matrix: list[list[int]]
    iterations: int
    exact_checks: int
    elapsed_seconds: float

    @property
    def identity(self) -> tuple[int, int, int]:
        return (self.matrix_class, self.delete_row, self.delete_column)

    def metadata(self, matrix_path: str | None = None) -> dict[str, object]:
        payload = matrix_bytes(self.matrix)
        record: dict[str, object] = {
            "class": self.matrix_class,
            "delete_row": self.delete_row,
            "delete_column": self.delete_column,
            "seed": self.seed,
            "absolute_determinant": str(self.score),
            "iterations": self.iterations,
            "exact_checks": self.exact_checks,
            "elapsed_seconds": self.elapsed_seconds,
            "matrix_sha256": hashlib.sha256(payload).hexdigest(),
        }
        if matrix_path is not None:
            record["matrix_path"] = matrix_path
        return record


def matrix_bytes(matrix: list[list[int]]) -> bytes:
    return "".join(" ".join(map(str, row)) + "\n" for row in matrix).encode(
        "ascii"
    )


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def design_catalog(path: Path) -> tuple[dict[int, list[str]], dict[int, int]]:
    """Parse Spence's 60 representatives and Hadamard-design counts."""

    try:
        text = path.read_text(encoding="ascii")
    except OSError as exc:
        raise ValueError(f"cannot read design catalog: {exc}") from exc
    headings = list(
        re.finditer(r"(?m)^[ \t]*(\d+)\.[ \t]+([01]{23})[ \t]*$", text)
    )
    cores: dict[int, list[str]] = {}
    counts: dict[int, int] = {}
    for position, heading in enumerate(headings):
        matrix_class = int(heading.group(1))
        end = (
            headings[position + 1].start()
            if position + 1 < len(headings)
            else len(text)
        )
        block = text[heading.end() : end]
        rows = [heading.group(2)] + re.findall(
            r"(?m)^[ \t]+([01]{23})[ \t]*$",
            block,
        )
        count_match = re.search(
            r"The number of non-isomorphic designs =[ \t]*(\d+)",
            block,
        )
        if matrix_class in cores or len(rows) != 23 or count_match is None:
            raise ValueError(
                f"invalid design-catalog block for class {matrix_class}"
            )
        cores[matrix_class] = rows
        counts[matrix_class] = int(count_match.group(1))
    expected = set(range(1, CLASS_COUNT + 1))
    if set(cores) != expected:
        raise ValueError("design catalog does not contain classes 1..60 exactly")
    return cores, counts


def normalized_core(matrix: list[list[int]]) -> list[str]:
    """Return the dephased (1,1)-minor using 1 for +1 and 0 for -1."""

    pivot = matrix[0][0]
    dephased = [
        [
            matrix[row][column]
            * matrix[row][0]
            * matrix[0][column]
            * pivot
            for column in range(HADAMARD_ORDER)
        ]
        for row in range(HADAMARD_ORDER)
    ]
    return [
        "".join("1" if value == 1 else "0" for value in row[1:])
        for row in dephased[1:]
    ]


def strict_matrix(path: Path, order: int) -> list[list[int]]:
    try:
        lines = path.read_text(encoding="ascii").splitlines()
    except OSError as exc:
        raise ValueError(f"cannot read result matrix: {exc}") from exc
    if len(lines) != order:
        raise ValueError(f"result has {len(lines)} rows, expected {order}")
    matrix: list[list[int]] = []
    for row_number, line in enumerate(lines, start=1):
        tokens = line.split()
        if len(tokens) != order:
            raise ValueError(
                f"result row {row_number} has {len(tokens)} entries, "
                f"expected {order}"
            )
        if any(token not in ("-1", "1") for token in tokens):
            raise ValueError(f"result row {row_number} has a non-sign entry")
        matrix.append([int(token) for token in tokens])
    return matrix


def exact_determinant(matrix: list[list[int]]) -> int:
    work = [row[:] for row in matrix]
    order = len(work)
    previous = 1
    sign = 1
    for column in range(order - 1):
        pivot_row = next(
            (
                row
                for row in range(column, order)
                if work[row][column] != 0
            ),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, order):
            for inner in range(column + 1, order):
                numerator = (
                    work[row][inner] * pivot
                    - work[row][column] * work[column][inner]
                )
                if column:
                    quotient, remainder = divmod(numerator, previous)
                    if remainder:
                        raise ArithmeticError("exact Bareiss division failed")
                    work[row][inner] = quotient
                else:
                    work[row][inner] = numerator
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def parse_classes(value: str) -> list[int]:
    classes: set[int] = set()
    for item in value.split(","):
        item = item.strip()
        if not item:
            raise argparse.ArgumentTypeError("empty class in --classes")
        try:
            matrix_class = int(item)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                f"invalid Hadamard class: {item}"
            ) from exc
        if not 1 <= matrix_class <= CLASS_COUNT:
            raise argparse.ArgumentTypeError(
                f"Hadamard class must be in 1..{CLASS_COUNT}"
            )
        classes.add(matrix_class)
    if not classes:
        raise argparse.ArgumentTypeError("--classes cannot be empty")
    return sorted(classes)


def positive_float(value: str) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number") from exc
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return number


def nonnegative_integer(value: str) -> int:
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if number < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return number


def positive_integer(value: str) -> int:
    number = nonnegative_integer(value)
    if number == 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def final_record(path: Path) -> dict[str, object]:
    try:
        records = [
            json.loads(line)
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot parse reactive log: {exc}") from exc
    if not records or records[-1].get("event") != "finished":
        raise ValueError("reactive search did not write a finished record")
    return records[-1]


def result_order(result: Result) -> tuple[int, int, int, int]:
    return (
        -result.score,
        result.matrix_class,
        result.delete_row,
        result.delete_column,
    )


def retain_top(results: list[Result], result: Result, count: int) -> None:
    results.append(result)
    results.sort(key=result_order)
    del results[count:]


def result_from_record(
    output_dir: Path,
    record: dict[str, object],
) -> Result:
    relative = Path(str(record["matrix_path"]))
    if relative.is_absolute():
        raise ValueError("resume matrix path must be relative")
    path = (output_dir / relative).resolve(strict=True)
    resolved_output = output_dir.resolve()
    if path != resolved_output and resolved_output not in path.parents:
        raise ValueError("resume matrix path escapes the output directory")
    matrix = strict_matrix(path, HADAMARD_ORDER - 1)
    payload_sha256 = hashlib.sha256(matrix_bytes(matrix)).hexdigest()
    if payload_sha256 != str(record["matrix_sha256"]):
        raise ValueError(f"resume matrix hash mismatch: {relative}")
    score = abs(exact_determinant(matrix))
    if score != int(str(record["absolute_determinant"])):
        raise ValueError(f"resume matrix determinant mismatch: {relative}")
    return Result(
        matrix_class=int(record["class"]),
        delete_row=int(record["delete_row"]),
        delete_column=int(record["delete_column"]),
        seed=int(record["seed"]),
        score=score,
        matrix=matrix,
        iterations=int(record["iterations"]),
        exact_checks=int(record["exact_checks"]),
        elapsed_seconds=float(record["elapsed_seconds"]),
    )


def load_resume(
    summary_path: Path,
    output_dir: Path,
    expected: dict[str, object],
    top_count: int,
) -> tuple[
    list[Result],
    dict[int, Result],
    dict[str, int],
    int,
    float,
    bool,
]:
    try:
        document = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot load resume summary: {exc}") from exc
    if not isinstance(document, dict):
        raise ValueError("resume summary must be a JSON object")
    for key, value in expected.items():
        if document.get(key) != value:
            raise ValueError(
                f"resume summary mismatch for {key}: "
                f"{document.get(key)!r} != {value!r}"
            )

    completed = int(document.get("tasks_completed", -1))
    selected = int(document["tasks_selected"])
    if not 0 <= completed <= selected:
        raise ValueError("resume tasks_completed is outside the task range")
    totals_record = document.get("totals")
    if not isinstance(totals_record, dict):
        raise ValueError("resume totals must be an object")
    totals = {
        "iterations": int(totals_record["iterations"]),
        "exact_checks": int(totals_record["exact_checks"]),
    }

    top_records = document.get("top")
    if not isinstance(top_records, list):
        raise ValueError("resume top must be an array")
    top = [
        result_from_record(output_dir, record)
        for record in top_records
        if isinstance(record, dict)
    ]
    if len(top) != len(top_records):
        raise ValueError("resume top contains a non-object record")
    top.sort(key=result_order)
    del top[top_count:]

    class_records = document.get("class_best")
    if not isinstance(class_records, dict):
        raise ValueError("resume class_best must be an object")
    class_best: dict[int, Result] = {}
    for key, record in class_records.items():
        if not isinstance(record, dict):
            raise ValueError("resume class_best contains a non-object record")
        result = result_from_record(output_dir, record)
        if str(result.matrix_class) != key:
            raise ValueError("resume class_best key does not match its record")
        class_best[result.matrix_class] = result

    elapsed = float(document.get("elapsed_wall_seconds", 0.0))
    if not math.isfinite(elapsed) or elapsed < 0:
        raise ValueError("resume elapsed_wall_seconds is invalid")
    return (
        top,
        class_best,
        totals,
        completed,
        elapsed,
        bool(document.get("complete", False)),
    )


def atomic_json(path: Path, document: dict[str, object]) -> None:
    payload = (
        json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode("utf-8")
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_retained(
    output_dir: Path,
    top: list[Result],
    class_best: dict[int, Result],
) -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    elite_dir = output_dir / "elites"
    elite_dir.mkdir(parents=True, exist_ok=True)
    top_records: list[dict[str, object]] = []
    for result in top:
        relative = (
            Path("elites")
            / (
                f"class-{result.matrix_class}-r{result.delete_row}"
                f"-c{result.delete_column}.matrix.txt"
            )
        )
        write_matrix(output_dir / relative, result.matrix)
        top_records.append(result.metadata(str(relative)))

    class_records: dict[str, dict[str, object]] = {}
    for matrix_class, result in sorted(class_best.items()):
        relative = Path(f"class-{matrix_class}-best.matrix.txt")
        write_matrix(output_dir / relative, result.matrix)
        class_records[str(matrix_class)] = result.metadata(str(relative))
    return top_records, class_records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "--reactive-binary",
        type=Path,
        default=DEFAULT_BINARY,
    )
    parser.add_argument(
        "--design-catalog",
        type=Path,
        default=DEFAULT_DESIGN_CATALOG,
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--classes",
        type=parse_classes,
        default=parse_classes("14,42,51"),
    )
    parser.add_argument("--seconds-per-minor", type=positive_float, default=0.25)
    parser.add_argument("--seed-base", type=nonnegative_integer, default=30000)
    parser.add_argument("--top-count", type=positive_integer, default=20)
    parser.add_argument("--shard-count", type=positive_integer, default=1)
    parser.add_argument("--shard-index", type=nonnegative_integer, default=0)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="resume the completed task prefix from output-dir/summary.json",
    )
    parser.add_argument(
        "--max-runs",
        type=positive_integer,
        help="optional prefix cap after sharding (for smoke tests)",
    )
    args = parser.parse_args()

    if args.shard_index >= args.shard_count:
        parser.error("--shard-index must be smaller than --shard-count")
    try:
        input_path = args.input.expanduser().resolve(strict=True)
        binary = args.reactive_binary.expanduser().resolve(strict=True)
        catalog_path = args.design_catalog.expanduser().resolve(strict=True)
        if not binary.is_file() or not os.access(binary, os.X_OK):
            raise ValueError("reactive binary is not an executable regular file")
        output_dir = args.output_dir.expanduser().absolute()
        output_dir.mkdir(parents=True, exist_ok=True)
        if not output_dir.is_dir():
            raise ValueError("output directory is not a directory")
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    tasks = [
        (matrix_class, delete_row, delete_column)
        for matrix_class in args.classes
        for delete_row in range(1, HADAMARD_ORDER + 1)
        for delete_column in range(1, HADAMARD_ORDER + 1)
    ]
    tasks = [
        task
        for index, task in enumerate(tasks)
        if index % args.shard_count == args.shard_index
    ]
    if args.max_runs is not None:
        tasks = tasks[: args.max_runs]

    matrices = {
        matrix_class: load_hadamard(input_path, matrix_class)
        for matrix_class in args.classes
    }
    catalog_cores, catalog_counts = design_catalog(catalog_path)
    for matrix_class, matrix in matrices.items():
        if normalized_core(matrix) != catalog_cores[matrix_class]:
            parser.error(
                f"input class {matrix_class} does not match the design catalog"
            )
    top: list[Result] = []
    class_best: dict[int, Result] = {}
    totals = {"iterations": 0, "exact_checks": 0}
    completed = 0
    previous_elapsed = 0.0
    if args.resume:
        summary_path = output_dir / "summary.json"
        expected = {
            "schema_version": 1,
            "method": "reactive-tabu-over-order-24-Hadamard-deletions",
            "classes": args.classes,
            "deletions_per_class": HADAMARD_ORDER * HADAMARD_ORDER,
            "tasks_selected": len(tasks),
            "seconds_per_minor": args.seconds_per_minor,
            "seed_base": args.seed_base,
            "shard_count": args.shard_count,
            "shard_index": args.shard_index,
            "input_sha256": file_sha256(input_path),
            "design_catalog_sha256": file_sha256(catalog_path),
            "reactive_binary_sha256": file_sha256(binary),
        }
        try:
            (
                top,
                class_best,
                totals,
                completed,
                previous_elapsed,
                already_complete,
            ) = load_resume(
                summary_path,
                output_dir,
                expected,
                args.top_count,
            )
        except ValueError as exc:
            parser.error(str(exc))
        if already_complete:
            print(f"campaign already complete: {completed}/{len(tasks)}")
            return 0
    started = time.monotonic()
    complete = False
    failure: str | None = None

    def checkpoint() -> None:
        top_records, class_records = write_retained(
            output_dir,
            top,
            class_best,
        )
        document: dict[str, object] = {
            "schema_version": 1,
            "method": "reactive-tabu-over-order-24-Hadamard-deletions",
            "complete": complete,
            "classes": args.classes,
            "deletions_per_class": HADAMARD_ORDER * HADAMARD_ORDER,
            "tasks_selected": len(tasks),
            "tasks_completed": completed,
            "seconds_per_minor": args.seconds_per_minor,
            "seed_base": args.seed_base,
            "top_count": args.top_count,
            "shard_count": args.shard_count,
            "shard_index": args.shard_index,
            "input_path": str(input_path),
            "input_sha256": file_sha256(input_path),
            "design_catalog_path": str(catalog_path),
            "design_catalog_sha256": file_sha256(catalog_path),
            "catalog_match_for_selected_classes": True,
            "inequivalent_designs_in_catalog": sum(catalog_counts.values()),
            "inequivalent_designs_in_selected_classes": sum(
                catalog_counts[matrix_class] for matrix_class in args.classes
            ),
            "reactive_binary_path": str(binary),
            "reactive_binary_sha256": file_sha256(binary),
            "reactive_source_sha256": (
                file_sha256(DEFAULT_SOURCE) if DEFAULT_SOURCE.is_file() else None
            ),
            "campaign_source_sha256": file_sha256(Path(__file__)),
            "seed_generator_source_sha256": file_sha256(
                Path(__file__).with_name("hadamard24_seed.py")
            ),
            "elapsed_wall_seconds": (
                previous_elapsed + time.monotonic() - started
            ),
            "totals": totals,
            "top": top_records,
            "class_best": class_records,
        }
        if failure is not None:
            document["failure"] = failure
        atomic_json(output_dir / "summary.json", document)

    try:
        with tempfile.TemporaryDirectory(
            dir=output_dir,
            prefix=".campaign-",
        ) as temporary_name:
            temporary = Path(temporary_name)
            seed_path = temporary / "seed.matrix.txt"
            result_path = temporary / "result.matrix.txt"
            log_path = temporary / "reactive.jsonl"
            for task_index, (
                matrix_class,
                delete_row,
                delete_column,
            ) in enumerate(tasks[completed:], start=completed):
                seed = (
                    args.seed_base
                    + (matrix_class - 1) * HADAMARD_ORDER * HADAMARD_ORDER
                    + (delete_row - 1) * HADAMARD_ORDER
                    + delete_column
                )
                start_matrix = delete_row_and_column(
                    matrices[matrix_class],
                    delete_row,
                    delete_column,
                )
                write_matrix(seed_path, start_matrix)
                result_path.unlink(missing_ok=True)
                log_path.unlink(missing_ok=True)
                command = [
                    str(binary),
                    "--start",
                    str(seed_path),
                    "--output",
                    str(result_path),
                    "--log",
                    str(log_path),
                    "--seed",
                    str(seed),
                    "--seconds",
                    str(args.seconds_per_minor),
                    "--heartbeat-seconds",
                    "0",
                ]
                try:
                    process = subprocess.run(
                        command,
                        check=False,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.PIPE,
                        text=True,
                        timeout=max(30.0, 3.0 * args.seconds_per_minor + 10.0),
                    )
                except subprocess.TimeoutExpired as exc:
                    raise RuntimeError(
                        f"class {matrix_class}, deletion "
                        f"({delete_row},{delete_column}) timed out"
                    ) from exc
                if process.returncode != 0:
                    raise RuntimeError(
                        f"class {matrix_class}, deletion "
                        f"({delete_row},{delete_column}) failed: "
                        f"{process.stderr.strip()}"
                    )
                final = final_record(log_path)
                matrix = strict_matrix(result_path, HADAMARD_ORDER - 1)
                score = abs(exact_determinant(matrix))
                logged_score = int(str(final["absolute_determinant"]))
                if score != logged_score:
                    raise ArithmeticError(
                        f"independent determinant mismatch: {score} != "
                        f"{logged_score}"
                    )
                result = Result(
                    matrix_class=matrix_class,
                    delete_row=delete_row,
                    delete_column=delete_column,
                    seed=seed,
                    score=score,
                    matrix=matrix,
                    iterations=int(final["iterations"]),
                    exact_checks=int(final["exact_checks"]),
                    elapsed_seconds=float(final["elapsed_seconds"]),
                )
                retain_top(top, result, args.top_count)
                incumbent = class_best.get(matrix_class)
                if incumbent is None or result_order(result) < result_order(
                    incumbent
                ):
                    class_best[matrix_class] = result
                totals["iterations"] += result.iterations
                totals["exact_checks"] += result.exact_checks
                completed = task_index + 1
                if completed % CHECKPOINT_INTERVAL == 0:
                    checkpoint()
    except (KeyboardInterrupt, Exception) as exc:
        failure = f"{type(exc).__name__}: {exc}"
        checkpoint()
        if isinstance(exc, KeyboardInterrupt):
            return 130
        raise

    complete = True
    checkpoint()
    if top:
        print(
            f"completed {completed}/{len(tasks)} deletion searches; "
            f"best={top[0].score} class={top[0].matrix_class} "
            f"row={top[0].delete_row} column={top[0].delete_column}"
        )
    else:
        print("no tasks selected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
