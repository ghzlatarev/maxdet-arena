#!/usr/bin/env python3
"""Exact odd alternating-cycle pilot for dephased order-23 matrices.

After sign dephasing, an order-23 sign matrix is equivalent to a 22-by-22
binary matrix B and

    det(A) = 2^22 det(B).

Treat B as a bipartite adjacency matrix.  Its determinant is the signed
imbalance of its perfect matchings.  Relative to a fixed perfect matching,
an alternating cycle on k matched pairs changes matching sign by
(-1)^(k-1).  This helper chooses a matching with the sign of det(B), permutes
it to the diagonal, and exactly tests absent odd cycles in two ways:

* completion: add the k cycle edges;
* exchange: remove the k diagonal matching edges and add the cycle edges.

Only the retained family winners are materialized.  Every retained winner is
recomputed with integer Bareiss elimination on both its binary core and its
23-by-23 sign matrix.  The rank-k determinant lemma is only an exact screening
accelerator; it is never the final audit.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
from itertools import combinations, permutations
import json
import os
from pathlib import Path
import sys
import tempfile
import time
from typing import Iterable, Sequence


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from maxdet.contract import load_contract, load_matrix  # noqa: E402
from maxdet.exact import (  # noqa: E402
    bareiss_determinant,
    matrix_text,
    normalize_signs,
)
from maxdet.receipt import canonical_json_bytes, verify_matrix  # noqa: E402


ORDER = 23
CORE_ORDER = ORDER - 1
SCALE = 1 << CORE_ORDER
METHOD = "dephased-bipartite-odd-matching-cycle-v1"

Matrix = list[list[int]]
Cycle = tuple[int, ...]
Update = tuple[int, int, int]


@dataclass(frozen=True)
class InputSpec:
    label: str
    path: Path


@dataclass
class FamilyBest:
    length: int
    mode: str
    cycles: int = 0
    candidates: int = 0
    strict_source_improvements: int = 0
    strict_frontier_wins: int = 0
    best_quotient: int | None = None
    best_cycle: Cycle | None = None

    def consider(
        self,
        quotient: int,
        cycle: Cycle,
        source_absolute_quotient: int,
        frontier: int,
    ) -> None:
        self.candidates += 1
        absolute = abs(quotient)
        if absolute > source_absolute_quotient:
            self.strict_source_improvements += 1
        if absolute * SCALE > frontier:
            self.strict_frontier_wins += 1
        if (
            self.best_quotient is None
            or absolute > abs(self.best_quotient)
            or (
                absolute == abs(self.best_quotient)
                and cycle < self.best_cycle
            )
        ):
            self.best_quotient = quotient
            self.best_cycle = cycle


def parse_input_spec(value: str) -> InputSpec:
    label, separator, raw_path = value.partition("=")
    if not separator or not label or not raw_path:
        raise argparse.ArgumentTypeError("input must have the form LABEL=PATH")
    if any(character not in "abcdefghijklmnopqrstuvwxyz0123456789-_" for character in label):
        raise argparse.ArgumentTypeError(
            "input label may contain only lowercase letters, digits, '-' and '_'"
        )
    return InputSpec(label=label, path=Path(raw_path))


def parse_cycle_lengths(value: str) -> tuple[int, ...]:
    try:
        lengths = tuple(int(token) for token in value.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "cycle lengths must be comma-separated integers"
        ) from error
    if not lengths or len(set(lengths)) != len(lengths):
        raise argparse.ArgumentTypeError("cycle lengths must be nonempty and unique")
    if any(length < 3 or length > CORE_ORDER or length % 2 == 0 for length in lengths):
        raise argparse.ArgumentTypeError(
            "cycle lengths must be odd integers from 3 through 21"
        )
    return lengths


def atomic_write_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(payload)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
        directory = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        temporary.unlink(missing_ok=True)


def permutation_sign(permutation: Sequence[int]) -> int:
    inversions = sum(
        permutation[first] > permutation[second]
        for first in range(len(permutation))
        for second in range(first + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def signed_perfect_matching(
    matrix: Matrix,
    target_sign: int,
    node_limit: int,
) -> tuple[tuple[int, ...], int]:
    """Find a perfect matching of requested parity by bounded exact DFS."""

    order = len(matrix)
    matching = [-1] * order
    used_columns = [False] * order
    visited_nodes = 0

    def visit(assigned: int) -> tuple[int, ...] | None:
        nonlocal visited_nodes
        visited_nodes += 1
        if visited_nodes > node_limit:
            raise ValueError(
                f"matching parity search exceeded {node_limit} DFS nodes"
            )
        if assigned == order:
            result = tuple(matching)
            return result if permutation_sign(result) == target_sign else None

        choices: list[tuple[int, int, list[int]]] = []
        for row in range(order):
            if matching[row] >= 0:
                continue
            available = [
                column
                for column, value in enumerate(matrix[row])
                if value and not used_columns[column]
            ]
            if not available:
                return None
            choices.append((len(available), row, available))
        _, row, available = min(choices)

        for column in available:
            matching[row] = column
            used_columns[column] = True
            result = visit(assigned + 1)
            if result is not None:
                return result
            used_columns[column] = False
            matching[row] = -1
        return None

    result = visit(0)
    if result is None:
        raise ValueError(
            "no perfect matching with determinant-majority sign was found"
        )
    return result, visited_nodes


def exact_adjugate(matrix: Matrix) -> Matrix:
    order = len(matrix)
    adjugate = [[0] * order for _ in range(order)]
    for row in range(order):
        for column in range(order):
            minor = [
                [
                    matrix[minor_row][minor_column]
                    for minor_column in range(order)
                    if minor_column != column
                ]
                for minor_row in range(order)
                if minor_row != row
            ]
            cofactor = bareiss_determinant(minor)
            if (row + column) % 2:
                cofactor = -cofactor
            adjugate[column][row] = cofactor
    return adjugate


def check_adjugate(matrix: Matrix, determinant: int, adjugate: Matrix) -> None:
    order = len(matrix)
    for row in range(order):
        for column in range(order):
            product = sum(
                matrix[row][inner] * adjugate[inner][column]
                for inner in range(order)
            )
            expected = determinant if row == column else 0
            if product != expected:
                raise ArithmeticError("exact adjugate identity failed")


def determinant_after_updates(
    determinant: int,
    adjugate: Matrix,
    updates: Sequence[Update],
) -> int:
    """Evaluate a binary entry update by the exact determinant lemma."""

    rank = len(updates)
    if rank == 0:
        return determinant
    lemma_matrix = [[0] * rank for _ in range(rank)]
    for outer, (_, column, _) in enumerate(updates):
        for inner, (row, _, delta) in enumerate(updates):
            lemma_matrix[outer][inner] = (
                (determinant if outer == inner else 0)
                + adjugate[column][row] * delta
            )
    numerator = bareiss_determinant(lemma_matrix)
    denominator = determinant ** (rank - 1)
    quotient, remainder = divmod(numerator, denominator)
    if remainder:
        raise ArithmeticError("determinant-lemma division was not exact")
    return quotient


def dephased_binary_core(sign_matrix: Matrix) -> tuple[Matrix, Matrix]:
    dephased = normalize_signs(sign_matrix)
    binary = [
        [(1 - dephased[row + 1][column + 1]) // 2 for column in range(CORE_ORDER)]
        for row in range(CORE_ORDER)
    ]
    return dephased, binary


def sign_matrix_from_matched_core(
    matched_core: Matrix,
    matching: Sequence[int],
) -> Matrix:
    binary = [[0] * CORE_ORDER for _ in range(CORE_ORDER)]
    for row in range(CORE_ORDER):
        for matched_column, original_column in enumerate(matching):
            binary[row][original_column] = matched_core[row][matched_column]
    sign_matrix = [[1] * ORDER for _ in range(ORDER)]
    for row in range(CORE_ORDER):
        for column in range(CORE_ORDER):
            sign_matrix[row + 1][column + 1] = 1 - 2 * binary[row][column]
    return sign_matrix


def absent_cycles(matrix: Matrix, length: int) -> Iterable[Cycle]:
    """Yield each directed simple cycle once, modulo cyclic rotation."""

    for vertex_set in combinations(range(len(matrix)), length):
        first = vertex_set[0]
        for tail in permutations(vertex_set[1:]):
            cycle = (first,) + tail
            if all(
                matrix[cycle[index]][cycle[(index + 1) % length]] == 0
                for index in range(length)
            ):
                yield cycle


def cycle_updates(cycle: Cycle, mode: str) -> tuple[Update, ...]:
    additions = tuple(
        (cycle[index], cycle[(index + 1) % len(cycle)], 1)
        for index in range(len(cycle))
    )
    if mode == "completion":
        return additions
    if mode == "exchange":
        removals = tuple((vertex, vertex, -1) for vertex in cycle)
        return removals + additions
    raise ValueError(f"unknown cycle mode: {mode}")


def apply_updates(matrix: Matrix, updates: Sequence[Update]) -> Matrix:
    result = [row[:] for row in matrix]
    for row, column, delta in updates:
        updated = result[row][column] + delta
        if updated not in (0, 1):
            raise ArithmeticError("cycle update left the binary domain")
        result[row][column] = updated
    return result


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPOSITORY_ROOT))
    except ValueError:
        return str(path.resolve())


def finalize_report(report_path: Path) -> None:
    """Bind externally generated arena receipts and independently reverify."""

    resolved_report = (
        report_path
        if report_path.is_absolute()
        else REPOSITORY_ROOT / report_path
    ).resolve(strict=True)
    report = json.loads(resolved_report.read_text(encoding="ascii"))
    if report.get("method") != METHOD:
        raise ValueError("report method does not match this helper")

    contract = load_contract(REPOSITORY_ROOT / "challenge.json")
    inventory: list[dict[str, object]] = []
    for input_report in report["inputs"]:
        for retained in input_report["retained"]:
            matrix_path = Path(retained["matrix"])
            if not matrix_path.is_absolute():
                matrix_path = REPOSITORY_ROOT / matrix_path
            matrix_path = matrix_path.resolve(strict=True)
            receipt_path = matrix_path.with_name(
                matrix_path.name.replace(".matrix.txt", ".receipt.json")
            )
            raw_receipt = receipt_path.read_bytes()
            receipt = json.loads(raw_receipt)
            if canonical_json_bytes(receipt) != raw_receipt:
                raise ValueError(f"{receipt_path}: receipt is not canonical JSON")

            independently_verified = verify_matrix(matrix_path, contract)
            if independently_verified.receipt != receipt:
                raise ValueError(
                    f"{receipt_path}: receipt differs from independent verification"
                )
            if (
                receipt["score"]["absolute_determinant"]
                != retained["absolute_determinant"]
                or receipt["matrix"]["raw_sha256"]
                != retained["matrix_sha256"]
            ):
                raise ValueError(
                    f"{receipt_path}: receipt does not bind retained result"
                )

            binding = {
                "path": display_path(receipt_path),
                "file_sha256": hashlib.sha256(raw_receipt).hexdigest(),
                "receipt_sha256": receipt["receipt_sha256"],
                "absolute_determinant": receipt["score"][
                    "absolute_determinant"
                ],
                "independent_reverification": "passed",
            }
            retained["arena_receipt"] = binding
            inventory.append(
                {
                    "input": input_report["label"],
                    "cycle_length": retained["cycle_length"],
                    "mode": retained["mode"],
                    "matrix": retained["matrix"],
                    "matrix_sha256": retained["matrix_sha256"],
                    "receipt": binding["path"],
                    "receipt_file_sha256": binding["file_sha256"],
                    "receipt_internal_sha256": binding["receipt_sha256"],
                    "absolute_determinant": binding["absolute_determinant"],
                }
            )

    inventory.sort(key=lambda item: str(item["matrix"]))
    inventory_payload = b"".join(
        (
            json.dumps(item, sort_keys=True, separators=(",", ":")) + "\n"
        ).encode("ascii")
        for item in inventory
    )
    report["arena_audit_status"] = "passed"
    report["arena_audited_retained_count"] = len(inventory)
    report["receipt_inventory_sha256"] = hashlib.sha256(
        inventory_payload
    ).hexdigest()
    atomic_write_bytes(
        resolved_report,
        (json.dumps(report, indent=2, sort_keys=True) + "\n").encode("ascii"),
    )
    print(
        json.dumps(
            {
                "report": display_path(resolved_report),
                "arena_audited_retained_count": len(inventory),
                "receipt_inventory_sha256": report[
                    "receipt_inventory_sha256"
                ],
            },
            sort_keys=True,
        )
    )


def run_input(
    spec: InputSpec,
    output_directory: Path,
    cycle_lengths: Sequence[int],
    frontier: int,
    matching_node_limit: int,
) -> dict[str, object]:
    contract = load_contract(REPOSITORY_ROOT / "challenge.json")
    input_path = (
        spec.path if spec.path.is_absolute() else REPOSITORY_ROOT / spec.path
    ).resolve(strict=True)
    raw_bytes, sign_matrix = load_matrix(input_path, contract)
    source_determinant = bareiss_determinant(sign_matrix)
    dephased, binary = dephased_binary_core(sign_matrix)
    binary_determinant = bareiss_determinant(binary)
    if source_determinant != SCALE * binary_determinant:
        raise ArithmeticError("dephasing determinant identity failed")
    if dephased[0] != [1] * ORDER or any(row[0] != 1 for row in dephased):
        raise ArithmeticError("sign dephasing failed")

    target_sign = 1 if binary_determinant > 0 else -1
    matching, matching_nodes = signed_perfect_matching(
        binary, target_sign, matching_node_limit
    )
    matched = [
        [binary[row][matching[column]] for column in range(CORE_ORDER)]
        for row in range(CORE_ORDER)
    ]
    if any(matched[index][index] != 1 for index in range(CORE_ORDER)):
        raise ArithmeticError("selected matching did not become the diagonal")
    matched_determinant = bareiss_determinant(matched)
    if matched_determinant != abs(binary_determinant):
        raise ArithmeticError("matching orientation did not make det(B) positive")

    adjugate = exact_adjugate(matched)
    check_adjugate(matched, matched_determinant, adjugate)
    families = {
        (length, mode): FamilyBest(length=length, mode=mode)
        for length in cycle_lengths
        for mode in ("completion", "exchange")
    }

    started = time.perf_counter()
    for length in cycle_lengths:
        for cycle in absent_cycles(matched, length):
            for mode in ("completion", "exchange"):
                family = families[length, mode]
                family.cycles += 1
                quotient = determinant_after_updates(
                    matched_determinant,
                    adjugate,
                    cycle_updates(cycle, mode),
                )
                family.consider(
                    quotient,
                    cycle,
                    matched_determinant,
                    frontier,
                )
    elapsed = time.perf_counter() - started

    retained: list[dict[str, object]] = []
    for key in sorted(families):
        family = families[key]
        if family.best_cycle is None or family.best_quotient is None:
            continue
        updates = cycle_updates(family.best_cycle, family.mode)
        candidate_matched = apply_updates(matched, updates)
        direct_core_determinant = bareiss_determinant(candidate_matched)
        if direct_core_determinant != family.best_quotient:
            raise ArithmeticError("retained core failed direct Bareiss audit")
        candidate = sign_matrix_from_matched_core(candidate_matched, matching)
        direct_sign_determinant = bareiss_determinant(candidate)
        expected_sign_determinant = (
            SCALE * permutation_sign(matching) * direct_core_determinant
        )
        if direct_sign_determinant != expected_sign_determinant:
            raise ArithmeticError("retained sign matrix failed Bareiss audit")

        artifact_name = (
            f"{spec.label}-k{family.length}-{family.mode}-best.matrix.txt"
        )
        artifact_path = output_directory / artifact_name
        atomic_write_bytes(
            artifact_path,
            matrix_text(candidate).encode("ascii"),
        )
        retained.append(
            {
                "cycle_length": family.length,
                "mode": family.mode,
                "cycle_one_based": [
                    vertex + 1 for vertex in family.best_cycle
                ],
                "binary_updates_one_based": [
                    {
                        "row": row + 1,
                        "column": column + 1,
                        "delta": delta,
                    }
                    for row, column, delta in updates
                ],
                "core_determinant": str(direct_core_determinant),
                "determinant": str(direct_sign_determinant),
                "absolute_determinant": str(abs(direct_sign_determinant)),
                "source_delta": str(
                    abs(direct_sign_determinant) - abs(source_determinant)
                ),
                "frontier_delta": str(
                    abs(direct_sign_determinant) - frontier
                ),
                "matrix": display_path(artifact_path),
                "matrix_sha256": sha256_file(artifact_path),
                "exact_bareiss_core_audit": "passed",
                "exact_bareiss_sign_matrix_audit": "passed",
                "arena_receipt": None,
            }
        )

    family_reports = []
    for key in sorted(families):
        family = families[key]
        family_reports.append(
            {
                "cycle_length": family.length,
                "mode": family.mode,
                "absent_cycles": family.cycles,
                "candidates": family.candidates,
                "strict_source_improvements": family.strict_source_improvements,
                "strict_frontier_wins": family.strict_frontier_wins,
                "best_core_determinant": (
                    None
                    if family.best_quotient is None
                    else str(family.best_quotient)
                ),
                "best_absolute_determinant": (
                    None
                    if family.best_quotient is None
                    else str(abs(family.best_quotient) * SCALE)
                ),
            }
        )

    return {
        "label": spec.label,
        "path": display_path(input_path),
        "matrix_sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "source_determinant": str(source_determinant),
        "source_absolute_determinant": str(abs(source_determinant)),
        "binary_core_determinant": str(binary_determinant),
        "matching_permutation_one_based": [
            column + 1 for column in matching
        ],
        "matching_sign": permutation_sign(matching),
        "matching_search_nodes": matching_nodes,
        "matched_core_determinant": str(matched_determinant),
        "adjugate_identity": "passed",
        "elapsed_seconds": round(elapsed, 6),
        "families": family_reports,
        "retained": retained,
    }


def self_test() -> None:
    matrices = [
        [[1, 0, 1], [1, 1, 0], [0, 1, 1]],
        [[1, 1, 0], [0, 1, 1], [1, 0, 0]],
    ]
    update_sets = [
        ((0, 1, 1),),
        ((0, 1, 1), (1, 1, -1)),
        ((0, 1, 1), (1, 2, 1), (2, 0, -1)),
    ]
    checks = 0
    for matrix in matrices:
        determinant = bareiss_determinant(matrix)
        if determinant == 0:
            continue
        adjugate = exact_adjugate(matrix)
        check_adjugate(matrix, determinant, adjugate)
        for updates in update_sets:
            try:
                updated = apply_updates(matrix, updates)
            except ArithmeticError:
                continue
            predicted = determinant_after_updates(
                determinant, adjugate, updates
            )
            observed = bareiss_determinant(updated)
            if predicted != observed:
                raise ArithmeticError("determinant-lemma self-test failed")
            checks += 1
    if checks < 2:
        raise ArithmeticError("insufficient determinant-lemma self-test coverage")
    print(f"self-test passed: {checks} exact update checks")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        action="append",
        type=parse_input_spec,
        default=[],
        help="repeatable LABEL=PATH input",
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--cycle-lengths",
        type=parse_cycle_lengths,
        default=(3, 5),
    )
    parser.add_argument("--matching-node-limit", type=int, default=1_000_000)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument(
        "--finalize-report",
        type=Path,
        help=(
            "bind existing .receipt.json files produced by ./arena verify "
            "and independently reverify every retained matrix"
        ),
    )
    args = parser.parse_args()

    if args.finalize_report is not None:
        finalize_report(args.finalize_report)
        if not args.input:
            return 0
    if args.self_test:
        self_test()
        if not args.input:
            return 0
    if not args.input:
        parser.error("at least one --input is required")
    if args.output_dir is None:
        parser.error("--output-dir is required")
    if args.matching_node_limit < CORE_ORDER:
        parser.error("--matching-node-limit is too small")
    labels = [spec.label for spec in args.input]
    if len(set(labels)) != len(labels):
        parser.error("input labels must be unique")

    output_directory = (
        args.output_dir
        if args.output_dir.is_absolute()
        else REPOSITORY_ROOT / args.output_dir
    )
    if output_directory.exists() and any(output_directory.iterdir()):
        parser.error("output directory must be absent or empty")
    output_directory.mkdir(parents=True, exist_ok=True)

    frontier_data = json.loads(
        (REPOSITORY_ROOT / "data/frontier.json").read_text(encoding="utf-8")
    )
    frontier = max(
        int(frontier_data["target_to_beat"]["absolute_determinant"]),
        int(frontier_data["arena_best"]["absolute_determinant"]),
    )
    started = time.perf_counter()
    inputs = [
        run_input(
            spec,
            output_directory,
            args.cycle_lengths,
            frontier,
            args.matching_node_limit,
        )
        for spec in args.input
    ]
    report = {
        "schema_version": 1,
        "method": METHOD,
        "representation": (
            "dephased 22x22 binary matrix as a bipartite graph"
        ),
        "theorem": (
            "relative to a perfect matching, an alternating cycle on k "
            "matched pairs changes permutation sign by (-1)^(k-1)"
        ),
        "search_consequence": (
            "complete or degree-preservingly exchange absent odd cycles "
            "relative to a determinant-majority-sign matching"
        ),
        "claim_boundary": (
            "Exact only for the enumerated absent cycles of lengths listed, "
            "at the fixed first-row/first-column dephasing and one "
            "deterministic majority-sign matching per input. It is not a "
            "complete neighborhood audit, optimality result, or novelty claim."
        ),
        "frontier_absolute_determinant": str(frontier),
        "cycle_lengths": list(args.cycle_lengths),
        "input_count": len(inputs),
        "elapsed_seconds": round(time.perf_counter() - started, 6),
        "inputs": inputs,
        "retained_count": sum(
            len(input_report["retained"]) for input_report in inputs
        ),
        "strict_frontier_win_count": sum(
            family["strict_frontier_wins"]
            for input_report in inputs
            for family in input_report["families"]
        ),
        "arena_audit_status": "pending external ./arena verify",
    }
    report_path = output_directory / "report.json"
    atomic_write_bytes(
        report_path,
        (json.dumps(report, indent=2, sort_keys=True) + "\n").encode("ascii"),
    )
    print(
        json.dumps(
            {
                "report": display_path(report_path),
                "retained_count": report["retained_count"],
                "strict_frontier_win_count": report[
                    "strict_frontier_win_count"
                ],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
