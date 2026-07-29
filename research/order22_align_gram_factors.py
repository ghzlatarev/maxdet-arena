#!/usr/bin/env python3
"""Align order-22 maximal-determinant factors to exact target row Grams.

The row Gram of a sign matrix is unchanged up to simultaneous signed row
permutation under H-equivalence.  A colored signed-double-cover graph gives
an exact certificate and an explicit signed permutation from each source
Gram to one supplied target Gram.

The output includes every oriented H-class factor, one representative for
each HT-class, exact alignment witnesses, and normalized column masks in the
same convention as ``order22_gram_shell.cpp``.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import sys
import tempfile
from typing import Any

import pynauty

from h_equivalence_audit import determinant, h_certificate, transpose


ORDER = 22
EXPECTED_ABS_DETERMINANT = 409_600_000_000_000
PINNED_PYNAUTY = "2.8.8.1"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = (
    ROOT
    / (
        "runs/direct-search/order22-plateau-harvest-wave3-20260729/"
        "closure/report.json"
    )
)
DEFAULT_OUTPUT = (
    ROOT
    / "runs/direct-search/order22-gram-aligned-factors-20260729"
)
DEFAULT_TARGETS = {
    "mendeley": (
        ROOT
        / "runs/direct-search/order22-border/mendeley-order22.matrix.txt"
    ),
    "gsds": (
        ROOT
        / (
            "runs/direct-search/order22-gsds-border-20260729/"
            "gsds-order22.matrix.txt"
        )
    ),
}
DEFAULT_SHELLS = {
    "mendeley": (
        ROOT
        / (
            "runs/direct-search/order22-factor-recovery-20260729/"
            "mendeley-shell.json"
        )
    ),
    "gsds": (
        ROOT
        / (
            "runs/direct-search/order22-factor-recovery-20260729/"
            "gsds-shell.json"
        )
    ),
}


Matrix = list[list[int]]
Gram = list[list[int]]


def sha256_bytes(contents: bytes) -> str:
    return hashlib.sha256(contents).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def display_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def read_matrix(path: Path) -> Matrix:
    try:
        matrix = [
            [int(token) for token in line.split()]
            for line in path.read_text(encoding="ascii").splitlines()
            if line.strip()
        ]
    except (OSError, UnicodeDecodeError, ValueError) as error:
        raise ValueError(f"cannot parse {path}: {error}") from error
    if len(matrix) != ORDER or any(
        len(row) != ORDER or any(value not in (-1, 1) for value in row)
        for row in matrix
    ):
        raise ValueError(f"{path}: expected exactly {ORDER}x{ORDER} signs")
    return matrix


def matrix_bytes(matrix: list[list[int]]) -> bytes:
    return "".join(
        " ".join(map(str, row)) + "\n" for row in matrix
    ).encode("ascii")


def gram(matrix: Matrix) -> Gram:
    return [
        [
            sum(
                matrix[row][index] * matrix[column][index]
                for index in range(ORDER)
            )
            for column in range(ORDER)
        ]
        for row in range(ORDER)
    ]


def validate_order22_gram(matrix: Gram, label: str) -> None:
    if len(matrix) != ORDER or any(len(row) != ORDER for row in matrix):
        raise ValueError(f"{label}: Gram is not {ORDER}x{ORDER}")
    for row in range(ORDER):
        if matrix[row][row] != ORDER:
            raise ValueError(f"{label}: Gram diagonal is not {ORDER}")
        for column in range(row):
            if matrix[row][column] != matrix[column][row]:
                raise ValueError(f"{label}: Gram is not symmetric")
            if matrix[row][column] not in (-2, 0, 2):
                raise ValueError(
                    f"{label}: unsupported off-diagonal Gram value "
                    f"{matrix[row][column]}"
                )


def signed_gram_graph(
    matrix: Gram,
) -> tuple[pynauty.Graph, dict[int, set[int]]]:
    """Encode a {-2,0,2} signed Gram up to signed row permutation.

    Vertices 2*i and 2*i+1 are the two signs of row i.  Vertex 2*n+i
    is a differently colored fiber gadget joining those two clones.
    Positive Gram edges join equal parities and negative edges cross them.
    """

    validate_order22_gram(matrix, "signed-graph input")
    adjacency = {vertex: set() for vertex in range(3 * ORDER)}

    def add_edge(left: int, right: int) -> None:
        adjacency[left].add(right)
        adjacency[right].add(left)

    for row in range(ORDER):
        add_edge(2 * ORDER + row, 2 * row)
        add_edge(2 * ORDER + row, 2 * row + 1)
    for row in range(ORDER):
        for column in range(row):
            value = matrix[row][column]
            if value == 2:
                add_edge(2 * row, 2 * column)
                add_edge(2 * row + 1, 2 * column + 1)
            elif value == -2:
                add_edge(2 * row, 2 * column + 1)
                add_edge(2 * row + 1, 2 * column)

    graph = pynauty.Graph(
        number_of_vertices=3 * ORDER,
        directed=False,
        adjacency_dict={
            vertex: sorted(neighbors)
            for vertex, neighbors in adjacency.items()
        },
        vertex_coloring=[
            set(range(2 * ORDER)),
            set(range(2 * ORDER, 3 * ORDER)),
        ],
    )
    return graph, adjacency


def graph_certificate_sha256(matrix: Gram) -> str:
    graph, _ = signed_gram_graph(matrix)
    return sha256_bytes(pynauty.certificate(graph))


def canonical_graph_isomorphism(
    source: Gram,
    target: Gram,
) -> dict[int, int] | None:
    source_graph, source_adjacency = signed_gram_graph(source)
    target_graph, target_adjacency = signed_gram_graph(target)
    if pynauty.certificate(source_graph) != pynauty.certificate(target_graph):
        return None

    source_label = pynauty.canon_label(source_graph)
    target_label = pynauty.canon_label(target_graph)
    if len(source_label) != 3 * ORDER or len(target_label) != 3 * ORDER:
        raise ArithmeticError("pynauty returned an invalid canonical labeling")
    mapping = {
        source_label[position]: target_label[position]
        for position in range(3 * ORDER)
    }
    if set(mapping) != set(range(3 * ORDER)) or set(mapping.values()) != set(
        range(3 * ORDER)
    ):
        raise ArithmeticError("canonical correspondence is not a bijection")
    for left in range(3 * ORDER):
        for right in range(left):
            source_edge = right in source_adjacency[left]
            target_edge = mapping[right] in target_adjacency[mapping[left]]
            if source_edge != target_edge:
                raise ArithmeticError(
                    "canonical correspondence failed exact edge validation"
                )
    return mapping


def align_factor(
    source: Matrix,
    target_gram: Gram,
) -> tuple[Matrix, list[int], list[int]] | None:
    """Return aligned factor, target-to-source row map, and row signs."""

    source_gram = gram(source)
    mapping = canonical_graph_isomorphism(source_gram, target_gram)
    if mapping is None:
        return None

    aligned: list[list[int] | None] = [None] * ORDER
    target_to_source: list[int | None] = [None] * ORDER
    row_signs: list[int | None] = [None] * ORDER
    for source_row in range(ORDER):
        mapped_fiber = mapping[2 * ORDER + source_row]
        if not 2 * ORDER <= mapped_fiber < 3 * ORDER:
            raise ArithmeticError("a row fiber mapped outside the fiber color")
        target_row = mapped_fiber - 2 * ORDER
        mapped_positive_clone = mapping[2 * source_row]
        mapped_negative_clone = mapping[2 * source_row + 1]
        if mapped_positive_clone // 2 != target_row:
            raise ArithmeticError("positive clone and fiber mappings disagree")
        if mapped_negative_clone // 2 != target_row:
            raise ArithmeticError("negative clone and fiber mappings disagree")
        if {
            mapped_positive_clone % 2,
            mapped_negative_clone % 2,
        } != {0, 1}:
            raise ArithmeticError("a clone pair did not map to a clone pair")
        sign = 1 if mapped_positive_clone % 2 == 0 else -1
        if aligned[target_row] is not None:
            raise ArithmeticError("two source rows mapped to one target row")
        aligned[target_row] = [
            sign * value for value in source[source_row]
        ]
        target_to_source[target_row] = source_row
        row_signs[target_row] = sign

    if any(row is None for row in aligned):
        raise ArithmeticError("signed row permutation was incomplete")
    if any(row is None for row in target_to_source):
        raise ArithmeticError("row permutation was incomplete")
    if any(sign is None for sign in row_signs):
        raise ArithmeticError("row signs were incomplete")
    result = [row for row in aligned if row is not None]
    permutation = [row for row in target_to_source if row is not None]
    signs = [sign for sign in row_signs if sign is not None]
    if gram(result) != target_gram:
        raise ArithmeticError("aligned factor does not reproduce target Gram")
    return result, permutation, signs


def normalized_column_masks(matrix: Matrix) -> list[int]:
    """Match ``normalized_factor_columns`` in order22_gram_shell.cpp."""

    masks = []
    for column in range(ORDER):
        switch_sign = matrix[0][column]
        mask = 0
        for row in range(ORDER):
            if matrix[row][column] * switch_sign == 1:
                mask |= 1 << row
        if not mask & 1:
            raise ArithmeticError("normalized column mask omitted row zero")
        masks.append(mask)
    if len(set(masks)) != ORDER:
        raise ArithmeticError("nonsingular factor has duplicate columns")
    return masks


def resolve_report_matrix(report_path: Path, recorded: str) -> Path:
    path = Path(recorded)
    if path.is_file():
        return path.resolve()
    fallback = report_path.parent / "classes" / path.name
    if fallback.is_file():
        return fallback.resolve()
    raise FileNotFoundError(f"cannot resolve report matrix {recorded}")


def write_bytes(path: Path, contents: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(contents)


def parse_named_paths(
    values: list[str],
    defaults: dict[str, Path],
    option: str,
) -> dict[str, Path]:
    if not values:
        return defaults.copy()
    paths: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"{option} must be NAME=PATH")
        name, raw_path = value.split("=", 1)
        if not name or any(
            character
            not in "abcdefghijklmnopqrstuvwxyz0123456789-"
            for character in name
        ):
            raise ValueError(
                "names must contain only lowercase letters, digits, or -"
            )
        if name in paths:
            raise ValueError(f"duplicate name for {option}: {name}")
        paths[name] = Path(raw_path).expanduser().resolve()
    return paths


def parse_shell_report(
    path: Path,
    target_factor: Matrix,
) -> tuple[list[int], dict[str, Any]]:
    report = json.loads(path.read_text(encoding="utf-8"))
    if (
        report.get("engine") != "order22-gram-shell-v1"
        or report.get("order") != ORDER
        or report.get("complete") is not True
        or report.get("assignments_completed") != 1 << (ORDER - 1)
        or report.get("factor_columns_in_shell") is not True
    ):
        raise ValueError(f"{path}: incomplete or unexpected shell report")
    raw_masks = report.get("shell_masks")
    if (
        not isinstance(raw_masks, list)
        or report.get("shell_size") != len(raw_masks)
        or any(
            type(mask) is not int
            or mask < 0
            or mask >= 1 << ORDER
            or mask & 1 == 0
            for mask in raw_masks
        )
        or len(set(raw_masks)) != len(raw_masks)
    ):
        raise ValueError(f"{path}: invalid shell masks")
    known = normalized_column_masks(target_factor)
    if sorted(report.get("known_factor_column_masks", [])) != sorted(known):
        raise ValueError(f"{path}: shell is not bound to the target factor")
    return raw_masks, report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plateau-report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument(
        "--target",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="exact target factor; defaults to Mendeley and GSDS",
    )
    parser.add_argument(
        "--shell",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="complete target shell report; defaults to Mendeley and GSDS",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--expect-h-classes", type=int, default=30)
    parser.add_argument("--expect-ht-classes", type=int, default=26)
    arguments = parser.parse_args()

    if getattr(pynauty, "__version__", None) != PINNED_PYNAUTY:
        raise RuntimeError(f"requires pynauty=={PINNED_PYNAUTY}")
    report_path = arguments.plateau_report.expanduser().resolve()
    output_dir = arguments.output_dir.expanduser().resolve()
    if output_dir.exists():
        raise FileExistsError(f"output directory already exists: {output_dir}")
    targets = parse_named_paths(arguments.target, DEFAULT_TARGETS, "--target")
    shells = parse_named_paths(arguments.shell, DEFAULT_SHELLS, "--shell")
    if len(targets) < 2:
        raise ValueError("provide at least two distinct target Grams")
    if set(shells) != set(targets):
        raise ValueError("--shell and --target names must match exactly")
    for name, path in targets.items():
        if not path.is_file():
            raise FileNotFoundError(f"missing target {name}: {path}")
    for name, path in shells.items():
        if not path.is_file():
            raise FileNotFoundError(f"missing shell {name}: {path}")

    report = json.loads(report_path.read_text(encoding="utf-8"))
    classes = report.get("classes")
    if not isinstance(classes, list):
        raise ValueError("plateau report has no classes list")
    if len(classes) != arguments.expect_h_classes:
        raise ValueError(
            f"expected {arguments.expect_h_classes} H classes, "
            f"found {len(classes)}"
        )
    ht_classes = {
        record["ht_certificate_sha256"] for record in classes
    }
    if len(ht_classes) != arguments.expect_ht_classes:
        raise ValueError(
            f"expected {arguments.expect_ht_classes} HT classes, "
            f"found {len(ht_classes)}"
        )

    target_data: dict[str, dict[str, Any]] = {}
    target_certificate_to_name: dict[str, str] = {}
    for name, path in sorted(targets.items()):
        factor = read_matrix(path)
        signed_determinant = determinant(factor)
        if abs(signed_determinant) != EXPECTED_ABS_DETERMINANT:
            raise ValueError(
                f"target {name} determinant is {signed_determinant}"
            )
        target_gram = gram(factor)
        validate_order22_gram(target_gram, f"target {name}")
        shell_masks, shell_report = parse_shell_report(
            shells[name], factor
        )
        certificate = graph_certificate_sha256(target_gram)
        if certificate in target_certificate_to_name:
            raise ValueError(
                f"targets {name} and "
                f"{target_certificate_to_name[certificate]} are equivalent"
            )
        target_certificate_to_name[certificate] = name
        gram_contents = matrix_bytes(target_gram)
        target_data[name] = {
            "factor": factor,
            "factor_path": path,
            "factor_sha256": sha256_file(path),
            "factor_determinant": signed_determinant,
            "gram": target_gram,
            "gram_bytes": gram_contents,
            "gram_sha256": sha256_bytes(gram_contents),
            "signed_gram_graph_certificate_sha256": certificate,
            "shell_path": shells[name],
            "shell_sha256": sha256_file(shells[name]),
            "shell_masks": shell_masks,
            "shell_size": len(shell_masks),
            "shell_engine": shell_report["engine"],
        }

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{output_dir.name}.",
            dir=output_dir.parent,
        )
    )
    try:
        for name, data in target_data.items():
            write_bytes(
                temporary / "targets" / f"{name}.gram.txt",
                data["gram_bytes"],
            )

        factor_records = []
        target_h_classes: dict[str, list[str]] = defaultdict(list)
        target_ht_classes: dict[str, set[str]] = defaultdict(set)
        ht_members: dict[str, list[str]] = defaultdict(list)
        ht_target: dict[str, str] = {}
        for record in sorted(
            classes,
            key=lambda item: item["h_certificate_sha256"],
        ):
            source_path = resolve_report_matrix(
                report_path, record["path"]
            )
            source = read_matrix(source_path)
            signed_determinant = determinant(source)
            if abs(signed_determinant) != EXPECTED_ABS_DETERMINANT:
                raise ValueError(
                    f"{source_path}: determinant is {signed_determinant}"
                )
            direct_h = sha256_bytes(h_certificate(source))
            transpose_h = sha256_bytes(h_certificate(transpose(source)))
            ht = min(direct_h, transpose_h)
            for key, checked in (
                ("h_certificate_sha256", direct_h),
                ("transpose_h_certificate_sha256", transpose_h),
                ("ht_certificate_sha256", ht),
            ):
                if record[key] != checked:
                    raise ArithmeticError(
                        f"{source_path}: report {key} failed verification"
                    )

            matches = []
            for name, data in target_data.items():
                alignment = align_factor(source, data["gram"])
                if alignment is not None:
                    matches.append((name, alignment))
            if len(matches) != 1:
                raise ArithmeticError(
                    f"{source_path}: expected exactly one target, "
                    f"matched {[name for name, _ in matches]}"
                )
            target_name, (aligned, permutation, signs) = matches[0]
            aligned_contents = matrix_bytes(aligned)
            aligned_sha256 = sha256_bytes(aligned_contents)
            masks = normalized_column_masks(aligned)
            if not set(masks).issubset(
                set(target_data[target_name]["shell_masks"])
            ):
                raise ArithmeticError(
                    f"{source_path}: aligned factor leaves target shell"
                )
            factor_stem = f"h-{direct_h[:12]}-to-{target_name}"
            factor_relative = Path("factors") / f"{factor_stem}.matrix.txt"
            masks_relative = (
                Path("factors") / f"{factor_stem}.support-masks.json"
            )
            write_bytes(
                temporary / factor_relative,
                aligned_contents,
            )
            masks_document = {
                "schema_version": 1,
                "factor": str(factor_relative),
                "factor_sha256": aligned_sha256,
                "normalization": (
                    "For each column, multiply by its first entry; bit r "
                    "(zero-based) is 1 iff the normalized row-r entry is +1."
                ),
                "bit_zero_is_always_one": True,
                "column_order": "matrix column order, zero-based in this list",
                "column_masks_decimal": masks,
                "column_masks_hex": [f"0x{mask:06x}" for mask in masks],
            }
            write_bytes(
                temporary / masks_relative,
                (
                    json.dumps(masks_document, indent=2, sort_keys=True)
                    + "\n"
                ).encode("utf-8"),
            )

            aligned_h = sha256_bytes(h_certificate(aligned))
            aligned_transpose_h = sha256_bytes(
                h_certificate(transpose(aligned))
            )
            if aligned_h != direct_h or aligned_transpose_h != transpose_h:
                raise ArithmeticError(
                    f"{source_path}: alignment changed an H certificate"
                )
            if gram(aligned) != target_data[target_name]["gram"]:
                raise ArithmeticError(
                    f"{source_path}: final exact Gram check failed"
                )

            factor_record = {
                "source_path": display_path(source_path),
                "source_sha256": sha256_file(source_path),
                "source_determinant": signed_determinant,
                "h_certificate_sha256": direct_h,
                "transpose_h_certificate_sha256": transpose_h,
                "ht_certificate_sha256": ht,
                "target": target_name,
                "aligned_factor_path": str(factor_relative),
                "aligned_factor_sha256": aligned_sha256,
                "aligned_factor_determinant": determinant(aligned),
                "support_masks_path": str(masks_relative),
                "normalized_column_masks_decimal": masks,
                "alignment_witness": {
                    "equation": (
                        "aligned[target_row] = row_sign * "
                        "source[source_row]"
                    ),
                    "indices": "target_to_source_row is 1-based",
                    "target_to_source_row": [
                        index + 1 for index in permutation
                    ],
                    "target_row_sign": signs,
                    "verified_exact_equation": (
                        f"aligned * aligned^T = targets/"
                        f"{target_name}.gram.txt"
                    ),
                },
            }
            factor_records.append(factor_record)
            target_h_classes[target_name].append(direct_h)
            target_ht_classes[target_name].add(ht)
            ht_members[ht].append(direct_h)
            previous_target = ht_target.setdefault(ht, target_name)
            if previous_target != target_name:
                raise ArithmeticError(
                    f"HT class {ht} split across target Grams"
                )

        ht_records = []
        for ht in sorted(ht_classes):
            members = sorted(ht_members[ht])
            if not members:
                raise ArithmeticError(f"HT class {ht} has no H factors")
            representative = ht if ht in members else members[0]
            factor_record = next(
                item
                for item in factor_records
                if item["h_certificate_sha256"] == representative
            )
            ht_records.append(
                {
                    "ht_certificate_sha256": ht,
                    "target": ht_target[ht],
                    "h_certificate_sha256_members": members,
                    "representative_h_certificate_sha256": representative,
                    "representative_aligned_factor_path": factor_record[
                        "aligned_factor_path"
                    ],
                    "representative_support_masks_path": factor_record[
                        "support_masks_path"
                    ],
                }
            )

        target_records = []
        for name, data in sorted(target_data.items()):
            known_factor_records = sorted(
                (
                    record
                    for record in factor_records
                    if record["target"] == name
                ),
                key=lambda record: record["h_certificate_sha256"],
            )
            known_list_relative = (
                Path("cpsat") / f"{name}-known-factors.txt"
            )
            known_list_contents = "".join(
                f"../{record['aligned_factor_path']}\n"
                for record in known_factor_records
            ).encode("utf-8")
            write_bytes(
                temporary / known_list_relative,
                known_list_contents,
            )
            cpsat_input_relative = Path("cpsat") / f"{name}-inputs.json"
            cpsat_input = {
                "schema_version": 1,
                "target": name,
                "solver": "research/order22_gram_factor_cpsat.py",
                "shell_report": display_path(data["shell_path"]),
                "target_factor": display_path(data["factor_path"]),
                "known_factor_list": str(known_list_relative),
                "known_factor_count": len(known_factor_records),
                "known_factors": [
                    {
                        "h_certificate_sha256": record[
                            "h_certificate_sha256"
                        ],
                        "ht_certificate_sha256": record[
                            "ht_certificate_sha256"
                        ],
                        "hint_factor": display_path(
                            output_dir / record["aligned_factor_path"]
                        ),
                        "exclude_factor": display_path(
                            output_dir / record["aligned_factor_path"]
                        ),
                        "normalized_column_masks_decimal": record[
                            "normalized_column_masks_decimal"
                        ],
                    }
                    for record in known_factor_records
                ],
                "usage": (
                    "Choose any known_factors[].hint_factor for "
                    "--hint-factor. Pass known_factor_list to "
                    "--exclude-factor-list to exclude every exact known "
                    "representative support."
                ),
                "command_template": shlex.join(
                    [
                        "/tmp/maxdet-order22-factor-venv/bin/python",
                        "research/order22_gram_factor_cpsat.py",
                        "--shell-report",
                        display_path(data["shell_path"]),
                        "--factor",
                        display_path(data["factor_path"]),
                        "--hint-factor",
                        "HINT_FACTOR",
                        "--exclude-factor-list",
                        display_path(output_dir / known_list_relative),
                        "--output",
                        "OUTPUT_MATRIX",
                        "--metadata",
                        "OUTPUT_METADATA",
                    ]
                ),
                "exclusion_boundary": (
                    "Each listed exclusion blocks one exact normalized "
                    "column support. It does not by itself block the full "
                    "signed-Gram-automorphism orbit of that support."
                ),
            }
            write_bytes(
                temporary / cpsat_input_relative,
                (
                    json.dumps(cpsat_input, indent=2, sort_keys=True)
                    + "\n"
                ).encode("utf-8"),
            )
            target_records.append(
                {
                    "name": name,
                    "source_factor_path": display_path(data["factor_path"]),
                    "source_factor_sha256": data["factor_sha256"],
                    "source_factor_determinant": data[
                        "factor_determinant"
                    ],
                    "exact_gram_path": f"targets/{name}.gram.txt",
                    "exact_gram_sha256": data["gram_sha256"],
                    "exact_gram_determinant": str(
                        data["factor_determinant"]
                        * data["factor_determinant"]
                    ),
                    "signed_gram_graph_certificate_sha256": data[
                        "signed_gram_graph_certificate_sha256"
                    ],
                    "complete_shell_report": display_path(
                        data["shell_path"]
                    ),
                    "complete_shell_report_sha256": data["shell_sha256"],
                    "complete_shell_size": data["shell_size"],
                    "aligned_factor_columns_checked_in_shell": (
                        len(known_factor_records) * ORDER
                    ),
                    "cpsat_inputs_path": str(cpsat_input_relative),
                    "cpsat_known_factor_list_path": str(
                        known_list_relative
                    ),
                    "h_class_count": len(target_h_classes[name]),
                    "h_certificate_sha256": sorted(
                        target_h_classes[name]
                    ),
                    "ht_class_count": len(target_ht_classes[name]),
                    "ht_certificate_sha256": sorted(
                        target_ht_classes[name]
                    ),
                }
            )

        script_path = Path(__file__).resolve()
        manifest = {
            "schema_version": 1,
            "claim": (
                "Every listed aligned factor is an exact order-22 sign "
                "factor of exactly one listed target row Gram. The signed "
                "row permutation witness and the equality A A^T = G were "
                "checked exactly."
            ),
            "claim_boundary": (
                "This enumerates factors recovered in the supplied seeded "
                "one-flip plateau report; it is not a complete enumeration "
                "of all maximal-determinant order-22 H or HT classes."
            ),
            "order": ORDER,
            "target_absolute_factor_determinant": (
                EXPECTED_ABS_DETERMINANT
            ),
            "plateau_report": display_path(report_path),
            "plateau_report_sha256": sha256_file(report_path),
            "h_class_count": len(factor_records),
            "ht_class_count": len(ht_records),
            "target_gram_class_count": len(target_records),
            "targets": target_records,
            "ht_classes": ht_records,
            "h_factors": factor_records,
            "support_mask_convention": {
                "compatible_with": (
                    "research/order22_gram_shell.cpp "
                    "normalized_factor_columns"
                ),
                "definition": (
                    "For each column, multiply by its row-0 sign; bit r is "
                    "1 iff the normalized entry at row r is +1."
                ),
                "bits": ORDER,
                "bit_numbering": "zero-based, least-significant bit is row 0",
            },
            "dependencies": {
                "python": "3",
                "pynauty": PINNED_PYNAUTY,
            },
            "generator": {
                "path": display_path(script_path),
                "sha256": sha256_file(script_path),
                "command": shlex.join(
                    [
                        "/tmp/maxdet-h-audit/bin/python",
                        "research/order22_align_gram_factors.py",
                        *sys.argv[1:],
                    ]
                ),
            },
        }
        write_bytes(
            temporary / "manifest.json",
            (
                json.dumps(manifest, indent=2, sort_keys=True) + "\n"
            ).encode("utf-8"),
        )
        os.replace(temporary, output_dir)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise

    print(
        f"aligned {len(factor_records)} H classes / "
        f"{len(ht_records)} HT classes to {len(target_records)} exact Grams"
    )
    print(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
