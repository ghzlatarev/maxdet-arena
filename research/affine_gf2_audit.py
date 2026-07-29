#!/usr/bin/env python3
"""Audit entry-flip cubes as exact affine spaces after sign dephasing.

An order-23 sign matrix has a canonical 22-by-22 dephased core.  Entry
flips act linearly on that core over GF(2):

* an interior entry toggles one core bit;
* a first-row entry toggles one core column;
* a first-column entry toggles one core row; and
* the top-left entry toggles every core bit.

This tool converts an input cube to that exact affine representation, computes
a canonical reduced-row-echelon certificate, and compares it with completed
cube reports retained from the principal-minor campaigns in this repository.
Before comparison it binds every report to its start/support files, campaign
summary, assignment count, and engine fingerprint.  It performs no determinant
search and writes no files.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import random
import sys
from typing import Any, Iterable, Optional, Sequence


ORDER = 23
CORE_ORDER = ORDER - 1
CORE_BITS = CORE_ORDER * CORE_ORDER
VECTOR_BYTES = (CORE_BITS + 7) // 8

SCHEMA_VERSION = 2
METHOD = "exact-dephased-gf2-affine-rref-v1"
AFFINE_TAG = b"maxdet-dephased-affine-gf2-v1\0"
INVENTORY_TAG = b"maxdet-dephased-cube-inventory-v1\0"
REMAINDER_TAG = b"maxdet-dephased-union-remainder-v1\0"
STATE_UNION_TAG = b"maxdet-dephased-state-union-v1\0"
SUPPORT_TAG = b"maxdet-entry-support-v1\0"
ENTRY_CUBE_TAG = b"maxdet-entry-cube-v1\0"
EVIDENCE_TAG = b"maxdet-affine-audit-evidence-inventory-v1\0"
ENGINE = "fast-principal-minor-entry-cube-v1"
ENGINE_DIMENSION = 27
ENGINE_ASSIGNMENTS = 1 << ENGINE_DIMENSION
DEEP_DIMENSION = 32
DEEP_LEAF_COUNT = 1 << (DEEP_DIMENSION - ENGINE_DIMENSION)
DEEP_ASSIGNMENTS = 1 << DEEP_DIMENSION

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BATCH_DIRECTORY = Path(
    "runs/direct-search/fast-principal-cube/"
    "batch200-mixed-h012-20260728"
)
DEFAULT_LNPS_DIRECTORY = Path(
    "runs/direct-search/fast-principal-cube/"
    "lnps-beam-h012-20260728"
)
DEFAULT_DEEP32_DIRECTORIES = (
    Path(
        "runs/direct-search/fast-principal-cube/"
        "partition32-h0-a0-20260728"
    ),
    Path(
        "runs/direct-search/fast-principal-cube/"
        "partition32-h2-bridge-20260728"
    ),
)
FAMILY_ORDER = {"batch200": 0, "lnps": 1, "deep32": 2}

Coordinate = tuple[int, int]
Matrix = tuple[tuple[int, ...], ...]
Basis = tuple[int, ...]


@dataclass(frozen=True)
class Cube:
    """One affine cube in the canonical dephased-core coordinate space."""

    family: str
    label: str
    start_path: Path
    start_sha256: str
    support: tuple[Coordinate, ...]
    support_set_certificate_sha256: str
    base: int
    basis: Basis
    affine_certificate_sha256: str

    @property
    def dimension(self) -> int:
        return len(self.support)

    @property
    def rank(self) -> int:
        return len(self.basis)


@dataclass(frozen=True)
class Intersection:
    """Exact pairwise affine-intersection result."""

    intersects: bool
    union_rank: int
    delta_remainder: int
    dimension: Optional[int] = None
    point: Optional[int] = None
    basis: Basis = ()
    affine_certificate_sha256: Optional[str] = None
    candidate_in_prior: bool = False
    prior_in_candidate: bool = False

    @property
    def state_count(self) -> int:
        return 0 if self.dimension is None else 1 << self.dimension


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def file_sha256(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def display_path(path: Path, root: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(root.resolve()))
    except ValueError:
        return str(resolved)


def resolve_path(root: Path, path: Path) -> Path:
    return path if path.is_absolute() else root / path


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return value


def read_matrix(path: Path) -> Matrix:
    try:
        rows = tuple(
            tuple(int(token) for token in line.split())
            for line in path.read_text(encoding="ascii").splitlines()
            if line.strip()
        )
    except (OSError, UnicodeDecodeError, ValueError) as error:
        raise ValueError(f"{path}: invalid ASCII integer matrix") from error
    if len(rows) != ORDER or any(len(row) != ORDER for row in rows):
        raise ValueError(f"{path}: expected exactly 23 rows of 23 entries")
    if any(value not in (-1, 1) for row in rows for value in row):
        raise ValueError(f"{path}: entries must all be -1 or +1")
    return rows


def parse_support(path: Path) -> tuple[Coordinate, ...]:
    coordinates: list[Coordinate] = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), 1
    ):
        stripped = line.partition("#")[0].strip()
        if not stripped:
            continue
        tokens = stripped.split()
        if len(tokens) != 2:
            raise ValueError(
                f"{path}:{line_number}: expected 'ROW COLUMN'"
            )
        try:
            row, column = (int(token) for token in tokens)
        except ValueError as error:
            raise ValueError(
                f"{path}:{line_number}: coordinates must be integers"
            ) from error
        if not (1 <= row <= ORDER and 1 <= column <= ORDER):
            raise ValueError(
                f"{path}:{line_number}: coordinate is outside 1..23"
            )
        coordinates.append((row, column))
    if not coordinates:
        raise ValueError(f"{path}: support is empty")
    if len(coordinates) != len(set(coordinates)):
        raise ValueError(f"{path}: support contains duplicate coordinates")
    return tuple(coordinates)


def parse_json_support(
    value: Any, source: str
) -> tuple[Coordinate, ...]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{source}: support must be a nonempty list")
    coordinates: list[Coordinate] = []
    for index, coordinate in enumerate(value):
        if (
            not isinstance(coordinate, list)
            or len(coordinate) != 2
            or not all(isinstance(item, int) for item in coordinate)
        ):
            raise ValueError(
                f"{source}: bad coordinate at support index {index}"
            )
        row, column = coordinate
        if not (1 <= row <= ORDER and 1 <= column <= ORDER):
            raise ValueError(
                f"{source}: coordinate at index {index} is outside 1..23"
            )
        coordinates.append((row, column))
    if len(coordinates) != len(set(coordinates)):
        raise ValueError(f"{source}: duplicate support coordinates")
    return tuple(coordinates)


def support_bytes(support: Sequence[Coordinate]) -> bytes:
    return "".join(
        f"{row} {column}\n" for row, column in sorted(support)
    ).encode("ascii")


def support_set_certificate(support: Sequence[Coordinate]) -> str:
    return sha256_bytes(SUPPORT_TAG + support_bytes(support))


def evidence_inventory_certificate(
    records: Sequence[dict[str, str]]
) -> str:
    serialized = json.dumps(
        sorted(
            records,
            key=lambda record: (
                record["id"],
                record.get("report_sha256", ""),
            ),
        ),
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return sha256_bytes(EVIDENCE_TAG + serialized)


def canonical_matrix_bytes(matrix: Matrix) -> bytes:
    return "".join(
        " ".join(str(value) for value in row) + "\n"
        for row in matrix
    ).encode("ascii")


def entry_cube_fingerprint(
    matrix: Matrix, support: Sequence[Coordinate]
) -> str:
    """Reproduce the retained engines' full 23x23 ternary fingerprint."""

    free = set(support)
    ternary = bytearray()
    for row in range(1, ORDER + 1):
        for column in range(1, ORDER + 1):
            if (row, column) in free:
                ternary.append(2)
            else:
                ternary.append(
                    1 if matrix[row - 1][column - 1] == 1 else 0
                )
    return sha256_bytes(ENTRY_CUBE_TAG + bytes(ternary))


def apply_mask(
    matrix: Matrix,
    support: Sequence[Coordinate],
    mask: int,
) -> Matrix:
    if mask < 0 or mask >> len(support):
        raise ValueError("entry-flip mask lies outside its support")
    result = [list(row) for row in matrix]
    for index, (row, column) in enumerate(support):
        if (mask >> index) & 1:
            result[row - 1][column - 1] *= -1
    return tuple(tuple(row) for row in result)


def dephased_core_bits(matrix: Matrix) -> int:
    """Return the canonical 22x22 dephased core, with -1 encoded as one."""

    result = 0
    pivot = matrix[0][0]
    for row in range(1, ORDER):
        for column in range(1, ORDER):
            dephased = (
                matrix[row][column]
                * matrix[row][0]
                * matrix[0][column]
                * pivot
            )
            if dephased == -1:
                result |= 1 << (
                    (row - 1) * CORE_ORDER + (column - 1)
                )
    return result


def direction_for_coordinate(coordinate: Coordinate) -> int:
    """Map one one-based full-matrix flip to its dephased GF(2) vector."""

    row, column = coordinate
    if not (1 <= row <= ORDER and 1 <= column <= ORDER):
        raise ValueError(f"coordinate is outside 1..23: {coordinate}")
    row -= 1
    column -= 1
    if row == 0 and column == 0:
        return (1 << CORE_BITS) - 1
    if row == 0:
        return sum(
            1 << (inner_row * CORE_ORDER + column - 1)
            for inner_row in range(CORE_ORDER)
        )
    if column == 0:
        return ((1 << CORE_ORDER) - 1) << (
            (row - 1) * CORE_ORDER
        )
    return 1 << ((row - 1) * CORE_ORDER + column - 1)


def rref_basis(vectors: Iterable[int]) -> Basis:
    """Return the unique descending-pivot reduced basis over GF(2)."""

    rows: dict[int, int] = {}
    for original in vectors:
        if original < 0 or original >> CORE_BITS:
            raise ValueError("GF(2) vector lies outside the 484-bit core")
        value = original
        for pivot in sorted(rows, reverse=True):
            if (value >> pivot) & 1:
                value ^= rows[pivot]
        if value == 0:
            continue
        pivot = value.bit_length() - 1
        for other in list(rows):
            if (rows[other] >> pivot) & 1:
                rows[other] ^= value
        rows[pivot] = value
    return tuple(rows[pivot] for pivot in sorted(rows, reverse=True))


def reduce_vector(value: int, basis: Basis) -> int:
    for vector in basis:
        pivot = vector.bit_length() - 1
        if (value >> pivot) & 1:
            value ^= vector
    return value


def affine_certificate(base: int, basis: Basis) -> str:
    reduced_base = reduce_vector(base, basis)
    payload = (
        AFFINE_TAG
        + len(basis).to_bytes(2, "big")
        + reduced_base.to_bytes(VECTOR_BYTES, "big")
        + b"".join(
            vector.to_bytes(VECTOR_BYTES, "big") for vector in basis
        )
    )
    return sha256_bytes(payload)


def remainder_certificate(remainder: int) -> str:
    return sha256_bytes(
        REMAINDER_TAG + remainder.to_bytes(VECTOR_BYTES, "big")
    )


def state_union_certificate(states: Iterable[int]) -> str:
    return sha256_bytes(
        STATE_UNION_TAG
        + b"".join(
            state.to_bytes(VECTOR_BYTES, "big")
            for state in sorted(states)
        )
    )


def solve_combination(
    generators: Sequence[int], target: int
) -> Optional[int]:
    """Return a generator-coefficient mask for target, or None."""

    rows: dict[int, tuple[int, int]] = {}
    for index, original in enumerate(generators):
        value = original
        coefficients = 1 << index
        for pivot in sorted(rows, reverse=True):
            if (value >> pivot) & 1:
                value ^= rows[pivot][0]
                coefficients ^= rows[pivot][1]
        if value == 0:
            continue
        pivot = value.bit_length() - 1
        for other in list(rows):
            if (rows[other][0] >> pivot) & 1:
                rows[other] = (
                    rows[other][0] ^ value,
                    rows[other][1] ^ coefficients,
                )
        rows[pivot] = (value, coefficients)

    value = target
    coefficients = 0
    for pivot in sorted(rows, reverse=True):
        if (value >> pivot) & 1:
            value ^= rows[pivot][0]
            coefficients ^= rows[pivot][1]
    return coefficients if value == 0 else None


def direction_intersection(first: Basis, second: Basis) -> Basis:
    """Return a canonical basis for span(first) intersect span(second)."""

    quotient_rows: dict[int, tuple[int, int]] = {}
    intersection_vectors: list[int] = []
    for index, original in enumerate(first):
        value = reduce_vector(original, second)
        coefficients = 1 << index
        for pivot in sorted(quotient_rows, reverse=True):
            if (value >> pivot) & 1:
                value ^= quotient_rows[pivot][0]
                coefficients ^= quotient_rows[pivot][1]
        if value == 0:
            intersection = 0
            for inner, vector in enumerate(first):
                if (coefficients >> inner) & 1:
                    intersection ^= vector
            if intersection:
                intersection_vectors.append(intersection)
            continue
        pivot = value.bit_length() - 1
        for other in list(quotient_rows):
            if (quotient_rows[other][0] >> pivot) & 1:
                quotient_rows[other] = (
                    quotient_rows[other][0] ^ value,
                    quotient_rows[other][1] ^ coefficients,
                )
        quotient_rows[pivot] = (value, coefficients)
    return rref_basis(intersection_vectors)


def intersect_affine(first: Cube, second: Cube) -> Intersection:
    union_basis = rref_basis((*first.basis, *second.basis))
    delta = first.base ^ second.base
    remainder = reduce_vector(delta, union_basis)
    if remainder:
        return Intersection(
            intersects=False,
            union_rank=len(union_basis),
            delta_remainder=remainder,
        )

    coefficients = solve_combination(
        (*first.basis, *second.basis), delta
    )
    if coefficients is None:
        raise ArithmeticError(
            "zero reduced remainder had no generator representation"
        )
    point = first.base
    for index, vector in enumerate(first.basis):
        if (coefficients >> index) & 1:
            point ^= vector

    intersection_basis = direction_intersection(
        first.basis, second.basis
    )
    expected_dimension = (
        first.rank + second.rank - len(union_basis)
    )
    if len(intersection_basis) != expected_dimension:
        raise ArithmeticError("affine intersection dimension mismatch")
    if (
        reduce_vector(point ^ first.base, first.basis)
        or reduce_vector(point ^ second.base, second.basis)
    ):
        raise ArithmeticError("computed intersection point is not shared")
    certificate = affine_certificate(point, intersection_basis)
    return Intersection(
        intersects=True,
        union_rank=len(union_basis),
        delta_remainder=0,
        dimension=expected_dimension,
        point=point,
        basis=intersection_basis,
        affine_certificate_sha256=certificate,
        candidate_in_prior=expected_dimension == first.rank,
        prior_in_candidate=expected_dimension == second.rank,
    )


def make_cube(
    *,
    family: str,
    label: str,
    start_path: Path,
    support: Sequence[Coordinate],
) -> Cube:
    matrix = read_matrix(start_path)
    canonical_support = tuple(support)
    basis = rref_basis(
        direction_for_coordinate(coordinate)
        for coordinate in canonical_support
    )
    base = dephased_core_bits(matrix)
    return Cube(
        family=family,
        label=label,
        start_path=start_path,
        start_sha256=file_sha256(start_path),
        support=canonical_support,
        support_set_certificate_sha256=support_set_certificate(
            canonical_support
        ),
        base=base,
        basis=basis,
        affine_certificate_sha256=affine_certificate(base, basis),
    )


def cube_payload(cube: Cube, root: Path) -> dict[str, Any]:
    return {
        "affine_certificate_sha256":
            cube.affine_certificate_sha256,
        "dephased_direction_rank": cube.rank,
        "dimension": cube.dimension,
        "family": cube.family,
        "full_dephased_rank": cube.rank == cube.dimension,
        "label": cube.label,
        "start": display_path(cube.start_path, root),
        "start_sha256": cube.start_sha256,
        "support": [
            [row, column] for row, column in cube.support
        ],
        "support_set_certificate_sha256":
            cube.support_set_certificate_sha256,
        "total_distinct_dephased_states": str(1 << cube.rank),
    }


def require_value(
    actual: Any, expected: Any, source: Path, field: str
) -> None:
    if actual != expected or type(actual) is not type(expected):
        raise ValueError(
            f"{source}: {field}={actual!r}, expected {expected!r}"
        )


def require_hash(
    actual: str, expected: Any, source: Path, field: str
) -> None:
    if (
        not isinstance(expected, str)
        or len(expected) != 64
        or any(character not in "0123456789abcdef" for character in expected)
        or actual != expected
    ):
        raise ValueError(
            f"{source}: {field} does not bind to the retained file"
        )


def require_reported_path(
    report: dict[str, Any],
    field: str,
    actual: Path,
    root: Path,
    source: Path,
) -> None:
    value = report.get(field)
    if not isinstance(value, str) or not value:
        raise ValueError(f"{source}: {field} is not a path")
    reported = Path(value)
    actual_resolved = actual.resolve()
    if reported.is_absolute() and reported.exists():
        if reported.resolve() != actual_resolved:
            raise ValueError(
                f"{source}: {field} resolves to a different retained file"
            )
        return
    try:
        suffix = actual_resolved.relative_to(root.resolve()).parts
    except ValueError:
        suffix = actual_resolved.parts
    if len(reported.parts) < len(suffix):
        raise ValueError(f"{source}: {field} has the wrong path suffix")
    if tuple(reported.parts[-len(suffix):]) != tuple(suffix):
        raise ValueError(f"{source}: {field} has the wrong path suffix")


def indexed_records(
    value: Any, key: str, source: Path
) -> dict[str, dict[str, Any]]:
    if not isinstance(value, list):
        raise ValueError(f"{source}: expected a record list")
    indexed: dict[str, dict[str, Any]] = {}
    for index, record in enumerate(value):
        if not isinstance(record, dict):
            raise ValueError(f"{source}: record {index} is not an object")
        identifier = record.get(key)
        if not isinstance(identifier, str) or not identifier:
            raise ValueError(
                f"{source}: record {index} has no string {key}"
            )
        if identifier in indexed:
            raise ValueError(
                f"{source}: duplicate {key} {identifier!r}"
            )
        indexed[identifier] = record
    return indexed


def require_engine_support_dimension(
    support: Sequence[Coordinate], source: Path
) -> None:
    if len(support) != ENGINE_DIMENSION:
        raise ValueError(
            f"{source}: support has {len(support)} coordinates, "
            f"expected {ENGINE_DIMENSION}"
        )


def validate_engine_report(
    *,
    root: Path,
    report_path: Path,
    start_path: Path,
    support_path: Path,
    expected_support: Optional[Sequence[Coordinate]] = None,
) -> tuple[
    dict[str, Any], Matrix, tuple[Coordinate, ...], str, str, str
]:
    """Bind one completed engine report to its exact start/support files."""

    report = read_json(report_path)
    require_value(report.get("engine"), ENGINE, report_path, "engine")
    require_value(report.get("complete"), True, report_path, "complete")
    require_value(
        report.get("all_assignments_bound_checked"),
        True,
        report_path,
        "all_assignments_bound_checked",
    )
    require_value(
        report.get("dimension"),
        ENGINE_DIMENSION,
        report_path,
        "dimension",
    )
    require_value(
        report.get("assignments"),
        ENGINE_ASSIGNMENTS,
        report_path,
        "assignments",
    )
    require_value(
        report.get("coordinate_indexing"),
        "one_based",
        report_path,
        "coordinate_indexing",
    )
    require_value(
        report.get("support_source"),
        "coordinate_file",
        report_path,
        "support_source",
    )
    if not start_path.is_file():
        raise ValueError(f"{report_path}: retained start matrix is missing")
    if not support_path.is_file():
        raise ValueError(f"{report_path}: retained support file is missing")
    require_reported_path(
        report, "start", start_path, root, report_path
    )
    require_reported_path(
        report, "coordinate_file", support_path, root, report_path
    )

    matrix = read_matrix(start_path)
    support = parse_support(support_path)
    require_engine_support_dimension(support, support_path)
    report_support = parse_json_support(
        report.get("support"), f"{report_path}:support"
    )
    if support != report_support:
        raise ValueError(
            f"{report_path}: report support order/content differs from file"
        )
    if (
        expected_support is not None
        and support != tuple(expected_support)
    ):
        raise ValueError(
            f"{report_path}: support differs from its manifest/partition"
        )

    start_raw_sha256 = file_sha256(start_path)
    start_parsed_sha256 = sha256_bytes(canonical_matrix_bytes(matrix))
    support_raw_sha256 = file_sha256(support_path)
    require_hash(
        start_raw_sha256,
        report.get("start_raw_sha256"),
        report_path,
        "start_raw_sha256",
    )
    require_hash(
        start_parsed_sha256,
        report.get("start_parsed_matrix_sha256"),
        report_path,
        "start_parsed_matrix_sha256",
    )
    require_hash(
        support_raw_sha256,
        report.get("coordinate_file_raw_sha256"),
        report_path,
        "coordinate_file_raw_sha256",
    )
    return (
        report,
        matrix,
        support,
        start_raw_sha256,
        support_raw_sha256,
        entry_cube_fingerprint(matrix, support),
    )


def load_batch_cubes(
    root: Path, relative_directory: Path, expected_count: int
) -> tuple[list[Cube], dict[str, Any]]:
    directory = resolve_path(root, relative_directory)
    manifest_path = directory / "manifest.json"
    aggregate_path = directory / "aggregate-report.json"
    manifest = read_json(manifest_path)
    aggregate = read_json(aggregate_path)
    manifest_runs = indexed_records(
        manifest.get("runs"), "id", manifest_path
    )
    aggregate_runs = indexed_records(
        aggregate.get("runs"), "id", aggregate_path
    )
    run_count = len(manifest_runs)
    if run_count != expected_count:
        raise ValueError(
            f"{manifest_path}: found {run_count} runs, "
            f"expected {expected_count}"
        )
    if set(aggregate_runs) != set(manifest_runs):
        raise ValueError(
            f"{aggregate_path}: manifest and aggregate run IDs differ"
        )
    for source, payload in (
        (manifest_path, manifest),
        (aggregate_path, aggregate),
    ):
        require_value(
            payload.get("engine"), ENGINE, source, "engine"
        )
        require_value(
            payload.get("dimension"),
            ENGINE_DIMENSION,
            source,
            "dimension",
        )
    require_value(
        aggregate.get("complete"), True, aggregate_path, "complete"
    )
    for field in ("completed_runs", "planned_runs"):
        require_value(
            aggregate.get(field), run_count, aggregate_path, field
        )
    require_value(
        manifest.get("planned_runs"),
        run_count,
        manifest_path,
        "planned_runs",
    )
    require_value(
        aggregate.get("total_assignments"),
        run_count * ENGINE_ASSIGNMENTS,
        aggregate_path,
        "total_assignments",
    )
    for source, payload in (
        (manifest_path, manifest),
        (aggregate_path, aggregate),
    ):
        require_value(
            payload.get("support_fingerprints_unique"),
            run_count,
            source,
            "support_fingerprints_unique",
        )
        require_value(
            payload.get("support_set_hashes_unique"),
            run_count,
            source,
            "support_set_hashes_unique",
        )

    cubes: list[Cube] = []
    evidence_records: list[dict[str, str]] = []
    seen_directories: set[Path] = set()
    seen_fingerprints: set[str] = set()
    for label, record in sorted(manifest_runs.items()):
        run_directory = (directory / label).resolve()
        if run_directory in seen_directories:
            raise ValueError(
                f"{manifest_path}: duplicate run directory {run_directory}"
            )
        seen_directories.add(run_directory)
        report_path = directory / label / "report.json"
        support = parse_json_support(
            record.get("support"), f"{manifest_path}:{label}"
        )
        start_value = record.get("start_path")
        if not isinstance(start_value, str) or not start_value:
            raise ValueError(
                f"{manifest_path}:{label}: start_path is not a path"
            )
        start_path = resolve_path(root, Path(start_value))
        support_path = directory / label / "support.coords.txt"
        (
            _report,
            _matrix,
            bound_support,
            start_raw_sha256,
            support_raw_sha256,
            fingerprint,
        ) = validate_engine_report(
            root=root,
            report_path=report_path,
            start_path=start_path,
            support_path=support_path,
            expected_support=support,
        )
        require_hash(
            start_raw_sha256,
            record.get("start_raw_sha256"),
            manifest_path,
            f"{label}.start_raw_sha256",
        )
        require_hash(
            support_raw_sha256,
            record.get("support_sha256"),
            manifest_path,
            f"{label}.support_sha256",
        )
        require_hash(
            fingerprint,
            record.get("fingerprint_sha256"),
            manifest_path,
            f"{label}.fingerprint_sha256",
        )
        aggregate_record = aggregate_runs[label]
        require_hash(
            support_raw_sha256,
            aggregate_record.get("support_sha256"),
            aggregate_path,
            f"{label}.support_sha256",
        )
        require_hash(
            fingerprint,
            aggregate_record.get("fingerprint_sha256"),
            aggregate_path,
            f"{label}.fingerprint_sha256",
        )
        require_value(
            aggregate_record.get("assignments"),
            ENGINE_ASSIGNMENTS,
            aggregate_path,
            f"{label}.assignments",
        )
        if fingerprint in seen_fingerprints:
            raise ValueError(
                f"{manifest_path}: duplicate affine fingerprint {fingerprint}"
            )
        seen_fingerprints.add(fingerprint)
        evidence_records.append(
            {
                "id": label,
                "report_sha256": file_sha256(report_path),
                "start_raw_sha256": start_raw_sha256,
                "support_raw_sha256": support_raw_sha256,
            }
        )
        cubes.append(
            make_cube(
                family="batch200",
                label=label,
                start_path=start_path,
                support=bound_support,
            )
        )
    return cubes, {
        "aggregate": display_path(aggregate_path, root),
        "aggregate_sha256": file_sha256(aggregate_path),
        "completed_all_bound_reports_verified": len(cubes),
        "engine_assignments_per_report": ENGINE_ASSIGNMENTS,
        "entry_cube_fingerprints_verified": len(seen_fingerprints),
        "evidence_inventory_certificate_sha256":
            evidence_inventory_certificate(evidence_records),
        "manifest": display_path(manifest_path, root),
        "manifest_sha256": file_sha256(manifest_path),
        "start_and_support_file_hashes_verified": len(cubes),
    }


def load_lnps_cubes(
    root: Path, relative_directory: Path, expected_count: int
) -> tuple[list[Cube], dict[str, Any]]:
    directory = resolve_path(root, relative_directory)
    aggregate_path = directory / "aggregate-report.json"
    aggregate = read_json(aggregate_path)
    require_value(
        aggregate.get("engine"), ENGINE, aggregate_path, "engine"
    )
    require_value(
        aggregate.get("dimension"),
        ENGINE_DIMENSION,
        aggregate_path,
        "dimension",
    )
    require_value(
        aggregate.get("completed_cubes"),
        expected_count,
        aggregate_path,
        "completed_cubes",
    )
    require_value(
        aggregate.get("total_assignment_visits"),
        expected_count * ENGINE_ASSIGNMENTS,
        aggregate_path,
        "total_assignment_visits",
    )

    report_paths = sorted(
        directory.glob("generation-*/cubes/*/report.json")
    )
    all_report_paths = {
        path.resolve() for path in directory.glob("**/report.json")
    }
    if all_report_paths != {path.resolve() for path in report_paths}:
        raise ValueError(
            f"{directory}: unexpected report path outside generation cubes"
        )
    summary_paths = {
        path.resolve()
        for path in directory.glob("**/lnps-summary.json")
    }
    expected_summary_paths = {
        (path.parent / "lnps-summary.json").resolve()
        for path in report_paths
    }
    if summary_paths != expected_summary_paths:
        raise ValueError(
            f"{directory}: LNPS report/summary inventory differs"
        )

    cubes: list[Cube] = []
    evidence_records: list[dict[str, str]] = []
    seen_ids: set[str] = set()
    seen_directories: set[Path] = set()
    seen_fingerprints: set[str] = set()
    generation_counts: Counter[int] = Counter()
    for report_path in report_paths:
        cube_directory = report_path.parent.resolve()
        if cube_directory in seen_directories:
            raise ValueError(
                f"{directory}: duplicate cube directory {cube_directory}"
            )
        seen_directories.add(cube_directory)
        generation_name = report_path.parents[2].name
        try:
            generation = int(generation_name.removeprefix("generation-"))
        except ValueError as error:
            raise ValueError(
                f"{report_path}: invalid generation directory"
            ) from error
        if generation_name != f"generation-{generation:02d}":
            raise ValueError(
                f"{report_path}: noncanonical generation directory"
            )
        generation_counts[generation] += 1

        summary_path = report_path.parent / "lnps-summary.json"
        summary = read_json(summary_path)
        identifier = summary.get("id")
        if (
            not isinstance(identifier, str)
            or identifier != report_path.parent.name
        ):
            raise ValueError(
                f"{summary_path}: id does not match its directory"
            )
        if identifier in seen_ids:
            raise ValueError(
                f"{directory}: duplicate LNPS cube ID {identifier!r}"
            )
        seen_ids.add(identifier)

        local_start = report_path.parent / "start.matrix.txt"
        support_path = report_path.parent / "support.coords.txt"
        (
            _report,
            _matrix,
            support,
            start_raw_sha256,
            support_raw_sha256,
            fingerprint,
        ) = validate_engine_report(
            root=root,
            report_path=report_path,
            start_path=local_start,
            support_path=support_path,
        )
        require_value(
            summary.get("assignments"),
            ENGINE_ASSIGNMENTS,
            summary_path,
            "assignments",
        )
        require_hash(
            support_raw_sha256,
            summary.get("support_sha256"),
            summary_path,
            "support_sha256",
        )
        require_hash(
            fingerprint,
            summary.get("fingerprint_sha256"),
            summary_path,
            "fingerprint_sha256",
        )
        if fingerprint in seen_fingerprints:
            raise ValueError(
                f"{directory}: duplicate LNPS fingerprint {fingerprint}"
            )
        seen_fingerprints.add(fingerprint)
        evidence_records.append(
            {
                "id": identifier,
                "report_sha256": file_sha256(report_path),
                "start_raw_sha256": start_raw_sha256,
                "summary_sha256": file_sha256(summary_path),
                "support_raw_sha256": support_raw_sha256,
            }
        )
        cubes.append(
            make_cube(
                family="lnps",
                label=display_path(report_path.parent, root),
                start_path=local_start,
                support=support,
            )
        )
    if len(cubes) != expected_count:
        raise ValueError(
            f"{directory}: found {len(cubes)} completed cubes, "
            f"expected {expected_count}"
        )

    generation_reports = aggregate.get("generation_reports")
    if not isinstance(generation_reports, list):
        raise ValueError(
            f"{aggregate_path}: generation_reports must be a list"
        )
    aggregate_generations: dict[int, dict[str, Any]] = {}
    for record in generation_reports:
        if not isinstance(record, dict):
            raise ValueError(
                f"{aggregate_path}: generation record is not an object"
            )
        generation = record.get("generation")
        if type(generation) is not int or generation < 0:
            raise ValueError(
                f"{aggregate_path}: invalid generation identifier"
            )
        if generation in aggregate_generations:
            raise ValueError(
                f"{aggregate_path}: duplicate generation {generation}"
            )
        aggregate_generations[generation] = record
    if set(aggregate_generations) != set(generation_counts):
        raise ValueError(
            f"{aggregate_path}: aggregate/directory generations differ"
        )
    for generation, count in generation_counts.items():
        require_value(
            aggregate_generations[generation].get("completed_cubes"),
            count,
            aggregate_path,
            f"generation-{generation:02d}.completed_cubes",
        )
    require_value(
        sum(
            int(record["completed_cubes"])
            for record in aggregate_generations.values()
        ),
        len(cubes),
        aggregate_path,
        "generation completed-cube sum",
    )
    require_value(
        aggregate.get("affine_fingerprints_evaluated"),
        len(seen_fingerprints),
        aggregate_path,
        "affine_fingerprints_evaluated",
    )
    return cubes, {
        "aggregate": display_path(aggregate_path, root),
        "aggregate_sha256": file_sha256(aggregate_path),
        "completed_all_bound_reports_verified": len(cubes),
        "engine_assignments_per_report": ENGINE_ASSIGNMENTS,
        "entry_cube_fingerprints_verified": len(seen_fingerprints),
        "evidence_inventory_certificate_sha256":
            evidence_inventory_certificate(evidence_records),
        "lnps_summaries_verified": len(cubes),
        "start_and_support_file_hashes_verified": len(cubes),
    }


def load_deep32_cubes(
    root: Path,
    relative_directories: Sequence[Path],
    expected_count: int,
) -> tuple[list[Cube], list[dict[str, Any]]]:
    if len(relative_directories) != expected_count:
        raise ValueError(
            f"received {len(relative_directories)} deep-32 directories, "
            f"expected {expected_count}"
        )
    resolved_directories = [
        resolve_path(root, directory).resolve()
        for directory in relative_directories
    ]
    if len(resolved_directories) != len(set(resolved_directories)):
        raise ValueError("duplicate deep-32 directory in inventory")
    labels = [directory.name for directory in resolved_directories]
    if len(labels) != len(set(labels)):
        raise ValueError("duplicate deep-32 directory label in inventory")

    cubes: list[Cube] = []
    sources: list[dict[str, Any]] = []
    for directory in resolved_directories:
        aggregate_path = directory / "aggregate-report.json"
        aggregate = read_json(aggregate_path)
        evidence_records: list[dict[str, str]] = []
        partition = aggregate.get("leaf_partition", {})
        for field, expected in (
            ("complete", True),
            ("completed_leaves", DEEP_LEAF_COUNT),
            ("planned_leaves", DEEP_LEAF_COUNT),
            ("dimension", DEEP_DIMENSION),
            ("engine_leaf_dimension", ENGINE_DIMENSION),
            ("total_assignments", DEEP_ASSIGNMENTS),
            ("leaf_affine_fingerprints_unique", DEEP_LEAF_COUNT),
        ):
            require_value(
                aggregate.get(field), expected, aggregate_path, field
            )
        if not isinstance(partition, dict):
            raise ValueError(
                f"{aggregate_path}: leaf_partition is not an object"
            )
        for field, expected in (
            ("disjoint", True),
            ("fixed_outer_bits", DEEP_DIMENSION - ENGINE_DIMENSION),
            ("leaf_assignments", ENGINE_ASSIGNMENTS),
            ("leaf_count", DEEP_LEAF_COUNT),
            ("total_assignments", DEEP_ASSIGNMENTS),
            ("union_complete", True),
        ):
            require_value(
                partition.get(field),
                expected,
                aggregate_path,
                f"leaf_partition.{field}",
            )
        support = parse_json_support(
            aggregate.get("full_support"),
            f"{aggregate_path}:full_support",
        )
        if len(support) != DEEP_DIMENSION:
            raise ValueError(
                f"{aggregate_path}: full_support is not 32 coordinates"
            )

        global_ties = indexed_records(
            aggregate.get("global_ties"),
            "global_mask_decimal",
            aggregate_path,
        )
        if "0" not in global_ties:
            raise ValueError(
                f"{aggregate_path}: expected one global mask-0 record"
            )
        mask_zero = global_ties["0"]
        artifact = mask_zero.get("artifact")
        if not isinstance(artifact, str) or not artifact:
            raise ValueError(
                f"{aggregate_path}: mask-0 artifact is not a path"
            )
        base = resolve_path(root, Path(artifact))
        expected_base = (
            directory / "global-frontier-ties" / "mask-0.matrix.txt"
        )
        if base.resolve() != expected_base.resolve():
            raise ValueError(
                f"{aggregate_path}: mask-0 artifact points outside partition"
            )
        if not base.is_file():
            raise ValueError(f"{directory}: mask-0 base matrix is missing")
        require_hash(
            file_sha256(base),
            mask_zero.get("raw_sha256"),
            aggregate_path,
            "global_ties.mask-0.raw_sha256",
        )
        evidence_records.append(
            {
                "id": "global-mask-0",
                "start_raw_sha256": file_sha256(base),
            }
        )
        base_matrix = read_matrix(base)
        require_hash(
            entry_cube_fingerprint(base_matrix, support),
            aggregate.get("full_affine_fingerprint_sha256"),
            aggregate_path,
            "full_affine_fingerprint_sha256",
        )

        aggregate_leaves = indexed_records(
            aggregate.get("leaf_reports"), "leaf_id", aggregate_path
        )
        if len(aggregate_leaves) != DEEP_LEAF_COUNT:
            raise ValueError(
                f"{aggregate_path}: expected 32 unique leaf reports"
            )
        report_paths = sorted(
            (directory / "leaves").glob("*/report.json")
        )
        all_report_paths = {
            path.resolve()
            for path in directory.glob("leaves/**/report.json")
        }
        if (
            len(report_paths) != DEEP_LEAF_COUNT
            or all_report_paths
            != {path.resolve() for path in report_paths}
        ):
            raise ValueError(
                f"{directory}: retained leaf report inventory is not 32"
            )
        report_ids = {path.parent.name for path in report_paths}
        if report_ids != set(aggregate_leaves):
            raise ValueError(
                f"{aggregate_path}: aggregate/directory leaf IDs differ"
            )

        inner_support = support[:ENGINE_DIMENSION]
        outer_support = support[ENGINE_DIMENSION:]
        seen_outer_masks: set[int] = set()
        seen_leaf_directories: set[Path] = set()
        seen_leaf_fingerprints: set[str] = set()
        for report_path in report_paths:
            leaf_directory = report_path.parent.resolve()
            if leaf_directory in seen_leaf_directories:
                raise ValueError(
                    f"{directory}: duplicate leaf directory {leaf_directory}"
                )
            seen_leaf_directories.add(leaf_directory)
            leaf_id = report_path.parent.name
            summary_path = report_path.parent / "partition-summary.json"
            summary = read_json(summary_path)
            if summary != aggregate_leaves[leaf_id]:
                raise ValueError(
                    f"{summary_path}: summary differs from aggregate record"
                )
            outer_value = summary.get("outer_mask_decimal")
            try:
                outer_mask = int(outer_value)
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"{summary_path}: invalid outer_mask_decimal"
                ) from error
            if (
                not isinstance(outer_value, str)
                or not 0 <= outer_mask < DEEP_LEAF_COUNT
            ):
                raise ValueError(
                    f"{summary_path}: outer mask is outside 0..31"
                )
            expected_leaf_id = (
                f"leaf-{outer_mask:02d}-outer-{outer_mask:05b}"
            )
            if (
                leaf_id != expected_leaf_id
                or summary.get("outer_mask_binary") != f"{outer_mask:05b}"
            ):
                raise ValueError(
                    f"{summary_path}: leaf ID/binary mask mismatch"
                )
            if outer_mask in seen_outer_masks:
                raise ValueError(
                    f"{aggregate_path}: duplicate outer mask {outer_mask}"
                )
            seen_outer_masks.add(outer_mask)
            reroot_value = summary.get("reroot_xor_mask_decimal")
            try:
                reroot_mask = int(reroot_value)
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"{summary_path}: invalid reroot mask"
                ) from error
            if (
                not isinstance(reroot_value, str)
                or not 0 <= reroot_mask < ENGINE_ASSIGNMENTS
            ):
                raise ValueError(
                    f"{summary_path}: reroot mask is outside 27 bits"
                )
            require_value(
                summary.get("rerooted"),
                bool(reroot_mask),
                summary_path,
                "rerooted",
            )
            require_value(
                summary.get("assignments"),
                ENGINE_ASSIGNMENTS,
                summary_path,
                "assignments",
            )

            local_start = report_path.parent / "start.matrix.txt"
            support_path = report_path.parent / "support.coords.txt"
            (
                _report,
                leaf_matrix,
                _bound_support,
                start_raw_sha256,
                support_raw_sha256,
                leaf_fingerprint,
            ) = validate_engine_report(
                root=root,
                report_path=report_path,
                start_path=local_start,
                support_path=support_path,
                expected_support=inner_support,
            )
            expected_leaf_matrix = apply_mask(
                apply_mask(base_matrix, outer_support, outer_mask),
                inner_support,
                reroot_mask,
            )
            if leaf_matrix != expected_leaf_matrix:
                raise ValueError(
                    f"{report_path}: leaf start is not its claimed slice"
                )
            require_hash(
                leaf_fingerprint,
                summary.get("leaf_affine_fingerprint_sha256"),
                summary_path,
                "leaf_affine_fingerprint_sha256",
            )
            if leaf_fingerprint in seen_leaf_fingerprints:
                raise ValueError(
                    f"{aggregate_path}: duplicate leaf fingerprint"
                )
            seen_leaf_fingerprints.add(leaf_fingerprint)
            evidence_records.append(
                {
                    "id": leaf_id,
                    "report_sha256": file_sha256(report_path),
                    "start_raw_sha256": start_raw_sha256,
                    "summary_sha256": file_sha256(summary_path),
                    "support_raw_sha256": support_raw_sha256,
                }
            )
        if seen_outer_masks != set(range(DEEP_LEAF_COUNT)):
            raise ValueError(
                f"{aggregate_path}: leaf masks do not cover 0..31"
            )

        cubes.append(
            make_cube(
                family="deep32",
                label=directory.name,
                start_path=base,
                support=support,
            )
        )
        sources.append(
            {
                "aggregate": display_path(aggregate_path, root),
                "aggregate_sha256": file_sha256(aggregate_path),
                "completed_all_bound_leaf_reports_verified":
                    len(seen_outer_masks),
                "disjoint_leaf_masks_verified": len(seen_outer_masks),
                "entry_cube_fingerprints_verified":
                    len(seen_leaf_fingerprints) + 1,
                "evidence_inventory_certificate_sha256":
                    evidence_inventory_certificate(evidence_records),
                "global_mask_zero_raw_sha256": file_sha256(base),
                "start_and_support_file_hashes_verified":
                    len(seen_outer_masks),
            }
        )
    return cubes, sources


def inventory_certificate(cubes: Sequence[Cube]) -> str:
    records = [
        {
            "affine_certificate_sha256":
                cube.affine_certificate_sha256,
            "family": cube.family,
            "label": cube.label,
            "rank": cube.rank,
            "start_sha256": cube.start_sha256,
            "support_set_certificate_sha256":
                cube.support_set_certificate_sha256,
        }
        for cube in sorted(
            cubes,
            key=lambda item: (
                FAMILY_ORDER[item.family], item.label
            ),
        )
    ]
    serialized = json.dumps(
        records, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return sha256_bytes(INVENTORY_TAG + serialized)


def require_unique_inventory(cubes: Sequence[Cube]) -> None:
    identifiers: set[tuple[str, str]] = set()
    affine_certificates: dict[str, tuple[str, str]] = {}
    for cube in cubes:
        identifier = (cube.family, cube.label)
        if identifier in identifiers:
            raise ValueError(
                "duplicate prior-cube inventory ID "
                f"{cube.family}:{cube.label}"
            )
        identifiers.add(identifier)
        prior_identifier = affine_certificates.get(
            cube.affine_certificate_sha256
        )
        if prior_identifier is not None:
            raise ValueError(
                "duplicate dephased affine cube in prior inventory: "
                f"{prior_identifier[0]}:{prior_identifier[1]} and "
                f"{cube.family}:{cube.label}"
            )
        affine_certificates[
            cube.affine_certificate_sha256
        ] = identifier


def enumerate_intersection_states(
    intersection: Intersection
) -> set[int]:
    if not intersection.intersects or intersection.point is None:
        return set()
    states: set[int] = set()
    for mask in range(1 << len(intersection.basis)):
        state = intersection.point
        for index, vector in enumerate(intersection.basis):
            if (mask >> index) & 1:
                state ^= vector
        states.add(state)
    if len(states) != intersection.state_count:
        raise ArithmeticError("intersection enumeration was not injective")
    return states


def dimension_distribution(
    comparisons: Sequence[dict[str, Any]]
) -> dict[str, int]:
    counts = Counter(
        str(record["intersection_dimension"])
        for record in comparisons
        if record["intersects"]
    )
    return {
        key: counts[key]
        for key in sorted(counts, key=lambda value: int(value))
    }


def comparison_record(
    prior: Cube, intersection: Intersection
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "candidate_in_prior": intersection.candidate_in_prior,
        "delta_union_remainder_sha256":
            remainder_certificate(intersection.delta_remainder),
        "delta_union_remainder_zero":
            intersection.delta_remainder == 0,
        "family": prior.family,
        "intersects": intersection.intersects,
        "label": prior.label,
        "prior_affine_certificate_sha256":
            prior.affine_certificate_sha256,
        "prior_dephased_direction_rank": prior.rank,
        "prior_in_candidate": intersection.prior_in_candidate,
        "union_direction_rank": intersection.union_rank,
    }
    if intersection.intersects:
        record.update(
            {
                "intersection_affine_certificate_sha256":
                    intersection.affine_certificate_sha256,
                "intersection_dimension": intersection.dimension,
                "intersection_state_count":
                    str(intersection.state_count),
            }
        )
    else:
        record.update(
            {
                "intersection_affine_certificate_sha256": None,
                "intersection_dimension": None,
                "intersection_state_count": "0",
            }
        )
    return record


def comparison_certificate(
    records: Sequence[dict[str, Any]]
) -> str:
    serialized = json.dumps(
        records, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return sha256_bytes(
        b"maxdet-dephased-comparison-records-v1\0" + serialized
    )


def audit(
    *,
    root: Path,
    start_path: Path,
    support: Sequence[Coordinate],
    label: str,
    batch_directory: Path,
    lnps_directory: Path,
    deep32_directories: Sequence[Path],
    expected_batch_count: int,
    expected_lnps_count: int,
    expected_deep32_count: int,
    maximum_enumerated_overlap_states: int,
    include_empty_comparisons: bool,
) -> dict[str, Any]:
    candidate = make_cube(
        family="candidate",
        label=label,
        start_path=start_path,
        support=support,
    )
    batch, batch_sources = load_batch_cubes(
        root, batch_directory, expected_batch_count
    )
    lnps, lnps_sources = load_lnps_cubes(
        root, lnps_directory, expected_lnps_count
    )
    deep32, deep32_sources = load_deep32_cubes(
        root, deep32_directories, expected_deep32_count
    )
    priors = sorted(
        [*batch, *lnps, *deep32],
        key=lambda cube: (FAMILY_ORDER[cube.family], cube.label),
    )
    require_unique_inventory(priors)

    records: list[dict[str, Any]] = []
    intersections: list[tuple[Cube, Intersection]] = []
    total_pairwise_overlap = 0
    for prior in priors:
        intersection = intersect_affine(candidate, prior)
        records.append(comparison_record(prior, intersection))
        if intersection.intersects:
            intersections.append((prior, intersection))
            total_pairwise_overlap += intersection.state_count

    if total_pairwise_overlap > maximum_enumerated_overlap_states:
        raise ValueError(
            "exact overlap union would enumerate "
            f"{total_pairwise_overlap} pairwise states, above "
            f"--maximum-enumerated-overlap-states="
            f"{maximum_enumerated_overlap_states}"
        )

    union_by_family = {
        family: set() for family in FAMILY_ORDER
    }
    union_all: set[int] = set()
    for prior, intersection in intersections:
        states = enumerate_intersection_states(intersection)
        union_by_family[prior.family].update(states)
        union_all.update(states)

    summaries: dict[str, Any] = {}
    for family in FAMILY_ORDER:
        family_records = [
            record for record in records
            if record["family"] == family
        ]
        nonempty = [
            record for record in family_records
            if record["intersects"]
        ]
        summaries[family] = {
            "comparison_certificate_sha256":
                comparison_certificate(family_records),
            "empty_intersection_count":
                len(family_records) - len(nonempty),
            "intersection_dimension_distribution":
                dimension_distribution(nonempty),
            "nonempty_intersection_count": len(nonempty),
            "pairwise_intersection_state_sum": str(
                sum(
                    int(record["intersection_state_count"])
                    for record in nonempty
                )
            ),
            "prior_cube_count": len(family_records),
            "unique_overlap_state_count":
                str(len(union_by_family[family])),
            "unique_overlap_state_union_sha256":
                state_union_certificate(union_by_family[family]),
        }

    containment_candidate = [
        {
            "family": record["family"],
            "label": record["label"],
        }
        for record in records
        if record["candidate_in_prior"]
    ]
    containment_prior = [
        {
            "family": record["family"],
            "label": record["label"],
        }
        for record in records
        if record["prior_in_candidate"]
    ]
    prior_rank_distribution = Counter(
        (cube.family, cube.rank) for cube in priors
    )
    total_candidate_states = 1 << candidate.rank
    unique_overlap = len(union_all)

    return {
        "schema_version": SCHEMA_VERSION,
        "method": METHOD,
        "coordinate_indexing": "one_based",
        "dephased_space": {
            "bit_count": CORE_BITS,
            "core_order": CORE_ORDER,
            "entry_flip_map": {
                "interior": "one core bit",
                "top_left": "all core bits",
                "first_row": "one core column",
                "first_column": "one core row",
            },
            "field": "GF(2)",
            "vector_serialization":
                "61-byte big-endian unsigned integer",
        },
        "candidate": cube_payload(candidate, root),
        "prior_corpus": {
            "count": len(priors),
            "counts_by_family": {
                family: sum(
                    cube.family == family for cube in priors
                )
                for family in FAMILY_ORDER
            },
            "dephased_affine_certificates_unique": len(priors),
            "inventory_certificate_sha256":
                inventory_certificate(priors),
            "rank_distribution": [
                {
                    "count": count,
                    "family": family,
                    "rank": rank,
                }
                for (family, rank), count in sorted(
                    prior_rank_distribution.items(),
                    key=lambda item: (
                        FAMILY_ORDER[item[0][0]], item[0][1]
                    ),
                )
            ],
            "sources": {
                "batch200": batch_sources,
                "lnps": lnps_sources,
                "deep32": deep32_sources,
            },
        },
        "comparisons": {
            "all_comparisons_certificate_sha256":
                comparison_certificate(records),
            "by_family": summaries,
            "candidate_contained_in_prior": containment_candidate,
            "empty_intersection_count":
                len(records) - len(intersections),
            "nonempty": [
                record for record in records if record["intersects"]
            ],
            "nonempty_intersection_count": len(intersections),
            "prior_contained_in_candidate": containment_prior,
            "records": records if include_empty_comparisons else None,
            "records_include_empty": include_empty_comparisons,
            "total_prior_cubes": len(records),
        },
        "coverage": {
            "candidate_new_dephased_state_count":
                str(total_candidate_states - unique_overlap),
            "candidate_total_dephased_state_count":
                str(total_candidate_states),
            "pairwise_intersection_state_sum":
                str(total_pairwise_overlap),
            "unique_prior_overlap_dephased_state_count":
                str(unique_overlap),
            "unique_prior_overlap_state_union_sha256":
                state_union_certificate(union_all),
            "warning":
                "Pairwise sums double-count states; unique union "
                "counts do not.",
        },
        "provenance": {
            "tool": display_path(Path(__file__), root),
            "tool_sha256": file_sha256(Path(__file__)),
        },
        "claim_boundary": [
            (
                "This is exact affine-set algebra after row/column sign "
                "dephasing, not determinant verification."
            ),
            (
                "A zero overlap applies only to the pinned local prior-cube "
                "inventory recorded in this report."
            ),
            (
                "The affine certificate does not quotient row/column "
                "permutations or transpose."
            ),
        ],
    }


def toy_cube(
    label: str, base: int, vectors: Sequence[int]
) -> Cube:
    basis = rref_basis(vectors)
    return Cube(
        family="toy",
        label=label,
        start_path=Path(label),
        start_sha256="",
        support=tuple((1, index + 1) for index in range(len(vectors))),
        support_set_certificate_sha256="",
        base=base,
        basis=basis,
        affine_certificate_sha256=affine_certificate(base, basis),
    )


def brute_states(base: int, vectors: Sequence[int]) -> set[int]:
    result: set[int] = set()
    for mask in range(1 << len(vectors)):
        value = base
        for index, vector in enumerate(vectors):
            if (mask >> index) & 1:
                value ^= vector
        result.add(value)
    return result


def run_self_test() -> dict[str, Any]:
    checks: list[str] = []

    all_positive: Matrix = tuple(
        tuple(1 for _ in range(ORDER)) for _ in range(ORDER)
    )
    base = dephased_core_bits(all_positive)
    if base != 0:
        raise AssertionError("all-positive dephased core was not zero")
    for row in range(1, ORDER + 1):
        for column in range(1, ORDER + 1):
            coordinate = (row, column)
            flipped = apply_mask(all_positive, (coordinate,), 1)
            observed = dephased_core_bits(flipped) ^ base
            if observed != direction_for_coordinate(coordinate):
                raise AssertionError(
                    f"boundary direction mismatch at {coordinate}"
                )
    checks.append("all 529 entry-flip dephasing directions")

    fingerprint_support = ((1, 1), (3, 4), (23, 23))
    fingerprint = entry_cube_fingerprint(
        all_positive, fingerprint_support
    )
    if entry_cube_fingerprint(
        apply_mask(all_positive, fingerprint_support, 0b101),
        tuple(reversed(fingerprint_support)),
    ) != fingerprint:
        raise AssertionError(
            "entry-cube fingerprint changed within its free support"
        )
    if entry_cube_fingerprint(
        apply_mask(all_positive, ((2, 2),), 1),
        fingerprint_support,
    ) == fingerprint:
        raise AssertionError(
            "entry-cube fingerprint ignored a fixed-entry change"
        )
    if apply_mask(
        apply_mask(all_positive, fingerprint_support, 0b111),
        fingerprint_support,
        0b111,
    ) != all_positive:
        raise AssertionError("entry-flip masks did not round-trip")
    checks.append("entry-cube fingerprint and mask semantics")

    duplicate_rejected = False
    try:
        indexed_records(
            [{"id": "same"}, {"id": "same"}],
            "id",
            Path("<self-test>"),
        )
    except ValueError:
        duplicate_rejected = True
    if not duplicate_rejected:
        raise AssertionError("duplicate evidence IDs were accepted")
    malformed_hash_rejected = False
    try:
        require_hash(
            "0" * 64,
            "not-a-sha256",
            Path("<self-test>"),
            "hash",
        )
    except ValueError:
        malformed_hash_rejected = True
    if not malformed_hash_rejected:
        raise AssertionError("malformed evidence hash was accepted")
    wrong_dimension_rejected = False
    try:
        require_engine_support_dimension(
            tuple((2, column) for column in range(1, ORDER + 1)),
            Path("<self-test>"),
        )
    except ValueError:
        wrong_dimension_rejected = True
    if not wrong_dimension_rejected:
        raise AssertionError("wrong engine-support dimension was accepted")
    duplicate_affine_rejected = False
    try:
        require_unique_inventory(
            (
                toy_cube("duplicate-affine-a", 0, (1, 2)),
                toy_cube("duplicate-affine-b", 3, (1, 2)),
            )
        )
    except ValueError:
        duplicate_affine_rejected = True
    if not duplicate_affine_rejected:
        raise AssertionError("duplicate dephased affine cubes were accepted")
    checks.append(
        "evidence hash, support-dimension, and duplicate guards"
    )

    first_row_directions = [
        direction_for_coordinate((1, column))
        for column in range(1, ORDER + 1)
    ]
    if len(rref_basis(first_row_directions)) != ORDER - 1:
        raise AssertionError("whole-row gauge relation has wrong rank")
    if (
        len(rref_basis(direction_for_coordinate((2, column))
                       for column in range(1, ORDER + 1)))
        != ORDER - 1
    ):
        raise AssertionError("nonfirst whole-row gauge relation is wrong")
    checks.append("row-sign gauge dependencies")

    first = toy_cube("first", 0, (1, 2))
    second = toy_cube("second", 4, (2, 4))
    shared = intersect_affine(first, second)
    if (
        not shared.intersects
        or shared.dimension != 1
        or enumerate_intersection_states(shared) != {0, 2}
    ):
        raise AssertionError("toy affine intersection is wrong")
    disjoint = intersect_affine(
        first, toy_cube("disjoint", 4, (1,))
    )
    if disjoint.intersects:
        raise AssertionError("toy disjoint affine spaces intersected")
    contained = intersect_affine(
        toy_cube("small", 0, (2,)), first
    )
    if not contained.candidate_in_prior:
        raise AssertionError("toy containment was not detected")
    checks.append("intersection and containment")

    if affine_certificate(0, rref_basis((1, 2))) != affine_certificate(
        3, rref_basis((2, 1, 3))
    ):
        raise AssertionError("affine certificate is not canonical")
    checks.append("canonical affine certificate")

    randomizer = random.Random(23)
    for trial in range(128):
        width = 7
        first_vectors = [
            randomizer.randrange(1 << width)
            for _ in range(randomizer.randrange(0, 6))
        ]
        second_vectors = [
            randomizer.randrange(1 << width)
            for _ in range(randomizer.randrange(0, 6))
        ]
        first_base = randomizer.randrange(1 << width)
        second_base = randomizer.randrange(1 << width)
        first_toy = toy_cube(
            f"random-first-{trial}", first_base, first_vectors
        )
        second_toy = toy_cube(
            f"random-second-{trial}", second_base, second_vectors
        )
        exact = (
            brute_states(first_base, first_vectors)
            & brute_states(second_base, second_vectors)
        )
        result = intersect_affine(first_toy, second_toy)
        observed = enumerate_intersection_states(result)
        if exact != observed:
            raise AssertionError(
                f"random affine differential failed at trial {trial}"
            )
    checks.append("128 brute-force differential intersections")

    return {
        "schema_version": SCHEMA_VERSION,
        "method": METHOD,
        "provenance": {
            "tool": str(Path(__file__).resolve()),
            "tool_sha256": file_sha256(Path(__file__)),
        },
        "self_test": {
            "checks": checks,
            "passed": True,
        },
    }


def positive_count(value: str, option: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            f"{option} must be an integer"
        ) from error
    if parsed <= 0:
        raise argparse.ArgumentTypeError(
            f"{option} must be positive"
        )
    return parsed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="run deterministic algebra and dephasing checks",
    )
    parser.add_argument("--start", type=Path)
    parser.add_argument(
        "--support",
        type=Path,
        help="one-based ROW COLUMN coordinate file",
    )
    parser.add_argument("--label", default="input-cube")
    parser.add_argument(
        "--repository-root", type=Path, default=REPOSITORY_ROOT
    )
    parser.add_argument(
        "--batch-directory",
        type=Path,
        default=DEFAULT_BATCH_DIRECTORY,
    )
    parser.add_argument(
        "--lnps-directory",
        type=Path,
        default=DEFAULT_LNPS_DIRECTORY,
    )
    parser.add_argument(
        "--deep32-directory",
        type=Path,
        action="append",
        dest="deep32_directories",
        help=(
            "completed deep-32 directory; repeat twice "
            "(defaults to the retained H0 and H2 partitions)"
        ),
    )
    parser.add_argument(
        "--expected-batch-count",
        type=lambda value: positive_count(
            value, "--expected-batch-count"
        ),
        default=200,
    )
    parser.add_argument(
        "--expected-lnps-count",
        type=lambda value: positive_count(
            value, "--expected-lnps-count"
        ),
        default=110,
    )
    parser.add_argument(
        "--expected-deep32-count",
        type=lambda value: positive_count(
            value, "--expected-deep32-count"
        ),
        default=2,
    )
    parser.add_argument(
        "--maximum-enumerated-overlap-states",
        type=lambda value: positive_count(
            value, "--maximum-enumerated-overlap-states"
        ),
        default=1_000_000,
        help=(
            "fail instead of materializing a larger exact overlap union"
        ),
    )
    parser.add_argument(
        "--include-empty-comparisons",
        action="store_true",
        help="include all pairwise records, not only nonempty intersections",
    )
    arguments = parser.parse_args()

    try:
        if arguments.self_test:
            if arguments.start is not None or arguments.support is not None:
                parser.error(
                    "--self-test cannot be combined with --start/--support"
                )
            result = run_self_test()
        else:
            if arguments.start is None or arguments.support is None:
                parser.error(
                    "audit mode requires both --start and --support"
                )
            root = arguments.repository_root.resolve()
            start = resolve_path(root, arguments.start)
            support_path = resolve_path(root, arguments.support)
            deep32_directories = (
                tuple(arguments.deep32_directories)
                if arguments.deep32_directories
                else DEFAULT_DEEP32_DIRECTORIES
            )
            result = audit(
                root=root,
                start_path=start,
                support=parse_support(support_path),
                label=arguments.label,
                batch_directory=arguments.batch_directory,
                lnps_directory=arguments.lnps_directory,
                deep32_directories=deep32_directories,
                expected_batch_count=arguments.expected_batch_count,
                expected_lnps_count=arguments.expected_lnps_count,
                expected_deep32_count=arguments.expected_deep32_count,
                maximum_enumerated_overlap_states=(
                    arguments.maximum_enumerated_overlap_states
                ),
                include_empty_comparisons=(
                    arguments.include_empty_comparisons
                ),
            )
    except (OSError, ValueError, ArithmeticError) as error:
        parser.error(str(error))

    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
