#!/usr/bin/env python3
"""Create a feasible CP-SAT hint inside a known signed-Gram support orbit.

Exact-support no-good cuts reject the listed factor bytes but not their
signed-row-automorphic images.  This tool walks the audited Gram automorphism
generators until it finds an image that contains requested symmetry-fixing
masks and differs from every excluded support.  The result is deliberately
H-equivalent to the input factor: it is a solver anchor, not a novelty claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import random
import time
from typing import Any

from order22_gram_factor_cpsat import (
    ORDER,
    InputError,
    atomic_write,
    bareiss_determinant,
    gram,
    matrix_bytes,
    normalized_column_masks,
    parse_mask,
    read_factor_list,
    read_json,
    read_matrix,
    reconstruct_factor,
    sha256,
)


def sha256_bytes(contents: bytes) -> str:
    return hashlib.sha256(contents).hexdigest()


def parse_generators(
    report: dict[str, Any],
    target_gram: list[list[int]],
) -> list[tuple[list[int], list[int]]]:
    group = report.get("group")
    raw_generators = (
        group.get("generators") if isinstance(group, dict) else None
    )
    if not isinstance(raw_generators, list) or not raw_generators:
        raise InputError("orbit report has no signed-Gram generators")
    result = []
    for generator_index, record in enumerate(raw_generators):
        permutation = (
            record.get("source_to_target_row_permutation")
            if isinstance(record, dict)
            else None
        )
        signs = record.get("source_row_signs") if isinstance(record, dict) else None
        if (
            not isinstance(permutation, list)
            or not isinstance(signs, list)
            or len(permutation) != ORDER
            or len(signs) != ORDER
            or set(permutation) != set(range(ORDER))
            or any(sign not in (-1, 1) for sign in signs)
        ):
            raise InputError(f"malformed generator {generator_index}")
        for first in range(ORDER):
            for second in range(ORDER):
                if (
                    target_gram[permutation[first]][permutation[second]]
                    != signs[first] * signs[second] * target_gram[first][second]
                ):
                    raise InputError(
                        f"generator {generator_index} fails signed Gram equality"
                    )
        result.append((permutation, signs))
    return result


def act_mask(
    mask: int, generator: tuple[list[int], list[int]]
) -> int:
    permutation, signs = generator
    source = [1 if mask >> row & 1 else -1 for row in range(ORDER)]
    target = [0] * ORDER
    for row in range(ORDER):
        target[permutation[row]] = signs[row] * source[row]
    if target[0] < 0:
        target = [-value for value in target]
    return sum((value == 1) << row for row, value in enumerate(target))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit-report", required=True, type=Path)
    parser.add_argument("--factor", required=True, type=Path)
    parser.add_argument(
        "--exclude-factor-list",
        action="append",
        default=[],
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--require-mask",
        action="append",
        default=[],
        required=True,
        type=parse_mask,
    )
    parser.add_argument("--seed", type=int, default=22)
    parser.add_argument("--max-steps", type=int, default=1_000_000)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    arguments = parser.parse_args()
    if arguments.max_steps <= 0:
        parser.error("--max-steps must be positive")
    if arguments.output.resolve() == arguments.metadata.resolve():
        parser.error("--output and --metadata must differ")
    if arguments.output.exists() or arguments.metadata.exists():
        parser.error("refusing to overwrite output or metadata")

    started = time.monotonic()
    report_path = arguments.orbit_report.expanduser().resolve()
    factor_path = arguments.factor.expanduser().resolve()
    report = read_json(report_path)
    if (
        report.get("engine") != "order22-signed-gram-shell-orbits-v1"
        or report.get("order") != ORDER
    ):
        raise InputError("orbit report has an unexpected schema")
    source = read_matrix(factor_path)
    target_gram = gram(source)
    source_determinant = bareiss_determinant(source)
    if source_determinant == 0:
        raise InputError("source factor is singular")
    report_factor = report.get("inputs", {}).get("factor", {})
    report_factor_path = Path(str(report_factor.get("path", "")))
    if (
        not report_factor_path.is_file()
        or report_factor.get("sha256") != sha256(report_factor_path)
        or gram(read_matrix(report_factor_path)) != target_gram
    ):
        raise InputError("orbit report is bound to a different target Gram")

    raw_orbits = report.get("shell_partition", {}).get("orbits")
    if not isinstance(raw_orbits, list) or not raw_orbits:
        raise InputError("orbit report has no shell partition")
    shell = {
        mask
        for record in raw_orbits
        for mask in (
            record.get("masks", []) if isinstance(record, dict) else []
        )
    }
    orbit_by_mask = {
        mask: orbit_index
        for orbit_index, record in enumerate(raw_orbits)
        for mask in record["masks"]
    }
    source_support = frozenset(normalized_column_masks(source))
    if not source_support <= shell:
        raise InputError("source factor leaves the audited shell")
    for mask in arguments.require_mask:
        if mask not in shell:
            raise InputError(f"required mask leaves the shell: {mask}")
        orbit_index = orbit_by_mask[mask]
        if not source_support.intersection(raw_orbits[orbit_index]["masks"]):
            raise InputError(
                f"source support has zero incidence in required orbit {orbit_index}"
            )

    exclusion_records = []
    excluded_supports = set()
    list_records = []
    for raw_list in arguments.exclude_factor_list:
        list_path = raw_list.expanduser().resolve()
        factors = read_factor_list(list_path)
        list_records.append(
            {
                "path": str(list_path),
                "sha256": sha256(list_path),
                "factor_count": len(factors),
            }
        )
        for path in factors:
            matrix = read_matrix(path)
            if gram(matrix) != target_gram:
                raise InputError(f"excluded factor has another Gram: {path}")
            support = frozenset(normalized_column_masks(matrix))
            excluded_supports.add(support)
            exclusion_records.append(
                {"path": str(path), "sha256": sha256(path)}
            )

    generators = parse_generators(report, target_gram)
    required = set(arguments.require_mask)
    generator_word: list[int] = []
    support = source_support
    generator_use_counts = [0] * len(generators)
    rng = random.Random(arguments.seed)
    found_step = None
    for step in range(1, arguments.max_steps + 1):
        generator_index = rng.randrange(len(generators))
        generator_word.append(generator_index)
        generator_use_counts[generator_index] += 1
        support = frozenset(
            act_mask(mask, generators[generator_index]) for mask in support
        )
        if len(support) != ORDER or not support <= shell:
            raise ArithmeticError("automorphism image is not a shell support")
        if required <= support and support not in excluded_supports:
            found_step = step
            break
    if found_step is None:
        raise RuntimeError(
            f"no admissible support found in {arguments.max_steps} steps"
        )

    selected_masks = sorted(support)
    factor = reconstruct_factor(selected_masks)
    if gram(factor) != target_gram:
        raise ArithmeticError("reconstructed hint fails exact Gram equality")
    determinant = bareiss_determinant(factor)
    if abs(determinant) != abs(source_determinant):
        raise ArithmeticError("reconstructed hint determinant changed")
    contents = matrix_bytes(factor)
    output_path = arguments.output.expanduser().resolve()
    metadata_path = arguments.metadata.expanduser().resolve()
    atomic_write(output_path, contents)
    metadata = {
        "engine": "order22-gram-support-orbit-hint-v1",
        "order": ORDER,
        "claim": (
            "The output is an exact signed-Gram-automorphic image of the input "
            "factor and a feasible CP-SAT anchor outside the listed exact supports."
        ),
        "claim_boundary": (
            "The output is deliberately H-equivalent to the source and is not "
            "a new factor-class claim."
        ),
        "inputs": {
            "orbit_report": {
                "path": str(report_path),
                "sha256": sha256(report_path),
            },
            "factor": {
                "path": str(factor_path),
                "sha256": sha256(factor_path),
            },
            "excluded_factor_lists": list_records,
            "excluded_factors": exclusion_records,
            "required_masks": sorted(required),
            "seed": arguments.seed,
            "max_steps": arguments.max_steps,
        },
        "walk": {
            "found_step": found_step,
            "generator_word": generator_word,
            "generator_use_counts": generator_use_counts,
        },
        "output": {
            "path": str(output_path),
            "sha256": sha256_bytes(contents),
            "determinant": determinant,
            "selected_masks": selected_masks,
            "exact_excluded_support_match": False,
        },
        "elapsed_seconds": time.monotonic() - started,
    }
    metadata_bytes = (
        json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    atomic_write(metadata_path, metadata_bytes)
    print(
        f"found_step={found_step} generators={len(generators)} "
        f"determinant={determinant} elapsed={metadata['elapsed_seconds']:.3f}s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
