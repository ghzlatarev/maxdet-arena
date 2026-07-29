#!/usr/bin/env python3
"""Emit radius-two core pairs hinged on exact neutral entry flips.

For a maximal order-22 core, a neutral entry is one whose exact cofactor is
zero, so flipping it preserves the determinant.  This helper emits every
unordered two-entry set containing at least one neutral entry.  Rebordering
those cores searches the complete one-entry shell around every immediate
neutral neighbor without enumerating the full radius-two sphere.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

from order22_border_structure import (
    bareiss_determinant,
    exact_adjugate,
    read_sign_matrix,
)


ORDER = 22


def parse_labeled_path(text: str) -> tuple[str, Path]:
    if "=" not in text:
        raise argparse.ArgumentTypeError("core must be LABEL=PATH")
    label, raw_path = text.split("=", 1)
    if not label or not raw_path:
        raise argparse.ArgumentTypeError("core must be LABEL=PATH")
    return label, Path(raw_path)


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_bytes(payload)
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--core",
        action="append",
        required=True,
        type=parse_labeled_path,
        metavar="LABEL=PATH",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    arguments = parser.parse_args()
    labels = [label for label, _ in arguments.core]
    if len(labels) != len(set(labels)):
        parser.error("core labels must be unique")

    output_dir = arguments.output_dir.resolve()
    manifest_path = output_dir / "neutral-hinge-pairs.json"
    if manifest_path.exists():
        raise FileExistsError(manifest_path)
    all_cells = [
        (row, column)
        for row in range(ORDER)
        for column in range(ORDER)
    ]
    records: list[dict[str, object]] = []
    for label, path in arguments.core:
        resolved = path.resolve()
        payload, core = read_sign_matrix(resolved, ORDER)
        determinant = bareiss_determinant(core)
        if not determinant:
            raise ArithmeticError(f"{resolved}: singular core")
        adjugate = exact_adjugate(core, determinant)
        neutral = [
            (row, column)
            for row, column in all_cells
            if adjugate[column][row] == 0
        ]
        neutral_set = set(neutral)
        pairs = [
            (first, second)
            for first_index, first in enumerate(all_cells)
            for second in all_cells[first_index + 1 :]
            if first in neutral_set or second in neutral_set
        ]
        lines = "".join(
            f"{first[0] + 1}\t{first[1] + 1}\t"
            f"{second[0] + 1}\t{second[1] + 1}\n"
            for first, second in pairs
        ).encode("ascii")
        tsv_path = output_dir / f"{label}.neutral-hinge-radius2.tsv"
        if tsv_path.exists():
            raise FileExistsError(tsv_path)
        atomic_write(tsv_path, lines)
        records.append(
            {
                "label": label,
                "core_path": str(resolved),
                "core_raw_sha256": hashlib.sha256(payload).hexdigest(),
                "determinant": str(determinant),
                "neutral_entries_1_based": [
                    [row + 1, column + 1] for row, column in neutral
                ],
                "neutral_entry_count": len(neutral),
                "pair_count": len(pairs),
                "radius2_tsv_path": str(tsv_path),
                "radius2_tsv_sha256": hashlib.sha256(lines).hexdigest(),
            }
        )
    manifest = {
        "engine": "order22-neutral-hinge-pairs-v1",
        "arithmetic": "exact-integer-and-rational-only",
        "definition": (
            "Every unordered two-entry core flip containing at least one "
            "entry whose exact cofactor is zero."
        ),
        "records": records,
    }
    atomic_write(
        manifest_path,
        (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode(
            "utf-8"
        ),
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
