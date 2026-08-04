#!/usr/bin/env python3
"""Encode an arena matrix for the MaxDetBounty23 Solidity contract."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from maxdet.contract import Contract, load_contract
from maxdet.receipt import VerifiedMatrix, verify_matrix


PUBLISHED_FRONTIER = 2_779_447_296_000_000
DETERMINANT_LATTICE_STEP = 1 << 22
MINIMUM_WINNING_DETERMINANT = PUBLISHED_FRONTIER + DETERMINANT_LATTICE_STEP
ENCODING = "uint32[23]-lsb-column-zero-plus-one-v1"


def encode_rows(matrix: list[list[int]]) -> list[int]:
    """Encode each row with bit j set exactly when column j is +1."""

    if len(matrix) != 23:
        raise ValueError(f"expected 23 rows, found {len(matrix)}")
    for row_index, row in enumerate(matrix):
        if len(row) != 23:
            raise ValueError(
                f"row {row_index + 1}: expected 23 entries, found {len(row)}"
            )
        if any(entry not in {-1, 1} for entry in row):
            raise ValueError(f"row {row_index + 1}: entries must be -1 or 1")

    return [
        sum(1 << column for column, entry in enumerate(row) if entry == 1)
        for row in matrix
    ]


def verify_and_encode(
    matrix_path: Path,
    contract: Contract,
) -> tuple[VerifiedMatrix, list[int]]:
    """Read once, then verify and encode the same parsed matrix."""

    verified = verify_matrix(matrix_path, contract)
    return verified, encode_rows(verified.matrix)


def cast_array(rows: list[int]) -> str:
    return "[" + ",".join(str(row) for row in rows) + "]"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Verify and encode a 23 x 23 matrix for MaxDetBounty23."
    )
    parser.add_argument(
        "matrix",
        nargs="?",
        type=Path,
        default=Path("candidate/matrix.txt"),
    )
    parser.add_argument(
        "--cast-array",
        action="store_true",
        help="print only the uint32[23] value accepted by cast",
    )
    parser.add_argument(
        "--require-winning",
        action="store_true",
        help="fail unless the exact score qualifies for the bounty",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    contract = load_contract(REPOSITORY_ROOT / "challenge.json")
    verified, rows = verify_and_encode(args.matrix, contract)
    qualifies = verified.abs_determinant >= MINIMUM_WINNING_DETERMINANT

    if args.require_winning and not qualifies:
        raise SystemExit(
            "matrix does not qualify: "
            f"{verified.abs_determinant} < {MINIMUM_WINNING_DETERMINANT}"
        )

    if args.cast_array:
        print(cast_array(rows))
        return 0

    print(
        json.dumps(
            {
                "encoding": ENCODING,
                "rows": rows,
                "cast_array": cast_array(rows),
                "absolute_determinant": str(verified.abs_determinant),
                "minimum_winning_determinant": str(
                    MINIMUM_WINNING_DETERMINANT
                ),
                "qualifies": qualifies,
                "matrix_sha256": verified.receipt["matrix"]["raw_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
