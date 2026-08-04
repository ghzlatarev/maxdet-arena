from __future__ import annotations

import random
import unittest
from pathlib import Path

from maxdet.contract import load_contract, load_matrix
from maxdet.exact import bareiss_determinant
from tools.encode_onchain_matrix import (
    DETERMINANT_LATTICE_STEP,
    MINIMUM_WINNING_DETERMINANT,
    PUBLISHED_FRONTIER,
    cast_array,
    encode_rows,
    verify_and_encode,
)


ROOT = Path(__file__).resolve().parent.parent


def decode_rows(rows: list[int]) -> list[list[int]]:
    return [
        [1 if encoded & (1 << column) else -1 for column in range(23)]
        for encoded in rows
    ]


def reduced_core(matrix: list[list[int]]) -> list[list[int]]:
    corner = matrix[0][0]
    return [
        [
            (
                1
                - matrix[row][column]
                * matrix[0][column]
                * matrix[row][0]
                * corner
            )
            // 2
            for column in range(1, 23)
        ]
        for row in range(1, 23)
    ]


class OnchainEncodingTests(unittest.TestCase):
    def setUp(self) -> None:
        contract = load_contract(ROOT / "challenge.json")
        _, self.reference = load_matrix(
            ROOT / "references" / "orrick-et-al-2003" / "matrix.txt",
            contract,
        )

    def test_reference_encoding_matches_solidity_vector(self) -> None:
        self.assertEqual(
            encode_rows(self.reference),
            [
                8_380_805,
                8_362_086,
                8_357_403,
                1_679_244,
                6_692_756,
                3_343_234,
                5_014_914,
                3_443_193,
                4_919_801,
                2_949_180,
                5_406_812,
                2_816_154,
                5_601_562,
                1_980_417,
                6_412_289,
                478_415,
                447_279,
                240_951,
                373_463,
                6_826_831,
                7_353_519,
                5_786_295,
                3_691_863,
            ],
        )

    def test_encoding_round_trip(self) -> None:
        rows = encode_rows(self.reference)
        self.assertEqual(decode_rows(rows), self.reference)
        self.assertTrue(cast_array(rows).startswith("[8380805,"))

    def test_encoder_rejects_invalid_direct_inputs(self) -> None:
        with self.assertRaisesRegex(ValueError, "expected 23 rows"):
            encode_rows([[1] * 23 for _ in range(22)])

        wrong_width = [[1] * 23 for _ in range(23)]
        wrong_width[8] = [1] * 22
        with self.assertRaisesRegex(ValueError, "row 9: expected 23 entries"):
            encode_rows(wrong_width)

        wrong_entry = [[1] * 23 for _ in range(23)]
        wrong_entry[4][12] = 0
        with self.assertRaisesRegex(ValueError, "row 5: entries must be -1 or 1"):
            encode_rows(wrong_entry)

    def test_verify_and_encode_reads_candidate_once(self) -> None:
        contract = load_contract(ROOT / "challenge.json")
        reference_bytes = (
            ROOT / "references" / "orrick-et-al-2003" / "matrix.txt"
        ).read_bytes()
        replacement_bytes = (
            (" ".join(["1"] * 23) + "\n") * 23
        ).encode("ascii")

        class ChangingPath:
            def __init__(self) -> None:
                self.read_count = 0

            def read_bytes(self) -> bytes:
                self.read_count += 1
                if self.read_count == 1:
                    return reference_bytes
                return replacement_bytes

        changing_path = ChangingPath()
        verified, rows = verify_and_encode(changing_path, contract)  # type: ignore[arg-type]

        self.assertEqual(changing_path.read_count, 1)
        self.assertEqual(rows, encode_rows(self.reference))
        self.assertEqual(verified.abs_determinant, PUBLISHED_FRONTIER)

    def test_reduced_determinant_identity_on_random_matrices(self) -> None:
        random_source = random.Random(23)
        for _ in range(32):
            matrix = [
                [random_source.choice((-1, 1)) for _ in range(23)]
                for _ in range(23)
            ]
            determinant = abs(bareiss_determinant(matrix))
            reduced = abs(bareiss_determinant(reduced_core(matrix)))
            self.assertEqual(determinant, reduced * (1 << 22))

    def test_published_reference_does_not_qualify(self) -> None:
        determinant = abs(bareiss_determinant(self.reference))
        self.assertEqual(determinant, PUBLISHED_FRONTIER)
        self.assertEqual(
            MINIMUM_WINNING_DETERMINANT,
            PUBLISHED_FRONTIER + DETERMINANT_LATTICE_STEP,
        )
        self.assertLess(determinant, MINIMUM_WINNING_DETERMINANT)


if __name__ == "__main__":
    unittest.main()
