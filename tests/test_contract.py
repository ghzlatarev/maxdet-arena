from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from maxdet.contract import load_contract, parse_matrix_bytes
from maxdet.errors import ContractError, MatrixFormatError

ROOT = Path(__file__).resolve().parent.parent


class ContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = load_contract(ROOT / "challenge.json")

    def valid_matrix_bytes(self) -> bytes:
        return (
            "\n".join(
                " ".join("-1" if row == column else "1" for column in range(23))
                for row in range(23)
            )
            + "\n"
        ).encode("ascii")

    def test_contract_is_pinned(self) -> None:
        self.assertEqual(self.contract.challenge_id, "maxdet-23-v1")
        self.assertEqual(self.contract.order, 23)
        self.assertEqual(
            self.contract.modular_primes,
            (998244353, 1000000007, 1000000009),
        )

    def test_accepts_ascii_whitespace_and_extra_trailing_newline(self) -> None:
        raw = b"\n\t" + self.valid_matrix_bytes().replace(b" ", b"\t") + b"\n "
        matrix = parse_matrix_bytes(raw, self.contract)
        self.assertEqual(len(matrix), 23)
        self.assertTrue(all(len(row) == 23 for row in matrix))

    def test_rejects_unicode_whitespace(self) -> None:
        raw = self.valid_matrix_bytes().replace(b" ", "\u00a0".encode("utf-8"), 1)
        with self.assertRaisesRegex(MatrixFormatError, "ASCII whitespace"):
            parse_matrix_bytes(raw, self.contract)

    def test_rejects_non_domain_entry(self) -> None:
        raw = self.valid_matrix_bytes().replace(b"-1", b"0", 1)
        with self.assertRaises(MatrixFormatError):
            parse_matrix_bytes(raw, self.contract)

    def test_rejects_extra_column(self) -> None:
        raw = self.valid_matrix_bytes().replace(b"\n", b" 1\n", 1)
        with self.assertRaisesRegex(MatrixFormatError, "expected 23 entries"):
            parse_matrix_bytes(raw, self.contract)

    def test_rejects_missing_row(self) -> None:
        raw = b"\n".join(self.valid_matrix_bytes().splitlines()[:-1])
        with self.assertRaisesRegex(MatrixFormatError, "expected 23 rows"):
            parse_matrix_bytes(raw, self.contract)

    def test_rejects_oversized_file_before_parsing(self) -> None:
        raw = self.valid_matrix_bytes() + b" " * 9000
        with self.assertRaisesRegex(MatrixFormatError, "byte limit"):
            parse_matrix_bytes(raw, self.contract)

    def test_contract_hash_uses_exact_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            altered = Path(temporary) / "challenge.json"
            altered.write_bytes((ROOT / "challenge.json").read_bytes() + b" ")
            altered_contract = load_contract(altered)
            self.assertNotEqual(altered_contract.sha256, self.contract.sha256)

    def test_semantic_contract_mutation_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            altered = Path(temporary) / "challenge.json"
            data = json.loads((ROOT / "challenge.json").read_text())
            data["objective"]["direction"] = "minimize"
            altered.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaisesRegex(ContractError, "immutable"):
                load_contract(altered)

    def test_duplicate_contract_key_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            altered = Path(temporary) / "challenge.json"
            raw = (ROOT / "challenge.json").read_bytes()
            altered.write_bytes(raw.replace(b'"schema_version": 1,', b'"schema_version": 1, "schema_version": 1,', 1))
            with self.assertRaisesRegex(ContractError, "duplicate"):
                load_contract(altered)


if __name__ == "__main__":
    unittest.main()
