from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from maxdet.contract import load_contract
from maxdet.exact import bareiss_determinant, matrix_text
from maxdet.receipt import canonical_json_bytes, verify_matrix

ROOT = Path(__file__).resolve().parent.parent


class ReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = load_contract(ROOT / "challenge.json")

    def test_genesis_score_and_checks(self) -> None:
        verified = verify_matrix(ROOT / "records/genesis/matrix.txt", self.contract)
        self.assertEqual(verified.abs_determinant, 21 * (2**22))
        self.assertEqual(verified.receipt["checks"]["gram_identity"], "passed")
        self.assertEqual(
            len(verified.receipt["checks"]["modular_determinants"]),
            3,
        )
        self.assertEqual(
            verified.receipt["checks"]["crt_unique_reconstruction"],
            "passed",
        )
        for witness in verified.receipt["checks"]["modular_determinants"]:
            self.assertIn("gram_determinant_mod_prime", witness)

    def test_receipt_is_deterministic(self) -> None:
        path = ROOT / "records/genesis/matrix.txt"
        first = canonical_json_bytes(verify_matrix(path, self.contract).receipt)
        second = canonical_json_bytes(verify_matrix(path, self.contract).receipt)
        self.assertEqual(first, second)

    def test_sign_equivalent_matrices_share_normalized_hash(self) -> None:
        original_path = ROOT / "records/genesis/matrix.txt"
        original = verify_matrix(original_path, self.contract)
        changed = [row[:] for row in original.matrix]
        changed[4] = [-value for value in changed[4]]
        for row in changed:
            row[7] *= -1

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "matrix.txt"
            path.write_text(matrix_text(changed), encoding="ascii")
            equivalent = verify_matrix(path, self.contract)

        self.assertEqual(original.abs_determinant, equivalent.abs_determinant)
        self.assertEqual(
            original.receipt["matrix"]["sign_normalized_sha256"],
            equivalent.receipt["matrix"]["sign_normalized_sha256"],
        )

    def test_exact_score_is_not_floating_point(self) -> None:
        verified = verify_matrix(ROOT / "records/genesis/matrix.txt", self.contract)
        score = verified.receipt["score"]
        self.assertIsInstance(score["absolute_determinant"], str)
        self.assertEqual(
            int(score["determinant_squared"]),
            bareiss_determinant(verified.matrix) ** 2,
        )


if __name__ == "__main__":
    unittest.main()
