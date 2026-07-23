from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path

from maxdet.contract import load_contract
from maxdet.errors import SubmissionError
from maxdet.receipt import canonical_json_bytes, verify_matrix
from maxdet.submission import discover_submission_directories, verify_submission

ROOT = Path(__file__).resolve().parent.parent


class SubmissionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = load_contract(ROOT / "challenge.json")

    def make_submission(self, root: Path) -> Path:
        directory = root / "entry"
        directory.mkdir()
        matrix = (ROOT / "records/genesis/matrix.txt").read_bytes()
        (directory / "matrix.txt").write_bytes(matrix)
        metadata = {
            "schema_version": 1,
            "submission_id": "test-001",
            "handle": "researcher",
            "method": "unit test",
            "parent_receipt_sha256": None,
            "artifact_license": "CC0-1.0",
        }
        (directory / "metadata.json").write_text(
            json.dumps(metadata, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        receipt = verify_matrix(directory / "matrix.txt", self.contract).receipt
        (directory / "receipt.json").write_bytes(canonical_json_bytes(receipt))
        return directory

    def test_complete_submission(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = self.make_submission(Path(temporary))
            result = verify_submission(directory, self.contract)
            self.assertEqual(result["absolute_determinant"], str(21 * 2**22))

    def test_stale_receipt_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = self.make_submission(Path(temporary))
            (directory / "receipt.json").write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(SubmissionError, "stale or altered"):
                verify_submission(directory, self.contract)

    def test_unknown_file_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = self.make_submission(Path(temporary))
            (directory / "solver.py").write_text("print('never run')\n")
            with self.assertRaisesRegex(SubmissionError, "unexpected"):
                verify_submission(directory, self.contract)

    @unittest.skipUnless(hasattr(os, "symlink"), "symlinks unavailable")
    def test_symlinked_matrix_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directory = self.make_submission(root)
            (directory / "matrix.txt").unlink()
            os.symlink(ROOT / "records/genesis/matrix.txt", directory / "matrix.txt")
            with self.assertRaisesRegex(SubmissionError, "not a link"):
                verify_submission(directory, self.contract)

    def test_invalid_handle_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = self.make_submission(Path(temporary))
            metadata_path = directory / "metadata.json"
            metadata = json.loads(metadata_path.read_text())
            metadata["handle"] = "../../escape"
            metadata_path.write_text(json.dumps(metadata) + "\n")
            with self.assertRaisesRegex(SubmissionError, "safe lowercase"):
                verify_submission(directory, self.contract)

    def test_boolean_runtime_is_not_an_integer(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = self.make_submission(Path(temporary))
            metadata_path = directory / "metadata.json"
            metadata = json.loads(metadata_path.read_text())
            metadata["runtime_seconds"] = True
            metadata_path.write_text(json.dumps(metadata) + "\n")
            with self.assertRaisesRegex(SubmissionError, "non-negative integer"):
                verify_submission(directory, self.contract)

    def test_duplicate_metadata_key_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = self.make_submission(Path(temporary))
            metadata_path = directory / "metadata.json"
            raw = metadata_path.read_bytes()
            metadata_path.write_bytes(
                raw.replace(
                    b'"schema_version": 1',
                    b'"schema_version": 1, "schema_version": 1',
                    1,
                )
            )
            with self.assertRaisesRegex(SubmissionError, "duplicate"):
                verify_submission(directory, self.contract)

    def test_non_utf8_metadata_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = self.make_submission(Path(temporary))
            metadata_path = directory / "metadata.json"
            metadata_path.write_bytes(metadata_path.read_bytes() + b"\xff")
            with self.assertRaisesRegex(SubmissionError, "UTF-8"):
                verify_submission(directory, self.contract)

    def test_discovery_rejects_incomplete_tree(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "submissions"
            (root / "researcher" / "result-001").mkdir(parents=True)
            with self.assertRaisesRegex(SubmissionError, "missing submission files"):
                directory = discover_submission_directories(root)[0]
                verify_submission(directory, self.contract)

    @unittest.skipUnless(hasattr(os, "symlink"), "symlinks unavailable")
    def test_discovery_rejects_symlinked_handle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            submissions = temporary_root / "submissions"
            elsewhere = temporary_root / "elsewhere"
            submissions.mkdir()
            elsewhere.mkdir()
            os.symlink(elsewhere, submissions / "researcher")
            with self.assertRaisesRegex(SubmissionError, "real directory"):
                discover_submission_directories(submissions)


if __name__ == "__main__":
    unittest.main()
