from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from maxdet.cli import command_prepare, command_verify
from maxdet.contract import load_contract
from maxdet.errors import SubmissionError
from maxdet.receipt import verify_matrix
from maxdet.submission import verify_submission

ROOT = Path(__file__).resolve().parent.parent


class VerifyCommandTests(unittest.TestCase):
    def test_json_output_rejects_a_matrix_like_path(self) -> None:
        arguments = argparse.Namespace(
            matrix=Path("candidate/matrix.txt"),
            contract=None,
            json_path=Path("references/seed.matrix.txt"),
            quiet=True,
        )
        with self.assertRaisesRegex(ValueError, "must end in .json"):
            command_verify(arguments)

    def test_json_output_rejects_the_input_path(self) -> None:
        arguments = argparse.Namespace(
            matrix=Path("candidate/receipt.json"),
            contract=None,
            json_path=Path("candidate/receipt.json"),
            quiet=True,
        )
        with self.assertRaisesRegex(ValueError, "must differ"):
            command_verify(arguments)


class PrepareCommandTests(unittest.TestCase):
    def make_root(self, temporary: str) -> Path:
        root = Path(temporary)
        (root / "candidate").mkdir()
        (root / "data").mkdir()
        (root / "challenge.json").write_bytes((ROOT / "challenge.json").read_bytes())
        (root / "candidate" / "matrix.txt").write_bytes(
            (ROOT / "references/orrick-et-al-2003/matrix.txt").read_bytes()
        )
        contract = load_contract(root / "challenge.json")
        floor = verify_matrix(ROOT / "records/genesis/matrix.txt", contract)
        frontier = {
            "schema_version": 1,
            "challenge_id": contract.challenge_id,
            "status": "private-dogfooding",
            "target_to_beat": {
                "absolute_determinant": str(floor.abs_determinant),
                "label": "test floor",
                "receipt_sha256": floor.receipt["receipt_sha256"],
                "source": "records/genesis",
            },
            "arena_best": {
                "absolute_determinant": str(floor.abs_determinant),
                "label": "test best",
                "receipt_sha256": floor.receipt["receipt_sha256"],
                "source": "records/genesis",
            },
            "claim_boundary": "test only",
            "updated": "2026-07-23",
        }
        (root / "data" / "frontier.json").write_text(
            json.dumps(frontier),
            encoding="utf-8",
        )
        return root

    def args(
        self,
        *,
        method: str = "unit-test search",
        parent: str | None = None,
        agent: str | None = None,
        runtime_seconds: int | None = None,
        seed: int | None = None,
    ) -> argparse.Namespace:
        return argparse.Namespace(
            submission_id="result-001",
            handle="researcher",
            method=method,
            parent=parent,
            agent=agent,
            runtime_seconds=runtime_seconds,
            seed=seed,
            notes=None,
        )

    def test_prepare_creates_a_complete_verified_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = self.make_root(temporary)
            with patch("maxdet.cli.repository_root", return_value=root):
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(
                        command_prepare(
                            self.args(
                                agent="codex",
                                runtime_seconds=12,
                                seed=2301,
                            )
                        ),
                        0,
                    )
            bundle = root / "submissions" / "researcher" / "result-001"
            contract = load_contract(root / "challenge.json")
            result = verify_submission(bundle, contract)
            self.assertEqual(result["absolute_determinant"], "2779447296000000")
            metadata = json.loads((bundle / "metadata.json").read_text())
            self.assertEqual(metadata["agent"], "codex")
            self.assertEqual(metadata["runtime_seconds"], 12)
            self.assertEqual(metadata["seed"], 2301)
            self.assertFalse(any((root / "runs").iterdir()))

    def test_invalid_metadata_leaves_no_partial_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = self.make_root(temporary)
            with (
                patch("maxdet.cli.repository_root", return_value=root),
                self.assertRaisesRegex(SubmissionError, "not be blank"),
            ):
                command_prepare(self.args(method="   "))
            self.assertFalse((root / "submissions").exists())

    def test_non_improvement_is_rejected_before_creating_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = self.make_root(temporary)
            frontier_path = root / "data" / "frontier.json"
            frontier = json.loads(frontier_path.read_text())
            reference = verify_matrix(
                ROOT / "references/orrick-et-al-2003/matrix.txt",
                load_contract(root / "challenge.json"),
            )
            frontier["target_to_beat"]["absolute_determinant"] = str(
                reference.abs_determinant
            )
            frontier["target_to_beat"]["receipt_sha256"] = reference.receipt[
                "receipt_sha256"
            ]
            frontier_path.write_text(json.dumps(frontier), encoding="utf-8")
            with (
                patch("maxdet.cli.repository_root", return_value=root),
                self.assertRaisesRegex(ValueError, "does not beat"),
            ):
                command_prepare(self.args())
            self.assertFalse((root / "submissions").exists())

    def test_unknown_parent_is_rejected_before_creating_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = self.make_root(temporary)
            with (
                patch("maxdet.cli.repository_root", return_value=root),
                self.assertRaisesRegex(ValueError, "trusted artifacts"),
            ):
                command_prepare(self.args(parent="a" * 64))
            self.assertFalse((root / "submissions").exists())

    @unittest.skipUnless(hasattr(os, "symlink"), "symlinks unavailable")
    def test_prepare_rejects_a_symlinked_submissions_root(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = self.make_root(temporary)
            elsewhere = root / "elsewhere"
            elsewhere.mkdir()
            os.symlink(elsewhere, root / "submissions")
            with (
                patch("maxdet.cli.repository_root", return_value=root),
                self.assertRaisesRegex(SubmissionError, "real directory"),
            ):
                command_prepare(self.args())
            self.assertEqual(list(elsewhere.iterdir()), [])
            self.assertFalse((root / "runs").exists())


if __name__ == "__main__":
    unittest.main()
