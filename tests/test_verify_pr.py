from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from maxdet.contract import load_contract
from maxdet.exact import matrix_text
from maxdet.receipt import canonical_json_bytes, verify_matrix

ROOT = Path(__file__).resolve().parent.parent


class PullRequestVerifierTests(unittest.TestCase):
    def git(self, root: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["git", "-C", str(root), *arguments],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def clone(self, source: Path, target: Path) -> None:
        subprocess.run(
            ["git", "clone", "--quiet", "--no-hardlinks", str(source), str(target)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.git(target, "config", "user.name", "Arena Test")
        self.git(target, "config", "user.email", "arena-test@example.invalid")

    def commit_all(self, root: Path, message: str) -> None:
        self.git(root, "add", "-A")
        self.git(root, "commit", "--quiet", "-m", message)

    def add_bundle(
        self,
        root: Path,
        *,
        matrix_source: Path,
        parent: str | None,
    ) -> Path:
        bundle = root / "submissions" / "researcher" / "result-001"
        bundle.mkdir(parents=True)
        shutil.copyfile(matrix_source, bundle / "matrix.txt")
        contract = load_contract(root / "challenge.json")
        verified = verify_matrix(bundle / "matrix.txt", contract)
        metadata = {
            "schema_version": 1,
            "submission_id": "result-001",
            "handle": "researcher",
            "method": "pull-request integration test",
            "parent_receipt_sha256": parent,
            "artifact_license": "CC0-1.0",
        }
        (bundle / "metadata.json").write_text(
            json.dumps(metadata, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        (bundle / "receipt.json").write_bytes(
            canonical_json_bytes(verified.receipt)
        )
        return bundle

    def run_verifier(
        self,
        trusted: Path,
        untrusted: Path,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(trusted / "tools" / "verify_pr.py"),
                "--untrusted-root",
                str(untrusted),
            ],
            cwd=trusted,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def test_accepts_one_data_only_strict_improvement(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            trusted = temporary_root / "trusted"
            untrusted = temporary_root / "untrusted"
            self.clone(ROOT, trusted)

            contract = load_contract(trusted / "challenge.json")
            genesis = verify_matrix(
                trusted / "records" / "genesis" / "matrix.txt",
                contract,
            )
            frontier_path = trusted / "data" / "frontier.json"
            frontier = json.loads(frontier_path.read_text(encoding="utf-8"))
            frontier["target_to_beat"] = {
                "absolute_determinant": str(genesis.abs_determinant),
                "label": "integration-test floor",
                "receipt_sha256": genesis.receipt["receipt_sha256"],
                "source": "records/genesis",
            }
            frontier["arena_best"] = {
                "absolute_determinant": str(genesis.abs_determinant),
                "label": "integration-test checkpoint",
                "receipt_sha256": genesis.receipt["receipt_sha256"],
                "source": "records/genesis",
            }
            frontier_path.write_text(
                json.dumps(frontier, indent=2) + "\n",
                encoding="utf-8",
            )
            self.commit_all(trusted, "Lower floor for integration test")
            self.clone(trusted, untrusted)
            novel_matrix = [
                [int(token) for token in line.split()]
                for line in (
                    untrusted / "records" / "prelaunch-smoke" / "matrix.txt"
                ).read_text(encoding="ascii").splitlines()
            ]
            novel_matrix[0][0] *= -1
            novel_path = temporary_root / "novel-matrix.txt"
            novel_path.write_text(matrix_text(novel_matrix), encoding="ascii")
            self.add_bundle(
                untrusted,
                matrix_source=novel_path,
                parent=genesis.receipt["receipt_sha256"],
            )
            self.commit_all(untrusted, "Add improving matrix")

            result = self.run_verifier(trusted, untrusted)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("PULL REQUEST SUBMISSION VERIFIED", result.stdout)
            self.assertGreater(
                verify_matrix(novel_path, contract).abs_determinant,
                genesis.abs_determinant,
            )

    def test_rejects_a_tie_with_the_current_frontier(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            trusted = temporary_root / "trusted"
            untrusted = temporary_root / "untrusted"
            self.clone(ROOT, trusted)
            self.clone(trusted, untrusted)
            reference_receipt = json.loads(
                (
                    untrusted
                    / "references"
                    / "orrick-et-al-2003"
                    / "receipt.json"
                ).read_text(encoding="utf-8")
            )
            self.add_bundle(
                untrusted,
                matrix_source=(
                    untrusted
                    / "references"
                    / "orrick-et-al-2003"
                    / "matrix.txt"
                ),
                parent=reference_receipt["receipt_sha256"],
            )
            self.commit_all(untrusted, "Add tied matrix")

            result = self.run_verifier(trusted, untrusted)
            self.assertEqual(result.returncode, 2)
            self.assertIn("does not beat current frontier", result.stderr)

    def test_rejects_any_non_submission_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            trusted = temporary_root / "trusted"
            untrusted = temporary_root / "untrusted"
            self.clone(ROOT, trusted)
            self.clone(trusted, untrusted)
            with (untrusted / "README.md").open("a", encoding="utf-8") as stream:
                stream.write("\nuntrusted edit\n")
            self.commit_all(untrusted, "Edit trusted code")

            result = self.run_verifier(trusted, untrusted)
            self.assertEqual(result.returncode, 2)
            self.assertIn("only one submission bundle", result.stderr)

    @unittest.skipUnless(hasattr(os, "symlink"), "symlinks unavailable")
    def test_rejects_a_symlinked_submission_blob(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            trusted = temporary_root / "trusted"
            untrusted = temporary_root / "untrusted"
            self.clone(ROOT, trusted)
            self.clone(trusted, untrusted)
            reference = (
                untrusted
                / "references"
                / "orrick-et-al-2003"
                / "matrix.txt"
            )
            reference_receipt = json.loads(
                reference.with_name("receipt.json").read_text(encoding="utf-8")
            )
            bundle = self.add_bundle(
                untrusted,
                matrix_source=reference,
                parent=reference_receipt["receipt_sha256"],
            )
            (bundle / "matrix.txt").unlink()
            os.symlink(
                "../../../references/orrick-et-al-2003/matrix.txt",
                bundle / "matrix.txt",
            )
            self.commit_all(untrusted, "Add symlinked matrix")

            result = self.run_verifier(trusted, untrusted)
            self.assertEqual(result.returncode, 2)
            self.assertIn("regular 0644 blob", result.stderr)


if __name__ == "__main__":
    unittest.main()
