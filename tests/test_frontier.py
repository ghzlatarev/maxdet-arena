from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from maxdet.contract import load_contract
from maxdet.frontier import effective_frontier
from maxdet.receipt import canonical_json_bytes, verify_matrix

ROOT = Path(__file__).resolve().parent.parent


class FrontierTests(unittest.TestCase):
    def make_root(self, temporary: str) -> tuple[Path, Path]:
        root = Path(temporary)
        (root / "data").mkdir()
        (root / "submissions" / "researcher" / "result-001").mkdir(parents=True)
        (root / "challenge.json").write_bytes((ROOT / "challenge.json").read_bytes())
        contract = load_contract(root / "challenge.json")

        genesis = verify_matrix(ROOT / "records/genesis/matrix.txt", contract)
        reference = verify_matrix(
            ROOT / "references/orrick-et-al-2003/matrix.txt",
            contract,
        )
        frontier = {
            "schema_version": 1,
            "challenge_id": contract.challenge_id,
            "status": "private-dogfooding",
            "target_to_beat": {
                "absolute_determinant": str(genesis.abs_determinant),
                "label": "test floor",
                "receipt_sha256": genesis.receipt["receipt_sha256"],
                "source": "records/genesis",
            },
            "arena_best": {
                "absolute_determinant": str(genesis.abs_determinant),
                "label": "test best",
                "receipt_sha256": genesis.receipt["receipt_sha256"],
                "source": "records/genesis",
            },
            "claim_boundary": "test only",
            "updated": "2026-07-23",
        }
        (root / "data" / "frontier.json").write_text(
            json.dumps(frontier),
            encoding="utf-8",
        )

        bundle = root / "submissions" / "researcher" / "result-001"
        bundle.joinpath("matrix.txt").write_bytes(
            (ROOT / "references/orrick-et-al-2003/matrix.txt").read_bytes()
        )
        bundle.joinpath("metadata.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "submission_id": "result-001",
                    "handle": "researcher",
                    "method": "frontier unit test",
                    "parent_receipt_sha256": genesis.receipt["receipt_sha256"],
                    "artifact_license": "CC0-1.0",
                }
            )
            + "\n",
            encoding="utf-8",
        )
        bundle.joinpath("receipt.json").write_bytes(
            canonical_json_bytes(reference.receipt)
        )
        return root, bundle

    def test_effective_frontier_includes_verified_submissions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, _ = self.make_root(temporary)
            contract = load_contract(root / "challenge.json")
            frontier = effective_frontier(root, contract)
            self.assertEqual(frontier.absolute_determinant, 2779447296000000)
            self.assertEqual(
                frontier.source,
                "submissions/researcher/result-001",
            )

    def test_bundle_can_be_excluded_for_eligibility_check(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, bundle = self.make_root(temporary)
            contract = load_contract(root / "challenge.json")
            frontier = effective_frontier(
                root,
                contract,
                exclude_directory=bundle,
            )
            self.assertEqual(frontier.absolute_determinant, 21 * 2**22)


if __name__ == "__main__":
    unittest.main()
