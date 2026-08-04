from __future__ import annotations

import hashlib
import re
import unittest
from pathlib import Path

from maxdet.json_tools import loads_strict_json


ROOT = Path(__file__).resolve().parent.parent
DEPLOYMENT_PATH = ROOT / "contracts" / "deployments" / "sepolia.json"
ADDRESS_PATTERN = re.compile(r"^0x[0-9a-fA-F]{40}$")
HASH_PATTERN = re.compile(r"^0x[0-9a-fA-F]{64}$")


class BountyDeploymentRegistryTests(unittest.TestCase):
    def test_registry_is_pinned_to_the_compiled_source(self) -> None:
        deployment = loads_strict_json(DEPLOYMENT_PATH.read_bytes())
        self.assertIsInstance(deployment, dict)
        source_path = ROOT / deployment["source"]

        self.assertEqual(deployment["schema_version"], 1)
        self.assertEqual(deployment["network"], "Sepolia")
        self.assertEqual(deployment["chain_id"], 11_155_111)
        self.assertEqual(deployment["asset"], "Sepolia ETH")
        self.assertEqual(deployment["contract_name"], "MaxDetBounty23")
        self.assertEqual(
            deployment["minimum_winning_determinant"],
            "2779447300194304",
        )
        self.assertEqual(
            deployment["source_sha256"],
            hashlib.sha256(source_path.read_bytes()).hexdigest(),
        )
        self.assertRegex(deployment["expected_runtime_codehash"], HASH_PATTERN)
        self.assertEqual(
            deployment["compiler"],
            {
                "version": "0.8.36",
                "evm_version": "osaka",
                "optimizer": True,
                "optimizer_runs": 200,
                "via_ir": True,
            },
        )

    def test_live_fields_are_all_present_or_all_absent(self) -> None:
        deployment = loads_strict_json(DEPLOYMENT_PATH.read_bytes())
        status = deployment["status"]
        address = deployment["contract_address"]
        transaction = deployment["deployment_transaction"]
        block = deployment["deployment_block"]

        self.assertIn(status, {"awaiting-deployment", "live"})
        if status == "awaiting-deployment":
            self.assertIsNone(address)
            self.assertIsNone(transaction)
            self.assertIsNone(block)
        else:
            self.assertRegex(address, ADDRESS_PATTERN)
            self.assertRegex(transaction, HASH_PATTERN)
            self.assertIsInstance(block, int)
            self.assertGreater(block, 0)


if __name__ == "__main__":
    unittest.main()
