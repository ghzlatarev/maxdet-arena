from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
PUBLISHER = ROOT / "tools" / "update_search_progress.py"
FRONTIER = 2_779_447_296_000_000


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def write_jsonl(path: Path, values: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(f"{json.dumps(value)}\n" for value in values),
        encoding="utf-8",
    )


def campaign_fixture(arm_ids: list[str]) -> dict[str, object]:
    return {
        "schema_version": 1,
        "campaign_id": "test-gram-campaign",
        "claim_boundary": (
            "Live Gram sketches are invariant descriptors, not canonical "
            "classes."
        ),
        "known_core_frontier": 662_671_875,
        "strict_core_target": 662_671_876,
        "frontier_absolute_determinant": str(FRONTIER),
        "first_strict_absolute_determinant": "2779447300194304",
        "quotient_gate": 620_000_000,
        "gate_absolute_determinant": "2600468480000000",
        "budget_seconds_per_arm": 120,
        "seed_basins": [
            {
                "id": "seed-a",
                "label": "Seed A",
                "path": "/private/source/matrix.txt",
                "absolute_determinant": str(FRONTIER),
                "ht_gram_basin_key_sha256": "a" * 64,
            },
            {
                "id": "seed-b",
                "label": "Seed B",
                "path": "/private/source/other.txt",
                "absolute_determinant": "2700000000000000",
                "ht_gram_basin_key_sha256": "b" * 64,
            },
        ],
        "arms": [
            {
                "id": arm_id,
                "random_seed": 37_000 + index,
                "kick_flips": 12 * (index + 1),
                "epoch_moves": 10_000,
            }
            for index, arm_id in enumerate(arm_ids)
        ],
    }


class SearchProgressPublisherTests(unittest.TestCase):
    def run_publisher(
        self, run_root: Path, *extra_arguments: str
    ) -> dict[str, object]:
        output = run_root / "published.json"
        subprocess.run(
            [
                sys.executable,
                str(PUBLISHER),
                "--run-root",
                str(run_root),
                "--output",
                str(output),
                *extra_arguments,
            ],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        return json.loads(output.read_text(encoding="utf-8"))

    def test_campaign_schema_is_compact_exact_and_sanitized(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_root = Path(temporary)
            write_json(
                run_root / "campaign.json",
                campaign_fixture(["arm-k12"]),
            )
            events: list[dict[str, object]] = []
            for index in range(100):
                events.append(
                    {
                        "event": "epoch",
                        "engine": "test-hopper",
                        "elapsed_seconds": float(index),
                        "epochs_completed": index,
                        "archive_size": min(256, index + 8),
                        "sketch_discoveries": index * 3,
                        "search_added_sketch_discoveries": index * 2,
                        "best_core_quotient": 662_671_875,
                        "best_absolute_determinant": str(FRONTIER),
                        "strict_target_states": 0,
                    }
                )
            events.append(
                {
                    **events[-1],
                    "event": "finished",
                    "elapsed_seconds": 120.0,
                }
            )
            write_jsonl(run_root / "arm-k12" / "run.jsonl", events)
            write_json(
                run_root / "arm-k12" / "summary.json",
                {
                    "complete": True,
                    "engine": "test-hopper",
                    "elapsed_seconds": 120,
                    "best_core_quotient": 662_671_875,
                    "best_absolute_determinant": str(FRONTIER),
                    "archive_size": 107,
                    "archive_capacity": 256,
                    "sketch_discoveries": 297,
                    "search_added_sketch_discoveries": 198,
                    "epochs_completed": 99,
                    "strict_target_states": 0,
                },
            )

            result = self.run_publisher(
                run_root, "--max-history-points", "7"
            )

            self.assertEqual(result["schema_version"], 2)
            self.assertEqual(result["campaign"]["status"], "complete")
            self.assertEqual(result["best_absolute_determinant"], str(FRONTIER))
            self.assertEqual(result["totals"]["archive_cells"], 107)
            self.assertEqual(result["totals"]["sketch_discoveries"], 198)
            self.assertEqual(len(result["seed_basins"]), 2)
            self.assertNotIn("path", result["seed_basins"][0])
            self.assertEqual(
                result["seed_basins"][0]["gram_basin_key_sha256"], "a" * 64
            )
            arm = result["arms"][0]
            self.assertEqual(arm["label"], "k=12")
            self.assertEqual(arm["seed_basin_ids"], ["seed-a", "seed-b"])
            self.assertEqual(len(arm["history"]), 7)
            self.assertEqual(arm["history"][0]["elapsed_seconds"], 0.0)
            self.assertEqual(arm["history"][-1]["elapsed_seconds"], 120.0)
            self.assertEqual(arm["history"][-1]["sketch_discoveries"], 198)
            self.assertEqual(arm["archive_size"], 107)
            self.assertIsInstance(result["source_updated_at"], str)

    def test_source_mtime_marks_unfinished_arm_stale(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_root = Path(temporary)
            write_json(
                run_root / "campaign.json",
                campaign_fixture(["arm-k12"]),
            )
            log = run_root / "arm-k12" / "run.jsonl"
            write_jsonl(
                log,
                [
                    {
                        "event": "start",
                        "engine": "test-hopper",
                        "elapsed_seconds": 2.0,
                        "archive_size": 8,
                        "sketch_discoveries": 8,
                        "epochs_completed": 0,
                        "best_absolute_determinant": str(FRONTIER),
                    }
                ],
            )
            old = time.time() - 120
            os.utime(log, (old, old))

            result = self.run_publisher(
                run_root, "--stale-seconds", "5"
            )

            self.assertEqual(result["campaign"]["status"], "stale")
            self.assertEqual(result["arms"][0]["status"], "stale")
            self.assertNotEqual(
                result["source_updated_at"], result["updated_at"]
            )

    def test_exact_coronal_manifests_are_sanitized_and_deduplicated(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_root = Path(temporary)
            campaign = campaign_fixture(["arm-k12"])
            campaign["parent_policy"] = "exact_coronal_pareto"
            write_json(run_root / "campaign.json", campaign)
            write_jsonl(
                run_root / "arm-k12" / "run.jsonl",
                [
                    {
                        "event": "finished",
                        "elapsed_seconds": 120,
                        "best_absolute_determinant": str(FRONTIER),
                    }
                ],
            )
            leaked_path = str(run_root / "private" / "frontier.matrix.txt")
            write_json(
                run_root / "arm-k12" / "archive" / "manifest.json",
                {
                    "schema_version": 1,
                    "elites": [
                        {
                            "rank": 0,
                            "path": leaked_path,
                            "origin": "seed",
                            "core_quotient": 662_671_875,
                            "absolute_determinant": str(FRONTIER),
                            "coronal_points": [
                                {
                                    "det_m": "267647730468750000",
                                    "kappa_numerator": "6455",
                                    "kappa_denominator": "2736",
                                },
                                {
                                    "det_m": "267647730468750000",
                                    "kappa_numerator": "6455",
                                    "kappa_denominator": "2736",
                                },
                            ],
                        },
                        {
                            "rank": 1,
                            "path": "/private/seed-b.matrix.txt",
                            "origin": "seed",
                            "core_quotient": 643_730_957,
                            "absolute_determinant": "2700000000000000",
                            "coronal_points": [
                                {
                                    "det_m": "300000000000000000",
                                    "kappa_numerator": "5",
                                    "kappa_denominator": "2",
                                }
                            ],
                        },
                    ],
                },
            )

            result = self.run_publisher(run_root)

            projection = result["coronal_projection"]
            self.assertEqual(projection["kind"], "retained_exact_coronal")
            self.assertEqual(projection["x_preference"], "larger")
            self.assertEqual(projection["y_preference"], "smaller")
            points = projection["points"]
            self.assertEqual(len(points), 2)
            frontier = points[0]
            self.assertEqual(frontier["arm_id"], "arm-k12")
            self.assertEqual(frontier["rank"], 0)
            self.assertEqual(frontier["orientation"], 0)
            self.assertEqual(frontier["core_quotient"], "662671875")
            self.assertEqual(
                frontier["absolute_determinant"], str(FRONTIER)
            )
            self.assertEqual(frontier["det_m"], "267647730468750000")
            self.assertEqual(frontier["kappa_numerator"], "6455")
            self.assertEqual(frontier["kappa_denominator"], "2736")
            self.assertAlmostEqual(
                frontier["kappa_decimal"], 6455 / 2736
            )
            self.assertEqual(frontier["identity"], "frontier")
            self.assertEqual(frontier["seed_id"], "seed-a")
            self.assertIs(frontier["pareto_front"], True)
            self.assertIs(points[1]["pareto_front"], True)
            serialized = json.dumps(result)
            self.assertNotIn(leaked_path, serialized)
            self.assertNotIn("/private/seed-b.matrix.txt", serialized)
            self.assertNotIn('"path"', serialized)

    def test_missing_or_malformed_coronal_manifest_keeps_time_score_fallback(
        self,
    ) -> None:
        for mode in ("missing", "malformed"):
            with self.subTest(mode=mode):
                with tempfile.TemporaryDirectory() as temporary:
                    run_root = Path(temporary)
                    campaign = campaign_fixture(["arm-k12"])
                    campaign["parent_policy"] = "exact_coronal_pareto"
                    write_json(run_root / "campaign.json", campaign)
                    write_jsonl(
                        run_root / "arm-k12" / "run.jsonl",
                        [
                            {
                                "event": "finished",
                                "elapsed_seconds": 120,
                                "best_absolute_determinant": str(FRONTIER),
                            }
                        ],
                    )
                    if mode == "malformed":
                        manifest = (
                            run_root
                            / "arm-k12"
                            / "archive"
                            / "manifest.json"
                        )
                        manifest.parent.mkdir(parents=True, exist_ok=True)
                        manifest.write_text(
                            '{"elites": [this is incomplete',
                            encoding="utf-8",
                        )

                    result = self.run_publisher(run_root)

                    self.assertEqual(result["schema_version"], 2)
                    self.assertNotIn("coronal_projection", result)
                    self.assertEqual(result["arms"][0]["id"], "arm-k12")

    def test_legacy_run_root_keeps_schema_one(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_root = Path(temporary)
            write_jsonl(
                run_root / "arm0" / "run.jsonl",
                [
                    {
                        "event": "finished",
                        "engine": "legacy-test",
                        "elapsed_seconds": 10,
                        "seconds": 10,
                        "best_absolute_determinant": str(FRONTIER),
                    }
                ],
            )

            result = self.run_publisher(run_root, "--arms", "1")

            self.assertEqual(result["schema_version"], 1)
            self.assertEqual(result["arms"][0]["status"], "complete")
            self.assertNotIn("campaign", result)


if __name__ == "__main__":
    unittest.main()
