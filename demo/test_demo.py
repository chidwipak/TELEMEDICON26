"""Regression checks for the live fold-3 example and known-reference arithmetic."""
from __future__ import annotations

import csv
import sys
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

sys.path.insert(0, str(Path(__file__).resolve().parent))
from app import app
from checkpoint_verify import RESULTS_ROOT


class DemoRegressionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.client_context = TestClient(app)
        cls.client = cls.client_context.__enter__()
        cls.example = cls.client.get("/api/example").json()
        cls.prediction_response = cls.client.post("/api/predict", json={
            "candidate_smiles": cls.example["candidate_smiles"],
            "reference_smiles": cls.example["reference_smiles"],
            "reference_activity": cls.example["reference_activity"],
            "unit": cls.example["unit"],
        })
        cls.prediction_response.raise_for_status()
        cls.prediction = cls.prediction_response.json()

    @classmethod
    def tearDownClass(cls) -> None:
        cls.client_context.__exit__(None, None, None)

    def test_archive_parity_and_example_provenance(self) -> None:
        with (RESULTS_ROOT / "fold_3/pair_predictions.csv").open(newline="") as handle:
            archive = next(row for row in csv.DictReader(handle) if row["source_pair_id"] == "10765")
        result = self.prediction
        self.assertEqual(result["fold"], 3)
        self.assertEqual(result["example_pair_id"], 10765)
        self.assertEqual(result["mmp_status"], "supplied Dablander MMP")
        self.assertEqual(archive["setting"], "test_train")
        self.assertEqual(archive["a_id"], "1655")
        self.assertEqual(archive["b_id"], "1614")
        self.assertEqual(archive["class"], "AC")
        for model in ("cliffforge", "chemprop"):
            self.assertAlmostEqual(result[model]["mean_delta"], float(archive[model]), delta=1e-4)
            for seed in range(3):
                self.assertAlmostEqual(
                    result[model]["seed_deltas"][seed],
                    float(archive[f"{model}_seed_{seed}"]),
                    delta=1e-4,
                )
        self.assertTrue(result["cliffforge"]["predicted_large_difference"])
        self.assertFalse(result["chemprop"]["predicted_large_difference"])
        self.assertNotIn("true_delta", result)
        self.assertNotIn("candidate_measured_pIC50", result)
        truth = self.client.get("/api/example/reveal").json()
        self.assertAlmostEqual(truth["true_delta"], float(archive["true_delta"]), places=6)

    def test_graph_visuals_match_input_graphs(self) -> None:
        for molecule in self.prediction["molecules"].values():
            self.assertEqual(molecule["directed_edges"], 2 * molecule["bonds"])
            self.assertEqual(molecule["graph_svg"].count("<circle"), molecule["atoms"])
            self.assertEqual(molecule["graph_svg"].count("<line"), molecule["bonds"])
            self.assertIn("<svg", molecule["svg"])

    def test_chemprop_uses_measured_reference(self) -> None:
        changed = self.client.post("/api/predict", json={
            "candidate_smiles": self.example["candidate_smiles"],
            "reference_smiles": self.example["reference_smiles"],
            "reference_activity": self.example["reference_pIC50"] + 0.5,
            "unit": "pIC50",
        })
        changed.raise_for_status()
        changed_result = changed.json()
        original = self.prediction
        self.assertIsNone(changed_result["example_pair_id"])
        for model in ("cliffforge", "chemprop"):
            before = original[model]
            after = changed_result[model]
            if model == "chemprop":
                for first, second in zip(before["seed_candidate_pIC50"], after["seed_candidate_pIC50"]):
                    self.assertAlmostEqual(first, second, places=5)
                for first, second in zip(before["seed_deltas"], after["seed_deltas"]):
                    self.assertAlmostEqual(second, first - 0.5, places=5)
            else:
                self.assertAlmostEqual(before["mean_delta"], after["mean_delta"], places=5)

    def test_invalid_smiles_is_rejected(self) -> None:
        bad = self.client.post("/api/predict", json={
            "candidate_smiles": "not-a-smiles",
            "reference_smiles": self.example["reference_smiles"],
            "reference_activity": 85.874,
            "unit": "ic50_um",
        })
        self.assertEqual(bad.status_code, 422)
        self.assertIn("could not parse SMILES", bad.json()["detail"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
