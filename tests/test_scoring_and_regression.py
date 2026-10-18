"""Tests that need the externally supplied relative-model weights and/or the curated Moonshot table.

They are skipped (not failed) when those resources are absent: see models/README.md and data/README.md.
Set CLIFFFORGE_SKIP_SLOW=1 to skip the ~40 s end-to-end demo regression test.
"""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from _common import HAVE_DATA, HAVE_WEIGHTS, ROOT, SKIP_SLOW


class TestFeaturization(unittest.TestCase):
    def test_graph_shapes_match_checkpoint_convention(self):
        from rdkit import Chem
        from cliffforge.features import NODE_DIM, EDGE_DIM, mol_to_whole_graph
        g = mol_to_whole_graph(Chem.MolFromSmiles("CC(=O)Nc1ccccc1"))
        self.assertEqual(g["x"].shape, (10, NODE_DIM))
        self.assertEqual(g["edge_attr"].shape[1], EDGE_DIM)
        self.assertEqual((NODE_DIM, EDGE_DIM), (48, 8))


@unittest.skipUnless(HAVE_WEIGHTS, "relative-model weights not present (models/README.md)")
class TestRelativeScorer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import torch
        from cliffforge.relmodel import RelScorer
        cls.scorer = RelScorer(device=torch.device("cpu"))

    def test_ensemble_size(self):
        self.assertEqual(len(self.scorer.models), 3)

    def test_pair_delta_is_antisymmetric_and_zero_on_identical_molecules(self):
        a = ["CNC(=O)CN1CC(C(=O)Nc2cncc3ccc(NS(C)(=O)=O)cc23)c2cc(Cl)ccc2C1=O"]
        b = ["CNC(=O)CN1CC(C(=O)Nc2cncc3ccc(NS(C)(=O)=O)cc23)c2cc(F)ccc2C1=O"]
        ab, _ = self.scorer.pair_delta(a, b)
        ba, _ = self.scorer.pair_delta(b, a)
        aa, _ = self.scorer.pair_delta(a, a)
        self.assertAlmostEqual(float(ab[0]), -float(ba[0]), places=5)
        self.assertAlmostEqual(float(aa[0]), 0.0, places=5)

    def test_scores_are_deterministic_and_finite(self):
        a = ["Cc1ccc(NC(=O)Cc2cccc(Cl)c2)cn1"] * 2
        b = ["Cc1ccc(NC(=O)Cc2cccc(F)c2)cn1", "Cc1ccc(NC(=O)Cc2cccc(Br)c2)cn1"]
        m1, s1 = self.scorer.pair_delta(a, b)
        m2, _ = self.scorer.pair_delta(a, b)
        self.assertTrue(np.isfinite(m1).all() and np.isfinite(s1).all())
        np.testing.assert_allclose(m1, m2)


@unittest.skipUnless(HAVE_DATA and HAVE_WEIGHTS and not SKIP_SLOW, "needs data + weights (and CLIFFFORGE_SKIP_SLOW unset)")
class TestDemoRegression(unittest.TestCase):
    def test_demo_reproduces_reference_funnel(self):
        from cliffforge import demo
        out = Path(tempfile.mkdtemp())
        demo.main(["--auto-seed", "--out-dir", str(out)])
        got = json.load(open(out / "funnel.json"))
        ref = json.load(open(ROOT / "results" / "metrics" / "funnel.json"))
        for k in ("generated_edits", "valid_rdkit", "unique_non_parent", "pass_all_filters", "scored_by_model"):
            self.assertEqual(got[k], ref[k], k)
        self.assertEqual(got["seed"]["all_csv_row"], ref["seed"]["all_csv_row"])


if __name__ == "__main__":
    unittest.main()
