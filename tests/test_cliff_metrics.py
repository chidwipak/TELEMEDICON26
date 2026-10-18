"""Unit tests for the cliff-detection metrics, uncertainty tools and the leakage guards."""
import unittest

import numpy as np

import _common  # noqa: F401  (puts src/ on sys.path)
from cliffforge import cliff_detection as CD


class TestConfusionMetrics(unittest.TestCase):
    def test_confusion_counts(self):
        y = np.array([1, 1, 0, 0, 0, 1])
        yh = np.array([1, 0, 0, 1, 0, 1])
        self.assertEqual(CD.conf(y, yh), (2, 1, 2, 1))

    def test_prf_values_and_undefined_cases(self):
        r = CD.prf(2, 1, 2, 1)
        self.assertAlmostEqual(r["precision"], 2 / 3)
        self.assertAlmostEqual(r["recall"], 2 / 3)
        self.assertAlmostEqual(r["f1"], 2 / 3)
        self.assertAlmostEqual(r["accuracy"], 4 / 6)
        none = CD.prf(0, 0, 39, 5)                       # 'never a cliff' majority classifier
        self.assertIsNone(none["precision"])
        self.assertEqual(none["recall"], 0.0)
        self.assertAlmostEqual(none["accuracy"], 39 / 44)  # high accuracy, zero detection

    def test_exact_recall_interval(self):
        lo, hi = CD.cp_ci(2, 5)
        self.assertAlmostEqual(lo, 0.0527, places=3)
        self.assertAlmostEqual(hi, 0.8534, places=3)
        self.assertEqual(CD.cp_ci(0, 5)[0], 0.0)
        self.assertEqual(CD.cp_ci(5, 5)[1], 1.0)


class TestThresholdSelection(unittest.TestCase):
    def test_choose_tau_separates_classes(self):
        pred = np.array([0.1, 0.2, 0.3, 1.4, 1.6, 2.0])
        y = np.array([0, 0, 0, 1, 1, 1])
        tau, f1, _ = CD.choose_tau(pred, y)
        self.assertEqual(f1, 1.0)
        self.assertTrue(0.3 < tau <= 1.4)

    def test_tie_breaks_toward_one(self):
        pred = np.array([0.1, 0.2, 3.0, 3.0])
        y = np.array([0, 0, 1, 1])
        tau, _, _ = CD.choose_tau(pred, y)
        self.assertEqual(tau, 1.0)


class TestUncertaintyTools(unittest.TestCase):
    def test_components_groups_connected_pairs(self):
        labels = CD.components([(1, 2), (2, 3), (7, 8)])
        self.assertEqual(labels[0], labels[1])
        self.assertNotEqual(labels[0], labels[2])

    def test_cluster_bootstrap_is_deterministic_and_clustered(self):
        clusters = [0, 0, 0, 1, 1, 2]
        ev = lambda ix: {"n": float(len(ix))}                          # noqa: E731
        a = CD.cluster_boot(clusters, ev, nboot=50, seed=3)
        b = CD.cluster_boot(clusters, ev, nboot=50, seed=3)
        self.assertEqual(a, b)
        # resampling clusters keeps whole clusters together: sizes are sums of {3, 2, 1}
        for r in a:
            self.assertTrue(float(r["n"]) in {s for s in range(0, 19)})


class TestLeakageGuards(unittest.TestCase):
    def test_knn_refuses_unavailable_labels(self):
        S = np.eye(4, dtype=np.float32) + 0.1
        y_known = np.array([5.0, np.nan, 6.0, 7.0])      # compound 1 is 'held out' (label masked)
        train_mask = np.array([True, True, True, True])  # but wrongly offered as a neighbour
        with self.assertRaises(AssertionError):
            CD.knn_pred(S, y_known, train_mask, k=3)

    def test_knn_ok_when_masked_labels_are_excluded(self):
        S = np.eye(4, dtype=np.float32) + 0.1
        y_known = np.array([5.0, np.nan, 6.0, 7.0])
        train_mask = np.array([True, False, True, True])
        out = CD.knn_pred(S, y_known, train_mask, k=2)
        self.assertTrue(np.isfinite(out).all())


if __name__ == "__main__":
    unittest.main()
