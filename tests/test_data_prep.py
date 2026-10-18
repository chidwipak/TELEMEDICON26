"""data_prep on a SYNTHETIC raw table: schema, censored-row removal, replicate averaging and the within-series split."""
import csv
import math
import tempfile
import unittest
from pathlib import Path

import _common  # noqa: F401
from cliffforge import data_prep


class TestDataPrep(unittest.TestCase):
    def test_curation_and_split(self):
        subs = ["F", "Cl", "Br", "C", "O", "N", "CC", "OC", "C#N", "S", "CO", "CCC"]   # 12 benzene analogues = one series
        rows = [(f"c1ccc({s})cc1", 10.0 ** (k % 4 - 1)) for k, s in enumerate(subs)]
        rows += [("c1ccc(F)cc1", 2.0)]                                       # replicate of the first compound
        rows += [("C1CCNCC1", 5.0), ("C1CCOCC1", 7.0)]                       # small 'series' (always train)
        rows += [("c1ccc(Cl)cc1C", ">99.5"), ("c1ccc(Cl)cc1CC", "0")]       # censored / non-positive -> dropped
        tmp = Path(tempfile.mkdtemp())
        raw = tmp / "raw.csv"
        with open(raw, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["suspected_SMILES", "Canonical_PostEra_ID", "IC50_(µM)"])
            for i, (smi, ic50) in enumerate(rows):
                w.writerow([smi, f"TOY-{i}", ic50])
        out = tmp / "out"
        data_prep.main(["--raw", str(raw), "--out-dir", str(out)])
        got = list(csv.DictReader(open(out / "all.csv")))
        self.assertEqual(list(got[0].keys()), ["smiles", "pIC50", "series", "series_size", "split", "n_rep", "rep_sd"])
        smiles = [r["smiles"] for r in got]
        self.assertEqual(len(smiles), len(set(smiles)))                      # duplicates collapsed
        self.assertNotIn("Cc1ccccc1Cl", smiles)                               # censored rows removed
        by = {r["smiles"]: r for r in got}
        rep = by["Fc1ccccc1"]
        self.assertEqual(int(rep["n_rep"]), 2)                                # replicate measurements merged
        self.assertAlmostEqual(float(rep["pIC50"]), 6 - (math.log10(0.1) + math.log10(2.0)) / 2, places=3)
        series = [r for r in got if int(r["series_size"]) >= 8]
        self.assertEqual(len(series), 12)
        n_test = sum(r["split"] == "test" for r in series)
        self.assertGreaterEqual(n_test, 3)
        self.assertTrue(all(r["split"] == "train" for r in got if int(r["series_size"]) < 8))


if __name__ == "__main__":
    unittest.main()
