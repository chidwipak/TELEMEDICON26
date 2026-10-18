"""Synthetic checks for training on the published, molecule-split MMP rows."""
import csv
import tempfile
import unittest
from pathlib import Path

import _common  # noqa: F401
from cliffforge.pairs import PairSet, load_mmp_train_pairs
from cliffforge.relmodel import load_mmp_train_items


PAIR_COLUMNS = [
    "smiles_1", "smiles_2", "transformation", "constant_part",
    "pIC50_1", "pIC50_2", "delta_pIC50", "abs_delta_pIC50",
    "pair_class", "cliff_label", "split",
]


def pair_row(first, second, y_first, y_second, pair_class="half_cliff", split="train"):
    delta = y_first - y_second
    return {
        "smiles_1": first, "smiles_2": second,
        "transformation": "toy transformation", "constant_part": "C",
        "pIC50_1": y_first, "pIC50_2": y_second,
        "delta_pIC50": delta, "abs_delta_pIC50": abs(delta),
        "pair_class": pair_class,
        "cliff_label": "" if pair_class == "half_cliff" else int(pair_class == "cliff"),
        "split": split,
    }


class TestPublishedMmpPairs(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.pair_csv = self.root / "mmp_train.csv"
        self.items = [
            {"smiles": "CCO", "y": 6.0, "split": "train", "data": "graph-CCO"},
            {"smiles": "CCN", "y": 5.0, "split": "train", "data": "graph-CCN"},
            {"smiles": "CCC", "y": 4.5, "split": "train", "data": "graph-CCC"},
        ]

    def write_pairs(self, rows):
        with self.pair_csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=PAIR_COLUMNS)
            writer.writeheader()
            writer.writerows(rows)

    def test_published_order_targets_and_half_cliffs_are_preserved(self):
        self.write_pairs([
            pair_row("CCN", "CCO", 5.0, 6.0, "non_cliff"),
            pair_row("CCO", "CCC", 6.0, 4.5, "half_cliff"),
        ])
        pairs = load_mmp_train_pairs(self.items, self.pair_csv)

        # Exactly the two published rows, in their original orientation and order.
        self.assertEqual(pairs, [(1, 0, -1.0), (0, 2, 1.5)])
        dataset = PairSet(self.items, pairs)
        self.assertEqual(dataset[0], ("graph-CCN", "graph-CCO", -1.0))
        self.assertEqual(dataset[1][2], 1.5)

    def test_pair_endpoint_must_be_in_training_molecules(self):
        self.write_pairs([pair_row("CCO", "N#N", 6.0, 5.0)])
        with self.assertRaises(ValueError):
            load_mmp_train_pairs(self.items, self.pair_csv)

    def test_test_mixed_and_unknown_split_rows_are_rejected(self):
        for split in ("test", "mixed", "unknown"):
            with self.subTest(split=split):
                self.write_pairs([pair_row("CCO", "CCN", 6.0, 5.0, split=split)])
                with self.assertRaises(ValueError):
                    load_mmp_train_pairs(self.items, self.pair_csv)

    def test_pair_activity_delta_must_agree_with_ordered_molecules(self):
        row = pair_row("CCN", "CCO", 5.0, 6.0)
        row["delta_pIC50"] = 1.0  # reversed sign would train the wrong target
        self.write_pairs([row])
        with self.assertRaises(ValueError):
            load_mmp_train_pairs(self.items, self.pair_csv)

    def test_duplicate_train_smiles_are_rejected(self):
        train_csv = self.root / "train.csv"
        with train_csv.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["smiles", "pIC50"])
            writer.writerows([("CCO", 6.0), ("CCO", 6.0)])
        with self.assertRaises(ValueError):
            load_mmp_train_items(train_csv)


if __name__ == "__main__":
    unittest.main()
