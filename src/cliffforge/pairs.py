"""Pair construction/collation for relative-affinity training (extracted unchanged from moonshot_ddg.py)."""
from collections import defaultdict

import numpy as np
import torch
from torch.utils.data import Dataset
from torch_geometric.data import Batch

MAX_PAIRS = 600


def make_pairs(items, seed):
    by = defaultdict(list)
    for i, it in enumerate(items):
        if it["split"] == "train":
            by[it["series"]].append(i)
    rng = np.random.RandomState(seed)
    pairs = []
    for s, idx in by.items():
        if len(idx) < 2:
            continue
        loc = [(i, j) for i in idx for j in idx if i != j]
        if len(loc) > MAX_PAIRS:
            sel = rng.choice(len(loc), MAX_PAIRS, replace=False)
            loc = [loc[k] for k in sel]
        pairs.extend(loc)
    return pairs


class PairSet(Dataset):
    def __init__(self, items, pairs):
        self.items = items; self.pairs = pairs
    def __len__(self): return len(self.pairs)
    def __getitem__(self, k):
        i, j = self.pairs[k]
        return self.items[i]["data"], self.items[j]["data"], \
            float(self.items[i]["y"] - self.items[j]["y"])


def collate(b):
    return (Batch.from_data_list([x[0] for x in b]),
            Batch.from_data_list([x[1] for x in b]),
            torch.tensor([x[2] for x in b], dtype=torch.float))
