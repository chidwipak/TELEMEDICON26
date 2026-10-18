"""
cliffforge.relmodel
======================
Target-specific (SARS-CoV-2 Mpro, COVID Moonshot pIC50) RELATIVE-affinity model for
the CliffForge proof of concept.

This is the same recipe as moonshot_ddg.py (internship research repository) (whole-molecule Siamese D-MPNN,
SmoothL1 on pIC50 differences, within-series TRAIN pairs only, 70 epochs, hidden 128,
depth 4, lr 5e-4) but it (a) keeps the SMILES so predictions are auditable and
(b) SAVES the checkpoint (the internship script saved only aggregate metrics).

Only compounds with split == "train" in data/moonshot/all.csv are ever used for
training. Test labels are used only for evaluation in cliffforge_cliffs.py.

CLI (train one seed):
    python -m cliffforge.relmodel train --seed 0 [--epochs 70] [--data-dir D] [--model-dir M]
Library:
    RelScorer().pair_delta(smiles_a, smiles_b) -> mean, sd over seeds of predicted
        pIC50(a) - pIC50(b)   (antisymmetrised: (f(a,b) - f(b,a)) / 2)
    RelScorer().anchored(smiles, ref_smiles, ref_y) -> reference-anchored pIC50
"""
import argparse
import csv
import glob
import os
import sys
import time
from collections import defaultdict

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from torch_geometric.data import Batch, Data
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
from . import paths as P
from .features import mol_to_whole_graph, NODE_DIM, EDGE_DIM
from .networks import SiameseAVES
from .pairs import PairSet, collate, make_pairs
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def smiles_to_data(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    g = mol_to_whole_graph(mol)
    if g is None:
        return None
    return Data(x=torch.from_numpy(g["x"]).float(),
                edge_index=torch.from_numpy(g["edge_index"]).long(),
                edge_attr=torch.from_numpy(g["edge_attr"]).float())


def load_items():
    items = []
    with open(P.processed_csv(), newline="") as f:
        for d in csv.DictReader(f):
            data = smiles_to_data(d["smiles"])
            if data is None:
                continue
            items.append({"data": data, "y": float(d["pIC50"]), "series": d["series"],
                          "split": d["split"], "smiles": d["smiles"]})
    return items


def train(seed, epochs=70, model_dir=None):
    model_dir = str(model_dir or P.model_dir())
    items = load_items()
    torch.manual_seed(seed)
    np.random.seed(seed)
    pairs = make_pairs(items, seed)
    assert all(items[i]["split"] == "train" and items[j]["split"] == "train" for i, j in pairs)
    loader = DataLoader(PairSet(items, pairs), batch_size=32, shuffle=True, collate_fn=collate)
    model = SiameseAVES(NODE_DIM, EDGE_DIM, 128, 4, 0.1).to(DEVICE)
    opt = torch.optim.Adam(model.parameters(), lr=5e-4, weight_decay=1e-5)
    lossf = nn.SmoothL1Loss()
    print(f"seed={seed} device={DEVICE} molecules={len(items)} train_pairs={len(pairs)} "
          f"batches/epoch={len(loader)}", flush=True)
    t0 = time.time()
    for ep in range(epochs):
        model.train()
        tot, n = 0.0, 0
        for c, r, y in loader:
            c, r, y = c.to(DEVICE), r.to(DEVICE), y.to(DEVICE)
            opt.zero_grad()
            loss = lossf(model(c, r), y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            tot += float(loss) * len(y)
            n += len(y)
        if ep % 5 == 0 or ep == epochs - 1:
            print(f"  ep{ep:03d} loss={tot / n:.4f} t={time.time() - t0:.0f}s", flush=True)
    os.makedirs(model_dir, exist_ok=True)
    path = os.path.join(model_dir, f"rel_siamese_s{seed}.pt")
    torch.save({"state_dict": model.state_dict(), "seed": seed, "epochs": epochs,
                "n_train_pairs": len(pairs), "recipe": "internship moonshot_ddg.py (whole-molecule Siamese D-MPNN)",
                "train_split": "data/moonshot/all.csv split==train only"}, path)
    print(f"saved {path}", flush=True)


class RelScorer:
    """Ensemble of saved Siamese models; all inference on CPU or GPU as available."""

    def __init__(self, model_dir=None, device=None):
        model_dir = str(model_dir or P.model_dir())
        self.device = device or DEVICE
        self.models = []
        self.paths = sorted(glob.glob(os.path.join(model_dir, "rel_siamese_s*.pt")))
        for p in self.paths:
            m = SiameseAVES(NODE_DIM, EDGE_DIM, 128, 4, 0.1).to(self.device)
            m.load_state_dict(torch.load(p, map_location=self.device)["state_dict"])
            m.eval()
            self.models.append(m)
        if not self.models:
            raise FileNotFoundError(f"no rel_siamese_s*.pt in {model_dir}")
        self._cache = {}

    def _graph(self, smi):
        if smi not in self._cache:
            self._cache[smi] = smiles_to_data(smi)
        return self._cache[smi]

    @torch.no_grad()
    def _raw(self, model, ga, gb, bs=128):
        out = []
        for i in range(0, len(ga), bs):
            a = Batch.from_data_list(ga[i:i + bs]).to(self.device)
            b = Batch.from_data_list(gb[i:i + bs]).to(self.device)
            out.append(model(a, b).cpu().numpy())
        return np.concatenate(out) if out else np.zeros(0)

    def pair_delta_all(self, smiles_a, smiles_b):
        """[n_models, n_pairs] antisymmetrised predicted pIC50(a) - pIC50(b)."""
        ga = [self._graph(s) for s in smiles_a]
        gb = [self._graph(s) for s in smiles_b]
        return np.stack([(self._raw(m, ga, gb) - self._raw(m, gb, ga)) / 2.0 for m in self.models])

    def pair_delta(self, smiles_a, smiles_b):
        d = self.pair_delta_all(smiles_a, smiles_b)
        return d.mean(0), d.std(0)

    def anchored(self, smiles, ref_smiles, ref_y):
        """Reference-anchored pIC50 (mean over refs of ref_y + f(query, ref)), per model."""
        ref_y = np.asarray(ref_y, dtype=float)
        gr = [self._graph(s) for s in ref_smiles]
        res = np.zeros((len(self.models), len(smiles)))
        for k, s in enumerate(smiles):
            gq = [self._graph(s)] * len(gr)
            for mi, m in enumerate(self.models):
                res[mi, k] = float(np.mean(ref_y + self._raw(m, gq, gr)))
        return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["train"])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--epochs", type=int, default=70)
    P.add_path_args(ap)
    a = ap.parse_args()
    P.apply_path_args(a)
    train(a.seed, a.epochs, a.model_dir)
