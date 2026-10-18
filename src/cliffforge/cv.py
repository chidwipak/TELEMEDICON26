"""
cliffforge.cv -- out-of-fold VALIDATION predictions for cliff-threshold selection.
=====================================================================================
Uses ONLY the 961-compound TRAIN partition of data/moonshot/all.csv (the 141 test
compounds are never loaded for labels, training or prediction here).

5-fold random CV over the train partition (KFold shuffle, random_state=0). For each fold F:
  * relative Siamese model: 3 seeds (0,1,2), same recipe as cliffforge.relmodel
    (whole-molecule Siamese D-MPNN, 70 epochs, same-series train pairs of the fold-train set)
  * ChemProp absolute model: same command as the existing models/chemprop_moonshot
    (chemprop_train regression, 60 epochs, defaults otherwise, internal random 0.8/0.1/0.1)
  trained on the other 4 folds, then used to predict the held-out fold.

Commands
  python -m cliffforge.cv rel --fold F --seed S
  python -m cliffforge.cv chemprop --fold F
Outputs: outputs/cv/fold{F}/{rel_siamese_s{S}.pt, chemprop/, chemprop_preds.csv}
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd

from . import paths as P
from . import chemprop_io
N_FOLDS = 5
CV_SEED = 0


def cv_dir():
    return P.output_dir() / "cv"


def fold_assignment():
    """fold id (0..4) for every TRAIN-partition row of all.csv; -1 for test-partition rows."""
    from sklearn.model_selection import KFold
    df = pd.read_csv(P.processed_csv())
    tr = np.where(df.split == "train")[0]
    fold = np.full(len(df), -1, dtype=int)
    for f, (_, val) in enumerate(KFold(N_FOLDS, shuffle=True, random_state=CV_SEED).split(tr)):
        fold[tr[val]] = f
    return df, fold


def train_rel(fold_id, seed, epochs=70):
    import torch
    import torch.nn as nn
    from torch.utils.data import DataLoader
    torch.set_num_threads(1)
    from .relmodel import load_items, DEVICE, NODE_DIM, EDGE_DIM, SiameseAVES
    from .pairs import PairSet, collate, make_pairs
    df, fold = fold_assignment()
    items = load_items()
    assert len(items) == len(df)
    for i, it in enumerate(items):  # fold-train = other train folds; test partition + val fold excluded
        it["split"] = "train" if (fold[i] >= 0 and fold[i] != fold_id) else "heldout"
    torch.manual_seed(seed)
    np.random.seed(seed)
    pairs = make_pairs(items, seed)
    assert all(items[i]["split"] == "train" and items[j]["split"] == "train" for i, j in pairs)
    loader = DataLoader(PairSet(items, pairs), batch_size=32, shuffle=True, collate_fn=collate)
    model = SiameseAVES(NODE_DIM, EDGE_DIM, 128, 4, 0.1).to(DEVICE)
    opt = torch.optim.Adam(model.parameters(), lr=5e-4, weight_decay=1e-5)
    lossf = nn.SmoothL1Loss()
    print(f"fold={fold_id} seed={seed} train_mols={sum(it['split'] == 'train' for it in items)} pairs={len(pairs)}", flush=True)
    for ep in range(epochs):
        model.train()
        for c, r, y in loader:
            c, r, y = c.to(DEVICE), r.to(DEVICE), y.to(DEVICE)
            opt.zero_grad()
            loss = lossf(model(c, r), y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
        if ep % 10 == 0 or ep == epochs - 1:
            print(f"  ep{ep:03d} loss={float(loss):.4f}", flush=True)
    d = os.path.join(str(cv_dir()), f"fold{fold_id}")
    os.makedirs(d, exist_ok=True)
    torch.save({"state_dict": model.state_dict(), "seed": seed, "fold": fold_id, "epochs": epochs},
               os.path.join(d, f"rel_siamese_s{seed}.pt"))
    print("saved", d, flush=True)


def train_chemprop(fold_id, epochs=60):
    df, fold = fold_assignment()
    d = os.path.join(str(cv_dir()), f"fold{fold_id}")
    os.makedirs(d, exist_ok=True)
    tr = df[(fold >= 0) & (fold != fold_id)][["smiles", "pIC50"]]
    tr.to_csv(os.path.join(d, "train.csv"), index=False)
    df[fold >= 0][["smiles"]].to_csv(os.path.join(d, "all_train_smiles.csv"), index=False)
    chemprop_io.train_fold(os.path.join(d, "train.csv"), os.path.join(d, "chemprop"), epochs=epochs)
    chemprop_io.predict_csv(os.path.join(d, "all_train_smiles.csv"), os.path.join(d, "chemprop"),
                            os.path.join(d, "chemprop_preds.csv"))
    print("chemprop fold done", fold_id, flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["rel", "chemprop"])
    ap.add_argument("--fold", type=int, required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--epochs", type=int, default=None, help="override epochs (default: 70 relative / 60 ChemProp); for smoke tests")
    P.add_path_args(ap)
    a = ap.parse_args()
    P.apply_path_args(a)
    if a.cmd == "rel":
        train_rel(a.fold, a.seed, a.epochs or 70)
    else:
        train_chemprop(a.fold, a.epochs or 60)
