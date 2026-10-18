"""Thin wrappers around the ChemProp 1.x command-line tools (CPU inference with an existing checkpoint)."""
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from . import paths as P


def predict_all(smiles, checkpoint_dir=None, cache=None, threads=4):
    """Predict pIC50 for `smiles` with the ChemProp checkpoint(s) in `checkpoint_dir` (default models/chemprop_moonshot).

    Results are cached in outputs/chemprop/all_preds.csv and re-used when the SMILES list is identical."""
    smiles = list(smiles)
    cache = Path(cache) if cache else P.output_dir() / "chemprop" / "all_preds.csv"
    if cache.exists():
        old = pd.read_csv(cache)
        if list(old["smiles"]) == smiles:
            return old["pIC50"].values.astype(float)
    ckpt = Path(checkpoint_dir) if checkpoint_dir else P.chemprop_dir()
    if not ckpt.exists():
        raise FileNotFoundError(f"ChemProp checkpoint directory not found: {ckpt} (see models/README.md)")
    cache.parent.mkdir(parents=True, exist_ok=True)
    inp = cache.parent / "input_smiles.csv"
    pd.DataFrame({"smiles": smiles}).to_csv(inp, index=False)
    env = dict(os.environ, CUDA_VISIBLE_DEVICES="", OMP_NUM_THREADS=str(threads))
    subprocess.run([P.bin_path("chemprop_predict"), "--test_path", str(inp), "--checkpoint_dir", str(ckpt),
                    "--preds_path", str(cache), "--no_cuda"], check=True, env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out = pd.read_csv(cache)
    assert list(out["smiles"]) == smiles, "ChemProp returned a different molecule order"
    return out["pIC50"].values.astype(float)


def train_fold(train_csv, save_dir, epochs=60, threads=2):
    """Same ChemProp command as the baseline checkpoint (regression, 60 epochs, defaults otherwise)."""
    env = dict(os.environ, CUDA_VISIBLE_DEVICES="", OMP_NUM_THREADS=str(threads))
    subprocess.run([P.bin_path("chemprop_train"), "--data_path", str(train_csv), "--dataset_type", "regression",
                    "--smiles_columns", "smiles", "--target_columns", "pIC50", "--save_dir", str(save_dir),
                    "--epochs", str(epochs), "--num_folds", "1", "--metric", "rmse", "--extra_metrics", "r2",
                    "--no_cuda", "--num_workers", "0", "--quiet"], check=True, env=env)


def predict_csv(test_csv, checkpoint_dir, preds_csv, threads=2):
    env = dict(os.environ, CUDA_VISIBLE_DEVICES="", OMP_NUM_THREADS=str(threads))
    subprocess.run([P.bin_path("chemprop_predict"), "--test_path", str(test_csv), "--checkpoint_dir", str(checkpoint_dir),
                    "--preds_path", str(preds_csv), "--no_cuda"], check=True, env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
