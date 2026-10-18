"""
moonshot_prep.py
================
Curate the Moonshot noncovalent Mpro set into an ML table with reproducible
congeneric series (Bemis-Murcko) and a within-series split, and estimate the
assay-noise ceiling from replicate measurements.

Outputs:
  data/moonshot/all.csv      (smiles,pIC50,series,series_size,split,n_rep,rep_sd)
  data/moonshot/train.csv    (smiles,pIC50)   <- for real ChemProp
  data/moonshot/test.csv     (smiles,pIC50)
  moonshot_prep.json         (noise ceiling, split stats; written next to all.csv)
"""
import csv
import json
import math
import os
from collections import defaultdict

import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from . import paths as P

OUT = None      # set in main()
RESULTS = None  # set in main()
NONCOV = None   # set in main()
MIN_SERIES = 8
TEST_FRAC = 0.35
SEED = 17


def canon2d(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None
    flat = Chem.MolToSmiles(m, isomericSmiles=False)
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(flat))
        scaf_smi = Chem.MolToSmiles(scaf)
    except Exception:
        scaf_smi = None
    return flat, (scaf_smi or "NONE")


def main(argv=None):
    import argparse
    global OUT, RESULTS, NONCOV
    ap = argparse.ArgumentParser(description="Curate the Moonshot non-covalent Mpro table (series + within-series split)")
    P.add_path_args(ap)
    ap.add_argument("--raw", help="path to cdd_noncovalent_dates_2023_10_18_filt.csv (default <data-dir>/raw/...)")
    ap.add_argument("--out-dir", help="where to write all.csv/train.csv/test.csv (default <data-dir>/moonshot)")
    args = ap.parse_args(argv)
    P.apply_path_args(args)
    NONCOV = args.raw or str(P.raw_csv())
    OUT = args.out_dir or str(P.processed_dir())
    RESULTS = OUT
    if not os.path.exists(NONCOV):
        raise SystemExit(f"raw Moonshot file not found: {NONCOV}  (see data/README.md)")
    os.makedirs(OUT, exist_ok=True)
    # load & group replicate measurements by 2D structure
    by_struct = defaultdict(list)
    struct_scaf = {}
    with open(NONCOV, newline="") as f:
        r = csv.DictReader(f)
        for d in r:
            smi = (d.get("suspected_SMILES") or "").strip()
            ic50 = d.get("IC50_(µM)") or d.get("IC50_(\u00b5M)")
            try:
                ic50 = float(ic50)
            except (TypeError, ValueError):
                continue
            if ic50 <= 0:
                continue
            flat, scaf = canon2d(smi)
            if flat is None:
                continue
            by_struct[flat].append(6.0 - math.log10(ic50))
            struct_scaf[flat] = scaf

    # dedup -> mean pIC50, replicate sd
    records = []
    rep_sds = []
    for flat, vals in by_struct.items():
        v = np.array(vals)
        records.append({"smiles": flat, "pIC50": float(v.mean()),
                        "scaf": struct_scaf[flat], "n_rep": len(v),
                        "rep_sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0})
        if len(v) > 1:
            rep_sds.append(v.std(ddof=1))

    # assign series by scaffold size
    scaf_members = defaultdict(list)
    for rec in records:
        scaf_members[rec["scaf"]].append(rec)
    sizes = {s: len(m) for s, m in scaf_members.items()}

    rng = np.random.RandomState(SEED)
    n_series_used = 0
    for s, members in scaf_members.items():
        size = len(members)
        if size >= MIN_SERIES:
            n_series_used += 1
            idx = rng.permutation(size)
            n_test = max(3, int(round(size * TEST_FRAC)))
            n_test = min(n_test, size - 3)
            test_ids = set(idx[:n_test].tolist())
            for i, rec in enumerate(members):
                rec["series"] = s
                rec["series_size"] = size
                rec["split"] = "test" if i in test_ids else "train"
        else:
            for rec in members:
                rec["series"] = s if size > 1 else "singleton"
                rec["series_size"] = size
                rec["split"] = "train"  # background; not evaluated within-series

    # write all.csv
    with open(os.path.join(OUT, "all.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["smiles", "pIC50", "series", "series_size", "split", "n_rep", "rep_sd"])
        for rec in records:
            w.writerow([rec["smiles"], f"{rec['pIC50']:.4f}", rec["series"],
                        rec["series_size"], rec["split"], rec["n_rep"], f"{rec['rep_sd']:.4f}"])

    # ChemProp train/test
    for split in ("train", "test"):
        with open(os.path.join(OUT, f"{split}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["smiles", "pIC50"])
            for rec in records:
                if rec["split"] == split:
                    w.writerow([rec["smiles"], f"{rec['pIC50']:.4f}"])

    n_train = sum(1 for r in records if r["split"] == "train")
    n_test = sum(1 for r in records if r["split"] == "test")
    mean_rep_sd = float(np.mean(rep_sds)) if rep_sds else 0.0
    # reliability-based correlation ceiling (Spearman/Pearson attenuation)
    all_p = np.array([r["pIC50"] for r in records])
    sig_var = float(all_p.var())
    noise_var = mean_rep_sd ** 2
    reliability = sig_var / (sig_var + noise_var) if (sig_var + noise_var) > 0 else float("nan")
    ceiling = math.sqrt(max(reliability, 0))

    rep = {
        "n_unique_2d": len(records), "n_train": n_train, "n_test": n_test,
        "n_series_used(ge%d)" % MIN_SERIES: n_series_used,
        "series_sizes_top": sorted(sizes.values(), reverse=True)[:15],
        "n_replicated_compounds": len(rep_sds),
        "mean_replicate_sd_pIC50": round(mean_rep_sd, 3),
        "global_pIC50_sd": round(float(all_p.std()), 3),
        "reliability": round(reliability, 3),
        "global_corr_ceiling_pearson": round(ceiling, 3),
    }
    json.dump(rep, open(os.path.join(RESULTS, "moonshot_prep.json"), "w"), indent=2)
    print("Curation:", json.dumps(rep, indent=2))
    print(f"\nWrote {OUT}/all.csv, train.csv ({n_train}), test.csv ({n_test})")
    print(f"Series used (>= {MIN_SERIES}): {n_series_used}")
    print(f"Assay noise: mean replicate SD = {mean_rep_sd:.3f} pIC50 over "
          f"{len(rep_sds)} replicated compounds; global-corr ceiling ~ {ceiling:.3f}")


if __name__ == "__main__":
    main()
