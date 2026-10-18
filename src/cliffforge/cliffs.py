"""
cliffforge_cliffs.py  --  CliffForge PoC, Experiment A
=====================================================
Activity-cliff analysis on REAL experimental COVID Moonshot SARS-CoV-2 Mpro pIC50 data
and an honest evaluation of whether affinity models get the DIRECTION of close-analogue
differences right.

Data / split (unchanged from the internship work): data/moonshot/all.csv
    1102 unique 2D structures, split column train=961 / test=141 (within-series split).
Cliff definition (fixed BEFORE any model was scored):
    PRIMARY      Morgan(r=2, 2048) Tanimoto >= 0.80  AND  |delta pIC50| >= 1.0
    SENSITIVITY  Tanimoto in {0.6, 0.7}, same |delta| >= 1.0  (reported, never swapped in)
Pair partitions (leakage control):
    test_test   both members in the held-out test partition  -> STRICT held-out evaluation
    test_train  one member in train                          -> EXPLORATORY (the train member's
                label was seen by the models; the pair itself was not)
    train_train both in train                                -> descriptive only, not scored
Predictors compared on each pair (predicted delta = pred(a) - pred(b)):
    chemprop_abs   existing models/chemprop_moonshot (absolute pIC50 D-MPNN, ChemProp 1.6.1)
    rel_siamese    3-seed ensemble of the whole-molecule Siamese relative model
                   (cliffforge.relmodel), direct pairwise, antisymmetrised
    rel_anchored   same models, reference-anchored absolute estimate (test_test only)
    knn_sim        similarity-weighted 5-NN over TRAIN labels (non-learning baseline)
    heavy_atoms    "bigger molecule is more potent" property baseline
    (random = 0.5 direction accuracy)
No FEP+ checkpoint is used for Moonshot (different target, different endpoint / units).
"""
import json
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from scipy.stats import binomtest, spearmanr

RDLogger.DisableLog("rdApp.*")
from . import paths as P
from . import chemprop_io

OUT = None  # set in main(): <output-dir>/cliff_pairs
SEED = 0
SIM_PRIMARY, DELTA_PRIMARY = 0.80, 1.0
SIM_SENS = (0.6, 0.7, 0.8)
MIN_SIM_STORE = 0.6
MIN_TEST_SERIES = 5


def morgan(smi):
    return AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), 2, 2048)


def knn_predict(S, y, is_train, k=5):
    """Similarity-weighted kNN over TRAIN labels; leave-one-out for train molecules."""
    n = len(y)
    pred = np.zeros(n)
    tr_idx = np.where(is_train)[0]
    for i in range(n):
        sims = S[i, tr_idx].copy()
        if is_train[i]:
            sims[tr_idx == i] = -1
        top = np.argsort(-sims)[:k]
        w = np.clip(sims[top], 1e-6, None)
        pred[i] = float(np.sum(w * y[tr_idx[top]]) / np.sum(w))
    return pred


def direction_acc(true_d, pred_d):
    """Fraction of pairs whose predicted sign matches the measured sign; exact ties = 0.5."""
    s_t, s_p = np.sign(true_d), np.sign(pred_d)
    ok = np.where(s_p == 0, 0.5, (s_t == s_p).astype(float))
    return ok


def boot_ci(x, n=2000, seed=SEED):
    rng = np.random.RandomState(seed)
    x = np.asarray(x, dtype=float)
    if len(x) < 2:
        return [None, None]
    m = [x[rng.randint(0, len(x), len(x))].mean() for _ in range(n)]
    return [float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))]


def summarize(true_d, preds, label):
    """preds: dict model -> predicted delta array (NaN allowed)."""
    out = {"n_pairs": int(len(true_d))}
    for name, pd_ in preds.items():
        ok_mask = np.isfinite(pd_)
        if ok_mask.sum() == 0:
            out[name] = None
            continue
        ok = direction_acc(true_d[ok_mask], pd_[ok_mask])
        rec = {"n_scored": int(ok_mask.sum()), "direction_accuracy": float(ok.mean()),
               "n_correct": float(ok.sum()), "ci95_pair_bootstrap": boot_ci(ok)}
        # exact binomial vs chance (ties counted as 0.5 -> round for the test only)
        rec["binom_p_vs_0.5"] = float(binomtest(int(round(ok.sum())), len(ok), 0.5).pvalue)
        if ok_mask.sum() >= 5 and np.std(pd_[ok_mask]) > 0:
            rec["spearman_pred_vs_true_delta"] = float(spearmanr(true_d[ok_mask], pd_[ok_mask])[0])
        out[name] = rec
    # paired exact (McNemar) test, relative Siamese vs ChemProp absolute, non-tied predictions only
    if "rel_siamese" in preds and "chemprop_abs" in preds:
        a_, b_ = preds["rel_siamese"], preds["chemprop_abs"]
        m = np.isfinite(a_) & np.isfinite(b_) & (a_ != 0) & (b_ != 0)
        ra = np.sign(a_[m]) == np.sign(true_d[m])
        rb = np.sign(b_[m]) == np.sign(true_d[m])
        only_a, only_b = int((ra & ~rb).sum()), int((~ra & rb).sum())
        out["mcnemar_rel_vs_chemprop"] = {
            "n_compared": int(m.sum()), "rel_only_correct": only_a, "chemprop_only_correct": only_b,
            "p_two_sided": float(binomtest(only_a, only_a + only_b, 0.5).pvalue) if only_a + only_b else None}
    return out


def within_series(df, pred_col, min_test=MIN_TEST_SERIES):
    """Mean within-series Spearman on the TEST partition (same definition as moonshot_cp_eval.py)."""
    rows = []
    t = df[df.split == "test"]
    for s, g in t.groupby("series"):
        if len(g) >= min_test and np.isfinite(g[pred_col]).all():
            rows.append(spearmanr(g["pIC50"], g[pred_col])[0])
    return {"n_series": len(rows), "within_series_spearman_mean": float(np.mean(rows)) if rows else None,
            "per_series": [float(x) for x in rows]}


def main(argv=None):
    import argparse
    global OUT
    ap = argparse.ArgumentParser(description="Descriptive cliff-pair analysis and within-series ranking (Mpro)")
    P.add_path_args(ap)
    ap.add_argument("--saved-chemprop-preds", help="optional CSV (smiles,pIC50) of earlier ChemProp test predictions to cross-check")
    args = ap.parse_args(argv)
    P.apply_path_args(args)
    OUT = str(P.output_dir() / "cliff_pairs")
    os.makedirs(OUT, exist_ok=True)
    df = pd.read_csv(P.processed_csv())
    n = len(df)
    y = df["pIC50"].values
    is_train = (df["split"] == "train").values
    is_test = ~is_train
    ha = np.array([Chem.MolFromSmiles(s).GetNumHeavyAtoms() for s in df.smiles])
    fps = [morgan(s) for s in df.smiles]
    S = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        S[i] = DataStructs.BulkTanimotoSimilarity(fps[i], fps)

    # ---------------- per-molecule predictions ----------------
    print("ChemProp (existing checkpoint) on all molecules ...", flush=True)
    cp_all = chemprop_io.predict_all(list(df.smiles))
    repro = None
    if args.saved_chemprop_preds:
        saved = pd.read_csv(args.saved_chemprop_preds).set_index("smiles")["pIC50"]
        te_idx = np.where(is_test)[0]
        saved_vec = np.array([saved.get(df.smiles[i], np.nan) for i in te_idx])
        repro = {"n_test_compared": int(np.isfinite(saved_vec).sum()),
                 "max_abs_diff_vs_saved_preds": float(np.nanmax(np.abs(saved_vec - cp_all[te_idx])))}
        print("  reproducibility vs saved ChemProp predictions:", repro, flush=True)
    knn = knn_predict(S, y, is_train)

    from .relmodel import RelScorer
    import torch
    scorer = RelScorer(device=torch.device("cpu"))
    n_models = len(scorer.models)
    print(f"Relative model ensemble: {n_models} checkpoints", flush=True)

    # ---------------- candidate pairs (sim >= 0.6) ----------------
    iu, ju = np.triu_indices(n, 1)
    keep = S[iu, ju] >= MIN_SIM_STORE
    iu, ju = iu[keep], ju[keep]
    part = np.where(is_test[iu] & is_test[ju], "test_test",
                    np.where(is_train[iu] & is_train[ju], "train_train", "test_train"))
    # score only pairs with >=1 test member (train_train is descriptive only)
    sc = part != "train_train"
    I, J = iu[sc], ju[sc]
    print(f"pairs sim>=0.6: {len(iu)}  (scored: {len(I)})", flush=True)
    rel_all = scorer.pair_delta_all(list(df.smiles[I]), list(df.smiles[J]))  # [models, pairs]

    # reference-anchored estimates for TEST molecules (same-series TRAIN refs)
    anch = np.full((n_models, n), np.nan)
    by_ser = defaultdict(list)
    for i in np.where(is_train)[0]:
        by_ser[df.series[i]].append(i)
    for s, g in df[is_test].groupby("series"):
        refs = by_ser.get(s, [])
        if not refs:
            continue
        q = list(g.index)
        anch[:, q] = scorer.anchored(list(df.smiles[q]), list(df.smiles[refs]), y[refs])
    anch_mean = anch.mean(0)
    anch_per_model = [anch[m] for m in range(n_models)]

    # ---------------- assemble the per-pair table ----------------
    rows = pd.DataFrame({
        "idx_a": I, "idx_b": J, "smiles_a": df.smiles[I].values, "smiles_b": df.smiles[J].values,
        "series_a": df.series[I].values, "series_b": df.series[J].values,
        "same_series": (df.series[I].values == df.series[J].values),
        "split_a": df.split[I].values, "split_b": df.split[J].values, "partition": part[sc],
        "tanimoto": S[I, J], "pIC50_a_measured": y[I], "pIC50_b_measured": y[J],
        "true_delta_a_minus_b": y[I] - y[J], "abs_true_delta": np.abs(y[I] - y[J]),
        "pred_delta_chemprop_abs": cp_all[I] - cp_all[J],
        "pred_delta_rel_siamese_mean": rel_all.mean(0), "pred_delta_rel_siamese_sd": rel_all.std(0),
        **{f"pred_delta_rel_seed{m}": rel_all[m] for m in range(n_models)},
        "pred_delta_rel_anchored": anch_mean[I] - anch_mean[J],
        "pred_delta_knn_sim": knn[I] - knn[J],
        "pred_delta_heavy_atoms": (ha[I] - ha[J]).astype(float),
    })
    rows["is_cliff_primary"] = (rows.tanimoto >= SIM_PRIMARY) & (rows.abs_true_delta >= DELTA_PRIMARY)
    rows["pair_class"] = np.select(
        [(rows.tanimoto >= SIM_PRIMARY) & (rows.abs_true_delta >= DELTA_PRIMARY),
         (rows.tanimoto >= SIM_PRIMARY)],
        ["activity_cliff", "similar_non_cliff"], default="sensitivity_only(sim<0.8)")
    rows["eligible_for_strict_heldout"] = rows.partition == "test_test"
    rows.to_csv(os.path.join(OUT, "activity_cliff_pairs.csv"), index=False, float_format="%.4f")

    # ---------------- dataset-level cliff counts ----------------
    dy = np.abs(y[:, None] - y[None, :])
    counts = {}
    for sim_t in (0.6, 0.7, 0.8, 0.9):
        for d_t in (0.5, 1.0, 1.5):
            m = (S[iu, ju] >= sim_t) & (dy[iu, ju] >= d_t)
            counts[f"sim>={sim_t}|d>={d_t}"] = {
                "all": int(m.sum()), "test_test": int((m & (part == "test_test")).sum()),
                "test_train": int((m & (part == "test_train")).sum()),
                "train_train": int((m & (part == "train_train")).sum())}
    n_sim_pairs = {sim_t: int((S[iu, ju] >= sim_t).sum()) for sim_t in (0.6, 0.7, 0.8)}

    # ---------------- model comparison ----------------
    models = {"chemprop_abs": "pred_delta_chemprop_abs", "rel_siamese": "pred_delta_rel_siamese_mean",
              "knn_sim": "pred_delta_knn_sim", "heavy_atoms": "pred_delta_heavy_atoms",
              **{f"rel_single_seed{m}": f"pred_delta_rel_seed{m}" for m in range(n_models)}}
    results = {}

    def evaluate(mask, with_anchored):
        sub = rows[mask]
        t = sub.true_delta_a_minus_b.values
        pr = {k: sub[c].values.astype(float) for k, c in models.items()}
        if with_anchored:
            pr["rel_anchored"] = sub["pred_delta_rel_anchored"].values.astype(float)
        return summarize(t, pr, "")

    for part_name, with_anch in (("test_test", True), ("test_train", False)):
        for sim_t in SIM_SENS:
            key = f"{part_name}|cliff(sim>={sim_t},|d|>={DELTA_PRIMARY})"
            m = (rows.partition == part_name) & (rows.tanimoto >= sim_t) & (rows.abs_true_delta >= DELTA_PRIMARY)
            results[key] = evaluate(m, with_anch)
        # primary-similarity strata by |delta| (cliff vs non-cliff vs within-noise)
        for lo, hi, nm in ((0.0, 0.4, "within_assay_noise(|d|<0.4)"), (0.4, 1.0, "moderate(0.4<=|d|<1.0)"),
                           (1.0, 99, "cliff(|d|>=1.0)")):
            m = (rows.partition == part_name) & (rows.tanimoto >= SIM_PRIMARY) & \
                (rows.abs_true_delta >= lo) & (rows.abs_true_delta < hi) & (rows.abs_true_delta > 0)
            results[f"{part_name}|sim>=0.8|{nm}"] = evaluate(m, with_anch)
        m = (rows.partition == part_name) & (rows.tanimoto >= SIM_PRIMARY) & (rows.abs_true_delta > 0)
        results[f"{part_name}|sim>=0.8|all"] = evaluate(m, with_anch)
        m = (rows.partition == part_name) & (rows.tanimoto >= 0.6) & (rows.abs_true_delta >= 0.5)
        results[f"{part_name}|sim>=0.6|d>=0.5(exploratory-large-n)"] = evaluate(m, with_anch)
    # everything touching the test partition at the primary cliff definition (exploratory pool)
    m = (rows.partition != "train_train") & rows.is_cliff_primary
    results["any_test_member|PRIMARY_cliff(pooled,exploratory)"] = evaluate(m, False)

    # ---------------- within-series check (reproduce internship numbers) ----------------
    df_ws = df.copy()
    df_ws["cp"] = cp_all
    df_ws["rel_anch_mean"] = anch_mean
    ws = {"chemprop_abs": within_series(df_ws, "cp"),
          "rel_anchored_ensemble": within_series(df_ws, "rel_anch_mean")}
    for m in range(n_models):
        df_ws[f"rel_anch_s{m}"] = anch_per_model[m]
        ws[f"rel_anchored_seed{m}"] = within_series(df_ws, f"rel_anch_s{m}")
    from scipy.stats import wilcoxon
    a_, b_ = ws["rel_anchored_ensemble"]["per_series"], ws["chemprop_abs"]["per_series"]
    ws["paired_wilcoxon_rel_ensemble_vs_chemprop"] = {
        "n_series": len(a_), "mean_diff": float(np.mean(np.array(a_) - np.array(b_))),
        "p": float(wilcoxon(a_, b_).pvalue) if len(a_) >= 6 else None}
    ws["reference_reported_in_repo"] = {"chemprop_abs (internship report)": 0.192,
                                       "whole_ddg 3 seeds (internship report)": 0.555}

    # per-test-molecule global metrics (context)
    t_idx = np.where(is_test)[0]
    glob = {"chemprop_abs_global_spearman_test": float(spearmanr(y[t_idx], cp_all[t_idx])[0]),
            "rel_anchored_global_spearman_test": float(spearmanr(y[t_idx][np.isfinite(anch_mean[t_idx])],
                                                                 anch_mean[t_idx][np.isfinite(anch_mean[t_idx])])[0])}

    summary = {
        "dataset": "data/moonshot/all.csv (COVID Moonshot noncovalent Mpro, pIC50 = 6 - log10(IC50 uM), "
                   "replicates averaged; 2D, stereochemistry removed)",
        "n_molecules": int(n), "n_train": int(is_train.sum()), "n_test": int(is_test.sum()),
        "n_test_series_with_>=5": ws["chemprop_abs"]["n_series"],
        "thresholds": {"primary": {"tanimoto": SIM_PRIMARY, "abs_delta_pIC50": DELTA_PRIMARY},
                       "sensitivity_tanimoto": list(SIM_SENS), "fingerprint": "Morgan r=2, 2048 bits",
                       "assay_noise_sd_pIC50": 0.396},
        "n_pairs_by_similarity": n_sim_pairs,
        "cliff_counts_grid": counts,
        "n_scored_pairs_sim>=0.6_with_test_member": int(len(rows)),
        "chemprop_reproducibility_check": repro,
        "relative_model": {"n_checkpoints": n_models, "paths": [os.path.relpath(p, str(P.REPO_ROOT)) for p in scorer.paths],
                           "prediction": "antisymmetrised pairwise delta, mean over checkpoints"},
        "pair_evaluation": results,
        "within_series_test": ws, "global_test": glob,
        "seed": SEED,
        "notes": ["test_test = strict held-out (both members in test); test_train = exploratory (train member's "
                  "label was used for training both learned models).",
                  "Pair-bootstrap CIs treat pairs as independent; pairs share molecules, so CIs are optimistic.",
                  "Stereochemistry was removed in data prep: enantiomer cliffs cannot be detected."],
    }
    json.dump(summary, open(os.path.join(OUT, "activity_cliff_metrics.json"), "w"), indent=2)
    print(json.dumps({k: summary[k] for k in ("n_pairs_by_similarity", "chemprop_reproducibility_check", "global_test")}, indent=1))
    for k, v in results.items():
        line = f"{k}: n={v['chemprop_abs']['n_scored'] if v.get('chemprop_abs') else 0}"
        for m in ("chemprop_abs", "rel_siamese", "rel_anchored", "knn_sim", "heavy_atoms"):
            if v.get(m):
                line += f" | {m}={v[m]['direction_accuracy']:.2f}"
        print(line)
    print("within-series:", {k: (v.get("within_series_spearman_mean") if isinstance(v, dict) else v) for k, v in ws.items()})


if __name__ == "__main__":
    main()
