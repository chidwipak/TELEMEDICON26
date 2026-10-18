"""
cliffforge_cliff_detection.py -- CliffForge PoC: DIRECT activity-cliff DETECTION evaluation
==========================================================================================
Question: does the relative-affinity model detect activity cliffs better than the existing
ChemProp Mpro model?  (Not: does it rank well, nor: how many candidates did we generate.)

DECLARED BEFORE ANY DETECTION METRIC WAS COMPUTED
-------------------------------------------------
Data / split   data/moonshot/all.csv, 961 train / 141 test (unchanged). pIC50 = 6 - log10 IC50[uM].
Cliff (true)   Morgan(r=2, 2048) Tanimoto >= 0.80  AND  |pIC50(a) - pIC50(b)| >= 1.0  (PRIMARY).
               Sensitivity only: Tanimoto >= 0.70 / 0.60 with the same |delta| (never promoted).
Predicted cliff  |predicted delta| >= tau.
               tau_val  = F1-maximising tau on OUT-OF-FOLD validation pairs from 5-fold CV inside the TRAIN
                          partition (cliffforge.cv; train-train pairs with Tanimoto >= 0.80), chosen per
                          model, tie -> closest to 1.0.  Test labels are never used for tau.
               tau_fixed = 1.0 (transparent reference).
Setting A      STRICT held-out pairs: both compounds in the test partition.  Predictions use structures only.
Setting B      LEAD-OPTIMISATION: one test compound (query) + one TRAIN reference.  The query's activity is
               never used; the reference's measured activity MAY be used.  Exploratory; NOT both-unseen.
Predicted delta (signed, query - reference or a - b)
   relative model   f(query, ref) ensemble mean (3 seeds), antisymmetrised.  Reference activity not needed.
   ChemProp (A)     chemprop(a) - chemprop(b)
   ChemProp (B)     chemprop(query) - measured_pIC50(ref)   <- uses the available reference activity
                    (the most favourable fair protocol for the absolute baseline); the pred-pred variant
                    chemprop(query) - chemprop(ref) is reported as secondary.
   5-NN baseline    similarity-weighted 5-NN over TRAIN labels, same protocols as ChemProp.
   majority class   "never a cliff"; chance AP = prevalence.
Uncertainty    cluster bootstrap (2000 resamples): A -> connected components of the pair graph;
               B -> by query (test) molecule (series-level variant also reported). Paired differences use the
               same resamples.  Exact Clopper-Pearson CI for recall.
DECISION RULE (outcome chosen mechanically from the metrics below)
   Setting A has far fewer than 10 positives -> descriptive only.
   (A) Preliminary evidence of better detection  iff ALL:
        1. B-primary: AP(rel) > AP(ChemProp) and the 95% cluster-bootstrap CI of delta-AP excludes 0;
        2. B-primary: F1(rel) > F1(ChemProp) at BOTH tau_val and tau=1.0;
        3. A-primary: AP(rel) >= AP(ChemProp) and F1(rel) >= F1(ChemProp) at both thresholds (point estimates).
   (B) No demonstrated improvement  iff B-primary delta-AP <= 0, or delta-F1 <= 0 at both thresholds.
   (C) Inconclusive otherwise.
Leakage guard: predictions are generated with y_known, a copy of the activities in which every held-out
compound is NaN; any use of a held-out label would produce NaN and trip an assertion.
"""
import json
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from scipy.stats import beta
from sklearn.metrics import average_precision_score, roc_auc_score

RDLogger.DisableLog("rdApp.*")
from . import paths as PATHS
from . import chemprop_io

OUT = None   # set in main(): <output-dir>/cliff_detection
CV = None    # set in main(): <output-dir>/cv (fold models from cliffforge.cv; optional)
SIM_PRIMARY, DELTA = 0.80, 1.0
SIMS = (0.8, 0.7, 0.6)
TAU_FIXED = 1.0
TAU_GRID = np.round(np.arange(0.05, 2.0001, 0.05), 2)
NBOOT = 2000  # overridable with --n-boot
MODELS = ["chemprop_abs", "rel_siamese", "knn_sim"]


def morgan(s):
    return AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, 2048)


def knn_pred(S, y_known, train_mask, k=5):
    """5-NN over labelled TRAIN molecules (train_mask); a molecule never uses its own label."""
    n = len(y_known)
    out = np.zeros(n)
    idx = np.where(train_mask)[0]
    for i in range(n):
        s = S[i, idx].copy()
        s[idx == i] = -1
        top = np.argsort(-s)[:k]
        w = np.clip(s[top], 1e-6, None)
        v = y_known[idx[top]]
        assert np.isfinite(v).all(), "kNN touched an unavailable label"
        out[i] = float(np.sum(w * v) / np.sum(w))
    return out


# ------------------------------------------------------------------ metrics
def conf(y, yh):
    return (int(((y == 1) & (yh == 1)).sum()), int(((y == 0) & (yh == 1)).sum()),
            int(((y == 0) & (yh == 0)).sum()), int(((y == 1) & (yh == 0)).sum()))


def prf(tp, fp, tn, fn):
    P = tp / (tp + fp) if tp + fp else None
    R = tp / (tp + fn) if tp + fn else None
    F = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else None
    spec = tn / (tn + fp) if tn + fp else None
    n = tp + fp + tn + fn
    den = np.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {"tp": tp, "fp": fp, "tn": tn, "fn": fn, "precision": P, "recall": R, "f1": F, "specificity": spec,
            "accuracy": (tp + tn) / n if n else None,
            "balanced_accuracy": (R + spec) / 2 if (R is not None and spec is not None) else None,
            "mcc": float((tp * tn - fp * fn) / den) if den > 0 else None}


def cp_ci(k, n, a=0.05):
    if n == 0:
        return [None, None]
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a / 2, k + 1, n - k)
    return [float(lo), float(hi)]


def scalar_metrics(y, score, yh):
    """f1/precision/recall + AP/ROC (NaN when undefined)."""
    tp, fp, tn, fn = conf(y, yh)
    r = {"precision": tp / (tp + fp) if tp + fp else np.nan, "recall": tp / (tp + fn) if tp + fn else np.nan,
         "f1": 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else np.nan}
    if y.sum() > 0 and y.sum() < len(y):
        r["ap"] = float(average_precision_score(y, score))
        r["roc_auc"] = float(roc_auc_score(y, score))
    else:
        r["ap"] = r["roc_auc"] = np.nan
    return r


def components(pairs_ab):
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in pairs_ab:
        parent[find(a)] = find(b)
    return [find(a) for a, _ in pairs_ab]


def cluster_boot(cluster_ids, evalfn, nboot=None, seed=0):
    nboot = nboot or NBOOT
    """evalfn(pair_index_array) -> dict of scalars. Returns list of dicts over resamples."""
    rng = np.random.RandomState(seed)
    cl = defaultdict(list)
    for i, c in enumerate(cluster_ids):
        cl[c].append(i)
    keys = list(cl)
    res = []
    for _ in range(nboot):
        pick = rng.randint(0, len(keys), len(keys))
        idx = np.array([j for p in pick for j in cl[keys[p]]], dtype=int)
        res.append(evalfn(idx))
    return res


def ci(vals):
    v = np.array([x for x in vals if x is not None and np.isfinite(x)])
    if len(v) < 20:
        return {"ci95": [None, None], "n_valid_resamples": int(len(v))}
    return {"ci95": [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))], "n_valid_resamples": int(len(v))}


# ------------------------------------------------------------------ predictions
def predict_pairs(df, S, ymask, ch_pred, rel_scorer, pairs, kind, knn):
    """pairs: list of (h, r) = (query/first, reference/second) row indices.
    ymask: activities available at prediction time (held-out compounds are NaN).
    kind[i] in {'AA','BQ'}: 'AA' both unseen (a,b), 'BQ' query + labelled reference."""
    H = np.array([p[0] for p in pairs])
    R = np.array([p[1] for p in pairs])
    rel_all = rel_scorer.pair_delta_all(list(df.smiles[H]), list(df.smiles[R]))  # [models, pairs]
    out = {"rel_siamese": rel_all.mean(0), "rel_sd": rel_all.std(0)}
    cp_pp = ch_pred[H] - ch_pred[R]
    kn_pp = knn[H] - knn[R]
    cp_main, kn_main = cp_pp.copy(), kn_pp.copy()
    isB = np.array([k == "BQ" for k in kind])
    yr = ymask[R]
    assert np.isfinite(yr[isB]).all(), "reference activity unavailable for a B pair"
    cp_main[isB] = ch_pred[H][isB] - yr[isB]
    kn_main[isB] = knn[H][isB] - yr[isB]
    out.update({"chemprop_abs": cp_main, "chemprop_predpred": cp_pp, "knn_sim": kn_main, "knn_predpred": kn_pp})
    for k, v in out.items():
        assert np.isfinite(v).all(), f"non-finite prediction for {k}: held-out label leak or failure"
    return out


def choose_tau(pred, ytrue):
    """F1-maximising tau on validation pairs; tie -> closest to 1.0."""
    best, best_f, curve = None, -1, []
    for t in TAU_GRID:
        tp, fp, tn, fn = conf(ytrue, (np.abs(pred) >= t).astype(int))
        f = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else 0.0
        curve.append((float(t), float(f)))
        if f > best_f + 1e-12 or (abs(f - best_f) <= 1e-12 and best is not None and abs(t - 1.0) < abs(best - 1.0)):
            best, best_f = float(t), f
    return best, best_f, curve


def main(argv=None):
    import argparse
    global OUT, CV, NBOOT
    ap = argparse.ArgumentParser(description="Direct activity-cliff detection evaluation (settings A and B)")
    PATHS.add_path_args(ap)
    ap.add_argument("--cv-dir", help="folder with fold0..fold4 models from cliffforge.cv (default <output-dir>/cv); "
                                     "if absent, reference validation thresholds are loaded from results/metrics/validation_thresholds.json")
    ap.add_argument("--n-boot", type=int, default=NBOOT, help="cluster-bootstrap resamples (default 2000; use 200 for a quick run)")
    args = ap.parse_args(argv)
    PATHS.apply_path_args(args)
    OUT = str(PATHS.output_dir() / "cliff_detection")
    CV = args.cv_dir or str(PATHS.output_dir() / "cv")
    NBOOT = args.n_boot
    os.makedirs(OUT, exist_ok=True)
    from .cv import fold_assignment
    from .relmodel import RelScorer
    import torch

    df, fold = fold_assignment()
    n = len(df)
    y = df.pIC50.values.astype(float)
    is_test = (df.split == "test").values
    is_train = ~is_test
    fps = [morgan(s) for s in df.smiles]
    S = np.array([DataStructs.BulkTanimotoSimilarity(f, fps) for f in fps], dtype=np.float32)

    # ---------- source-data verification of eligible pairs and counts ----------
    iu, ju = np.triu_indices(n, 1)
    sim = S[iu, ju]
    dy = np.abs(y[iu] - y[ju])
    ptype = np.where(is_test[iu] & is_test[ju], "A_test_test",
                     np.where(is_train[iu] & is_train[ju], "train_train", "B_test_train"))
    counts = {}
    for t in SIMS:
        for pt in ("A_test_test", "B_test_train", "train_train"):
            m = (sim >= t) & (ptype == pt)
            counts[f"sim>={t}|{pt}"] = {"pairs": int(m.sum()), "cliffs": int((m & (dy >= DELTA)).sum()),
                                        "tied_d0": int((m & (dy == 0)).sum())}
    near = ((sim >= 0.8) & (ptype != "train_train") & (np.abs(dy - 1.0) <= 0.4)).sum()

    # ---------- validation: out-of-fold predictions, train partition only ----------
    have_cv = all(os.path.exists(os.path.join(CV, f"fold{f}", "chemprop_preds.csv")) and
                  len(__import__("glob").glob(os.path.join(CV, f"fold{f}", "rel_siamese_s*.pt"))) == 3 for f in range(5))
    if have_cv:
        ch_oof = {}
        for f in range(5):
            p = pd.read_csv(os.path.join(CV, f"fold{f}", "chemprop_preds.csv"))
            tr_rows = np.where(fold >= 0)[0]
            assert len(p) == len(tr_rows)
            vec = np.full(n, np.nan)
            vec[tr_rows] = p["pIC50"].values
            ch_oof[f] = vec
        val_rows = []  # (h, r, kind, fold_of_model)
        mv = (sim >= SIM_PRIMARY) & (ptype == "train_train")
        for i, j in zip(iu[mv], ju[mv]):
            fi, fj = fold[i], fold[j]
            if fi == fj:
                val_rows.append((i, j, "AA", fi))
            else:
                h, r = (i, j) if fi < fj else (j, i)
                val_rows.append((h, r, "BQ", fold[h]))
        vpred = {k: np.zeros(len(val_rows)) for k in ["rel_siamese", "rel_sd", "chemprop_abs", "chemprop_predpred", "knn_sim", "knn_predpred"]}
        for f in range(5):
            sel = [k for k, v in enumerate(val_rows) if v[3] == f]
            if not sel:
                continue
            ymask = y.copy()
            ymask[fold == f] = np.nan      # fold-F molecules unavailable to fold-F models
            ymask[is_test] = np.nan
            knn_f = knn_pred(S, ymask, (fold >= 0) & (fold != f))
            sc = RelScorer(model_dir=os.path.join(CV, f"fold{f}"), device=torch.device("cpu"))
            assert len(sc.models) == 3
            o = predict_pairs(df, S, ymask, ch_oof[f], sc, [(val_rows[k][0], val_rows[k][1]) for k in sel],
                              [val_rows[k][2] for k in sel], knn_f)
            for key in vpred:
                vpred[key][sel] = o[key]
        vtrue = np.array([abs(y[h] - y[r]) >= DELTA for h, r, _, _ in val_rows]).astype(int)
        vd = np.array([y[h] - y[r] for h, r, _, _ in val_rows])
        taus, val_summary = {}, {}
        for m in MODELS:
            t, f, curve = choose_tau(vpred[m], vtrue)
            taus[m] = t
            slope = float(np.sum(vd * vpred[m]) / np.sum(vpred[m] ** 2))
            val_summary[m] = {"tau_val": t, "val_f1_at_tau": f, "f1_curve": curve,
                              "val_ap": float(average_precision_score(vtrue, np.abs(vpred[m]))),
                              "val_roc_auc": float(roc_auc_score(vtrue, np.abs(vpred[m]))),
                              "calibration_slope_true_on_pred_delta": slope,
                              "val_rmse_delta": float(np.sqrt(np.mean((vd - vpred[m]) ** 2)))}
        # POST-HOC SENSITIVITY (added after the primary test results were seen; primary tau/verdict unchanged):
        # tau chosen on the validation pair TYPE that matches the test setting (val-val for A, val-train for B).
        taus_kind = {}
        for kd, nm in (("AA", "A_test_test"), ("BQ", "B_test_train")):
            mk = np.array([v[2] == kd for v in val_rows])
            taus_kind[nm] = {m: choose_tau(vpred[m][mk], vtrue[mk])[0] for m in MODELS}
            taus_kind[nm]["n_val_pairs"] = int(mk.sum())
            taus_kind[nm]["n_val_cliffs"] = int(vtrue[mk].sum())
        pd.DataFrame({"idx_heldout_or_a": [v[0] for v in val_rows], "idx_ref_or_b": [v[1] for v in val_rows],
                      "setting": ["AA_val_val" if v[2] == "AA" else "BQ_val_train" for v in val_rows],
                      "cv_fold_of_model": [int(v[3]) for v in val_rows],
                      "tanimoto": S[[v[0] for v in val_rows], [v[1] for v in val_rows]],
                      "pred_delta_chemprop": vpred["chemprop_abs"], "pred_delta_rel": vpred["rel_siamese"],
                      "pred_delta_knn": vpred["knn_sim"],
                      "posthoc_true_delta": vd, "posthoc_true_cliff": vtrue}).to_csv(
            os.path.join(OUT, "activity_cliff_validation_predictions.csv"), index=False, float_format="%.4f")
        print("validation pairs:", len(val_rows), "positives:", int(vtrue.sum()), "taus:", taus, flush=True)
    else:
        ref = json.load(open(PATHS.reference_thresholds()))
        taus, taus_kind = ref["tau_val"], ref["tau_val_setting_specific"]
        val_summary, val_rows, vtrue = None, [], np.zeros(0, dtype=int)
        print("CV fold models not found -> using REFERENCE validation thresholds from", PATHS.reference_thresholds(), taus, flush=True)

    # ---------- test predictions (final models; test activities never used) ----------
    ymask = y.copy()
    ymask[is_test] = np.nan
    ch_pred = chemprop_io.predict_all(list(df.smiles))
    knn_t = knn_pred(S, ymask, is_train)
    scorer = RelScorer(device=torch.device("cpu"))
    assert len(scorer.models) == 3
    m6 = (sim >= 0.6) & (ptype != "train_train")
    pairs, kind, ptypes = [], [], []
    for i, j, pt in zip(iu[m6], ju[m6], ptype[m6]):
        if pt == "A_test_test":
            pairs.append((i, j)); kind.append("AA")
        else:                                    # B: query = the test member, reference = the train member
            h, r = (i, j) if is_test[i] else (j, i)
            pairs.append((h, r)); kind.append("BQ")
        ptypes.append(pt)
    P = predict_pairs(df, S, ymask, ch_pred, scorer, pairs, kind, knn_t)

    # ---- post-hoc ground truth (only now) ----
    H = np.array([p[0] for p in pairs]); Rr = np.array([p[1] for p in pairs])
    tsim = S[H, Rr]
    true_d = y[H] - y[Rr]
    ptypes = np.array(ptypes)
    def labels_for(setting, t):
        return (ptypes == setting) & (tsim >= t)

    results, boot_store = {}, {}
    for setting in ("A_test_test", "B_test_train"):
        for t in SIMS:
            mk = labels_for(setting, t)
            idx = np.where(mk)[0]
            ytrue = (np.abs(true_d[idx]) >= DELTA).astype(int)
            prev = float(ytrue.mean()) if len(ytrue) else None
            block = {"n_pairs": int(len(idx)), "n_cliffs": int(ytrue.sum()), "prevalence": prev,
                     "n_distinct_molecules_involved": int(len(set(H[idx]) | set(Rr[idx]))),
                     "n_distinct_query_molecules": int(len(set(H[idx]))) if setting == "B_test_train" else None,
                     "majority_class_baseline": {**prf(0, 0, int(len(idx) - ytrue.sum()), int(ytrue.sum())),
                                                 "note": "predicts 'no cliff' for every pair: accuracy = 1 - prevalence, recall 0, "
                                                         "precision undefined; chance AP = prevalence, ROC-AUC 0.5"},
                     "models": {}}
            proto = {"chemprop_abs": "chemprop_abs", "rel_siamese": "rel_siamese", "knn_sim": "knn_sim"}
            variants = dict(proto)
            if setting == "B_test_train":
                variants.update({"chemprop_predpred(secondary)": "chemprop_predpred", "knn_predpred(secondary)": "knn_predpred"})
            for name, key in variants.items():
                s = P[key][idx]
                base = name.split("(")[0]
                mm = {}
                tau_model = {"chemprop_abs": "chemprop_abs", "chemprop_predpred": "chemprop_abs",
                             "rel_siamese": "rel_siamese", "knn_sim": "knn_sim", "knn_predpred": "knn_sim"}[base]
                for rule in ("tau_val", "tau_fixed_1.0"):
                    tau = taus[tau_model] if rule == "tau_val" else TAU_FIXED
                    yh = (np.abs(s) >= tau).astype(int)
                    r = prf(*conf(ytrue, yh))
                    r["tau"] = tau
                    r["predicted_cliffs"] = int(yh.sum())
                    r["recall_exact_ci95"] = cp_ci(r["tp"], int(ytrue.sum()))
                    mm[rule] = r
                if ytrue.sum() > 0 and ytrue.sum() < len(ytrue):
                    mm["average_precision"] = float(average_precision_score(ytrue, np.abs(s)))
                    mm["roc_auc"] = float(roc_auc_score(ytrue, np.abs(s)))
                else:
                    mm["average_precision"] = mm["roc_auc"] = None
                sd = np.sign(true_d[idx]) != 0
                mm["direction_accuracy_secondary"] = {"n": int(sd.sum()), "value": float((np.sign(s[sd]) == np.sign(true_d[idx][sd])).mean())} if sd.sum() else None
                block["models"][name] = mm
            block["posthoc_setting_specific_tau"] = {
                "note": "POST-HOC sensitivity (protocol deviation disclosed): tau chosen on the matching validation pair type only; "
                        "does not affect the declared verdict",
                "thresholds": {m: taus_kind[setting][m] for m in MODELS},
                "models": {m: {**prf(*conf(ytrue, (np.abs(P[m][idx]) >= taus_kind[setting][m]).astype(int))),
                               "tau": taus_kind[setting][m]} for m in MODELS}}
            results[f"{setting}|sim>={t}"] = block

            # ---- cluster bootstrap (paired) ----
            if True:
                if setting == "A_test_test":
                    clus = components([(int(H[i]), int(Rr[i])) for i in idx])
                    cname = "connected_components_of_pair_graph"
                else:
                    clus = [int(H[i]) for i in idx]
                    cname = "query_test_molecule"
                series_clus = [df.series[int(H[i])] for i in idx] if setting == "B_test_train" else None
                scores = {m: np.abs(P[m][idx]) for m in MODELS}
                yh = {(m, rule): (np.abs(P[m][idx]) >= (taus[m] if rule == "tau_val" else TAU_FIXED)).astype(int)
                      for m in MODELS for rule in ("tau_val", "tau_fixed_1.0")}

                def ev(ix, scores=scores, yh=yh, ytrue=ytrue):
                    out = {}
                    yt = ytrue[ix]
                    has_both = 0 < yt.sum() < len(yt)
                    for m in MODELS:
                        out[f"{m}|ap"] = float(average_precision_score(yt, scores[m][ix])) if has_both else np.nan
                        out[f"{m}|roc"] = float(roc_auc_score(yt, scores[m][ix])) if has_both else np.nan
                        for rule in ("tau_val", "tau_fixed_1.0"):
                            h_ = yh[(m, rule)][ix]
                            tp = int(((yt == 1) & (h_ == 1)).sum()); fp = int(((yt == 0) & (h_ == 1)).sum())
                            fn = int(((yt == 1) & (h_ == 0)).sum())
                            out[f"{m}|precision|{rule}"] = tp / (tp + fp) if tp + fp else np.nan
                            out[f"{m}|recall|{rule}"] = tp / (tp + fn) if tp + fn else np.nan
                            out[f"{m}|f1|{rule}"] = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else np.nan
                    return out
                bs = cluster_boot(clus, ev)
                keys = bs[0].keys()
                vals = {k: [b[k] for b in bs] for k in keys}
                unc = {"cluster_unit": cname, "n_clusters": len(set(clus)), "n_resamples": NBOOT,
                       "resamples_without_positives": int(sum(1 for b in bs if not np.isfinite(b[f'{MODELS[0]}|ap']))),
                       "metrics": {k: ci(v) for k, v in vals.items()}, "paired_differences_rel_minus_chemprop": {}}
                for k in ("ap", "roc"):
                    d = [(a - b) if (np.isfinite(a) and np.isfinite(b)) else None
                         for a, b in zip(vals[f"rel_siamese|{k}"], vals[f"chemprop_abs|{k}"])]
                    unc["paired_differences_rel_minus_chemprop"][k] = ci(d)
                for rule in ("tau_val", "tau_fixed_1.0"):
                    for k in ("f1", "precision", "recall"):
                        d = [(a - b) if (np.isfinite(a) and np.isfinite(b)) else None
                             for a, b in zip(vals[f"rel_siamese|{k}|{rule}"], vals[f"chemprop_abs|{k}|{rule}"])]
                        unc["paired_differences_rel_minus_chemprop"][f"{k}|{rule}"] = ci(d)
                if series_clus is not None:
                    bs2 = cluster_boot(series_clus, ev, seed=1)
                    d2 = [(b["rel_siamese|ap"] - b["chemprop_abs|ap"]) if (np.isfinite(b["rel_siamese|ap"]) and np.isfinite(b["chemprop_abs|ap"])) else None for b in bs2]
                    unc["series_cluster_variant_delta_ap"] = {"n_clusters": len(set(series_clus)), **ci(d2)}
                block["uncertainty"] = unc

            # matched-pair discordance on TRUE cliffs (who finds the cliffs the other misses)
            for rule in ("tau_val", "tau_fixed_1.0"):
                yh_r = (np.abs(P["rel_siamese"][idx]) >= (taus["rel_siamese"] if rule == "tau_val" else TAU_FIXED)).astype(int)
                yh_c = (np.abs(P["chemprop_abs"][idx]) >= (taus["chemprop_abs"] if rule == "tau_val" else TAU_FIXED)).astype(int)
                pos = ytrue == 1
                block.setdefault("matched_pair_discordance_on_true_cliffs", {})[rule] = {
                    "rel_detects_chemprop_misses": int((pos & (yh_r == 1) & (yh_c == 0)).sum()),
                    "chemprop_detects_rel_misses": int((pos & (yh_c == 1) & (yh_r == 0)).sum()),
                    "both_detect": int((pos & (yh_r == 1) & (yh_c == 1)).sum()),
                    "neither_detects": int((pos & (yh_r == 0) & (yh_c == 0)).sum())}

    # ---------- decision rule ----------
    def g(setting, model, rule, key):
        return results[f"{setting}|sim>=0.8"]["models"][model][rule][key]
    A, Bk = "A_test_test", "B_test_train"
    apA = {m: results[f"{A}|sim>=0.8"]["models"][m]["average_precision"] for m in MODELS}
    apB = {m: results[f"{Bk}|sim>=0.8"]["models"][m]["average_precision"] for m in MODELS}
    unc_B = results[f"{Bk}|sim>=0.8"]["uncertainty"]["paired_differences_rel_minus_chemprop"]
    f1 = lambda s, m, r: (results[f"{s}|sim>=0.8"]["models"][m][r]["f1"] or 0.0)
    dAP = apB["rel_siamese"] - apB["chemprop_abs"]
    dF1 = {r: f1(Bk, "rel_siamese", r) - f1(Bk, "chemprop_abs", r) for r in ("tau_val", "tau_fixed_1.0")}
    c1 = dAP > 0 and unc_B["ap"]["ci95"][0] is not None and unc_B["ap"]["ci95"][0] > 0
    c2 = all(v > 0 for v in dF1.values())
    c3 = (apA["rel_siamese"] >= apA["chemprop_abs"]) and all(f1(A, "rel_siamese", r) >= f1(A, "chemprop_abs", r) for r in ("tau_val", "tau_fixed_1.0"))
    if c1 and c2 and c3:
        verdict = "A"
    elif dAP <= 0 or all(v <= 0 for v in dF1.values()):
        verdict = "B"
    else:
        verdict = "C"
    names = {"A": "Preliminary evidence of better activity-cliff detection under the evaluated conditions.",
             "B": "No demonstrated improvement over the baseline.",
             "C": "Inconclusive: insufficient positive cliffs or evaluation data to make a reliable comparison."}
    decision = {"verdict": verdict, "statement": names[verdict],
                "criteria": {"B_primary_delta_AP_rel_minus_chemprop": dAP, "B_primary_delta_AP_ci95_cluster": unc_B["ap"]["ci95"],
                             "criterion_1_AP_gain_with_CI_excluding_0": bool(c1),
                             "B_primary_delta_F1_by_rule": dF1, "criterion_2_F1_gain_at_both_thresholds": bool(c2),
                             "A_primary_AP": apA, "criterion_3_A_not_worse": bool(c3),
                             "A_primary_n_positive_cliffs": results[f"{A}|sim>=0.8"]["n_cliffs"],
                             "A_is_descriptive_only": results[f"{A}|sim>=0.8"]["n_cliffs"] < 10}}

    # ---------- endpoint-compatibility checks ----------
    te = np.where(is_test)[0]
    resid = ch_pred[te] - y[te]
    endpoint = {"chemprop_abs": {"trained_on": "data/moonshot train partition pIC50 (models/chemprop_moonshot/args.json)",
                                 "test_rmse_pIC50": float(np.sqrt(np.mean(resid ** 2))), "test_mean_bias": float(resid.mean())},
                "rel_siamese": {"trained_on": "pIC50 differences of same-series train pairs, data/moonshot/all.csv split==train",
                                "output_unit": "predicted pIC50(query) - pIC50(reference)"},
                "fep_models_used": False,
                "validation_calibration": ({m: {"slope_true_on_pred": val_summary[m]["calibration_slope_true_on_pred_delta"],
                                                "rmse_delta": val_summary[m]["val_rmse_delta"]} for m in MODELS}
                                           if val_summary else "not recomputed (reference thresholds used)")}

    # ---------- per-pair predictions CSV ----------
    rows = pd.DataFrame({
        "setting": np.where(ptypes == "A_test_test", "A_strict_test_test", "B_leadopt_test_train"),
        "idx_first_or_query": H, "idx_second_or_reference": Rr,
        "smiles_first_or_query": df.smiles.values[H], "smiles_second_or_reference": df.smiles.values[Rr],
        "split_first_or_query": df.split.values[H], "split_second_or_reference": df.split.values[Rr],
        "series_first_or_query": df.series.values[H],
        "tanimoto": tsim,
        "pred_delta_chemprop": P["chemprop_abs"], "pred_delta_chemprop_predpred": P["chemprop_predpred"],
        "pred_delta_rel_siamese": P["rel_siamese"], "pred_delta_rel_sd": P["rel_sd"], "pred_delta_knn": P["knn_sim"]})
    for m, col in (("chemprop_abs", "chemprop"), ("rel_siamese", "rel_siamese"), ("knn_sim", "knn")):
        rows[f"tau_val_{col}"] = taus[m]
        rows[f"predcliff_tauval_{col}"] = (np.abs(P[m]) >= taus[m]).astype(int)
        rows[f"predcliff_tau1.0_{col}"] = (np.abs(P[m]) >= TAU_FIXED).astype(int)
    rows["posthoc_measured_pIC50_first_or_query"] = y[H]
    rows["posthoc_measured_pIC50_second_or_reference"] = y[Rr]
    rows["posthoc_true_delta"] = true_d
    rows["posthoc_true_cliff_abs_delta_ge_1"] = (np.abs(true_d) >= DELTA).astype(int)
    rows["posthoc_true_cliff_PRIMARY"] = ((tsim >= SIM_PRIMARY) & (np.abs(true_d) >= DELTA)).astype(int)
    rows.to_csv(os.path.join(OUT, "activity_cliff_predictions.csv"), index=False, float_format="%.4f")

    metrics = {
        "declared_protocol": __doc__.split("DECLARED BEFORE ANY DETECTION METRIC WAS COMPUTED")[1].split("Leakage guard")[0].strip(),
        "source_data_counts_recomputed": counts,
        "pairs_within_0.4_of_cliff_boundary_(|d| in [0.6,1.4], test-involved sim>=0.8)": int(near),
        "duplicates": {"n_molecules": int(n), "n_unique_smiles": int(df.smiles.nunique()),
                       "compounds_with_replicate_measurements": int((df.n_rep > 1).sum()),
                       "assay_noise_sd_pIC50": 0.396},
        "models": {"chemprop_abs": "models/chemprop_moonshot/fold_0/model_0/model.pt (ChemProp 1.6.1 absolute-pIC50 baseline)",
                   "rel_siamese": "models/rel_siamese_s{0,1,2}.pt (train-partition only)",
                   "knn_sim": "similarity-weighted 5-NN over train labels (non-learning)"},
        "validation": {"scheme": "5-fold CV inside the train partition (random_state=0); fold models retrained with the "
                                 "same recipes; pairs = train-train Tanimoto>=0.80",
                       "n_validation_pairs": len(val_rows), "n_validation_cliffs": int(vtrue.sum()),
                       "n_val_val_pairs": int(sum(1 for v in val_rows if v[2] == "AA")),
                       "n_val_train_pairs": int(sum(1 for v in val_rows if v[2] == "BQ")),
                       "thresholds": taus, "posthoc_setting_specific_thresholds": taus_kind, "per_model": val_summary},
        "tau_fixed": TAU_FIXED,
        "endpoint_compatibility": endpoint,
        "results": results, "decision_rule_outcome": decision,
        "leakage_controls": ["held-out activities masked (NaN) in all inference code paths; assertions on finite predictions",
                             "tau chosen only on out-of-fold train-partition validation pairs",
                             "test partition never used in CV folds (asserted: no test rows with fold>=0)",
                             "same split as prior work"],
    }
    if not have_cv:
        metrics["validation"] = {"scheme": "REFERENCE thresholds loaded from results/metrics/validation_thresholds.json "
                                           "(5-fold CV not re-run in this invocation)",
                                 "thresholds": taus, "posthoc_setting_specific_thresholds": taus_kind}
    json.dump(metrics, open(os.path.join(OUT, "activity_cliff_detection_metrics.json"), "w"), indent=2, default=lambda o: bool(o) if isinstance(o, np.bool_) else (float(o) if isinstance(o, np.floating) else int(o) if isinstance(o, np.integer) else str(o)))
    print("VERDICT:", verdict, names[verdict])
    print(json.dumps(decision["criteria"], indent=1, default=str))
    for k, b in results.items():
        if "sim>=0.8" in k:
            print(k, "n=", b["n_pairs"], "cliffs=", b["n_cliffs"])
            for m, mm in b["models"].items():
                print("   ", m, "tauval", mm["tau_val"]["tau"], {x: mm["tau_val"][x] for x in ("tp", "fp", "tn", "fn")}, "| tau1.0", {x: mm["tau_fixed_1.0"][x] for x in ("tp", "fp", "tn", "fn")}, "| AP", None if mm["average_precision"] is None else round(mm["average_precision"], 3))


if __name__ == "__main__":
    main()
