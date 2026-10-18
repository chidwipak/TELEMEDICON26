"""
cliffforge.metrics -- shared evaluation functions for pair-level model comparison.
===================================================================================
Used by both CliffForge_MMP_Colab.ipynb and CliffForge_ChemProp_Colab.ipynb to
ensure identical metric computation across all models.

All functions take numpy arrays and return plain dicts so results can be directly
serialised to JSON for the combined summary.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Optional

import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# Direction accuracy
# ─────────────────────────────────────────────────────────────────────────────

def direction_accuracy(
    true_delta: np.ndarray,
    pred_delta: np.ndarray,
    n_bootstrap: int = 2000,
    seed: int = 0,
) -> dict:
    """
    Fraction of pairs where sign(pred) == sign(true).

    Parameters
    ----------
    true_delta : 1-D array  -- measured pIC50(a) - pIC50(b)
    pred_delta : 1-D array  -- predicted delta for the same pairs
    n_bootstrap : int       -- pair-bootstrap resamples for CI
    seed : int

    Returns
    -------
    dict with keys: n, direction_accuracy, n_correct, ci95, binom_p_vs_0.5
    """
    from scipy.stats import binomtest

    mask = np.isfinite(pred_delta) & (true_delta != 0) & (pred_delta != 0)
    if mask.sum() == 0:
        return {"n": 0, "direction_accuracy": None, "n_correct": 0,
                "ci95": [None, None], "binom_p_vs_0.5": None}

    correct = np.sign(true_delta[mask]) == np.sign(pred_delta[mask])
    n = int(mask.sum())
    n_ok = int(correct.sum())

    rng = np.random.RandomState(seed)
    bs = [correct[rng.randint(0, n, n)].mean() for _ in range(n_bootstrap)]
    p = float(binomtest(n_ok, n, 0.5).pvalue)

    return {
        "n": n,
        "direction_accuracy": float(correct.mean()),
        "n_correct": n_ok,
        "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
        "binom_p_vs_0.5": p,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Binary cliff detection
# ─────────────────────────────────────────────────────────────────────────────

def cliff_detection_metrics(
    is_cliff: np.ndarray,
    pred_abs_delta: np.ndarray,
    tau: float = 1.0,
) -> dict:
    """
    Binary cliff detection at a fixed threshold tau.

    Parameters
    ----------
    is_cliff      : boolean 1-D array   -- ground-truth cliff label
    pred_abs_delta: float 1-D array     -- |predicted delta| (higher = more likely cliff)
    tau           : float               -- threshold: predict cliff if pred_abs_delta >= tau

    Returns
    -------
    dict with keys: n, n_cliffs, prevalence, tau, tp, fp, tn, fn,
                    precision, recall, f1, ap, roc_auc
    """
    from sklearn.metrics import average_precision_score, roc_auc_score

    y = is_cliff.astype(int)
    scores = pred_abs_delta
    y_pred = (scores >= tau).astype(int)

    tp = int(((y == 1) & (y_pred == 1)).sum())
    fp = int(((y == 0) & (y_pred == 1)).sum())
    tn = int(((y == 0) & (y_pred == 0)).sum())
    fn = int(((y == 1) & (y_pred == 0)).sum())

    P = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    R = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    F = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else 0.0

    out = {
        "n": len(y), "n_cliffs": int(y.sum()),
        "prevalence": float(y.mean()),
        "tau": tau,
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "precision": P, "recall": R, "f1": F,
    }
    if 0 < y.sum() < len(y):
        out["ap"]      = float(average_precision_score(y, scores))
        out["roc_auc"] = float(roc_auc_score(y, scores))
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Cluster-bootstrap AP difference
# ─────────────────────────────────────────────────────────────────────────────

def bootstrap_ap_diff(
    is_cliff: np.ndarray,
    scores_a: np.ndarray,
    scores_b: np.ndarray,
    cluster_ids: np.ndarray,
    n_bootstrap: int = 2000,
    seed: int = 0,
) -> dict:
    """
    Cluster-bootstrap 95% CI for AP(model_a) - AP(model_b).

    Parameters
    ----------
    is_cliff   : boolean array -- cliff labels
    scores_a   : float array  -- |pred_delta| from model A (higher = more likely cliff)
    scores_b   : float array  -- |pred_delta| from model B
    cluster_ids: array        -- cluster label per pair (resampling unit)
    n_bootstrap: int
    seed       : int

    Returns
    -------
    dict: ci95_lo, ci95_hi, n_valid_resamples, point_estimate
    """
    from sklearn.metrics import average_precision_score

    y = is_cliff.astype(int)
    cl = defaultdict(list)
    for i, c in enumerate(cluster_ids):
        cl[c].append(i)
    keys = list(cl)

    rng = np.random.RandomState(seed)
    diffs = []
    for _ in range(n_bootstrap):
        pick = rng.randint(0, len(keys), len(keys))
        idx = np.array([j for p in pick for j in cl[keys[p]]], dtype=int)
        yt = y[idx]
        if yt.sum() == 0 or yt.sum() == len(yt):
            continue
        diffs.append(
            average_precision_score(yt, scores_a[idx]) -
            average_precision_score(yt, scores_b[idx])
        )

    if len(diffs) < 20:
        return {"ci95_lo": None, "ci95_hi": None,
                "n_valid_resamples": len(diffs), "point_estimate": None}

    point = float(
        average_precision_score(y[y.sum() > 0], scores_a[y.sum() > 0]) -  # type: ignore
        average_precision_score(y[y.sum() > 0], scores_b[y.sum() > 0])    # type: ignore
    ) if 0 < y.sum() < len(y) else None

    return {
        "ci95_lo": float(np.percentile(diffs, 2.5)),
        "ci95_hi": float(np.percentile(diffs, 97.5)),
        "n_valid_resamples": len(diffs),
        "point_estimate": float(np.mean(diffs)),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Within-series Spearman ρ
# ─────────────────────────────────────────────────────────────────────────────

def within_series_spearman(
    true_y: np.ndarray,
    pred_y: np.ndarray,
    series_ids: np.ndarray,
    min_series_size: int = 5,
) -> dict:
    """
    Mean Spearman ρ across chemical series.

    Parameters
    ----------
    true_y          : measured pIC50 per molecule
    pred_y          : predicted pIC50 per molecule (absolute scale)
    series_ids      : series label per molecule
    min_series_size : only series with >= this many members are included

    Returns
    -------
    dict: n_series, mean_rho, per_series (list of floats)
    """
    from scipy.stats import spearmanr

    per_series = []
    unique = np.unique(series_ids)
    for s in unique:
        mask = series_ids == s
        if mask.sum() < min_series_size:
            continue
        yt = true_y[mask]
        yp = pred_y[mask]
        if not np.isfinite(yp).all():
            continue
        rho, _ = spearmanr(yt, yp)
        per_series.append(float(rho))

    return {
        "n_series": len(per_series),
        "mean_rho": float(np.mean(per_series)) if per_series else None,
        "per_series": per_series,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Combine results from multiple models into a unified dict
# ─────────────────────────────────────────────────────────────────────────────

def combine_summaries(base_summary: dict, new_model_name: str, new_model_results: dict) -> dict:
    """
    Merge a new model's results into an existing summary dict.

    new_model_results keys expected:
        cliff_detection   -- dict keyed by 'test_test', 'test_train'
        direction_accuracy -- dict keyed by subset labels
        within_series_spearman -- dict

    Returns a new summary dict with the new model added.
    """
    import copy
    out = copy.deepcopy(base_summary)

    # Cliff detection
    for setting, metrics in new_model_results.get("cliff_detection", {}).items():
        key = f"{setting}|{new_model_name}"
        out.setdefault("cliff_detection", {})[key] = metrics

    # Direction accuracy
    for subset, metrics in new_model_results.get("direction_accuracy", {}).items():
        out.setdefault("direction_accuracy", {}).setdefault(subset, {})[new_model_name] = metrics

    # Within-series Spearman
    out.setdefault("within_series_spearman", {}).setdefault("models", {})[new_model_name] = \
        new_model_results.get("within_series_spearman", {})

    return out
