"""
cliffforge_recall.py  --  CliffForge PoC, Experiment D (generation <-> held-out outcomes)
=======================================================================================
Retrospective, leakage-controlled test of the GENERATE -> SCORE loop.

* Transformations are mined from the TRAIN partition only (TransformDB).
* Every held-out TEST compound is used as a seed. The generator proposes analogues; a
  proposal that is *exactly* another measured TEST compound is a "recovered" analogue, so
  the true delta pIC50(candidate - seed) is known and was never used for mining/training.
* Recall denominator: all test-test single-edit matched pairs (shared scaffold, sidechains
  <= 12 heavy atoms, scaffold >= 60 % of the molecule), i.e. the pairs a perfect library
  could reach.
* Direction accuracy of: relative Siamese ensemble, ChemProp absolute, and the raw
  transformation PRECEDENT prior (sign of the mean observed delta in train; model-free).
Run:  python -m cliffforge.recall
"""
import json
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from rdkit import Chem

from . import generate as G
from . import paths as P
from . import chemprop_io

OUT = None  # set in main(): <output-dir>/recall


def acc(true, pred):
    m = np.isfinite(pred) & (pred != 0) & (true != 0)
    if m.sum() == 0:
        return None
    ok = (np.sign(true[m]) == np.sign(pred[m]))
    rng = np.random.RandomState(0)
    bs = [ok[rng.randint(0, len(ok), len(ok))].mean() for _ in range(2000)]
    return {"n": int(m.sum()), "direction_accuracy": float(ok.mean()), "n_correct": int(ok.sum()),
            "ci95_pair_bootstrap": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
            "binom_p_vs_0.5": float(binomtest(int(ok.sum()), len(ok), 0.5).pvalue)}


def main(argv=None):
    import argparse
    global OUT
    ap = argparse.ArgumentParser(description="Generate->score retrospective check on held-out measured analogues")
    P.add_path_args(ap)
    P.apply_path_args(ap.parse_args(argv))
    OUT = str(P.output_dir() / "recall")
    os.makedirs(OUT, exist_ok=True)
    tr, te, al = G.load_split()
    tdb = G.TransformDB(tr)
    known = {s: sp for s, sp in zip(al.smiles, al.split)}
    y = dict(zip(al.smiles, al.pIC50))
    row_of = {s: i for i, s in enumerate(al.smiles)}
    test_smiles = list(te.smiles)

    # ---- denominator: test-test single-edit matched pairs ----
    idx = defaultdict(lambda: defaultdict(set))
    for s in test_smiles:
        m = Chem.MolFromSmiles(s)
        for scaf, side, n, env in G.cuts(m, G.PARAMS["max_side_heavy_atoms"], G.PARAMS["min_scaffold_fraction"]):
            idx[(scaf, env)][side].add(s)
    mmp_pairs = set()
    for sides in idx.values():
        items = list(sides.items())
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                for a in items[i][1]:
                    for b in items[j][1]:
                        if a != b:
                            mmp_pairs.add(tuple(sorted((a, b))))

    # ---- generate from every test seed; keep hits that are measured TEST compounds ----
    rows = []
    for k, s in enumerate(test_smiles):
        df = G.generate(s, tdb, known, seed_id=f"t{k}", min_precedent=1)
        if not len(df):
            continue
        hits = df[(df.known_in_moonshot == "test") & (df.valid_rdkit == True) & (df.unique_candidate == True)]  # noqa: E712
        for r in hits.itertuples():
            rows.append({"seed_smiles": s, "cand_smiles": r.candidate_smiles, "transformation": r.transformation,
                         "n_precedent_pairs": int(r.n_precedent_pairs),
                         "precedent_mean_delta": float(r.precedent_mean_delta_pIC50),
                         "other_filters_pass": (r.filter_status == "pass") or (r.filter_reason in ("",))})
    hits = pd.DataFrame(rows)
    if hits.empty:
        json.dump({"n_recovered_pairs": 0}, open(os.path.join(OUT, "generation_recall.json"), "w"))
        print("no recovered pairs")
        return
    hits["pair_key"] = [tuple(sorted(p)) for p in zip(hits.seed_smiles, hits.cand_smiles)]
    # one record per unordered pair: keep the transformation with the most precedent
    hits = hits.sort_values("n_precedent_pairs", ascending=False).drop_duplicates("pair_key")
    hits["true_delta"] = [y[c] - y[s] for s, c in zip(hits.seed_smiles, hits.cand_smiles)]
    hits["same_series"] = [al.series[row_of[s]] == al.series[row_of[c]] for s, c in zip(hits.seed_smiles, hits.cand_smiles)]

    # ---- predictions ----
    from .relmodel import RelScorer
    import torch
    sc = RelScorer(device=torch.device("cpu"))
    d_all = sc.pair_delta_all(list(hits.cand_smiles), list(hits.seed_smiles))
    hits["pred_delta_rel_mean"] = d_all.mean(0)
    cp = pd.Series(chemprop_io.predict_all(list(al.smiles)), index=list(al.smiles))
    hits["pred_delta_chemprop"] = [cp[c] - cp[s] for s, c in zip(hits.seed_smiles, hits.cand_smiles)]
    hits["pred_delta_precedent_prior"] = hits.precedent_mean_delta

    t = hits.true_delta.values
    subsets = {
        "all_recovered": np.ones(len(hits), bool),
        "precedent>=2 (default generator filter)": (hits.n_precedent_pairs >= 2).values,
        "|true delta|>=0.4 (above assay noise)": (np.abs(t) >= 0.4),
        "|true delta|>=1.0 (cliff-size)": (np.abs(t) >= 1.0),
        "same_series": hits.same_series.values,
    }
    res = {}
    for nm, mask in subsets.items():
        res[nm] = {"n_pairs": int(mask.sum()),
                   "rel_siamese": acc(t[mask], hits.pred_delta_rel_mean.values[mask]),
                   "chemprop_abs": acc(t[mask], hits.pred_delta_chemprop.values[mask]),
                   "precedent_prior": acc(t[mask], hits.pred_delta_precedent_prior.values[mask])}
        ra_, rb_ = hits.pred_delta_rel_mean.values[mask], hits.pred_delta_chemprop.values[mask]
        tt_ = t[mask]
        mm = (tt_ != 0) & (ra_ != 0) & (rb_ != 0)
        oa = int(((np.sign(ra_[mm]) == np.sign(tt_[mm])) & (np.sign(rb_[mm]) != np.sign(tt_[mm]))).sum())
        ob = int(((np.sign(ra_[mm]) != np.sign(tt_[mm])) & (np.sign(rb_[mm]) == np.sign(tt_[mm]))).sum())
        res[nm]["mcnemar_rel_vs_chemprop"] = {"rel_only_correct": oa, "chemprop_only_correct": ob,
                                              "p_two_sided": float(binomtest(oa, oa + ob, 0.5).pvalue) if oa + ob else None}
    summary = {
        "n_test_seeds": len(test_smiles),
        "n_test_test_single_edit_pairs_(reachable_universe)": len(mmp_pairs),
        "n_recovered_pairs_(min_precedent=1)": int(len(hits)),
        "recall_of_reachable_pairs": float(len(hits) / max(len(mmp_pairs), 1)),
        "n_recovered_with_precedent>=2": int((hits.n_precedent_pairs >= 2).sum()),
        "transformation_db": tdb.summary(),
        "direction_accuracy": res,
        "caveats": ["Recovered pairs are a biased sample: only transformations seen in train can recover a pair.",
                    "Pairs share molecules (not independent); bootstrap CIs are optimistic.",
                    "Test labels were used only to score recovered pairs, never for mining/training."],
    }
    hits.drop(columns=["pair_key"]).to_csv(os.path.join(OUT, "generation_recall_pairs.csv"), index=False, float_format="%.4f")
    json.dump(summary, open(os.path.join(OUT, "generation_recall.json"), "w"), indent=2)
    print(json.dumps({k: v for k, v in summary.items() if k != "direction_accuracy"}, indent=1))
    for nm, v in res.items():
        print(nm, "n=", v["n_pairs"], {m: (round(v[m]["direction_accuracy"], 3) if v[m] else None)
                                       for m in ("rel_siamese", "chemprop_abs", "precedent_prior")})


if __name__ == "__main__":
    main()
