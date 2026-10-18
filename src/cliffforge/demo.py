"""
CliffForge proof-of-concept demo: seed ligand -> constrained analogue generation ->
filters -> relative-affinity scoring + activity-cliff alerts -> ranked shortlist.

PROTOTYPE. Generation is example-based (matched-molecular-pair transformation transfer
mined from the Moonshot TRAIN partition), NOT a neural generator. Scores are 2D
model predictions of delta pIC50 vs the seed, NOT measured potencies and NOT docking.
Pocket-aware docking was NOT run (no docking software / Mpro structure in this repo).

Examples
  # headline demo: held-out Moonshot test compound #44 as the seed
  python -m cliffforge.demo --seed-index 44
  # any SMILES (generation + filters only needs RDKit; scoring needs the checkpoints)
  python -m cliffforge.demo --seed-smiles "<SMILES>" \
      --seed-id mylig --out-dir /tmp/out
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs

from . import generate as G
from . import paths as P

HIGH_SD = 0.3          # ensemble disagreement flagged as "high uncertainty"
CLIFF_SIM, CLIFF_DELTA = 0.8, 1.0


def postera_id(flat):
    """Map a 2D SMILES back to the raw Moonshot record ID (provenance)."""
    raw = str(P.raw_csv())
    ids = []
    if not os.path.exists(raw):      # provenance lookup is optional (needs the raw Moonshot file)
        return ids
    for smi, cid in pd.read_csv(raw, usecols=["suspected_SMILES", "Canonical_PostEra_ID"]).itertuples(index=False):
        m = Chem.MolFromSmiles(str(smi)) if isinstance(smi, str) else None
        if m is not None and Chem.MolToSmiles(m, isomericSmiles=False) == flat:
            ids.append(cid)
    return ids


def score(df, seed_smiles, train_df, model_dir):
    """Add relative-affinity scores (delta pIC50 vs seed) for filter-passing candidates."""
    cols = ["predicted_relative_activity", "predicted_rel_sd", "max_train_tanimoto", "cliff_alert",
            "score_type", "score_source", "uncertainty_or_caveat"]
    for c in cols:
        df[c] = np.nan if c in ("predicted_relative_activity", "predicted_rel_sd", "max_train_tanimoto") else ""
    ok = df.filter_status == "pass"
    try:
        from .relmodel import RelScorer
        import torch
        scorer = RelScorer(model_dir=model_dir, device=torch.device("cpu"))
    except FileNotFoundError:
        df["score_type"] = "not_run"
        df["score_source"] = "no relative-model checkpoints found; affinity scoring not performed"
        return df, False, None
    if ok.sum() == 0:
        return df, True, scorer
    cand = df.loc[ok, "candidate_smiles"].tolist()
    mean, sd = scorer.pair_delta(cand, [seed_smiles] * len(cand))
    tr_fps = [G.fp(Chem.MolFromSmiles(s)) for s in train_df.smiles]
    maxsim = [max(DataStructs.BulkTanimotoSimilarity(G.fp(Chem.MolFromSmiles(s)), tr_fps)) for s in cand]
    df.loc[ok, "predicted_relative_activity"] = mean
    df.loc[ok, "predicted_rel_sd"] = sd
    df.loc[ok, "max_train_tanimoto"] = maxsim
    df.loc[ok, "score_type"] = "model-predicted delta pIC50 (candidate - seed); NOT measured"
    df.loc[ok, "score_source"] = f"Siamese D-MPNN ensemble x{len(scorer.models)} (cliffforge.relmodel)"
    alert, cav = [], []
    for i in df.index[ok]:
        r = df.loc[i]
        a = "predicted_cliff" if (r.parent_tanimoto >= CLIFF_SIM and abs(r.predicted_relative_activity) >= CLIFF_DELTA) else ""
        c = [f"ensemble SD {r.predicted_rel_sd:.2f}" + (" (HIGH)" if r.predicted_rel_sd >= HIGH_SD else ""),
             f"{int(r.n_precedent_pairs)} train precedent pair(s)",
             "2D model: no pocket / stereo / docking"]
        if r.known_in_moonshot != "no":
            c.append(f"already measured in Moonshot ({r.known_in_moonshot} partition)")
        alert.append(a)
        cav.append("; ".join(c))
    df.loc[ok, "cliff_alert"] = alert
    df.loc[ok, "uncertainty_or_caveat"] = cav
    return df, True, scorer


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--seed-index", type=int, help="row index in data/moonshot/all.csv (NOT the test-set-local index)")
    g.add_argument("--seed-smiles", help="any seed SMILES")
    g.add_argument("--auto-seed", action="store_true",
                   help="pre-declared rule: held-out test compound, pIC50>=7, most filter-passing candidates")
    ap.add_argument("--seed-id", default=None)
    ap.add_argument("--out-dir", default=None, help="default: <output-dir>/demo")
    P.add_path_args(ap)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--no-score", action="store_true")
    a = ap.parse_args(argv)
    P.apply_path_args(a)
    a.out_dir = a.out_dir or str(P.output_dir() / "demo")
    a.model_dir = a.model_dir or str(P.model_dir())
    os.makedirs(a.out_dir, exist_ok=True)

    train, test, alldf = G.load_split()
    seed_label = None
    seed_part = "external"
    tdb = G.TransformDB(train)
    print("transformation DB (train partition only):", tdb.summary())
    if a.auto_seed:
        a.seed_index, table = G.pick_seed(tdb, alldf)
        table.to_csv(os.path.join(a.out_dir, "seed_selection_table.csv"), index=False)
        print(f"auto-seed rule picked all.csv row {a.seed_index}; candidates table -> seed_selection_table.csv")
    if a.seed_index is not None:
        r = alldf.iloc[a.seed_index]
        seed_smiles, seed_part, seed_label = r.smiles, r.split, float(r.pIC50)
        seed_id = a.seed_id or f"moonshot_idx{a.seed_index}"
        ids = sorted(set(postera_id(seed_smiles)))
    else:
        seed_smiles = a.seed_smiles
        seed_id = a.seed_id or "seed"
        flat = G.flat_smiles(seed_smiles)
        hit = alldf[alldf.smiles == flat]
        if len(hit):
            seed_part, seed_label = hit.iloc[0].split, float(hit.iloc[0].pIC50)
        ids = []
    print(f"seed {seed_id}: {seed_smiles}  (partition={seed_part})")
    known = {s: sp for s, sp in zip(alldf.smiles, alldf.split)}
    df = G.generate(seed_smiles, tdb, known, seed_id=seed_id)
    df.insert(df.columns.get_loc("parent_id") + 1, "target", "SARS-CoV-2 Mpro (COVID Moonshot)")
    scored, scorer = False, None
    if len(df) and not a.no_score:
        df, scored, scorer = score(df, seed_smiles, train, a.model_dir)
    elif len(df):
        df["score_type"] = "not_run"
    fun = G.funnel(df)
    fun["scored_by_model"] = int(df["predicted_relative_activity"].notna().sum()) if "predicted_relative_activity" in df else 0
    print("funnel:", fun)

    df.to_csv(os.path.join(a.out_dir, "generated_candidates.csv"), index=False, float_format="%.4f")
    ranked = df[df.filter_status == "pass"].copy()
    if scored and ranked.predicted_relative_activity.notna().any():
        ranked["priority_score_lcb"] = ranked.predicted_relative_activity - ranked.predicted_rel_sd
        ranked = ranked.sort_values("priority_score_lcb", ascending=False)
        ranked.insert(0, "rank", range(1, len(ranked) + 1))
        ranked["reason_retained"] = [
            f"{r.transformation} ({int(r.n_precedent_pairs)} train pairs, observed mean dpIC50 "
            f"{r.precedent_mean_delta_pIC50:+.2f}); Tanimoto to seed {r.parent_tanimoto:.2f}; edit "
            f"{int(r.edit_size_heavy_atoms)} heavy atoms; passed validity/PAINS/SA/element filters"
            for r in ranked.itertuples()]
    ranked.to_csv(os.path.join(a.out_dir, "ranked_candidates.csv"), index=False, float_format="%.4f")

    # ---- retrospective check: ANY valid unique candidate (filtered or not) that is a MEASURED held-out analogue ----
    retro = None
    if seed_label is not None and seed_part == "test" and scorer is not None:
        meas = alldf.set_index("smiles")
        h = df[(df.known_in_moonshot == "test") & (df.valid_rdkit == True) & (df.unique_candidate == True)]  # noqa: E712
        if len(h):
            mean, sd = scorer.pair_delta(h.candidate_smiles.tolist(), [seed_smiles] * len(h))
            rt = pd.DataFrame({"candidate_id": h.candidate_id.values, "candidate_smiles": h.candidate_smiles.values,
                               "transformation": h.transformation.values, "filter_status": h.filter_status.values,
                               "filter_reason": h.filter_reason.values,
                               "measured_delta_vs_seed": [float(meas.loc[c, "pIC50"]) - seed_label for c in h.candidate_smiles],
                               "predicted_delta_vs_seed": mean, "predicted_sd": sd})
            rt.to_csv(os.path.join(a.out_dir, "retrospective_heldout_check.csv"), index=False, float_format="%.4f")
            nz = rt[rt.measured_delta_vs_seed != 0]
            retro = {"n_recovered_heldout_measured_analogues": int(len(rt)),
                     "n_pass_filters": int((rt.filter_status == "pass").sum()),
                     "direction_accuracy": float((np.sign(nz.measured_delta_vs_seed) == np.sign(nz.predicted_delta_vs_seed)).mean()) if len(nz) else None}
            print("retrospective check:", retro)
    fun["retrospective"] = retro
    fun["seed"] = {"id": seed_id, "smiles": seed_smiles, "partition": seed_part,
                   "measured_pIC50": seed_label, "postera_ids": ids,
                   "all_csv_row": a.seed_index}
    fun["params"] = G.PARAMS
    json.dump(fun, open(os.path.join(a.out_dir, "funnel.json"), "w"), indent=2)

    # ---- shortlist.md ----
    L = [f"# CliffForge PoC shortlist — seed `{seed_id}`", "",
         "**PROTOTYPE OUTPUT. Candidates are model-prioritized, not validated. Scores are 2D-model "
         "predicted ΔpIC50 vs the seed (not measured potency, not docking, not binding free energy).**", "",
         f"- Seed SMILES (2D): `{seed_smiles}`",
         f"- Seed partition: `{seed_part}`; measured pIC50: {seed_label if seed_label is None else round(seed_label, 2)}; "
         f"PostEra ID(s): {', '.join(ids) if ids else 'n/a'}",
         "- Target: SARS-CoV-2 Mpro (COVID Moonshot noncovalent set). Pocket-aware docking: **not run**.",
         f"- Funnel: {fun['generated_edits']} attempted edits → {fun['valid_rdkit']} valid → "
         f"{fun['unique_non_parent']} unique non-parent → {fun['pass_all_filters']} pass all filters → "
         f"{fun['scored_by_model']} model-scored", ""]
    if scored and len(ranked):
        L += [f"## Top {min(a.top, len(ranked))} by lower-confidence-bound (predicted Δ − ensemble SD)", "",
              "| # | candidate SMILES | edit | Tan. | pred ΔpIC50 vs seed (±SD) | cliff alert | caveat |",
              "|--:|---|---|--:|--:|---|---|"]
        for r in ranked.head(a.top).itertuples():
            L.append(f"| {r.rank} | `{r.candidate_smiles}` | `{r.transformation}` | {r.parent_tanimoto:.2f} | "
                     f"{r.predicted_relative_activity:+.2f} ± {r.predicted_rel_sd:.2f} | {r.cliff_alert or '-'} | "
                     f"{r.uncertainty_or_caveat} |")
    else:
        L += ["Affinity scoring not run (no model checkpoints / `--no-score`): list below is chemistry-only.", ""]
        for r in ranked.head(a.top).itertuples():
            L.append(f"- `{r.candidate_smiles}`  `{r.transformation}`  Tanimoto {r.parent_tanimoto:.2f}")
    if retro:
        L += ["", "## Retrospective check on measured held-out analogues", "",
              f"{retro['n_recovered_heldout_measured_analogues']} generated candidates ({retro['n_pass_filters']} of them pass the filters) "
              "are real Moonshot compounds from the held-out test partition (never used to mine transformations or to train). "
              f"Direction accuracy of predicted vs measured Δ(candidate − seed): {retro['direction_accuracy']}. "
              "See `retrospective_heldout_check.csv` (tiny n; illustrative only)."]
    open(os.path.join(a.out_dir, "shortlist.md"), "w").write("\n".join(L) + "\n")
    print("wrote:", ", ".join(os.path.join(a.out_dir, f) for f in
                              ("generated_candidates.csv", "ranked_candidates.csv", "shortlist.md", "funnel.json")))


if __name__ == "__main__":
    main()
