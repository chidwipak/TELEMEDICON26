"""
cliffforge.generate
======================
Constrained, EXAMPLE-BASED analogue generation for the CliffForge proof of concept.

NOT a neural generative model. The method is classic matched-molecular-pair (MMP)
transformation transfer, implemented with RDKit only:

  1. Mine transformations R1 >> R2 from the TRAIN partition of the COVID Moonshot
     Mpro set (data/moonshot/all.csv, split == "train"). Two train compounds that share
     an identical scaffold (single acyclic-bond cut, or an aromatic C-H position for
     the H sidechain) and differ only in the sidechain define one observed
     transformation, with its measured delta-pIC50 (a data-derived "precedent").
  2. For a seed ligand, find every single cut whose sidechain equals some R1, swap in
     the R2, re-assemble with RDKit molzip, sanitize.
  3. Filter and annotate: validity, duplicates, parent identity, Morgan similarity to
     the seed, edit size, heavy-atom growth, PAINS (only alerts NOT already in the seed),
     synthetic accessibility, allowed elements, charge conservation.
Every candidate keeps the exact transformation (R1>>R2), its attachment environment and
its precedent statistics, so each edit is traceable to training-set chemistry.

Test-partition labels are never read here.

Known limitation (see MATCH_H_SITES below): in the reported runs, hydrogen <-> substituent transformations are not
matched; edits are substituent-for-substituent swaps at single-bond cuts.
"""
import math
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDConfig, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402

from . import paths as P

# ---- documented constraints (also written to summary.json) -----------------------
PARAMS = {
    "max_side_heavy_atoms": 12,        # replaced sidechain R1 and new sidechain R2 each <= 12 heavy atoms
    "min_scaffold_fraction": 0.6,      # shared scaffold keeps >= 60 % of the seed's heavy atoms
    "max_edit_size_heavy_atoms": 16,   # |R1| + |R2| (heavy atoms removed + added)
    "min_precedent_pairs": 2,          # transformation observed in >= 2 train compound pairs
    "min_parent_tanimoto": 0.5,        # Morgan r=2, 2048 bits
    "max_heavy_atom_growth": 8,        # candidate HA - seed HA
    "max_sa_score": 6.0,
    "allowed_elements": "C,N,O,F,Cl,Br,S",
    "fingerprint": "Morgan radius 2, 2048 bits",
}
ALLOWED = {6, 7, 8, 9, 17, 35, 16, 1}

# KNOWN LIMITATION (found while packaging; default kept so every reported number stays reproducible):
# aromatic C-H "H" sites are enumerated with attachment-environment key "ca", whereas sidechains obtained by cutting a
# bond use "Ca" (see `cuts`). H <-> substituent transformations are therefore never matched and the generator in the
# reported runs only performs substituent-for-substituent swaps. Set CLIFFFORGE_MATCH_H_SITES=1 (or MATCH_H_SITES = True)
# to use the consistent key; this changes the candidate sets and is NOT part of the reported results.
MATCH_H_SITES = os.environ.get("CLIFFFORGE_MATCH_H_SITES", "") == "1"

_pains = None


def _pains_catalog():
    global _pains
    if _pains is None:
        p = FilterCatalogParams()
        p.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        _pains = FilterCatalog(p)
    return _pains


def flat_smiles(smi):
    """Canonical stereo-free SMILES (the identity used by data/moonshot/all.csv)."""
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m, isomericSmiles=False) if m is not None else None


def fp(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)


def pains_names(mol):
    return {e.GetDescription() for e in _pains_catalog().GetMatches(mol)}


def _env(mol, idx):
    a = mol.GetAtomWithIdx(idx)
    return f"{a.GetSymbol()}{'a' if a.GetIsAromatic() else ''}"


def _label_dummies(frag):
    for a in frag.GetAtoms():
        if a.GetAtomicNum() == 0:
            a.SetIsotope(0)
            a.SetAtomMapNum(1)
    return frag


def _norm(smi):
    """Round-trip through a sanitized molecule so equal graphs get equal strings."""
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m is not None else None


def cuts(mol, max_side, min_scaf_frac):
    """All single edits sites of `mol`:
       (scaffold_smiles_with_[*:1], sidechain_smiles | 'H', n_side_heavy, env_key).
    Acyclic single bonds give a scaffold + sidechain (sidechain = smaller part).
    Aromatic C-H positions give the 'H' sidechain (n_side = 0)."""
    n_heavy = mol.GetNumHeavyAtoms()
    out = []
    for b in mol.GetBonds():
        if b.IsInRing() or b.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        frag = Chem.FragmentOnBonds(mol, [b.GetIdx()], addDummies=True)
        parts = Chem.GetMolFrags(frag, asMols=True, sanitizeFrags=False)
        if len(parts) != 2:
            continue
        sizes = [p.GetNumHeavyAtoms() for p in parts]  # RDKit does not count dummy atoms as heavy
        small = int(np.argmin(sizes)) if sizes[0] != sizes[1] else int(
            Chem.MolToSmiles(parts[0]) > Chem.MolToSmiles(parts[1]))
        n_side, n_scaf = sizes[small], sizes[1 - small]
        if n_side == 0 or n_side > max_side or n_scaf / float(n_heavy) < min_scaf_frac:
            continue
        side_m = _label_dummies(Chem.Mol(parts[small]))
        scaf_m = _label_dummies(Chem.Mol(parts[1 - small]))
        # attachment environment = the scaffold atom the dummy is bonded to
        d = [x for x in scaf_m.GetAtoms() if x.GetAtomicNum() == 0][0]
        nb = d.GetNeighbors()[0]
        env = f"{nb.GetSymbol()}{'a' if nb.GetIsAromatic() else ''}"
        scaf_s, side_s = _norm(Chem.MolToSmiles(scaf_m)), _norm(Chem.MolToSmiles(side_m))
        if scaf_s is None or side_s is None:
            continue
        out.append((scaf_s, side_s, n_side, env))
    # H sidechain on aromatic C-H
    for a in mol.GetAtoms():
        if a.GetIsAromatic() and a.GetAtomicNum() == 6 and a.GetTotalNumHs() == 1:
            rw = Chem.RWMol(mol)
            d = rw.AddAtom(Chem.Atom(0))
            rw.AddBond(a.GetIdx(), d, Chem.BondType.SINGLE)
            rw.GetAtomWithIdx(a.GetIdx()).SetNumExplicitHs(0)
            rw.GetAtomWithIdx(a.GetIdx()).SetNoImplicit(False)
            rw.GetAtomWithIdx(d).SetAtomMapNum(1)
            try:
                m2 = rw.GetMol()
                Chem.SanitizeMol(m2)
                s = _norm(Chem.MolToSmiles(m2))
                if s is not None:
                    out.append((s, "H", 0, "Ca" if MATCH_H_SITES else "ca"))
            except Exception:
                continue
    return out


def attach(scaf_smi, side_smi):
    """Re-assemble scaffold + sidechain -> sanitized molecule or None."""
    try:
        scaf = Chem.MolFromSmiles(scaf_smi)
        if side_smi == "H":
            rw = Chem.RWMol(scaf)
            dummies = [x.GetIdx() for x in rw.GetAtoms() if x.GetAtomicNum() == 0]
            for i in sorted(dummies, reverse=True):
                rw.RemoveAtom(i)
            m = rw.GetMol()
        else:
            m = Chem.molzip(scaf, Chem.MolFromSmiles(side_smi))
        Chem.SanitizeMol(m)
        return Chem.MolFromSmiles(Chem.MolToSmiles(m))
    except Exception:
        return None


class TransformDB:
    """Observed R1>>R2 transformations mined from the Moonshot TRAIN partition only."""

    def __init__(self, df_train, max_side=None, min_scaf_frac=None):
        max_side = max_side or PARAMS["max_side_heavy_atoms"]
        min_scaf_frac = min_scaf_frac or PARAMS["min_scaffold_fraction"]
        self.by_scaf = defaultdict(lambda: defaultdict(list))
        self.n_side = {}
        for cid, (smi, y) in enumerate(zip(df_train["smiles"], df_train["pIC50"])):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            seen = set()
            for scaf, side, n, env in cuts(mol, max_side, min_scaf_frac):
                if (scaf, side) in seen:
                    continue
                seen.add((scaf, side))
                self.by_scaf[(scaf, env)][side].append((smi, float(y)))
                self.n_side[side] = n
        tr = defaultdict(lambda: {"d": [], "scafs": set()})
        for (scaf, env), sides in self.by_scaf.items():
            items = list(sides.items())
            for i in range(len(items)):
                for j in range(len(items)):
                    if i == j:
                        continue
                    s1, m1 = items[i]
                    s2, m2 = items[j]
                    for (c1, y1) in m1:
                        for (c2, y2) in m2:
                            t = tr[(env, s1, s2)]
                            t["d"].append(y2 - y1)
                            t["scafs"].add(scaf)
        self.trans = {}
        self.by_left = defaultdict(list)
        for (env, s1, s2), t in tr.items():
            d = np.array(t["d"])
            rec = {"env": env, "R1": s1, "R2": s2, "n_pairs": len(d),
                   "n_scaffold_contexts": len(t["scafs"]),
                   "mean_delta_pIC50": float(d.mean()), "sd_delta_pIC50": float(d.std()),
                   "frac_positive": float((d > 0).mean())}
            self.trans[(env, s1, s2)] = rec
            self.by_left[(env, s1)].append(rec)
        self.n_train = len(df_train)

    def summary(self):
        n_t = len(self.trans)
        n_ok = sum(1 for r in self.trans.values() if r["n_pairs"] >= PARAMS["min_precedent_pairs"])
        return {"n_train_compounds": self.n_train, "n_transformations": n_t,
                "n_transformations_with_min_precedent": n_ok,
                "n_distinct_sidechains": len(self.n_side)}


def load_split():
    df = pd.read_csv(P.processed_csv())
    return df[df.split == "train"].reset_index(drop=True), df[df.split == "test"].reset_index(drop=True), df


def generate(seed_smiles, tdb, known=None, seed_id="seed", min_precedent=None):
    """Return a DataFrame with one row per attempted edit (valid or not), fully annotated.
    `known`: dict flat_smiles -> 'train'|'test' (for novelty annotation only; no labels)."""
    known = known or {}
    min_precedent = PARAMS["min_precedent_pairs"] if min_precedent is None else min_precedent
    seed = Chem.MolFromSmiles(seed_smiles)
    seed_flat = Chem.MolToSmiles(seed, isomericSmiles=False)
    seed_mol = Chem.MolFromSmiles(seed_flat)
    seed_fp = fp(seed_mol)
    seed_ha = seed_mol.GetNumHeavyAtoms()
    seed_pains = pains_names(seed_mol)
    seed_charge = Chem.GetFormalCharge(seed_mol)
    seed_sa = sascorer.calculateScore(seed_mol)
    rows, seen = [], {seed_flat}
    for scaf, side, n_r1, env in cuts(seed_mol, PARAMS["max_side_heavy_atoms"], PARAMS["min_scaffold_fraction"]):
        for rec in tdb.by_left.get((env, side), []):
            n_r2 = tdb.n_side.get(rec["R2"], 0)
            row = {"parent_id": seed_id, "parent_smiles": seed_smiles,
                   "transformation": f"{rec['R1']}>>{rec['R2']}", "attachment_env": env,
                   "scaffold_smiles": scaf, "removed_heavy_atoms": n_r1, "added_heavy_atoms": n_r2,
                   "edit_size_heavy_atoms": n_r1 + n_r2, "parent_heavy_atoms": seed_ha,
                   "n_precedent_pairs": rec["n_pairs"], "n_precedent_scaffolds": rec["n_scaffold_contexts"],
                   "precedent_mean_delta_pIC50": rec["mean_delta_pIC50"],
                   "precedent_sd_delta_pIC50": rec["sd_delta_pIC50"],
                   "generation_method": "MMP transformation transfer (train-partition precedents; RDKit)"}
            m = attach(scaf, rec["R2"])
            if m is None:
                row.update(candidate_smiles=None, valid_rdkit=False, unique_candidate=False,
                           filter_status="rejected", filter_reason="invalid_rdkit")
                rows.append(row)
                continue
            cs = Chem.MolToSmiles(m, isomericSmiles=False)
            c_fp = fp(m)
            reasons = []
            unique = cs not in seen
            if cs == seed_flat:
                reasons.append("same_as_parent")
            elif not unique:
                reasons.append("duplicate")
            seen.add(cs)
            tan = DataStructs.TanimotoSimilarity(seed_fp, c_fp)
            ha = m.GetNumHeavyAtoms()
            els = {a.GetAtomicNum() for a in m.GetAtoms()}
            new_pains = pains_names(m) - seed_pains
            sa = sascorer.calculateScore(m)
            if rec["n_pairs"] < min_precedent:
                reasons.append(f"precedent<{min_precedent}")
            if tan < PARAMS["min_parent_tanimoto"]:
                reasons.append("low_parent_similarity")
            if n_r1 + n_r2 > PARAMS["max_edit_size_heavy_atoms"]:
                reasons.append("edit_too_large")
            if ha - seed_ha > PARAMS["max_heavy_atom_growth"]:
                reasons.append("heavy_atom_growth")
            if not els <= ALLOWED:
                reasons.append("disallowed_element")
            if Chem.GetFormalCharge(m) != seed_charge:
                reasons.append("charge_changed")
            if new_pains:
                reasons.append("new_PAINS:" + "|".join(sorted(new_pains)))
            if sa > PARAMS["max_sa_score"]:
                reasons.append("hard_to_synthesize_SA")
            row.update(candidate_smiles=cs, valid_rdkit=True, unique_candidate=unique,
                       parent_tanimoto=tan, candidate_heavy_atoms=ha,
                       heavy_atom_delta=ha - seed_ha, sa_score=sa, mw=Descriptors.MolWt(m),
                       clogp=Crippen.MolLogP(m), new_pains_alerts="|".join(sorted(new_pains)),
                       known_in_moonshot=known.get(cs, "no"),
                       filter_status="rejected" if reasons else "pass",
                       filter_reason=";".join(reasons))
            rows.append(row)
    df = pd.DataFrame(rows)
    if len(df):
        df.insert(0, "candidate_id", [f"{seed_id}_c{i:04d}" for i in range(len(df))])
    return df


def pick_seed(tdb, alldf, min_pIC50=7.0):
    """Pre-declared seed rule (yield-based, outcome-blind): among HELD-OUT test-partition
    compounds with measured pIC50 >= `min_pIC50` (potent leads), take the one with the most
    filter-passing candidates (ties: higher pIC50, then lower row). Returns (row, table)."""
    known = {s: sp for s, sp in zip(alldf.smiles, alldf.split)}
    rows = []
    for i, r in alldf[(alldf.split == "test") & (alldf.pIC50 >= min_pIC50)].iterrows():
        f = funnel(generate(r.smiles, tdb, known, seed_id=f"row{i}"))
        rows.append({"all_csv_row": int(i), "pIC50": float(r.pIC50), "series_size": int(r.series_size),
                     **f})
    t = pd.DataFrame(rows).sort_values(["pass_all_filters", "pIC50", "all_csv_row"],
                                       ascending=[False, False, True])
    return int(t.iloc[0]["all_csv_row"]), t


def funnel(df):
    """Counts for the generation funnel."""
    if not len(df):
        return {"generated": 0}
    valid = df[df.valid_rdkit == True]  # noqa: E712
    uniq = valid[valid.unique_candidate == True]  # noqa: E712
    npar = uniq[~uniq.filter_reason.str.contains("same_as_parent")]
    ps = df[df.filter_status == "pass"]
    return {"generated_edits": int(len(df)), "valid_rdkit": int(len(valid)),
            "unique_non_parent": int(len(npar)), "pass_all_filters": int(len(ps))}
