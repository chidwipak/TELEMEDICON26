"""Whole-molecule 2D graph featurization (ChemProp-style) used by the relative-affinity model.

Extracted unchanged from the internship research code (featurize.py); the fragment/ghost-anchor channels of that
research are not used by CliffForge, so the corresponding node dimensions stay zero. Keep NODE_DIM/EDGE_DIM
unchanged: they define the shapes stored in the model checkpoints.
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdPartialCharges

ATOM_LIST = [5, 6, 7, 8, 9, 15, 16, 17, 35, 53]  # B C N O F P S Cl Br I
DEGREE_LIST = [0, 1, 2, 3, 4, 5]
CHARGE_LIST = [-2, -1, 0, 1, 2]
NUMH_LIST = [0, 1, 2, 3, 4]
HYB_LIST = [
    Chem.HybridizationType.SP,
    Chem.HybridizationType.SP2,
    Chem.HybridizationType.SP3,
    Chem.HybridizationType.SP3D,
    Chem.HybridizationType.SP3D2,
]

N_GHOST_DESC = 5  # inductive, resonance, electronegativity, polarizability, excluded-volume


def _onehot(value, choices):
    vec = [0.0] * (len(choices) + 1)
    try:
        vec[choices.index(value)] = 1.0
    except ValueError:
        vec[-1] = 1.0
    return vec


def atom_base_features(atom: Chem.Atom) -> list:
    """Standard atom features shared by every node (real atoms AND ghosts)."""
    feats = []
    feats += _onehot(atom.GetAtomicNum(), ATOM_LIST)          # 11
    feats += _onehot(atom.GetTotalDegree(), DEGREE_LIST)      # 7
    feats += _onehot(atom.GetFormalCharge(), CHARGE_LIST)     # 6
    feats += _onehot(atom.GetTotalNumHs(), NUMH_LIST)         # 6
    feats += _onehot(atom.GetHybridization(), HYB_LIST)       # 6
    feats.append(1.0 if atom.GetIsAromatic() else 0.0)        # 1
    feats.append(1.0 if atom.IsInRing() else 0.0)             # 1
    gc = atom.GetDoubleProp("_GasteigerCharge") if atom.HasProp("_GasteigerCharge") else 0.0
    if not np.isfinite(gc):
        gc = 0.0
    feats.append(float(gc))                                   # 1
    return feats   # length 39


ATOM_BASE_DIM = 39
# node layout: base(39) + is_ghost(1) + ghost_desc(5) + local_xyz(3)
NODE_DIM = ATOM_BASE_DIM + 1 + N_GHOST_DESC + 3


def bond_features(bond: Chem.Bond, length: float = 0.0) -> list:
    """Edge features. The 5th flag marks a *virtual* ghost-anchor edge."""
    if bond is None:  # virtual ghost edge
        return [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, float(length)]
    btype = bond.GetBondType()
    return [
        1.0 if btype == Chem.BondType.SINGLE else 0.0,
        1.0 if btype == Chem.BondType.DOUBLE else 0.0,
        1.0 if btype == Chem.BondType.TRIPLE else 0.0,
        1.0 if btype == Chem.BondType.AROMATIC else 0.0,
        0.0,  # virtual flag
        1.0 if bond.GetIsConjugated() else 0.0,
        1.0 if bond.IsInRing() else 0.0,
        float(length),
    ]   # length 8


EDGE_DIM = 8


def compute_gasteiger(mol: Chem.Mol) -> None:
    try:
        rdPartialCharges.ComputeGasteigerCharges(mol)
    except Exception:
        for a in mol.GetAtoms():
            a.SetDoubleProp("_GasteigerCharge", 0.0)


def mol_to_whole_graph(mol: Chem.Mol) -> dict | None:
    """Featurize a whole molecule into a 2D topological graph (ChemProp-style).
    Ghost / 3D channels stay zero for the baseline."""
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    mol = Chem.Mol(mol)
    compute_gasteiger(mol)
    n = mol.GetNumAtoms()
    x = np.zeros((n, NODE_DIM), dtype=np.float32)
    for atom in mol.GetAtoms():
        x[atom.GetIdx(), :ATOM_BASE_DIM] = atom_base_features(atom)
    src, dst, eattr = [], [], []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bf = bond_features(bond)
        src += [i, j]
        dst += [j, i]
        eattr += [bf, bf]
    if len(src) == 0:
        src, dst, eattr = [0], [0], [bond_features(None)]
    return {
        "x": x,
        "edge_index": np.array([src, dst], dtype=np.int64),
        "edge_attr": np.array(eattr, dtype=np.float32),
    }
