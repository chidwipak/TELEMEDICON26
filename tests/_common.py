"""Shared helpers for the test-suite (no external data or weights required)."""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DATA_CSV = Path(os.environ.get("CLIFFFORGE_DATA", ROOT / "data")) / "moonshot" / "all.csv"
MODEL_DIR = Path(os.environ.get("CLIFFFORGE_MODELS", ROOT / "models"))
HAVE_DATA = DATA_CSV.exists()
HAVE_WEIGHTS = all((MODEL_DIR / f"rel_siamese_s{s}.pt").exists() for s in (0, 1, 2))
SKIP_SLOW = os.environ.get("CLIFFFORGE_SKIP_SLOW", "") == "1"

# Synthetic TOY fixtures used only to unit-test code paths. They are NOT experimental data and are never used for results.
TOY_FS_PRECEDENT = [  # para-F -> para-Cl observed in three different scaffold contexts
    ("NC(=O)c1ccc(F)cc1", 5.0), ("NC(=O)c1ccc(Cl)cc1", 5.6),
    ("CNC(=O)c1ccc(F)cc1", 5.1), ("CNC(=O)c1ccc(Cl)cc1", 5.8),
    ("CCNC(=O)c1ccc(F)cc1", 5.2), ("CCNC(=O)c1ccc(Cl)cc1", 5.9),
]
TOY_H_PRECEDENT = [  # H -> Cl at an aromatic C-H, observed in two scaffold contexts
    ("NC(=O)c1ccccc1", 4.0), ("NC(=O)c1ccc(Cl)cc1", 4.9),
    ("CNC(=O)c1ccccc1", 4.2), ("CNC(=O)c1ccc(Cl)cc1", 5.0),
]

REQUIRED_COLUMNS = [
    "candidate_id", "parent_id", "parent_smiles", "candidate_smiles", "transformation", "valid_rdkit",
    "unique_candidate", "parent_tanimoto", "parent_heavy_atoms", "candidate_heavy_atoms",
    "edit_size_heavy_atoms", "removed_heavy_atoms", "added_heavy_atoms", "heavy_atom_delta",
    "generation_method", "filter_status", "filter_reason",
]


def toy_db(rows):
    import pandas as pd
    from cliffforge import generate as G
    return G.TransformDB(pd.DataFrame(rows, columns=["smiles", "pIC50"]))
