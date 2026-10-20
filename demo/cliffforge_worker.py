"""
cliffforge_worker.py
=====================
CliffForge inference for the demo. Loads three fold-3 SiameseAVES checkpoints
(one per seed) onto CPU at startup. Provides antisymmetrised pair-delta inference.

Architecture (from protocol.json + relmodel.py):
  - NODE_DIM = 48, EDGE_DIM = 8, hidden = 128, depth = 4, dropout = 0.1
  - Three seeds: 0, 1, 2

Formula (from notebook fold inference):
  Δ_CF,s = (f_s(candidate, reference) - f_s(reference, candidate)) / 2
  mean = mean over seeds
  pop_sd = population standard deviation over seeds (ddof=0)
  predicted_candidate_pIC50 = reference_measured_pIC50 + mean(Δ_CF)
"""
from __future__ import annotations

import sys
import os

# Add repo src to path so we can import cliffforge.{features,networks}
_REPO_SRC = os.path.join(os.path.dirname(__file__), "..", "src")
if _REPO_SRC not in sys.path:
    sys.path.insert(0, _REPO_SRC)

import numpy as np
import torch
from rdkit import Chem, RDLogger
from torch_geometric.data import Batch, Data

from cliffforge.features import mol_to_whole_graph, NODE_DIM, EDGE_DIM
from cliffforge.networks import SiameseAVES

RDLogger.DisableLog("rdApp.*")

DEVICE = torch.device("cpu")  # CPU-only for local demo

# Verify node/edge dimensions match protocol
assert NODE_DIM == 48, f"Expected NODE_DIM=48, got {NODE_DIM}"
assert EDGE_DIM == 8, f"Expected EDGE_DIM=8, got {EDGE_DIM}"


def _smiles_to_data(smi: str) -> Data | None:
    """Convert SMILES to a torch_geometric Data object (same as relmodel.smiles_to_data)."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    g = mol_to_whole_graph(mol)
    if g is None:
        return None
    return Data(
        x=torch.from_numpy(g["x"]).float(),
        edge_index=torch.from_numpy(g["edge_index"]).long(),
        edge_attr=torch.from_numpy(g["edge_attr"]).float(),
    )


class CliffForgeWorker:
    """
    Loads three SiameseAVES checkpoints for fold 3 and performs
    antisymmetrised pair-delta inference.
    """

    def __init__(self, checkpoint_paths: dict[str, str]):
        """
        checkpoint_paths: dict with keys 'cf_seed_0', 'cf_seed_1', 'cf_seed_2'
          mapping to absolute paths of final.pt files.
        """
        self.models: list[SiameseAVES] = []
        self.checkpoint_info: list[dict] = []

        for seed in range(3):
            key = f"cf_seed_{seed}"
            path = checkpoint_paths[key]
            model = SiameseAVES(NODE_DIM, EDGE_DIM, hidden=128, depth=4, dropout=0.1)
            # Use weights_only=False for checkpoint compatibility; scoped to trusted verified files
            state = torch.load(path, map_location=DEVICE, weights_only=False)
            # Checkpoints saved with state_dict at top level or under 'state_dict' key
            if "state_dict" in state:
                model.load_state_dict(state["state_dict"])
            else:
                # Saved directly as state_dict
                model.load_state_dict(state)
            model.eval()
            self.models.append(model)
            self.checkpoint_info.append({"seed": seed, "path": path})

        self._cache: dict[str, Data | None] = {}

    def _graph(self, smi: str) -> Data | None:
        if smi not in self._cache:
            self._cache[smi] = _smiles_to_data(smi)
        return self._cache[smi]

    @torch.no_grad()
    def _forward_one(self, model: SiameseAVES, smi_a: str, smi_b: str) -> float:
        """Raw model output f(a, b): scalar."""
        ga = self._graph(smi_a)
        gb = self._graph(smi_b)
        if ga is None or gb is None:
            raise ValueError(f"Could not featurize SMILES: {smi_a!r} or {smi_b!r}")
        ba = Batch.from_data_list([ga]).to(DEVICE)
        bb = Batch.from_data_list([gb]).to(DEVICE)
        out = model(ba, bb)
        return float(out.squeeze().cpu())

    def predict_pair(
        self,
        candidate_smiles: str,
        reference_smiles: str,
        reference_pIC50: float,
    ) -> dict:
        """
        Compute antisymmetrised Δ_CF for each seed, plus derived quantities.

        Returns dict with:
          seed_deltas: [Δ_CF,0, Δ_CF,1, Δ_CF,2]
          mean_delta: float
          pop_sd: float   (population SD, ddof=0)
          predicted_candidate_pIC50: reference_pIC50 + mean_delta
          raw_forward: [f_s(cand,ref) for s in seeds]
          raw_backward: [f_s(ref,cand) for s in seeds]
        """
        seed_deltas = []
        raw_fwd = []
        raw_bwd = []

        for model in self.models:
            fwd = self._forward_one(model, candidate_smiles, reference_smiles)
            bwd = self._forward_one(model, reference_smiles, candidate_smiles)
            delta = (fwd - bwd) / 2.0
            raw_fwd.append(fwd)
            raw_bwd.append(bwd)
            seed_deltas.append(delta)

        arr = np.array(seed_deltas, dtype=np.float64)
        mean_delta = float(arr.mean())
        pop_sd = float(arr.std(ddof=0))  # population SD
        predicted_candidate_pIC50 = reference_pIC50 + mean_delta

        return {
            "seed_deltas": seed_deltas,
            "mean_delta": mean_delta,
            "pop_sd": pop_sd,
            "predicted_candidate_pIC50": predicted_candidate_pIC50,
            "raw_forward": raw_fwd,
            "raw_backward": raw_bwd,
        }
