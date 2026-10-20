"""
chemprop_worker.py
===================
ChemProp 1.6.1 inference for the demo. Loads three fold-3 model.pt checkpoints
at startup. Predicts absolute pIC50 for a candidate SMILES, one seed at a time,
so individual seed predictions are returned (not just an ensemble mean).

Environment note:
  ChemProp 1.6.1 requires NumPy 1.x. This local deployment was checked with
  NumPy 1.26.4; its live predictions are compared against the Colab archive.

Security note:
  TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1 is set only within this module's scope,
  scoped to the trusted, hash-verified ChemProp checkpoints.

Formula (from notebook fold inference):
  Δ_CP,s = predicted_candidate_pIC50_s - reference_measured_pIC50
  mean_delta = mean over seeds
  (ChemProp uses the measured reference value, same as the five-fold notebook.)
"""
from __future__ import annotations

import os
import logging
import tempfile
import threading
from pathlib import Path

import numpy as np

from chemprop.train.make_predictions import make_predictions, load_model
from chemprop.args import PredictArgs

# Silence chemprop verbose logging
logging.getLogger("chemprop").setLevel(logging.WARNING)
logging.getLogger("hyperopt").setLevel(logging.WARNING)


class ChemPropWorker:
    """
    Loads three ChemProp 1.6.1 checkpoints for fold 3 individually.
    Runs per-seed predictions to get individual pIC50 values.
    """

    def __init__(self, checkpoint_paths: dict[str, str]):
        """
        checkpoint_paths: dict with keys 'cp_seed_0', 'cp_seed_1', 'cp_seed_2'
          mapping to absolute paths of model.pt files.
        """
        self.model_objects_list: list[tuple] = []
        self._lock = threading.Lock()

        for seed in range(3):
            key = f"cp_seed_{seed}"
            path = checkpoint_paths[key]

            # Build minimal PredictArgs for this single checkpoint
            args = PredictArgs()
            args.checkpoint_paths = [path]
            args.num_workers = 0
            args.batch_size = 50
            args.no_cuda = True
            args.gpu = None
            # Uncertainty-related (not needed)
            args.uncertainty_method = None
            args.calibration_method = None
            args.evaluation_methods = None
            args.calibration_path = None

            # load_model handles update_prediction_args which syncs train_args
            old_load_setting = os.environ.get("TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD")
            os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
            try:
                model_objects = load_model(args, generator=False)
            finally:
                if old_load_setting is None:
                    os.environ.pop("TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD", None)
                else:
                    os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = old_load_setting
            # model_objects = (args, train_args, [models], [scalers], num_tasks, task_names)
            self.model_objects_list.append(model_objects)

    def predict_single_seed(self, smiles: str, seed_idx: int) -> float:
        """
        Run ChemProp inference for one SMILES with one seed's model.
        Returns the prediction in pIC50 units (after inverse scaling).
        """
        model_objects = self.model_objects_list[seed_idx]
        args = model_objects[0]

        # ChemProp 1.6.1 writes a prediction CSV even when SMILES are passed in
        # memory. Give each call its own temporary path to prevent collisions.
        with tempfile.TemporaryDirectory(prefix="cliffforge_cp_") as tmp:
            args.preds_path = str(Path(tmp) / "predictions.csv")
            preds = make_predictions(
                args=args,
                smiles=[[smiles]],
                model_objects=model_objects,
            )
        # preds is List[List[Optional[float]]], shape [n_smiles, n_tasks]
        value = preds[0][0]
        if value is None:
            raise ValueError(f"ChemProp returned None for SMILES: {smiles!r}")
        return float(value)

    def predict_candidate(
        self,
        candidate_smiles: str,
        reference_measured_pIC50: float,
    ) -> dict:
        """
        Predict candidate pIC50 with each seed. Compute per-seed delta vs reference.

        Returns dict with:
          seed_pIC50: [pIC50_seed0, pIC50_seed1, pIC50_seed2]
          seed_deltas: [Δ_CP,0, Δ_CP,1, Δ_CP,2]  (candidate_pred - reference_measured)
          mean_delta: float
          mean_pIC50: float
        """
        seed_pIC50 = []
        seed_deltas = []

        # ChemProp updates process-wide featurization parameters during predict.
        # Serialize requests so two browser clicks cannot change those settings
        # or an args.preds_path while another prediction is still running.
        with self._lock:
            for seed_idx in range(3):
                pred = self.predict_single_seed(candidate_smiles, seed_idx)
                delta = pred - reference_measured_pIC50
                seed_pIC50.append(pred)
                seed_deltas.append(delta)

        mean_delta = float(np.mean(seed_deltas))
        mean_pIC50 = float(np.mean(seed_pIC50))

        return {
            "seed_pIC50": seed_pIC50,
            "seed_deltas": seed_deltas,
            "mean_delta": mean_delta,
            "mean_pIC50": mean_pIC50,
            "pop_sd": float(np.std(seed_deltas, ddof=0)),
        }
