"""
validate_inputs.py
==================
Input validation for the /api/predict endpoint.

Validates:
  - Candidate and reference SMILES (parseable, not single-atom/no-bond,
    preserves stereochemistry, max length)
  - Reference activity value (finite, positive for IC50 µM, not excessive)
  - Unit choice (pIC50 or ic50_um)

pIC50 / IC50 µM conversion:
  pIC50 = 6 - log10(IC50_µM)
  IC50_µM = 10^(6 - pIC50)

Positive Δ = candidate is MORE potent than reference.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

MAX_SMILES_LEN = 512
MIN_IC50_UM = 1e-6       # ~12 pIC50 max
MAX_IC50_UM = 1e6        # 0 pIC50 min — anything larger is outside the model's domain
MIN_PICSO = -6.0         # extremely high IC50; allowed but suspicious
MAX_PICSO = 14.0         # unrealistically potent


class ValidationError(ValueError):
    """Raised when an input fails validation."""


@dataclass
class ValidatedInputs:
    candidate_smiles: str       # Canonical, with stereo preserved
    reference_smiles: str       # Canonical, with stereo preserved
    reference_pIC50: float      # Always in pIC50 units internally
    unit: str                   # "pIC50" or "ic50_um"


def _parse_smiles(smi: str, label: str) -> str:
    """
    Parse SMILES, returning RDKit canonical SMILES (stereo preserved).
    Rejects None, empty, too long, unparseable, single-atom/no-bond molecules.
    """
    if not smi or not isinstance(smi, str):
        raise ValidationError(f"{label}: SMILES is empty or not a string")
    smi = smi.strip()
    if len(smi) > MAX_SMILES_LEN:
        raise ValidationError(
            f"{label}: SMILES too long ({len(smi)} chars, max {MAX_SMILES_LEN})"
        )

    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValidationError(f"{label}: could not parse SMILES: {smi!r}")

    n_atoms = mol.GetNumAtoms()
    n_bonds = mol.GetNumBonds()

    if n_atoms == 0:
        raise ValidationError(f"{label}: molecule has no atoms")

    if n_atoms == 1 and n_bonds == 0:
        raise ValidationError(
            f"{label}: single-atom structure (no bonds); "
            "the CliffForge encoder cannot safely process this"
        )

    if n_bonds == 0:
        raise ValidationError(
            f"{label}: structure has no bonds; "
            "the CliffForge encoder requires at least one bond"
        )

    # Return canonical SMILES with stereochemistry preserved
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    return canonical


def _validate_activity(value: float, unit: str) -> float:
    """
    Validate the reference activity and return it as pIC50.
    """
    if not math.isfinite(value):
        raise ValidationError(
            f"Reference activity ({value}) is not finite"
        )

    if unit == "ic50_um":
        if value <= 0:
            raise ValidationError(
                f"IC50 µM value must be positive, got {value}"
            )
        if value > MAX_IC50_UM:
            raise ValidationError(
                f"IC50 µM value {value} is unrealistically large (max {MAX_IC50_UM})"
            )
        pIC50 = 6.0 - math.log10(value)
    elif unit == "pIC50":
        pIC50 = value
    else:
        raise ValidationError(f"Unknown unit: {unit!r}. Must be 'pIC50' or 'ic50_um'")

    if not math.isfinite(pIC50):
        raise ValidationError(
            f"Converted pIC50 is not finite (got {pIC50} from {value} {unit})"
        )
    if pIC50 < MIN_PICSO or pIC50 > MAX_PICSO:
        raise ValidationError(
            f"pIC50 value {pIC50:.3f} is outside the valid range "
            f"[{MIN_PICSO}, {MAX_PICSO}]"
        )

    return pIC50


def validate_predict_inputs(
    candidate_smiles: str,
    reference_smiles: str,
    reference_activity: float,
    unit: str,
) -> ValidatedInputs:
    """
    Full validation for the /api/predict endpoint.
    Raises ValidationError with a human-readable message on failure.
    Returns ValidatedInputs with canonical SMILES and pIC50.
    """
    unit = (unit or "").strip().lower().replace("-", "_")
    if unit not in ("pic50", "ic50_um"):
        # Accept both "pIC50" / "pic50" and "ic50_um" / "IC50_uM" etc.
        if unit in ("pic50", "pic_50"):
            unit = "pIC50"
        elif unit in ("ic50_um", "ic50um", "ic50"):
            unit = "ic50_um"
        else:
            raise ValidationError(
                f"Unknown unit {unit!r}. Use 'pIC50' or 'ic50_um'"
            )

    # Normalise
    if unit == "pic50":
        unit = "pIC50"

    cand_canonical = _parse_smiles(candidate_smiles, "Candidate")
    ref_canonical = _parse_smiles(reference_smiles, "Reference")

    ref_pIC50 = _validate_activity(reference_activity, unit)

    return ValidatedInputs(
        candidate_smiles=cand_canonical,
        reference_smiles=ref_canonical,
        reference_pIC50=ref_pIC50,
        unit=unit,
    )


def pIC50_to_IC50_uM(pIC50: float) -> float:
    """Convert pIC50 to IC50 in µM."""
    return 10.0 ** (6.0 - pIC50)


def IC50_uM_to_pIC50(ic50_um: float) -> float:
    """Convert IC50 µM to pIC50."""
    return 6.0 - math.log10(ic50_um)
