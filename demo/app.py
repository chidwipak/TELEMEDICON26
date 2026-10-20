"""Local, live CliffForge and ChemProp evaluator demo."""
from __future__ import annotations

import importlib.metadata
import logging
import sys
import time
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

DEMO_DIR = Path(__file__).resolve().parent
REPO_ROOT = DEMO_DIR.parent
for entry in (str(REPO_ROOT / "src"), str(DEMO_DIR)):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from checkpoint_verify import checkpoint_absolute_paths, verify_checkpoints
from chemprop_worker import ChemPropWorker
from cliffforge_worker import CliffForgeWorker
from validate_inputs import ValidationError, pIC50_to_IC50_uM, validate_predict_inputs

log = logging.getLogger("cliffforge.demo")
logging.basicConfig(level=logging.INFO)

EXAMPLE = {
    "source_pair_id": 10765,
    "fold": 3,
    "setting": "test_train",
    "candidate_molecule_id": 1655,
    "reference_molecule_id": 1614,
    "candidate_smiles": "O=C(Cc1cc(Cl)cc(OC2CC(=O)N2)c1)Nc1cncc2ccccc12",
    "reference_smiles": "O=C(Cc1cc(Cl)cc(OC2CCC(=O)N2)c1)Nc1cncc2ccccc12",
    "reference_activity": 85.8740204655224,
    "unit": "ic50_um",
    "reference_pIC50": 4.066138203747178,
    "candidate_measured_pIC50": 6.673440265206619,
    "candidate_measured_ic50_um": 0.212109311760226,
    "true_delta": 2.607302061459441,
}

app = FastAPI(title="CliffForge + ChemProp local demo")
app.mount("/static", StaticFiles(directory=str(DEMO_DIR / "static")), name="static")
state = {"status": "starting", "error": None, "checkpoints": [], "cf": None, "cp": None}


@app.on_event("startup")
def startup() -> None:
    try:
        state["checkpoints"] = verify_checkpoints()
        paths = checkpoint_absolute_paths()
        state["cf"] = CliffForgeWorker(paths)
        state["cp"] = ChemPropWorker(paths)
        state["status"] = "ready"
        log.info("Both models are ready: fold 3, three seeds each")
    except Exception as exc:
        state["status"] = "error"
        state["error"] = str(exc)
        log.exception("Demo startup failed")


@app.get("/")
def index() -> FileResponse:
    return FileResponse(DEMO_DIR / "static" / "index.html")


@app.get("/api/health")
def health() -> dict:
    result = {
        "status": state["status"],
        "error": state["error"],
        "fold": 3,
        "models_loaded": {
            "cliffforge": state["cf"] is not None,
            "chemprop": state["cp"] is not None,
        },
        "checkpoints": [
            {"path": row["relative_path"], "sha256": row["actual_sha256"]}
            for row in state["checkpoints"]
        ],
        "versions": {
            "chemprop": importlib.metadata.version("chemprop"),
            "rdkit": rdBase.rdkitVersion,
        },
    }
    if state["status"] != "ready":
        raise HTTPException(status_code=503, detail=result)
    return result


@app.get("/api/example")
def example() -> dict:
    """Input values only; measured candidate outcome is revealed separately."""
    return {key: EXAMPLE[key] for key in (
        "source_pair_id", "fold", "setting", "candidate_molecule_id",
        "reference_molecule_id", "candidate_smiles", "reference_smiles",
        "reference_activity", "unit", "reference_pIC50",
    )}


@app.get("/api/example/reveal")
def reveal_example() -> dict:
    return {
        "source_pair_id": EXAMPLE["source_pair_id"],
        "candidate_measured_pIC50": EXAMPLE["candidate_measured_pIC50"],
        "candidate_measured_ic50_um": EXAMPLE["candidate_measured_ic50_um"],
        "reference_measured_pIC50": EXAMPLE["reference_pIC50"],
        "true_delta": EXAMPLE["true_delta"],
        "measured_class": "AC",
        "selection_note": "This held-out, known-reference pair was selected because CliffForge correctly detected its cliff.",
    }


class PredictRequest(BaseModel):
    candidate_smiles: str
    reference_smiles: str
    reference_activity: float
    unit: str = "pIC50"


def molecule_view(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not draw the validated molecule")
    drawer = rdMolDraw2D.MolDraw2DSVG(600, 330)
    drawer.drawOptions().padding = 0.075
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    graph_mol = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(graph_mol)
    conformer = graph_mol.GetConformer()
    points = [conformer.GetAtomPosition(i) for i in range(graph_mol.GetNumAtoms())]
    x_min, x_max = min(p.x for p in points), max(p.x for p in points)
    y_min, y_max = min(p.y for p in points), max(p.y for p in points)
    scale = min(520 / max(x_max - x_min, 1), 250 / max(y_max - y_min, 1))
    coords = [
        (300 + (p.x - (x_min + x_max) / 2) * scale,
         165 - (p.y - (y_min + y_max) / 2) * scale)
        for p in points
    ]
    bonds = []
    for bond in graph_mol.GetBonds():
        x1, y1 = coords[bond.GetBeginAtomIdx()]
        x2, y2 = coords[bond.GetEndAtomIdx()]
        width = 2.5 if bond.GetBondTypeAsDouble() >= 2 else 1.5
        bonds.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#9bb4b1" stroke-width="{width}"/>')
    nodes = []
    for atom in graph_mol.GetAtoms():
        x, y = coords[atom.GetIdx()]
        symbol = atom.GetSymbol()
        color = {"O": "#d65056", "N": "#5467ca", "Cl": "#5ba663", "F": "#5ba663"}.get(symbol, "#087d72")
        nodes.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="white" stroke="{color}" stroke-width="2"/>')
        if symbol != "C":
            nodes.append(f'<text x="{x:.1f}" y="{y-10:.1f}" text-anchor="middle" font-size="11" font-family="Arial,sans-serif" fill="{color}">{symbol}</text>')
    graph_svg = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 600 330" role="img" aria-label="Molecular graph">' + ''.join(bonds + nodes) + '</svg>'
    return {
        "svg": drawer.GetDrawingText(),
        "graph_svg": graph_svg,
        "atoms": mol.GetNumAtoms(),
        "bonds": mol.GetNumBonds(),
        "directed_edges": 2 * mol.GetNumBonds(),
    }


def result_for(delta: float, seed_deltas: list[float], sd: float,
               candidate_pic50: float) -> dict:
    return {
        "seed_deltas": seed_deltas,
        "mean_delta": delta,
        "seed_sd": sd,
        "predicted_candidate_pIC50": candidate_pic50,
        "predicted_candidate_ic50_um": pIC50_to_IC50_uM(candidate_pic50),
        "potency_direction": "more_potent" if delta > 0 else "less_potent" if delta < 0 else "same",
        "predicted_large_difference": abs(delta) > 1.5,
    }


@app.post("/api/predict")
def predict(req: PredictRequest) -> dict:
    if state["status"] != "ready":
        raise HTTPException(status_code=503, detail="Models are not ready; check /api/health")
    start = time.perf_counter()
    try:
        valid = validate_predict_inputs(
            req.candidate_smiles, req.reference_smiles, req.reference_activity, req.unit
        )
    except ValidationError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    validation_ms = (time.perf_counter() - start) * 1000

    try:
        candidate_view = molecule_view(valid.candidate_smiles)
        reference_view = molecule_view(valid.reference_smiles)
        drawing_ms = (time.perf_counter() - start) * 1000 - validation_ms

        cf_start = time.perf_counter()
        cf = state["cf"].predict_pair(
            valid.candidate_smiles, valid.reference_smiles, valid.reference_pIC50
        )
        cf_ms = (time.perf_counter() - cf_start) * 1000

        cp_start = time.perf_counter()
        cp = state["cp"].predict_candidate(valid.candidate_smiles, valid.reference_pIC50)
        cp_ms = (time.perf_counter() - cp_start) * 1000
    except Exception as exc:
        log.exception("Inference failed")
        raise HTTPException(status_code=500, detail=f"Inference failed: {exc}") from exc

    cf_out = result_for(cf["mean_delta"], cf["seed_deltas"], cf["pop_sd"],
                        cf["predicted_candidate_pIC50"])
    cf_out["raw_candidate_reference"] = cf["raw_forward"]
    cf_out["raw_reference_candidate"] = cf["raw_backward"]
    cp_out = result_for(cp["mean_delta"], cp["seed_deltas"], cp["pop_sd"], cp["mean_pIC50"])
    cp_out["seed_candidate_pIC50"] = cp["seed_pIC50"]

    example_candidate = Chem.MolToSmiles(Chem.MolFromSmiles(EXAMPLE["candidate_smiles"]), isomericSmiles=True)
    example_reference = Chem.MolToSmiles(Chem.MolFromSmiles(EXAMPLE["reference_smiles"]), isomericSmiles=True)
    is_example = (
        valid.candidate_smiles == example_candidate
        and valid.reference_smiles == example_reference
        and abs(valid.reference_pIC50 - EXAMPLE["reference_pIC50"]) < 1e-6
    )
    return {
        "fold": 3,
        "candidate_smiles": valid.candidate_smiles,
        "reference_smiles": valid.reference_smiles,
        "reference_pIC50": valid.reference_pIC50,
        "reference_ic50_um": pIC50_to_IC50_uM(valid.reference_pIC50),
        "molecules": {"candidate": candidate_view, "reference": reference_view},
        "cliffforge": cf_out,
        "chemprop": cp_out,
        "example_pair_id": EXAMPLE["source_pair_id"] if is_example else None,
        "mmp_status": "supplied Dablander MMP" if is_example else "unverified for this manual pair",
        "timings_ms": {
            "validation": round(validation_ms, 2),
            "drawing": round(drawing_ms, 2),
            "cliffforge": round(cf_ms, 2),
            "chemprop": round(cp_ms, 2),
            "total": round((time.perf_counter() - start) * 1000, 2),
        },
        "trace": [
            "Validate SMILES and measured reference activity",
            "Draw 2D molecules and build molecular graphs",
            "CliffForge evaluates both pair orientations for each seed",
            "ChemProp predicts candidate pIC50 for each seed",
            "Subtract measured reference activity and compare results",
        ],
    }
