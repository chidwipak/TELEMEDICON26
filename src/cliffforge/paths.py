"""Repository-relative locations for data, model weights and outputs.

Nothing here assumes a particular machine. Defaults are relative to the repository root and can be
overridden with the environment variables CLIFFFORGE_DATA / CLIFFFORGE_MODELS / CLIFFFORGE_OUTPUT or with the
--data-dir / --model-dir / --output-dir flags that every command-line entry point accepts.
"""
import os
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_FILENAME = "cdd_noncovalent_dates_2023_10_18_filt.csv"
_CFG = {
    "data": Path(os.environ.get("CLIFFFORGE_DATA", REPO_ROOT / "data")),
    "models": Path(os.environ.get("CLIFFFORGE_MODELS", REPO_ROOT / "models")),
    "output": Path(os.environ.get("CLIFFFORGE_OUTPUT", REPO_ROOT / "outputs")),
}


def configure(data_dir=None, model_dir=None, output_dir=None):
    for key, val in (("data", data_dir), ("models", model_dir), ("output", output_dir)):
        if val:
            _CFG[key] = Path(val).expanduser().resolve()


def add_path_args(ap):
    ap.add_argument("--data-dir", help="directory holding raw/ and moonshot/ (default: <repo>/data)")
    ap.add_argument("--model-dir", help="directory holding rel_siamese_s*.pt and chemprop_moonshot/ (default: <repo>/models)")
    ap.add_argument("--output-dir", help="where to write new outputs (default: <repo>/outputs)")
    return ap


def apply_path_args(args):
    configure(getattr(args, "data_dir", None), getattr(args, "model_dir", None), getattr(args, "output_dir", None))


def data_dir():
    return _CFG["data"]


def model_dir():
    return _CFG["models"]


def output_dir():
    return _CFG["output"]


def raw_csv():
    return data_dir() / "raw" / RAW_FILENAME


def processed_dir():
    return data_dir() / "moonshot"


def processed_csv():
    return processed_dir() / "all.csv"


def chemprop_dir():
    return model_dir() / "chemprop_moonshot"


def reference_thresholds():
    return REPO_ROOT / "results" / "metrics" / "validation_thresholds.json"


def bin_path(name):
    """Console script installed next to the running interpreter (e.g. chemprop_predict), else on PATH."""
    cand = Path(sys.executable).parent / name
    if cand.exists():
        return str(cand)
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(f"cannot find '{name}'; install ChemProp (pip install chemprop==1.6.1)")
