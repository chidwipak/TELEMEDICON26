"""
checkpoint_verify.py
====================
Verify file sizes and SHA-256 hashes of the six fold-3 checkpoints
against sha256_manifest.csv. Raises RuntimeError on any mismatch.

Called at startup (from app.py) so /api/health can report details.
"""
from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path
from typing import Dict

# The run directory may be overridden when the demo is launched from another
# checkout. By default the complete metadata and the six local, ignored weights
# live together in this repository.
_REPO_ROOT = Path(__file__).resolve().parents[1]
_RUN_CANDIDATES = [
    _REPO_ROOT / "cliffforge_dablander_5fold",
    _REPO_ROOT / "results" / "dablander_mmp_5fold",
]
_DEFAULT_RUN = next(
    (path for path in _RUN_CANDIDATES if (path / "sha256_manifest.csv").is_file()),
    _RUN_CANDIDATES[0],
)
RESULTS_ROOT = Path(os.environ.get("CLIFFFORGE_RESULTS_DIR", str(_DEFAULT_RUN))).resolve()
MANIFEST_CSV = RESULTS_ROOT / "sha256_manifest.csv"

# The six checkpoints we need
REQUIRED_RELATIVE_PATHS = [
    "fold_3/cliffforge/seed_0/final.pt",
    "fold_3/cliffforge/seed_1/final.pt",
    "fold_3/cliffforge/seed_2/final.pt",
    "fold_3/chemprop/seed_0/fold_0/model_0/model.pt",
    "fold_3/chemprop/seed_1/fold_0/model_0/model.pt",
    "fold_3/chemprop/seed_2/fold_0/model_0/model.pt",
]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_manifest() -> Dict[str, dict]:
    """Return dict: relative_path -> {bytes: int, sha256: str}."""
    manifest: Dict[str, dict] = {}
    with open(MANIFEST_CSV, newline="") as f:
        for row in csv.DictReader(f):
            manifest[row["relative_path"]] = {
                "bytes": int(row["bytes"]),
                "sha256": row["sha256"].strip().lower(),
            }
    return manifest


def verify_checkpoints() -> list[dict]:
    """
    Verify each required checkpoint. Returns a list of dicts, one per file,
    with keys: relative_path, absolute_path, expected_bytes, actual_bytes,
    expected_sha256, actual_sha256, ok.

    Raises RuntimeError if any file is absent or mismatched.
    """
    manifest = load_manifest()
    results = []
    errors = []

    for rel in REQUIRED_RELATIVE_PATHS:
        abs_path = RESULTS_ROOT / rel
        expected = manifest.get(rel)
        if expected is None:
            errors.append(f"MANIFEST MISSING entry for {rel}")
            continue

        if not abs_path.exists():
            errors.append(f"FILE NOT FOUND: {abs_path}")
            results.append({
                "relative_path": rel,
                "absolute_path": str(abs_path),
                "expected_bytes": expected["bytes"],
                "actual_bytes": None,
                "expected_sha256": expected["sha256"],
                "actual_sha256": None,
                "ok": False,
                "error": "file not found",
            })
            continue

        actual_bytes = os.path.getsize(abs_path)
        actual_sha = _sha256(abs_path)
        ok = (actual_bytes == expected["bytes"] and
              actual_sha == expected["sha256"])

        entry = {
            "relative_path": rel,
            "absolute_path": str(abs_path),
            "expected_bytes": expected["bytes"],
            "actual_bytes": actual_bytes,
            "expected_sha256": expected["sha256"],
            "actual_sha256": actual_sha,
            "ok": ok,
        }
        if not ok:
            entry["error"] = (
                f"size {actual_bytes}!={expected['bytes']}" if actual_bytes != expected["bytes"]
                else f"sha256 mismatch"
            )
            errors.append(f"MISMATCH {rel}: {entry['error']}")
        results.append(entry)

    if errors:
        raise RuntimeError("Checkpoint verification failed:\n" + "\n".join(errors))

    return results


def checkpoint_absolute_paths() -> dict:
    """Return dict of logical name -> absolute path str."""
    root = RESULTS_ROOT
    return {
        "cf_seed_0": str(root / "fold_3/cliffforge/seed_0/final.pt"),
        "cf_seed_1": str(root / "fold_3/cliffforge/seed_1/final.pt"),
        "cf_seed_2": str(root / "fold_3/cliffforge/seed_2/final.pt"),
        "cp_seed_0": str(root / "fold_3/chemprop/seed_0/fold_0/model_0/model.pt"),
        "cp_seed_1": str(root / "fold_3/chemprop/seed_1/fold_0/model_0/model.pt"),
        "cp_seed_2": str(root / "fold_3/chemprop/seed_2/fold_0/model_0/model.pt"),
    }


if __name__ == "__main__":
    import json
    results = verify_checkpoints()
    print(json.dumps(results, indent=2))
    print("All 6 checkpoints verified OK.")
