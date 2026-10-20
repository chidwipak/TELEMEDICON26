#!/usr/bin/env bash
# Start the verified fold-3 evaluator demo on this computer.
set -euo pipefail

DEMO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$DEMO_DIR")"
PORT="${CLIFFFORGE_DEMO_PORT:-8765}"
PYTHON="${CLIFFFORGE_PYTHON:-python3}"
DEFAULT_RUN="$REPO_ROOT/cliffforge_dablander_5fold"
ALTERNATE_RUN="$REPO_ROOT/results/dablander_mmp_5fold"
ORIGINAL_RUN="/Volumes/Work/TELEMEDICON26/results/dablander_mmp_5fold"

if [[ -n "${CLIFFFORGE_RESULTS_DIR:-}" ]]; then
    RUN_DIR="$CLIFFFORGE_RESULTS_DIR"
elif [[ -f "$DEFAULT_RUN/sha256_manifest.csv" ]]; then
    RUN_DIR="$DEFAULT_RUN"
elif [[ -f "$ALTERNATE_RUN/sha256_manifest.csv" ]]; then
    RUN_DIR="$ALTERNATE_RUN"
else
    echo "Dablander five-fold results were not found. Set CLIFFFORGE_RESULTS_DIR." >&2
    exit 1
fi

# Managed Git worktrees omit ignored .pt files. Restore only missing weights
# from the complete local archive, then verify every checksum below.
for seed in 0 1 2; do
    for rel in "fold_3/cliffforge/seed_$seed/final.pt" "fold_3/chemprop/seed_$seed/fold_0/model_0/model.pt"; do
        if [[ ! -f "$RUN_DIR/$rel" && -f "$ORIGINAL_RUN/$rel" ]]; then
            mkdir -p "$(dirname "$RUN_DIR/$rel")"
            cp "$ORIGINAL_RUN/$rel" "$RUN_DIR/$rel"
        fi
    done
done
RUN_DIR="$(cd "$RUN_DIR" && pwd)"

export CLIFFFORGE_RESULTS_DIR="$RUN_DIR"
export PYTHONPATH="$REPO_ROOT/src:$DEMO_DIR${PYTHONPATH:+:$PYTHONPATH}"
export REPO_ROOT

"$PYTHON" - <<'PY'
import hashlib, importlib.metadata, json, os, pathlib
import numpy, torch, torch_geometric, rdkit, fastapi, uvicorn
root=pathlib.Path(os.environ["REPO_ROOT"])
run=pathlib.Path(os.environ["CLIFFFORGE_RESULTS_DIR"])
protocol=json.loads((run/"protocol.json").read_text())
for relative, expected in protocol["embedded_model_code_sha256"].items():
    actual=hashlib.sha256((root/relative).read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit(f"CliffForge source hash differs from the trained model: {relative}")
if importlib.metadata.version("chemprop") != "1.6.1":
    raise SystemExit("ChemProp 1.6.1 is required")
if int(numpy.__version__.split(".")[0]) != 1:
    raise SystemExit("This ChemProp 1.6.1 setup requires NumPy 1.x")
print(f"Runtime: Python {os.sys.version.split()[0]}, Torch {torch.__version__}, PyG {torch_geometric.__version__}, RDKit {rdkit.__version__}, ChemProp 1.6.1")
PY

"$PYTHON" "$DEMO_DIR/checkpoint_verify.py" > /dev/null
echo "Six fold-3 checkpoints and model code hashes verified."
if [[ "${1:-}" == "--verify-only" ]]; then
    exit 0
fi
echo "Serving local demo at http://127.0.0.1:$PORT"
cd "$DEMO_DIR"
exec "$PYTHON" -m uvicorn app:app --host 127.0.0.1 --port "$PORT" --log-level warning --no-access-log
