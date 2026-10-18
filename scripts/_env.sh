# Shared environment for the helper scripts (sourced, not executed).
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PYTHON:-python}"
# data / weights / outputs locations (override with CLIFFFORGE_DATA / CLIFFFORGE_MODELS / CLIFFFORGE_OUTPUT)
DATA_DIR="${CLIFFFORGE_DATA:-$ROOT/data}"
MODEL_DIR="${CLIFFFORGE_MODELS:-$ROOT/models}"
OUTPUT_DIR="${CLIFFFORGE_OUTPUT:-$ROOT/outputs}"
export CLIFFFORGE_DATA="$DATA_DIR" CLIFFFORGE_MODELS="$MODEL_DIR" CLIFFFORGE_OUTPUT="$OUTPUT_DIR"
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONDONTWRITEBYTECODE=1
cd "$ROOT"

need_file() {   # need_file <path> <hint>
  if [ ! -e "$1" ]; then echo "Missing prerequisite: $1" >&2; echo "  -> $2" >&2; exit 2; fi
}
