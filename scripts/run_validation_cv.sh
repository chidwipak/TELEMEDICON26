#!/usr/bin/env bash
# OPTIONAL (~20-45 min on CPU): 5-fold CV inside the train partition to choose cliff-detection thresholds without touching
# test labels. 5 folds x (3 relative seeds + 1 ChemProp). Output: outputs/cv/fold{0..4}/.
# PARALLEL=1 runs all folds at once (needs ~20 CPU threads); default is fold by fold.
set -euo pipefail
source "$(dirname "$0")/_env.sh"
need_file "$DATA_DIR/moonshot/all.csv" "run scripts/prepare_data.sh first"
export CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
run_fold() {
  local f="$1"
  for s in 0 1 2; do "$PYTHON" -m cliffforge.cv rel --fold "$f" --seed "$s" > "$OUTPUT_DIR/cv/logs/rel_f${f}_s${s}.log" 2>&1 & done
  "$PYTHON" -m cliffforge.cv chemprop --fold "$f" > "$OUTPUT_DIR/cv/logs/chemprop_f${f}.log" 2>&1 &
  wait
}
mkdir -p "$OUTPUT_DIR/cv/logs"
if [ "${PARALLEL:-0}" = "1" ]; then
  for f in 0 1 2 3 4; do run_fold "$f" & done; wait
else
  for f in 0 1 2 3 4; do echo "== fold $f"; run_fold "$f"; done
fi
echo "done: $OUTPUT_DIR/cv/"
