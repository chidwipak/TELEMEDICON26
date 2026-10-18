#!/usr/bin/env bash
# Activity-cliff evaluation: (1) ranking + descriptive cliff pairs, (2) direct cliff DETECTION, (3) generate->score recall.
#   bash scripts/run_evaluation.sh            # reference thresholds, 2000 bootstrap resamples (~10 min CPU)
#   bash scripts/run_evaluation.sh --fast     # 200 resamples (~3 min)
#   bash scripts/run_evaluation.sh --with-cv  # re-derive thresholds from fresh CV folds (run scripts/run_validation_cv.sh first)
# Outputs: outputs/cliff_pairs/, outputs/cliff_detection/, outputs/recall/
set -euo pipefail
source "$(dirname "$0")/_env.sh"
need_file "$DATA_DIR/moonshot/all.csv" "run scripts/prepare_data.sh"
need_file "$MODEL_DIR/rel_siamese_s0.pt" "see models/README.md"
need_file "$MODEL_DIR/chemprop_moonshot/fold_0/model_0/model.pt" "see models/README.md (scripts/train_chemprop_baseline.sh)"
export CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
NBOOT=2000; CVARGS=()
for a in "$@"; do
  case "$a" in
    --fast) NBOOT=200;;
    --with-cv) need_file "$OUTPUT_DIR/cv/fold0/chemprop_preds.csv" "run scripts/run_validation_cv.sh"; CVARGS=(--cv-dir "$OUTPUT_DIR/cv");;
  esac
done
echo "== 1/3 within-series ranking and descriptive cliff pairs";  "$PYTHON" -m cliffforge.cliffs
echo "== 2/3 direct activity-cliff detection (n_boot=$NBOOT)";    "$PYTHON" -m cliffforge.cliff_detection --n-boot "$NBOOT" "${CVARGS[@]}"
echo "== 3/3 generate->score retrospective check";                "$PYTHON" -m cliffforge.recall
echo "Compare with results/metrics/ (reference run). Done."
