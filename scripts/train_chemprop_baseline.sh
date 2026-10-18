#!/usr/bin/env bash
# Train the ChemProp absolute-pIC50 baseline on the TRAIN partition -> models/chemprop_moonshot/
# (same settings as the reference baseline; set CHEMPROP_DEVICE_ARGS="--gpu 0" to use a GPU).
set -euo pipefail
source "$(dirname "$0")/_env.sh"
need_file "$DATA_DIR/moonshot/train.csv" "run scripts/prepare_data.sh first"
CHEMPROP_TRAIN="$(dirname "$(command -v "$PYTHON")")/chemprop_train"
[ -x "$CHEMPROP_TRAIN" ] || CHEMPROP_TRAIN="$(command -v chemprop_train)"
# shellcheck disable=SC2086
"$CHEMPROP_TRAIN" --data_path "$DATA_DIR/moonshot/train.csv" --dataset_type regression \
  --smiles_columns smiles --target_columns pIC50 --save_dir "$MODEL_DIR/chemprop_moonshot" \
  --epochs 60 --num_folds 1 --metric rmse --extra_metrics r2 --num_workers 0 --quiet ${CHEMPROP_DEVICE_ARGS:---no_cuda}
