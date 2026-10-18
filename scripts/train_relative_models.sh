#!/usr/bin/env bash
# Train the 3-seed relative-affinity ensemble on the TRAIN partition -> models/rel_siamese_s{0,1,2}.pt (CPU ok, ~10 min/seed).
set -euo pipefail
source "$(dirname "$0")/_env.sh"
need_file "$DATA_DIR/moonshot/all.csv" "run scripts/prepare_data.sh first"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES-}"
for seed in 0 1 2; do
  echo "== training relative model, seed $seed"
  "$PYTHON" -m cliffforge.relmodel train --seed "$seed" "$@"
done
