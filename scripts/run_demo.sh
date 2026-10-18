#!/usr/bin/env bash
# End-to-end demo: seed ligand -> constrained analogue generation -> filters -> relative-affinity scoring -> ranked shortlist.
# Default seed: outcome-blind rule picking a held-out potent Mpro compound (reference demo, row 237).
#   bash scripts/run_demo.sh                       # reference demo   -> outputs/demo/
#   bash scripts/run_demo.sh --seed-smiles "<SMILES>" --seed-id mylig   # your own seed (replaces --auto-seed)
#   bash scripts/run_demo.sh --no-score            # generation + filtering only (no weights needed)
set -euo pipefail
source "$(dirname "$0")/_env.sh"
need_file "$DATA_DIR/moonshot/all.csv" "run scripts/prepare_data.sh (see data/README.md)"
case " $* " in *" --no-score "*) ;; *) need_file "$MODEL_DIR/rel_siamese_s0.pt" "see models/README.md (scripts/train_relative_models.sh)";; esac
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES-}"
ARGS=("$@")
case " $* " in *" --seed-smiles "*|*" --seed-index "*|*" --auto-seed "*) ;; *) ARGS=(--auto-seed "${ARGS[@]}");; esac
"$PYTHON" -m cliffforge.demo "${ARGS[@]}"
echo "Done. Output files are listed on the 'wrote:' line above (default: $OUTPUT_DIR/demo/)."
