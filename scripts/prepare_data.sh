#!/usr/bin/env bash
# Build data/moonshot/{all,train,test}.csv from the raw COVID Moonshot file (see data/README.md).
set -euo pipefail
source "$(dirname "$0")/_env.sh"
need_file "$DATA_DIR/raw/cdd_noncovalent_dates_2023_10_18_filt.csv" "download it as described in data/README.md"
"$PYTHON" -m cliffforge.data_prep "$@"
