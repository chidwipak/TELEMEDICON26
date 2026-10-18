#!/usr/bin/env bash
# Run the full submission validation (structure, secrets scan, results/README/slides consistency, unit tests).
set -euo pipefail
source "$(dirname "$0")/_env.sh"
"$PYTHON" scripts/validate_submission.py "$@"
