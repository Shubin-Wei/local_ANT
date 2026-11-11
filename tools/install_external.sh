#!/usr/bin/env bash
# Source-checkout convenience entry point. Installed packages use lant install-external.
set -euo pipefail
REPO_DIR=$(cd "$(dirname "$0")/.." && pwd)
export LANT_EXTERNAL=${LANT_EXTERNAL:-$REPO_DIR/external}
PYTHONPATH="$REPO_DIR${PYTHONPATH:+:$PYTHONPATH}" exec python -m lant install-external "$@"
