#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "${SCRIPT_DIR}/run_revision2_singer_rooted_bootstrap120.py" \
  -B 120 \
  -n 160 \
  -j 12 \
  --sfratios-jobs 1 \
  --ls-sims 1000
