#!/usr/bin/env bash
set -euo pipefail
source ./config.sh
mkdir -p "$OUTDIR"
python scripts/00_prepare_targets.py \
  --curated-tsv "$CURATED_TSV" \
  --outdir "$OUTDIR" \
  --keep-main Curated_main_set Tier3_rescued_by_extra_support
