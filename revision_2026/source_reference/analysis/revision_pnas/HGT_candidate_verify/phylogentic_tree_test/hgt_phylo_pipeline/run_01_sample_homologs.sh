#!/usr/bin/env bash
set -euo pipefail
source ./config.sh
mkdir -p "$OUTDIR/jobs" "$OUTDIR/logs"
python scripts/01_sample_homologs.py \
  --targets "$OUTDIR/targets_for_tree.tsv" \
  --genome-taxid-map "$GENOME_TAXID_MAP" \
  --protein-align-dir "$PROTEIN_ALIGN_DIR" \
  --protein-taxonomy "$PROTEIN_TAXONOMY" \
  --nt-align-dir "$NT_ALIGN_DIR" \
  --nt-taxonomy "$NT_TAXONOMY" \
  --outdir "$OUTDIR" \
  --max-self "$MAX_SELF" \
  --max-close "$MAX_CLOSE" \
  --max-out "$MAX_OUT" \
  --max-remote-euk "$MAX_REMOTE_EUK" \
  --max-total "$MAX_TOTAL" \
  --min-align-len "$MIN_ALIGN_LEN" \
  --min-pident "$MIN_PIDENT" \
  --fallback-top-hits "$FALLBACK_TOP_HITS" \
  --query-overlap "$QUERY_OVERLAP" \
  --query-coord-tolerance "$QUERY_COORD_TOLERANCE"
