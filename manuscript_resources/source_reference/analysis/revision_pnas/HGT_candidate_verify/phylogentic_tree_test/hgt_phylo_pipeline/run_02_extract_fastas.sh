#!/usr/bin/env bash
set -euo pipefail
source ./config.sh
mkdir -p "$OUTDIR/seq_cache" "$OUTDIR/logs"
python scripts/02_extract_fastas.py   --targets "$OUTDIR/targets_for_tree.tsv"   --jobs-dir "$OUTDIR/jobs"   --protein-db "$PROTEIN_DIAMOND_DB"   --nt-db "$NT_DIAMOND_DB"   --protein-db-fasta "$PROTEIN_DB_FASTA"   --nt-db-fasta "$NT_DB_FASTA"   --protein-query-fasta "$PROTEIN_QUERY_FASTA"   --nt-query-fasta "$NT_QUERY_FASTA"   --genome-fasta-dir "$GENOME_FASTA_DIR"   --outdir "$OUTDIR"   --flank-aa "$FLANK_AA"   --decompress-threads "$SEQKIT_THREADS"   --query-coord-tolerance "${QUERY_COORD_TOLERANCE:-2}"   --extract-sequences
