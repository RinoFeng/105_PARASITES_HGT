#!/usr/bin/env bash
set -euo pipefail

# ===== Input curated HGT table =====
CURATED_TSV="PROJECT_ROOT/analysis/revision_pnas/HGT_candidate_verify/phylogentic_tree_test/integrated_HGTdata_curated_with_contig_and_tree_priority.tsv"

# ===== Recipient genome -> taxid map =====
# Format: species_name<TAB>genome_fasta_file<TAB>taxid
# This is required so 01_sample_homologs.py can mimic HGT_verify_parasite.py usage:
#   python HGT_verify_parasite.py taxonomy.txt genome.acc2taxid.out "$taxid"
GENOME_TAXID_MAP="PROJECT_ROOT/105_name2accfil_withTaxid.txt"

# ===== Protein search results =====
PROTEIN_ALIGN_DIR="PROJECT_ROOT/align/protein_align/genome_divide"
PROTEIN_TAXONOMY="PROJECT_ROOT/align/protein_align/genome_divide/target_nt_taxid_lineage.txt"
PROTEIN_DIAMOND_DB="EXTERNAL_DATA_ROOT/DB/nr/diamond_nr/nr_rep_seq.dmnd"

# ===== Nucleotide/translated search results =====
NT_ALIGN_DIR="PROJECT_ROOT/align/second_align"
NT_TAXONOMY="PROJECT_ROOT/align/target_nt_taxid_lineage2.txt"
NT_DIAMOND_DB="EXTERNAL_DATA_ROOT/DB/core_nt_trans2prot_20241217/core_nt.prot.diamond_db.dmnd"

# ===== Query sequence source =====
# For protein HGT query sequences: usually extracted from nr db by accession.
# If some query accessions are absent from nr, provide a merged protein FASTA here; otherwise leave empty.
PROTEIN_QUERY_FASTA="PROJECT_ROOT/analysis/Results_article_new/Remove_EGT/rm_wrong_contig/final_prot_10k.prot_validated.fasta"

# For nucleotide HGT query regions, choose ONE option:
# 1) A FASTA containing HGT nucleotide regions with IDs exactly like Accession, e.g. CM000442.1:40422-40783
NT_QUERY_FASTA="PROJECT_ROOT/analysis/Results_article_new/Remove_EGT/rm_wrong_contig/final_nucl_10k.len_validated.fasta"
# 2) A directory containing genome FASTA files named exactly as Genome column, e.g. GCA_000002415.2_ASM241v2_genomic.fna
GENOME_FASTA_DIR="EXTERNAL_DATA_ROOT/Genome/parasite/188_parasites/Genome_files"

# ===== Output =====
OUTDIR="./hgt_phylo_out"
THREADS=4
PARALLEL_JOBS=8

# ===== Sampling parameters =====
MAX_SELF=12
MAX_CLOSE=30
MAX_OUT=70
MAX_REMOTE_EUK=10
MAX_TOTAL=130
MIN_ALIGN_LEN=80
MIN_PIDENT=20
# Diagnostic only; keep 0 for production. If >0, top N hits are sampled when all group classification fails.
FALLBACK_TOP_HITS=0
# For nucleotide HGT: if Accession does not exactly match .fna.out qseqid, match by same contig and interval overlap.
QUERY_OVERLAP=0.50
FLANK_AA=20
# For nucleotide HGT query ID matching: tolerate small 0/1-based coordinate shifts.
# Example: CM000444.2:175849-176494 can match CM000444.2:175848-176494 when this is >=1.
QUERY_COORD_TOLERANCE=2

# ===== Original FASTA files used to build DIAMOND databases =====
# Required by run_02_extract_fastas.sh. Many DIAMOND versions cannot extract sequences
# by accession list from .dmnd via getseq, so selected accessions are extracted from
# the original FASTA using seqkit grep or a Python streaming fallback.
# Set these to the FASTA files used when creating PROTEIN_DIAMOND_DB and NT_DIAMOND_DB.
PROTEIN_DB_FASTA="EXTERNAL_DATA_ROOT/DB/nr/nr_rep_seq.fasta"
NT_DB_FASTA="EXTERNAL_DATA_ROOT/DB/core_nt_trans2prot_20241217/core_nt.prot.faa.gz"
# auto: use seqkit if available; otherwise Python streaming. Options: auto, seqkit, python
EXTRACTOR="seqkit"
# Number of threads passed to seqkit grep. For very large gz FASTA, 16-32 is usually safer than 64 because I/O and decompression are often the bottleneck.
SEQKIT_THREADS=32
# 0: match FASTA first token exactly; 1: use seqkit grep -n to match the whole header/name when FASTA headers wrap accessions such as ref|XP_...|.
SEQKIT_MATCH_NAME=0
