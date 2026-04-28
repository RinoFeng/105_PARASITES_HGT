#!/bin/bash
#===============================================================================
# Script: step1_nucl_hgt.sh
# Description: Pipeline for nucleotide-based HGT identification in parasites.
#              Includes both the first-round screening and the mandatory 
#              second-round verification logic.
#===============================================================================

set -e
set -u

# ==============================================================================
# 0. Configuration and Environment Setup
# ==============================================================================
BASE_DIR=$(cd "$(dirname "$0")/.."; pwd)
SRC_DIR="${BASE_DIR}/scripts"
DATA_DIR="${BASE_DIR}/data"
GENOME_DIR="${DATA_DIR}/genomes"

# Directories for 1st and 2nd round alignments
ALIGN_R1_DIR="${DATA_DIR}/align/nucl_align_r1"
ALIGN_R2_DIR="${DATA_DIR}/align/nucl_align_r2"
OUT_DIR="${BASE_DIR}/results/nucl"

# Metadata files
NAME2GENOME="${DATA_DIR}/config/105_name2accfil_withTaxid.txt"
ACC2GENOME="${DATA_DIR}/config/all_accession_genome.txt"

# Scripts
NUCL2TAXID="${SRC_DIR}/verification/nucl2prot_acc2taxid"
DISTRIBUTE_PY="${SRC_DIR}/utils/distribute_blast.py"
HGT_VERIFY_PY="${SRC_DIR}/verification/HGT_verify_parasite.py"
HGT_FILTER_PY="${SRC_DIR}/verification/HGT_filter_entropy.prot.py"
EXTRACT_PY="${SRC_DIR}/utils/extract_HGT_nucl_position.py"
MERGE_PY="${SRC_DIR}/utils/merge_bed.py"

mkdir -p "${OUT_DIR}" "${ALIGN_R1_DIR}/genome_divided" "${ALIGN_R2_DIR}"

echo "[INFO] Starting Nucleotide-based HGT Pipeline..."

# ==============================================================================
# Step 1: First-Round Alignment Processing & Distribution
# ==============================================================================
echo "[INFO] Step 1: Distributing 1st-round BLAST/DIAMOND results..."
cd "${ALIGN_R1_DIR}"
for file in *.out; do
    if [ -f "$file" ]; then
        python "${DISTRIBUTE_PY}" "${ACC2GENOME}" "$file"
    fi
done

cd "${ALIGN_R1_DIR}/genome_divided"
echo "[INFO] Generating TaxID mappings and lineage info for 1st round..."
while read -r species genome taxid; do
    if [ -f "${genome}.out" ]; then
        ${NUCL2TAXID} 12305 "${genome}.out" > "${genome}.acc2taxid.out"
    fi
done < "${NAME2GENOME}"

cut -f4 *.acc2taxid.out | sort | uniq | taxonkit lineage | taxonkit reformat > target_nt_taxid_lineage.txt

# ==============================================================================
# Step 2: First-Round HGT Verification and merge the potential HGT region.
# ==============================================================================
# cause of the limitation of file size in github, only showing the top 10,000 alignments results.
echo "[INFO] Step 2: First-round HGT verification and sequence extraction..."
while read -r species genome taxid; do
    if [ ! -f "${ALIGN_R1_DIR}/genome_divided/${genome}.acc2taxid.out" ]; then
        continue
    fi

    # Verification and Filtration
    python "${HGT_VERIFY_PY}" target_nt_taxid_lineage.txt \
        "${ALIGN_R1_DIR}/genome_divided/${genome}.acc2taxid.out" "${taxid}" > "${OUT_DIR}/${genome}.HGT_verify_r1.out"
        
    python "${HGT_FILTER_PY}" "${OUT_DIR}/${genome}.HGT_verify_r1.out" \
        -o "${OUT_DIR}/${genome}.HGT_verify_filtered_r1.out"
    
    # Extract genomic coordinates and sequences for the 2nd round alignment
    sed -i '1d' "${OUT_DIR}/${genome}.HGT_verify_filtered_r1.out"
    python "${EXTRACT_PY}" "${ALIGN_R1_DIR}/genome_divided/${genome}.out" \
        "${OUT_DIR}/${genome}.HGT_verify_filtered_r1.out" > "${OUT_DIR}/${genome}.r1.bed"
        
    python "${MERGE_PY}" "${OUT_DIR}/${genome}.r1.bed" "${OUT_DIR}/${genome}.r1_merged.bed"
    
    seqtk subseq "${GENOME_DIR}/${genome}" "${OUT_DIR}/${genome}.r1_merged.bed" > "${OUT_DIR}/${genome}.r1_merged.fa"
done < "${NAME2GENOME}"

# ==============================================================================
# Step 3: Second-Round Alignment (Placeholder for HPC execution)
# ==============================================================================
echo "[INFO] Step 3: Second-round DIAMOND alignment..."
# Note: In actual production, this step is performed on an HPC cluster against the core_nt database.
# Here we assume the results are already placed in ${ALIGN_R2_DIR}
# Command example: diamond blastx -d core_nt -q all_r1_merged.fa -o 105_parasite.second_align.out ...

# ==============================================================================
# Step 4: Second-Round HGT Verification (Strict Filtration with entropy)
# ==============================================================================
echo "[INFO] Step 4: Processing secondary alignment results and finalizing candidates..."
cd "${ALIGN_R2_DIR}"

# Assume acc2taxid and target_nt_taxid_lineage.txt for R2 are generated similarly to Step 1
while read -r species genome taxid; do
    if [ -f "${genome}.acc2taxid.out" ] && [ -f "target_nt_taxid_lineage.txt" ]; then
        
        # Second round verification
        python "${HGT_VERIFY_PY}" target_nt_taxid_lineage.txt \
            "${genome}.acc2taxid.out" "${taxid}" > "${OUT_DIR}/${genome}.HGT_verify_r2.out"
            
        python "${HGT_FILTER_PY}" "${OUT_DIR}/${genome}.HGT_verify_r2.out" \
            -o "${OUT_DIR}/${genome}.HGT_verify_filtered_r2.out"
        
        # Extract final high-confidence HGT sequences
        sed '1d' "${OUT_DIR}/${genome}.HGT_verify_filtered_r2.out" | cut -f1 | \
            awk -F ":" '{print $1"\t"$2}' | awk -F "-" '{print $1"\t"$2}' > "${OUT_DIR}/${genome}.final.bed"
            
        seqtk subseq "${GENOME_DIR}/${genome}" "${OUT_DIR}/${genome}.final.bed" > "${OUT_DIR}/${genome}.final_HGT.fa"
    fi
done < "${NAME2GENOME}"

cat "${OUT_DIR}"/*.final_HGT.fa > "${OUT_DIR}/all_genome_HGT.nucl.fa"
echo "[SUCCESS] Nucleotide Pipeline completed!"