#!/bin/bash
#===============================================================================
# Script: step1_prot_hgt.sh
# Description: Pipeline for protein-based HGT identification in parasites.
#              Includes distribution of alignments, taxonomy lineage extraction,
#              verification, entropy filtration, and sequence clustering.
#===============================================================================

set -e
set -u

# ==============================================================================
# 0. Configuration and Environment Setup
# ==============================================================================
BASE_DIR=$(cd "$(dirname "$0")/.."; pwd)
SRC_DIR="${BASE_DIR}/scripts"
DATA_DIR="${BASE_DIR}/data"

ALIGN_DIR="${DATA_DIR}/align/protein_align"
OUT_DIR="${BASE_DIR}/results/protein_align/genome_divided"

# Metadata files
NAME2GENOME="${DATA_DIR}/config/105_name2accfil_withTaxid.txt"
ACC2GENOME="${DATA_DIR}/config/all_accession_genome_prot.txt"
PROTEIN_DB="${DATA_DIR}/genomes/all_parasite_protein_sequences.fasta" # Combined proteome

# Scripts
DISTRIBUTE_PY="${SRC_DIR}/utils/distribute_blast.py"
HGT_VERIFY_PROT="${SRC_DIR}/verification/HGT_verify_parasite.prot.py"
HGT_FILTER_PROT="${SRC_DIR}/verification/HGT_filter_entropy.prot.py"

mkdir -p "${ALIGN_DIR}" "${OUT_DIR}"

echo "======================================================================"
echo "[INFO] Starting Protein-based HGT Pipeline..."
echo "======================================================================"

# ==============================================================================
# Step 1: Distribute Alignment Results
# ==============================================================================
echo "[INFO] Step 1: Distributing protein alignment results to individual genomes..."
cd "${ALIGN_DIR}/genome_divided"

MAIN_ALIGNMENT="../157_parasite.taxonomy.blast.iden30.out"

# Check if distribution has already been done to save time
if [ -f "$MAIN_ALIGNMENT" ] && [ ! -f "distribution.done" ]; then
    python "${DISTRIBUTE_PY}" "${ACC2GENOME}" "$MAIN_ALIGNMENT"
    touch distribution.done
    echo "[SUCCESS] Alignment distribution completed."
else
    echo "[SKIP] Distribution already completed or main alignment file not found."
fi

# ==============================================================================
# Step 2: Extract Target Taxonomy Lineages
# ==============================================================================
echo "[INFO] Step 2: Extracting target taxonomy lineages..."

if [ ! -f "target_nt_taxid_lineage.txt" ] && [ -f "$MAIN_ALIGNMENT" ]; then
    cut -f4 "$MAIN_ALIGNMENT" | sort | uniq | taxonkit lineage | taxonkit reformat > target_nt_taxid_lineage_tmp.txt
    
    # Merge with nucleotide lineage if available to ensure comprehensive taxonomy coverage
    NUCL_LINEAGE="${DATA_DIR}/align/nucl_align/target_nt_taxid_lineage.txt"
    if [ -f "$NUCL_LINEAGE" ]; then
        cat target_nt_taxid_lineage_tmp.txt "$NUCL_LINEAGE" | sort | uniq > target_nt_taxid_lineage.txt
    else
        mv target_nt_taxid_lineage_tmp.txt target_nt_taxid_lineage.txt
    fi
    echo "[SUCCESS] Taxonomy lineage file generated."
fi

# ==============================================================================
# Step 3: Verification, Entropy Filtration, and Clustering
# ==============================================================================
echo "[INFO] Step 3: Verifying candidates, filtering by entropy, and clustering..."

while read -r species genome taxid; do
    ACC2TAXID_FILE="${genome}.acc2taxid.out"

    if [ -f "$ACC2TAXID_FILE" ]; then
        echo ">> Processing: ${genome} (TaxID: ${taxid})"

        # 1. HGT Verification based on Taxonomy
        python "${HGT_VERIFY_PROT}" target_nt_taxid_lineage.txt "$ACC2TAXID_FILE" "${taxid}" \
            > "${OUT_DIR}/${genome}.HGT_verify_new.out"
        
        # 2. Shannon Entropy Filtration
        python "${HGT_FILTER_PROT}" "${OUT_DIR}/${genome}.HGT_verify_new.out" \
            -o "${OUT_DIR}/${genome}.HGT_verify_filtered.entropy.out"
        
        # 3. Extract accessions (remove header)
        sed '1d' "${OUT_DIR}/${genome}.HGT_verify_filtered.entropy.out" | awk '{print $1}' \
            > "${OUT_DIR}/${genome}.tmp.acc"
        
        # 4. Extract sequences and remove redundancy using CD-HIT
        if [ -s "${OUT_DIR}/${genome}.tmp.acc" ] && [ -f "$PROTEIN_DB" ]; then
            seqtk subseq "$PROTEIN_DB" "${OUT_DIR}/${genome}.tmp.acc" \
                > "${OUT_DIR}/${genome}.HGT_verify_filtered.entropy.fa"
            
            # CD-HIT clustering at 90% identity
            cd-hit -c 0.9 -i "${OUT_DIR}/${genome}.HGT_verify_filtered.entropy.fa" \
                -o "${OUT_DIR}/${genome}.HGT_verify_filtered.entropy.clustered.fa" -M 0 -d 0 > /dev/null
        fi
    fi
done < "${NAME2GENOME}"

# ==============================================================================
# Step 4: Final Aggregation
# ==============================================================================
echo "[INFO] Step 4: Aggregating all clustered protein HGT sequences..."

cat "${OUT_DIR}"/*.entropy.clustered.fa > "${OUT_DIR}/105_parasite_prot_HGT.entropy.fa"

echo "======================================================================"
echo "[SUCCESS] Protein Pipeline completed successfully!"
echo "Final clustered sequences saved to: ${OUT_DIR}/105_parasite_prot_HGT.entropy.fa"
echo "======================================================================"