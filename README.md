# Identification of Horizontal Gene Transfer (HGT) in Parasitic Genomes

This repository contains the code, scripts, and analytical pipelines used for the comprehensive identification and refinement of Horizontal Gene Transfer (HGT) events across 105 pathogenic parasitic species. 

Our approach features a **Dual-Pipeline Detection Workflow** (Nucleotide Slicing & Protein Alignment) combined with a rigorous quality control framework, including Phylogenetic Entropy Filtering, Endosymbiotic Gene Transfer (EGT) exclusion, and a novel **Contig Voting Decontamination Strategy**.

---

## 📊 Project Overview

Parasitic genomes are often characterized by extensive gene loss, rapid divergence, and high risks of assembly contamination (e.g., from host or environmental microbiomes). To address this, we developed a highly robust pipeline to identify authentic HGT events:

1. **Phylogenetic Grouping:** Taxa are dynamically assigned into `Self-group` (same phylum), `Out-group` (potential donors: bacteria, archaea, plants, fungi, viruses, host class), and `Close-group` (other eukaryotes).
2. **Dual-Pipeline Alignment:** - *Genomic Slicing:* 1,500 bp sliding windows to capture HGTs across contig boundaries.
   - *Protein-Level:* Proteome-wide alignment against NCBI `core_nr`.
3. **Strict HGT Metrics:** Alien Index (AI) > 0, Outgroup Percentage > 90%, and Dominant Outgroup Count.
4. **Advanced Refinement:** - **Phylogenetic Entropy Filtering:** Distinguishes recent HGTs from ancient, widely conserved eukaryotic genes.
   - **EGT Exclusion:** Removes intracellular transfers from organellar progenitors (e.g., Rhodophyta, Cyanobacteria).
   - **Contig Voting Decontamination:** Assesses the integrity of HGT candidates by calculating the proportion of non-eukaryotic slices within the residing contig (>50% non-eukaryotic slices flags the contig as contaminated).

---

## 📁 Repository Structure

```text
HGT_parasite_Project/
├── data/                           # Input datasets and configurations
│   ├── 105parasite_info.txt        # Metadata for the 105 parasite species
│   ├── segments/                   # 1,500bp genomic slices generated for nucleotide pipeline
│   └── GCA_000227135.2_...fna      # Raw genome fasta files
├── requirements/                   # Environment dependencies (Conda/Pip)
├── results/                        # Pipeline outputs
│   ├── nucl_align_r1/              # Round 1: Nucleotide slice alignments
│   ├── nucl_align_r2/              # Round 2: Merged contiguous HGT slices alignment
│   ├── protein_align/              # Protein-level HGT identification & clustering
│   └── rmEGT_rmContaminatedcontig/ # Final QC results (EGT & contamination removed)
└── script/                         # Core pipeline scripts
    ├── 01_align/                   # Alignment execution scripts
    │   ├── 01_nucl_align.sh        # Nucleotide genomic slicing pipeline
    │   └── 01_prot_align.sh        # Protein alignment pipeline
    ├── 02_rmEGTand                 # EGT exclusion and Contig Voting script
    ├── perl_script/                # Utilities for sequence manipulation (e.g., sliding window)
    └── python_script/              # Custom Python scripts for AI calculation, Entropy, and Voting
```

---

## 🛠️ Dependencies

Ensure you have the following bioinformatics tools installed and accessible in your `$PATH`:

* **[DIAMOND](https://github.com/bbuchfink/diamond) (v2.1.7+):** For ultra-fast protein and translated DNA alignments (`blastp` and `blastx`).
* **[TaxonKit](https://bioinf.shenwei.me/taxonkit/):** For efficient NCBI taxonomy lineage resolution.
* **[CD-HIT](https://github.com/weizhongli/cdhit) (v4.8.1):** For clustering and redundancy reduction.
* **[BUSCO](https://busco.ezlab.org/) (v6.0.0+):** For genome completeness assessment.
* **Python 3.8+:** Required libraries include `pandas`, `numpy`, and `scipy` (for Shannon entropy calculations).

---

## 🚀 Usage / Pipeline Workflow

### Step 1: Pre-processing and Alignment
Run the dual-pipeline alignment scripts to perform DIAMOND searches and preliminary taxonomic assignments.
```bash
# 1.1 Run nucleotide slicing and alignment
bash script/01_align/01_nucl_align.sh

# 1.2 Run protein-level alignment
bash script/01_align/01_prot_align.sh
```

### Step 2: Preliminary HGT Verification
*Implemented within the step 1 scripts.* The pipeline calculates the Alien Index (AI), `Out_pct`, and Shannon entropy for each sequence, generating preliminary HGT candidate lists (`*HGT_verify.out`).

### Step 3: Quality Assurance (EGT & Contamination Removal)
Execute the refinement script to remove Endosymbiotic Gene Transfers (EGTs) and apply the **Contig Voting Strategy** to discard sequences residing on heavily contaminated contigs (>50% non-eukaryotic slices).
```bash
bash script/02_rmEGTand
```

The final high-confidence HGT candidates will be outputted to `results/rmEGT_rmContaminatedcontig/rm_wrong_contig/` as `final_nucl_10k.len_validated.txt` and `final_prot_10k.prot_validated.txt`.

---

## 📝 Citation


## 📄 License
This project is licensed under the MIT License - see the LICENSE file for details.
