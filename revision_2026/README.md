# Parasite HGT manuscript: 2026 analysis materials

This directory supplements the original `HGT_parasite_Project` code in the parent repository.
It contains a manuscript-to-analysis index, selected source data, provenance copies of analysis scripts,
and two tested downstream reproduction examples. Preparation date: 2026-10-22.

**Release status: preparation snapshot.** The two examples are independently runnable.
The `source_reference/` scripts document existing analyses and have not been packaged as a complete
raw-genome-to-final-catalogue workflow. External databases, complete phylogenies, alignments,
and some upstream input datasets are indexed but are not included here.

## Quick start

```bash
cd revision_2026
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-examples.txt
bash examples/run_examples.sh
```

Conda is optional; `environment-examples.yml` provides the tested Python/package versions.
Each example verifies its numerical output against independently saved manuscript-analysis results,
then writes editable-text SVG and PDF illustrations. The illustrations are example outputs,
not replacements for the manuscript's final manually formatted figures.

## Contents

| Location | Contents | Execution status |
|---|---|---|
| `examples/family_comparison/` | 48-species family/locus comparison | Tested against 9 saved contrasts |
| `examples/functional_contribution/` | Matched-family functional contributions and deletion summaries | Tested against 4 categories and 8 deletion summaries |
| `analyses/source_profiles/` | Existing R implementation of PERMANOVA/PERMDISP with input tables | Inputs included; not rerun in this preparation |
| `docs/final_analysis_manifest.tsv` | 23 analysis modules and their manuscript locations | Author-facing evidence index |
| `docs/artifact_manifest.tsv` | Project-relative scripts, inputs, outputs and original fingerprints | Indexed files verified to exist |
| `source_data/` | Selected manuscript source tables and literature-comparison summaries | Saved results, not a fresh full reanalysis |
| `source_reference/` | Selected earlier scripts in their original project-relative layout | Reference only; external paths require configuration |
| `docs/release_gaps.md` | Remaining upstream dependencies and manuscript consistency decisions | Resolve before claiming a complete reproduction release |

The quantity example starts from species counts. The functional example starts from already matched
family profiles. Neither example reruns sequence searches, assigns HGT origin, or infers new transfers.
The examples preserve the numerical definitions used for the archived results.

## Source-data versions and interpretation

The catalogue contains 8,989 entries. The Figure 2B measurement inputs located during this preparation
contain 8,985 entries (6,878 protein and 2,107 nucleotide records). These scopes must not be silently merged.
The main family-count panel has 48 species, 2,299 candidate loci and 1,973 species-by-candidate-family units.
The functional matched panel has 1,940 species-by-family observations. The conditional concentration
comparison has a more restricted scope, recorded in its source tables.

For literature comparisons, counts of matched literature proteins and counts of catalogue entries differ.
A sequence correspondence is not automatically a shared transfer event.
RNA-seq TPM and PRO-seq nascent-transcription signal are separate measurements.

## Provenance and access

The artifact index uses paths relative to the original project root to identify the source files.
Its `bundled_path` column specifies files actually included in this directory. An empty value means the
file remains external, even if it was found during local preparation.

Provenance scripts have local home-directory prefixes replaced by `PROJECT_ROOT`, `EXTERNAL_DATA_ROOT`,
or `HPC_DATA_ROOT`. These are literal documentary placeholders, not configured runtime variables.
Original fingerprints and published-file fingerprints are therefore recorded separately.
Do not execute reference scripts indiscriminately: some write outputs into their source project layout.
Use the tested examples for the first run.

This folder is intended to be added to the existing MIT-licensed repository without removing its older
files or replacing the repository license. Cite the manuscript and the eventual archived release DOI;
no DOI has been minted for this preparation snapshot.
