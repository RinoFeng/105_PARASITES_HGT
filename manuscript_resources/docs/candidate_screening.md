# Candidate screening: contigs and gene trees

This guide connects the initial search results to the curated HGT catalogue. It documents the source implementations included in this repository; it does not launch a new screen. Historical catalogue curation and the targeted viral reassessment have separate provenance. Parameters below are attributed to their implementations rather than assumed to have been applied in one rerun to all 8,989 entries.

## Contig screening

**Inputs:** candidate search summaries, genome sequences, candidate coordinates, window-level taxonomic evidence and assembly N50. Genomic placement is checked before contig context is interpreted.

1. **Initial contig filtering.** The saved [method table](../source_data/tables/S1/Method_parameters.tsv) records the original `initial_contig` rules, including the reported-window fraction and eukaryotic-anchor requirement. This stage is more specific than rejecting every contig with a majority of distant hits.
2. **Window composition in the viral reassessment.** [build_contig_evidence.py](../source_reference/analysis/revision_pnas/virus_supplement/scripts/build_contig_evidence.py) records eukaryotic, non-eukaryotic, conflicting and unassigned windows. The conservative fraction is `(non_eukaryotic + conflicting) / (non_eukaryotic + eukaryotic + conflicting)`; unassigned windows are omitted. A fraction above 0.5 is flagged. Missing composition evidence is not automatically counted as a passing result.
3. **Contig-context tiers.** The `tier()` function in [map_protein_candidates.py](../source_reference/analysis/revision_pnas/virus_supplement/scripts/map_protein_candidates.py) implements the table below. [harmonize_contig_support.py](../source_reference/analysis/revision_pnas/virus_supplement/scripts/harmonize_contig_support.py) evaluates the supporting features after both discovery routes are available.

| Check, in evaluation order | Implemented condition |
|---|---|
| Candidate-region fraction | Exclude if the candidate region occupies more than 70% of its contig |
| Tier 1 | Contig length ≥100 kb and distance to its nearest edge ≥10 kb |
| Tier 2 | Contig length ≥max(20 kb, min(50 kb, assembly N50)) and edge distance ≥5 kb |
| Tier 3 | Contig length ≥20 kb, edge distance ≥2 kb, and at least three of the five supporting features below |
| Otherwise | Exclude from this contig-context selection |

The five Tier 3 features are: support from both discovery routes; at least three close-group hits; a bacterial or archaeal best distant hit; candidate length ≥150; and entropy ≤2.5. Length is measured in amino acids for the protein route and nucleotides for the nucleotide route. Dual-route support uses at least 50% reciprocal genomic overlap on the same strand. Eligibility also requires the route-specific mapping/sequence checks and acceptable composition evidence; a tier alone is not final HGT acceptance. Bilateral flank evidence is recorded separately and checked in the final supplementary-candidate assessment.

**Outputs:** window votes, contig composition summaries, mapped candidate loci, contig tiers, supporting-feature flags and loci eligible for tree assessment. The standalone source scripts refer to the original project layout and additional local modules; they require complete input configuration before execution.

## Homolog sampling and tree construction

The historical [hgt_phylo_pipeline directory](../source_reference/analysis/revision_pnas/HGT_candidate_verify/phylogentic_tree_test/hgt_phylo_pipeline/) now contains the configuration, shell entry points and Python utilities. The stages are:

| Source entry | Input | Output or operation |
|---|---|---|
| `run_00_prepare_targets.sh` | Curated candidate table | Target identifiers for tree construction |
| `run_01_sample_homologs.sh` | Targets, search hits, taxonomy and recipient taxids | Selected homologs from self, close and distant groups |
| `run_02_extract_fastas.sh` | Selected identifiers, original reference FASTA and query sequences | Amino-acid FASTA for each candidate |
| `run_03_align_trim_tree.sh` | Per-candidate amino-acid FASTA | MAFFT alignment, trimAl alignment and IQ-TREE gene tree |

The included `config.sh` specifies sampling caps of 12 self-group, 30 close-group, 70 out-group and 10 remote-eukaryote references, with a maximum total of 130. Minimum alignment length is 80 and minimum identity is 20%. These are the saved configuration values, not evidence that every tree contains the maximum number of references. The selected-reference files are needed to establish the actual sampling for each tree.

The included tree-building script uses:

```bash
mafft --auto --thread "$threads" input.aa.fa > aln.fa
trimal -in aln.fa -out aln.trim.fa -automated1
iqtree2 -s aln.trim.fa -m MFP -B 1000 --alrt 1000 -T "$threads" -pre iqtree
```

The script skips a job with fewer than four sequences before or after trimming. Its configuration sets four threads per tree and up to eight parallel jobs. Required tools include MAFFT, trimAl, IQ-TREE 2 and Python with pandas; the extraction utilities use sequence databases and may use seqkit. Configure input paths and available resources before running. Local prefixes in the published source copies are literal placeholders, as explained in the [release README](../README.md).

## Gene-tree assessment

[classify_hgt_trees_v2_no_close_aware.py](../source_reference/analysis/revision_pnas/virus_supplement/scripts/classify_hgt_trees_v2_no_close_aware.py) reads candidate tree files and optional ORF-quality metadata. It examines query placement, patristic distances, nearby self/close/distant references, local mixed branches, branch support and sequence-quality flags.

Relevant defaults in this classifier include:

- SH-aLRT ≥80 and UFBoot ≥95 as strong thresholds when those values are available; unavailable support must be inspected separately.
- A terminal branch is flagged when its length is ≥2 substitutions per site and at least ten times the median of the other terminal branches, or when its absolute length is ≥5.
- Absent or sparse close-group sampling is recorded explicitly. A tree with only distant and self-group references is not equivalent to a tree with informative close-group alternatives.

The classifier writes per-tree evidence, category summaries, identifiers grouped by category and records for manual review. It summarizes tree evidence; its category labels alone are not the final catalogue-retention rule. Original retained assignments and subsequent review decisions must be interpreted with their provenance.

For the additional viral candidates, [finish_candidate_gates.py](../source_reference/analysis/revision_pnas/virus_supplement/reporting_cohort_20260915/final_closure/scripts/finish_candidate_gates.py) combines the original screening and exclusion checks, mapping and contig evidence, flank support, sequence consistency, tree categories, long-branch checks and deduplication. Its tree checks also require valid products and saved UFBoot convergence ≥0.99. This supplementary-candidate assessment should not be mistaken for a new uniform classification of the historical catalogue.

## Inspecting the released evidence

- [Final catalogue](../source_data/tables/S1/Catalogue_8989.tsv): retained records.
- [Gene-tree assessment](../source_data/tables/S1/Gene_tree_assessment.tsv): available per-record tree features and categories.
- [Additional tree support](../source_data/tables/S1/Additional_tree_support.tsv): additional tree evidence.
- [Method parameters](../source_data/tables/S1/Method_parameters.tsv): implementation-specific rules and scopes.
- [Analysis index](final_analysis_manifest.tsv), modules A01–A05 and A09: scripts, input and output correspondence.
- [Artifact index](artifact_manifest.tsv): original-file fingerprints and published locations.

Complete candidate alignments, trees and database snapshots still require a separate versioned data archive. Inclusion of a source script is distinct from a tested full reexecution. The two tested examples in this repository begin with downstream count or matched-family tables.
