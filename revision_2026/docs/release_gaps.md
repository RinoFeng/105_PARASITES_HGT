# Release scope and unresolved dependencies

1. **Measurement scope:** Figure 2B currently uses 8,985 measurement records, while one manuscript
sentence states 6,882 protein records. Decide between reporting the existing measured subset or adding
four verified measurement records and recalculating the affected figure/table. Protein and nucleotide
lengths must be tested in their own units; the older plotting script includes a mixed-unit length test.
2. **Discovery provenance:** keep the historical MMseqs workflow and the later DIAMOND analysis distinct.
The exact database snapshots and protein preprocessing stages (95%, 90%, 50%) require a consistent final
methods account. This package does not claim that the entire catalogue was newly screened under one command.
3. **Large external resources:** raw searches, reference genomes/proteomes, complete tree/sequence alignments,
OrthoFinder outputs and taxonomy snapshots are not bundled. Provide a versioned data archive and an index
linking candidate IDs to phylogenetic evidence. Document database dates and accessions.
4. **Source scripts:** `source_reference/` preserves selected implementations for inspection; it is not a
portable end-to-end workflow. Local directory prefixes were redacted without changing scientific thresholds.
Restore explicit configuration and complete input dependencies before running these scripts.
5. **PRO-seq:** the primary denominator uses all counted reference genes. A matrix restricted to the 27
candidate-associated genes cannot recreate it. Publish the full count/length matrix or the exact original
normalization factors. No new differential nascent-transcription model is provided here.
6. **Candidate trees:** do not equate the newer automatic classifier with the final historical judgement
for all 8,989 entries. Publish the original retained categories and any manual decisions with stable IDs.
7. **Optional analyses:** A15 is retained only if ecological/Leishmania comparisons remain cited in the final
manuscript. A23 is a recent deduplication sensitivity analysis not yet incorporated in the manuscript.
8. **Local-genome illustrations:** the final independent generation entry for the Plasmodium microsynteny
illustration is not located in this preparation. The PCR table is an experimental record, not an executable
analysis. Neither is presented as reproduced by these examples.
9. **Public release:** no files have been pushed, no data DOI has been assigned, and no release has been
published by this preparation. The manuscript must cite the final tested commit/release and data record.

The completed deliverables are the indexed source package and the two tested downstream examples.
Complete raw-data reexecution remains a separate release milestone.
