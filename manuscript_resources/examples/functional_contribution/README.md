# Matched-family functional contribution example

Run `python examples/functional_contribution/run.py` from `manuscript_resources/`.
The input is the frozen final `family_effects_augmented_same_reference.tsv.gz` table,
with 1,940 species-by-candidate-family observations across 48 species.
`OG` is the original gene-family identifier. Matching of background families has already been done;
this example does not reconstruct the matching or rerun annotation.

For each category, d is the candidate-family proportion minus its matched background proportion.
Family observations are equally weighted within species, species within clades, and the four clades
in the overall estimate. Signed, positive and negative contributions are w*d, w*max(d,0), and w*min(d,0).
Positive contributions are computed per species-family observation before summing by family across species.
Thus, opposite signs for the same family in different species are not cancelled before the positive sum.

Curves rank families by positive contribution and divide by the total positive contribution.
The script reproduces the number of positive families, n50/n80, top-five shares and net differences.
It separately removes the five families with the most candidate loci or the five largest positive
contributors. Each deletion uses its own retained-species baseline, at least ten remaining candidate
families per species and at least three species per retained clade. Matched controls remain fixed.

Expected summaries and deletion results are filtered copies of the saved 48-species manuscript results.
All four category summaries and eight deletion summaries are numerically checked.
The example does not recreate the Poisson-weighted interval calculation or the separate conditional
randomization inference; their source scripts and outputs are indexed in modules A13 and A14.
Annotations and their contributions do not measure pathway flux or gene essentiality.
