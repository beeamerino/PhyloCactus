# Evaluate DNA Barcoding Resolution and Taxonomic Classification Power

Not implemented. This function is a documented placeholder that signals
an error when called. It is retained so that the section keeps a stable
entry point in the reference index and in the tutorial series while its
analytical strategy is being defined.

## Usage

``` r
evaluate_dna_barcoding(
  fasta_dir,
  checklist_path,
  output_dir,
  distance_models = c("raw", "K80"),
  train_decipher = TRUE,
  confidence_threshold = 50,
  max_subsample = 20L
)
```

## Arguments

- fasta_dir:

  Character. Directory containing curated ingroup FASTA alignments.

- checklist_path:

  Character. Path to the accepted botanical checklist (Korotkova *et
  al*. 2021).

- output_dir:

  Character. Destination directory for barcoding summary tables.

- distance_models:

  Character vector. Substitution models for pairwise distance
  estimation.

- train_decipher:

  Logical. Train a taxonomic classifier?

- confidence_threshold:

  Numeric. Bootstrap confidence threshold for classification.

- max_subsample:

  Integer. Maximum number of sequences per species used in distance
  calculations.

## Value

None. The function signals an error.

## Details

The implementation withdrawn in v0.4.2 computed intra- and
inter-specific pairwise distance distributions, derived the DNA barcode
gap per locus, and trained a
[`DECIPHER::LearnTaxa`](https://rdrr.io/pkg/DECIPHER/man/LearnTaxa.html)
classifier against the botanical backbone of Korotkova *et al*. (2021).
It was withdrawn for two reasons. First, its accuracy estimates were
produced by evaluating the trained classifier on the same sequence set
used for training, which measures memorisation rather than diagnostic
performance and therefore reported values that could not be interpreted.
Second, the barcode gap thresholds it applied were fixed rather than
derived from the empirical distributions of the family, a decision that
cannot be made before the phylogenetic backbone and its visualisation
are complete.

The full rationale, the design decisions already taken, and the criteria
that must be met before the section is reopened are documented in
[`vignette("tutorial-5-cactus-phylogeny-barcoding")`](https://beeamerino.github.io/PhyloCactus/articles/tutorial-5-cactus-phylogeny-barcoding.md).
The withdrawn implementation remains recoverable from the repository
history.

## References

Korotkova, N., Aquino, D., Arias, S., Eggli, U., Franck, A.,
Gómez-Hinostrosa, C., Guerrero, P. C., Hernández, H. M., Kohlbecker, A.,
Köhler, M., Luther, K., Majure, L. C., Müller, A., Metzing, D.,
Nyffeler, R., Sánchez, D., Schlumpberger, B. O., & Berendsohn, W. G.
(2021). Cactaceae at Caryophyllales.org - a dynamic online species-level
taxonomic backbone for the family. *Willdenowia*, 51(2), 251-270.
[doi:10.3372/wi.51.51208](https://doi.org/10.3372/wi.51.51208)

## See also

[`infer_gene_trees()`](https://beeamerino.github.io/PhyloCactus/reference/infer_gene_trees.md)
for single-locus gene tree inference, which remains available.

## Examples

``` r
if (FALSE) { # \dontrun{
# Signals an error and points to the manifest:
evaluate_dna_barcoding("4_Cleaned/cleaned_markers_ingroup", "checklist.csv", "out")
} # }
```
