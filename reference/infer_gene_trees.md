# Infer Single-Locus Maximum-Likelihood Gene Trees and Assess Species Monophyly

Reconstructs individual gene trees for all curated marker alignments
using `RAxML-NG` (with automated fallback to `phangorn` or `ape`),
assesses reciprocal monophyly for all multi-accession species, and
prepares gene tree collections for coalescent analyses (e.g.,
`ASTRAL-III`).

## Usage

``` r
infer_gene_trees(
  fasta_dir,
  output_dir,
  include_outgroup = FALSE,
  method = c("auto", "raxml", "phangorn", "nj"),
  model = "GTR+G",
  threads = 2L,
  checklist_path = NULL
)
```

## Arguments

- fasta_dir:

  Character. Directory containing curated aligned FASTA sequence files
  (e.g., `4_Cleaned/cleaned_markers_ingroup` or
  `5_MAFFT_Cleaned/aligned_markers`).

- output_dir:

  Character. Root destination directory path to store Newick tree files
  and monophyly audit tables.

- include_outgroup:

  Logical. Include outgroup taxa in single-locus gene tree
  reconstruction? Defaults to `FALSE` to avoid Long-Branch Attraction
  (LBA) artifacts when evaluating species monophyly within the ingroup
  radiation. Set to `TRUE` when preparing unrooted gene trees for
  `ASTRAL-III`.

- method:

  Character. Phylogenetic inference method (`"auto"`, `"raxml"`,
  `"phangorn"`, or `"nj"`). Defaults to `"auto"`.

- model:

  Character. Nucleotide substitution model for maximum-likelihood
  search. Defaults to `"GTR+G"`.

- threads:

  Integer. Number of computational threads. Defaults to `2L`.

- checklist_path:

  Character. Path to accepted botanical checklist CSV (e.g.,
  `CactaceaeFullList_accepted.csv`).

## Value

A list containing the species monophyly summary table and the per-marker
tree paths.

## References

Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A.
(2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum
likelihood phylogenetic inference. *Bioinformatics*, 35(21), 4453–4455.
[doi:10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

Schliep, K. P. (2011). phangorn: phylogenetic analysis in R.
*Bioinformatics*, 27(4), 592–593.
[doi:10.1093/bioinformatics/btq706](https://doi.org/10.1093/bioinformatics/btq706)

## Examples

``` r
if (FALSE) { # \dontrun{
infer_gene_trees(
  fasta_dir = "4_Cleaned/cleaned_markers_ingroup",
  output_dir = "4_Gene_Trees",
  include_outgroup = FALSE
)
} # }
```
