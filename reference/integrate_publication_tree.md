# Render Final Publication Figures and Registry

Maps statistical support values (e.g., Transfer Bootstrap Expectation,
TBE) onto nodes of the final chronogram and maximum-likelihood
phylogeny. Nodes failing to meet the minimum support threshold
(\`collapse_cutoff\`) are systematically collapsed into soft polytomies
(analytical uncertainty) to produce conservative, publication-ready
figures.

## Usage

``` r
integrate_publication_tree(
  ml_support_tree_path,
  summary_chronogram_path,
  constraints_path,
  out_dir,
  collapse_cutoff = 0.7
)
```

## Arguments

- ml_support_tree_path:

  Character. Best ML support tree path.

- summary_chronogram_path:

  Character. Chronogram path with HPD annotations.

- constraints_path:

  Character. Taxonomy constraints CSV path.

- out_dir:

  Character. Publication figures directory.

- collapse_cutoff:

  Numeric. TBE threshold to collapse weak nodes into soft polytomies
  (0.0 to 1.0; default 0.70).

## Value

A data frame listing exported figure files and threshold parameters.
