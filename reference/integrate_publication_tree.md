# Render Final Publication Figures and Registry

Maps statistical support values (e.g., Felsenstein Bootstrap
Proportions, FBP) onto nodes of the final chronogram and
maximum-likelihood phylogeny. Nodes with support below `collapse_cutoff`
are collapsed into soft polytomies in the figures.

## Usage

``` r
integrate_publication_tree(
  ml_support_tree_path,
  summary_chronogram_path,
  constraints_path = NULL,
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

  Not used in this version; accepted so that earlier calls still run.

- out_dir:

  Character. Publication figures directory.

- collapse_cutoff:

  Numeric. Felsenstein Bootstrap Proportion (FBP) threshold below which
  weakly supported nodes are collapsed into soft polytomies (0.0 to 1.0;
  default 0.70; collapsing 840 of 1171 free nodes in the reference
  dataset).

## Value

A data frame listing exported figure files and threshold parameters.
