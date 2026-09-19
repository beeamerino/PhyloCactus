# Comparative Tree-Space Projection and Topological Validation Pipeline

Compares the focal tree produced by `PhyloCactus` with published
phylogenies (e.g., Amaral *et al.* 2022, Thompson *et al.* 2024, Zuntini
*et al.* 2024, de Vos *et al.* 2025). Pairwise distances are computed
with the Robinson-Foulds distance (Robinson & Foulds, 1981) and the
Different Phylogenetic Information distance (`TreeDist`; Smith, 2020),
both on unrooted trees, and projected into tree space by
multidimensional scaling (MDS).

## Usage

``` r
validate_phylogenies(
  trees_mapping_list,
  checklist_csv,
  constraints_map,
  out_dir,
  min_tips = 4L,
  pruning_strategy = c("common_set", "pairwise"),
  rooting_genus = "Leuenbergeria",
  tree_labels = NULL
)
```

## Arguments

- trees_mapping_list:

  Named list of file paths pointing to Newick/Nexus reference trees and
  focal tree.

- checklist_csv:

  Character. Path to accepted botanical checklist CSV or Excel file for
  taxon name standardization (Korotkova *et al.*, 2021).

- constraints_map:

  Character. Path to taxonomy mapping CSV file.

- out_dir:

  Character. Destination directory path to save validation output
  tables, tree-space plots, and logs.

- min_tips:

  Integer. Minimum number of tips for comparisons. Defaults to 4L.

- pruning_strategy:

  Character. Strategy for pruning trees before metric computation. Under
  `"common_set"`, the shared intersection across all trees is used for
  the multidimensional scaling (MDS) projection, whereas pairwise tables
  retain all shared tips between each pair. Defaults to `"common_set"`.

- rooting_genus:

  Character. Genus whose terminals root every tree after the family
  filter. Defaults to `"Leuenbergeria"`.

  This differs from the rooting used elsewhere in `PhyloCactus`. The
  comparison is restricted to **Cactaceae**, so the family filter
  removes the outgroup terminals used by the dating pipeline.
  *Leuenbergeria* is sister to the rest of the family in the focal tree.
  The distances and the tree-space projection are computed on unrooted
  trees and do not depend on this choice; the rooting affects only the
  exported trees. All terminals of the genus are used, not one, for the
  reasons given in
  [`resolve_rooting_outgroup()`](https://beeamerino.github.io/PhyloCactus/reference/resolve_rooting_outgroup.md).
  Trees lacking the genus are left unrooted and flagged in
  `SUPP_rooting_summary`.

- tree_labels:

  Named character vector or `NULL`. Legend labels for the trees in the
  tree-space figure, named by the entries of `trees_mapping_list`.
  Labels of the form `"Author et al. (year)"` are typeset with *et al.*
  in italics. Defaults to `NULL`, which labels the reference trees
  distributed with the package as `"Zuntini et al. (2024)"`,
  `"Thompson et al. (2024)"`, `"Amaral et al. (2022)"` and
  `"de Vos et al. (2025)"`, and the focal tree as `"This study"`.
  Supplied labels replace the defaults of the same name.

## Value

A named list containing file paths for generated tables (`tables`),
figures (`figures`), trees (`trees`), and the combined ggplot object
(`plot`).

## References

Korotkova, N., Aquino, D., Arias, S., Eggli, U., Franck, A.,
Gómez-Hinostrosa, C., Guerrero, P. C., Hernández, H. M., Kohlbecker, A.,
Köhler, M., Luther, K., Majure, L. C., Müller, A., Metzing, D.,
Nyffeler, R., Sánchez, D., Schlumpberger, B. O., & Berendsohn, W. G.
(2021). Cactaceae at Caryophyllales.org - a dynamic online species-level
taxonomic backbone for the family. *Willdenowia*, 51(2), 251-270.
[doi:10.3372/wi.51.51208](https://doi.org/10.3372/wi.51.51208)

Robinson, D. F., & Foulds, L. R. (1981). Comparison of phylogenetic
trees. *Mathematical Biosciences*, 53(1-2), 131-147.
[doi:10.1016/0025-5564(81)90043-2](https://doi.org/10.1016/0025-5564%2881%2990043-2)

Smith, M. R. (2020). Information theoretic generalized Robinson-Foulds
metrics for comparing phylogenetic trees. *Bioinformatics*, 36(20),
5007-5013.
[doi:10.1093/bioinformatics/btaa614](https://doi.org/10.1093/bioinformatics/btaa614)

## Examples

``` r
if (FALSE) { # \dontrun{
validate_phylogenies(
  trees_mapping_list = list(
    FocalTree = "bestTree.tree",
    Amaral = "Tree_80MD.tre"
  ),
  checklist_csv = "checklist.csv",
  constraints_map = "taxonomy.csv",
  out_dir = "9_Validation"
)
} # }
```
