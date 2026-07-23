# Map Transfer Bootstrap Expectation (TBE) Support Values onto Reference Phylogeny

Maps statistical clade support metrics derived from non-parametric
bootstrap replicates onto the best maximum-likelihood tree topology.
Implements Transfer Bootstrap Expectation (TBE; Lemoine \*et al.\*,
2018) as the primary support metric, which provides robust support
evaluation for large plant phylogenies without underestimating support
for minor clade position shifts.

## Usage

``` r
map_branch_supports(
  raxml_bin,
  best_tree,
  bootstraps_file,
  metric = "tbe",
  threads = 4,
  output_dir = dirname(best_tree),
  prefix = "cactus_support"
)
```

## Arguments

- raxml_bin:

  Character. System command or full path to executable \`RAxML-NG\`
  binary.

- best_tree:

  Character. Path to reference maximum-likelihood tree file.

- bootstraps_file:

  Character. Path to concatenated non-parametric bootstrap trees file.

- metric:

  Character. Bootstrap support metric: \`"tbe"\` (Transfer Bootstrap
  Expectation) or \`"fbp"\` (Felsenstein's Bootstrap Percentage).
  Defaults to \`"tbe"\`.

- threads:

  Integer. Number of CPU threads. Defaults to \`4\`.

- output_dir:

  Character. Output directory for annotated support tree. Defaults to
  \`dirname(best_tree)\`.

- prefix:

  Character. Output file prefix. Defaults to \`"cactus_support"\`.

## Value

Character path to the annotated support tree file (\`.raxml.support\`).

## References

Lemoine, F., Entfellner, J. B., Gascuel, O., & Gascuel, O. (2018).
Renewing Felsenstein’s phylogenetic bootstrap in the era of big data.
\*Nature\*, 556(7702), 452-456.
[doi:10.1038/s41586-018-0043-0](https://doi.org/10.1038/s41586-018-0043-0)

## Examples

``` r
if (FALSE) { # \dontrun{
map_branch_supports(
  raxml_bin = "raxml-ng",
  best_tree = "cactus_search.raxml.bestTree",
  bootstraps_file = "cactus_ALL_bootstraps.tree",
  metric = "tbe"
)
} # }
```
