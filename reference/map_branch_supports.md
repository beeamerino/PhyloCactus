# Map Bootstrap Support Values onto Reference Phylogeny

Maps clade support derived from non-parametric bootstrap replicates onto
the best maximum-likelihood tree topology. The replicates themselves
carry no metric: `RAxML-NG --bootstrap` writes plain topologies with
branch lengths, and the metric is chosen here, at the summarising step.
The same replicate file can therefore be summarised under both metrics
without recomputation.

## Usage

``` r
map_branch_supports(
  raxml_bin,
  best_tree,
  bootstraps_file,
  metric = c("fbp", "tbe"),
  threads = 4,
  output_dir = dirname(best_tree),
  prefix = "cactus_support"
)
```

## Arguments

- raxml_bin:

  Character. System command or full path to executable `RAxML-NG`
  binary.

- best_tree:

  Character. Path to reference maximum-likelihood tree file.

- bootstraps_file:

  Character. Path to concatenated non-parametric bootstrap trees file.

- metric:

  Character. Bootstrap support metric: `"fbp"` (Felsenstein's Bootstrap
  Percentage, the default and the comparable metric) or `"tbe"`
  (Transfer Bootstrap Expectation, secondary).

- threads:

  Integer. Number of CPU threads. Defaults to `4`.

- output_dir:

  Character. Output directory for annotated support tree. Defaults to
  `dirname(best_tree)`.

- prefix:

  Character. Output file prefix. Defaults to `"cactus_support"`.

## Value

Character path to the annotated support tree file (`.raxml.support`).

## Details

Defaults to Felsenstein's Bootstrap Percentage (FBP; Felsenstein, 1985),
which is the metric the Cactaceae and Caryophyllales dating literature
reports and the only one against which this tree can be compared.
Transfer Bootstrap Expectation (TBE; Lemoine *et al.*, 2018) is
available through `metric = "tbe"` and belongs in a clearly labelled
secondary column.

The two are not on a common scale and TBE must never be reported as
though it were a bootstrap percentage. TBE is bounded below by FBP and
its inflation grows with clade size. Measured on the 986 unconstrained
nodes of the 1023 terminal supermatrix tree (2026-09-06): median TBE
0.783 against median FBP 0.468, with the gap reaching 0.704 for clades
of 51 to 200 terminals. Reporting TBE would place 61 percent of nodes
above 0.70; FBP places 29 percent.

Support values are meaningless for any bipartition imposed through
`--tree-constraint`, because every replicate reproduces it by
construction. Those nodes must be reported as constrained rather than
supported. Use
[`classify_constrained_nodes()`](https://beeamerino.github.io/PhyloCactus/reference/classify_constrained_nodes.md)
to separate them before tabulating support.

## References

Felsenstein, J. (1985). Confidence limits on phylogenies: an approach
using the bootstrap. *Evolution*, 39(4), 783-791.
[doi:10.1111/j.1558-5646.1985.tb00420.x](https://doi.org/10.1111/j.1558-5646.1985.tb00420.x)

Lemoine, F., Domelevo Entfellner, J. B., Wilkinson, E., Correia, D.,
Davila Felipe, M., De Oliveira, T., & Gascuel, O. (2018). Renewing
Felsenstein's phylogenetic bootstrap in the era of big data. *Nature*,
556(7702), 452-456.
[doi:10.1038/s41586-018-0043-0](https://doi.org/10.1038/s41586-018-0043-0)

## Examples

``` r
if (FALSE) { # \dontrun{
# Primary metric, comparable with the published literature.
map_branch_supports(
  raxml_bin = "raxml-ng",
  best_tree = "cactus_search.raxml.bestTree",
  bootstraps_file = "cactus_ALL_bootstraps.tree",
  metric = "fbp",
  prefix = "cactus_support_fbp"
)

# Secondary metric, same replicates, no recomputation.
map_branch_supports(
  raxml_bin = "raxml-ng",
  best_tree = "cactus_search.raxml.bestTree",
  bootstraps_file = "cactus_ALL_bootstraps.tree",
  metric = "tbe",
  prefix = "cactus_support_tbe"
)
} # }
```
