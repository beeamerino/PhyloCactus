# Classify Every Internal Node as Constrained or Estimated

Compares the bipartitions of a reference topology against those of the
topological constraint it was inferred under, and reports for each
internal node whether the split was imposed or estimated. Support values
attached to imposed splits are returned unchanged but must not be
reported as support: when the bootstrap replicates were themselves run
with `--tree-constraint`, every replicate reproduces those splits and
the value is 1 by construction.

## Usage

``` r
classify_constrained_nodes(tree, constraint_tree)
```

## Arguments

- tree:

  An `ape` `phylo` object, or a path to a Newick file. Pass the
  annotated support tree (`.raxml.support`) to carry support values
  through into the result.

- constraint_tree:

  An `ape` `phylo` object, or a path to a Newick file, holding the
  topology passed to `RAxML-NG --tree-constraint`. Terminals absent from
  `tree` (collapsed as identical during parsing) are dropped before
  comparison.

## Value

A `data.frame` with one row per internal node: `node` (index in the
`ape` numbering), `n_tips` (size of the smaller side of the
bipartition), `constrained` (logical) and `support` (numeric, `NA` when
the tree carries no node labels).

## Examples

``` r
if (FALSE) { # \dontrun{
nodes <- classify_constrained_nodes(
  tree            = "7_Phylogenetics/cactus_support_fbp.raxml.support",
  constraint_tree = "7_Phylogenetics/cactus_constraints.tree"
)
table(nodes$constrained)
summary(nodes$support[!nodes$constrained])
} # }
```
