# Audit the Nodes a Calibration Table Addresses

For every row of a calibration table, resolves the node its terminals
address in the reference topology and reports whether that node was
imposed by the topological constraint and what support it carries. This
is the check that has to pass before a chronogram is interpreted: a
bound placed on a constrained node dates a node the analysis assumed,
and a bound placed on a weakly supported node propagates that
uncertainty into every age downstream of it without recording that it
did.

## Usage

``` r
audit_calibration_support(
  calibrations,
  constraints,
  tree,
  constraint_tree,
  support_trees = list(),
  active_only = TRUE
)
```

## Arguments

- calibrations:

  A `data.frame` in the format of
  `inst/extdata/calibrations_bounds.csv`, with columns `mrca`, `column`,
  `value`, `min`, `max` and `used_in_analysis`.

- constraints:

  A `data.frame` in the format of `inst/extdata/cactus_constraints.csv`,
  used to translate a taxon-set name into terminals.

- tree:

  An `ape` `phylo` object, or a path to a Newick file, holding the
  reference topology.

- constraint_tree:

  An `ape` `phylo` object, or a path to a Newick file, holding the
  topology passed to `RAxML-NG --tree-constraint`.

- support_trees:

  Named list of annotated support trees (`phylo` objects or paths), one
  per metric, for example
  `list(fbp = "cactus_support_fbp.raxml.support", tbe = "...")`. One
  column is added per element. Defaults to an empty list.

- active_only:

  Logical. Restrict to rows with `used_in_analysis == TRUE`. Defaults to
  `TRUE`.

## Value

A `data.frame` with one row per calibration: `mrca`, `n_tips`, `node`,
`constrained`, `is_root`, one support column per entry of
`support_trees`, and `min` and `max`.

## Examples

``` r
if (FALSE) { # \dontrun{
audit_calibration_support(
  calibrations    = read.csv(system.file("extdata", "calibrations_bounds.csv",
                                         package = "PhyloCactus")),
  constraints     = read.csv(system.file("extdata", "cactus_constraints.csv",
                                         package = "PhyloCactus")),
  tree            = "7_Phylogenetics/ml_search/cactus_search.raxml.bestTree",
  constraint_tree = "7_Phylogenetics/cactus_constraints.tree",
  support_trees   = list(fbp = "7_Phylogenetics/cactus_support_fbp.raxml.support",
                         tbe = "7_Phylogenetics/cactus_support_tbe.raxml.support")
)
} # }
```
