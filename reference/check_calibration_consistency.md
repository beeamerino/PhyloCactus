# Check That a Set of Calibration Bounds Is Internally Satisfiable

Bounds are declared per node, but the nodes are nested, and an
ultrametric tree forces every ancestor to be older than every one of its
descendants. Two bounds taken from analyses on different timescales
routinely violate that ordering, and `treePL` does not report the
violation: it returns a chronogram with the offending node pinned to a
bound, which reads as an estimate.

## Usage

``` r
check_calibration_consistency(calibs, tree, constraints, strict = TRUE)
```

## Arguments

- calibs:

  Data frame read from `calibrations_bounds.csv`. Rows with
  `used_in_analysis` not `TRUE` are ignored.

- tree:

  A rooted `phylo` object, typically the output of
  [`root_on_clade()`](https://beeamerino.github.io/PhyloCactus/reference/root_on_clade.md).

- constraints:

  Data frame of the taxonomic table (`cactus_constraints.csv`).

- strict:

  Logical. When `TRUE`, the default, an infeasible pair raises an error.

## Value

Invisibly, a data frame with one row per active calibration: resolved
node, terminals matched, and whether the set is monophyletic in `tree`.

## Details

- Infeasible pair:

  The ancestor's maximum is at or below the descendant's minimum. No
  assignment of ages satisfies both. Reported as an error.

- Non-binding ancestor floor:

  The ancestor's minimum is younger than the descendant's minimum, so
  the ancestor's own lower bound constrains nothing and the effective
  floor is inherited from the descendant. The signature of two sources
  on incompatible timescales, reported as a warning because the run is
  still valid.

A descendant maximum below the ancestor minimum is deliberately not
flagged: it only states that the descendant is necessarily younger,
which is the normal condition for nested bounds.

Two further conditions are checked because they silently reassign bounds
rather than break the run: a row resolving to fewer than two terminals,
which the configuration builder drops without reporting, and two rows
resolving to the same node, which makes the later `min`/`max` pair
overwrite the earlier one.

## Examples

``` r
if (FALSE) { # \dontrun{
calibs <- read.csv(system.file("extdata", "calibrations_bounds.csv", package = "PhyloCactus"))
constraints <- read.csv(system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus"))
rooted <- root_on_clade(ml_tree, resolve_rooting_outgroup(ml_tree$tip.label))
check_calibration_consistency(calibs, rooted, constraints)
} # }
```
