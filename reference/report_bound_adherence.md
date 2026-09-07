# Report Which Calibrated Nodes Came Back Sitting on a Bound

A node whose estimated age equals one of its own bounds was not
estimated. Penalized likelihood returned the constraint, and the number
carries the prior rather than the data. This is invisible in the output
chronogram, which looks like any other, so it has to be checked
explicitly before a date is reported or interpreted.

## Usage

``` r
report_bound_adherence(
  chronogram,
  calibs,
  constraints,
  tol = 0.05,
  bootstraps = NULL
)
```

## Arguments

- chronogram:

  An ultrametric, rooted `phylo` object, typically
  `BestTree_treePL.tree`.

- calibs:

  Data frame read from `calibrations_bounds.csv`.

- constraints:

  Data frame of the taxonomic table (`cactus_constraints.csv`).

- tol:

  Numeric. Absolute tolerance, in millions of years, within which an age
  counts as sitting on a bound. Defaults to `0.05`.

- bootstraps:

  Optional. The bootstrap chronograms, as a `multiPhylo`, a list of
  `phylo`, or a path to a Newick file with one tree per line
  (`bsTree_treePL.tree`). When supplied, the verdict is taken from the
  fraction of replicates sitting on a bound rather than from the single
  point estimate, and a node pinned in most replicates is reported even
  when the maximum-likelihood tree alone would have called it interior.
  Defaults to `NULL`.

## Value

Invisibly, a data frame with one row per calibrated node: its age in the
maximum-likelihood chronogram, its bounds, a `status` of `"at_min"`,
`"at_max"` or `"interior"`, and, when `bootstraps` is supplied, the
number of replicates evaluated, their median and 95% interval, and the
percentage of them sitting on a bound.

## Details

In the August 2026 run every one of the five calibrated nodes came back
on a bound, the root at its maximum and the rest at their minimum, and
the chronogram gave no sign of it.

**A single chronogram is not enough to answer this.** On 2026-09-02 the
maximum-likelihood tree placed `ACP_root` at 52.96 Ma, 0.41 Ma inside
its upper bound of 53.37, and this function reported it as interior.
Across the 100 bootstrap replicates the interval was 53.29 to 53.37 and
96 of them returned the bound exactly. The point estimate was one
realisation of a node whose age the data cannot identify, and it
happened to land just inside. Supplying `bootstraps` is what
distinguishes a node that was estimated from one that is unidentifiable
and collapsed onto its nearest constraint.

## Examples

``` r
if (FALSE) { # \dontrun{
chrono <- ape::read.tree(file.path("8_Dating", "BestTree_treePL.tree"))
calibs <- read.csv(system.file("extdata", "calibrations_bounds.csv", package = "PhyloCactus"))
constraints <- read.csv(system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus"))
report_bound_adherence(chrono, calibs, constraints,
                       bootstraps = file.path("8_Dating", "bsTree_treePL.tree"))
} # }
```
