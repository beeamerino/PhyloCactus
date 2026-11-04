# Prime, Cross-Validate and Date a Tree with treePL

Runs the three stages of the `treePL` protocol of Maurin (2020)
(priming, cross-validation and dating) from R, and returns what each
stage chose.

## Usage

``` r
run_treePL_cv(
  cfg_file,
  tree_file,
  label,
  n_prime = 10L,
  prime_rule = c("lowest", "modal"),
  cv_method = c("randomcv", "cv"),
  cvstart = 1000,
  cvstop = 1e-14,
  cv_nthreads = 1L,
  cv_outlier_factor = 5,
  treepl_bin = "treePL",
  work_dir = NULL,
  quiet = FALSE
)
```

## Arguments

- cfg_file:

  Path to the configuration carrying `numsites`, the calibrations and
  any `nthreads`/`seed` lines. The `treefile` line is written by this
  function.

- tree_file:

  Path to the rooted tree to date.

- label:

  Run label. Every file this function writes carries it.

- n_prime:

  Integer. Priming repeats. The default is 10.

  **`n_prime` and `prime_rule` interact, and not symmetrically.** Under
  `"modal"` the selection is stable in `n_prime`: more repeats sharpen
  an estimate of the most frequent combination. Under `"lowest"` it is
  not, because the minimum of a sample can only fall as the sample
  grows, so a hundred repeats will select lower parameters than ten of
  the same runs would. Maurin's instruction is to repeat the priming
  analysis a few times and take the lowest, which is what the default
  pairs with. Raising `n_prime` while keeping `"lowest"` is a change of
  analysis, not a refinement of one, and it should be recorded as such.

- prime_rule:

  `"lowest"` (Maurin) or `"modal"` (most frequent combination). See the
  note on `n_prime` above before changing either.

- cv_method:

  `"randomcv"` (Maurin's recommendation) or `"cv"` (leave-one-out).

- cvstart, cvstop:

  Ends of the smoothing grid, evaluated at one value per order of
  magnitude. The defaults span 1e+03 down to 1e-14. Maurin (2020)
  reports optima between 1e-06 and 1e-08 for trees whose branch lengths
  are rescaled as this package rescales them; the reference Cactaceae
  dataset has an interior minimum at 1e-05. Only a grid extending well
  beyond that range tells an interior minimum from an edge or a plateau.
  See the section on curve shapes below.

- cv_nthreads:

  Integer. Number of threads for the cross-validation stage. Defaults to
  `1L`. `treePL` evaluates cross-validation replicates in an OpenMP
  loop, inside which simulated annealing calls the non-reentrant C
  `rand()`; the chi-square sum is also accumulated through an OpenMP
  reduction, whose order may vary between runs. With more than one
  thread, two runs with the same `seed` return different
  cross-validation curves and can select different smoothing values.
  With one thread the curve is reproducible for a given `seed`, at a
  cost of about 4% in wall-clock time on the reference dataset.

- cv_outlier_factor:

  One number above 1, or `NULL`. A random subsample of the
  cross-validation that fails, or adds more than this many times the
  median chi-square increment at any smoothing value of 0.1 or below, is
  left out of the selection. Defaults to `5`; `NULL` selects on the sum
  over every subsample, as before 0.5.3. See the section on anomalous
  subsamples below.

- treepl_bin:

  Name or path of the `treePL` binary.

- work_dir:

  Directory to run in. Defaults to the current one. `treePL` writes
  beside its configuration, so every intermediate lands here.

- quiet:

  Suppress the per-stage progress messages.

## Value

Invisibly, a list with `smoothing`, `cv_table`, `cv_at_edge`,
`cv_shape`, `cv_failures` (per smoothing value, the subsamples evaluated
and how many failed; `NULL` when the log carries no per-subsample line),
`cv_smoothing_without_failures` (see the section on failed subsamples
below), `cv_smoothing_all_subsamples`, `cv_excluded_subsamples` and
`cv_table_selection` (see the section on anomalous subsamples below),
`prime_lines`, `primes`, `tree_file` (the chronogram written) and the
paths of the three configurations.

## Choices relative to the protocol

- **Priming selection.** Maurin instructs repeating the priming analysis
  and taking the lines with the *lowest* `opt` and `optad`. `prime_rule`
  defaults to this rule; `"modal"` takes the most frequent combination
  across repeats instead.

- **Cross-validation method.** Maurin recommends `randomcv`, random
  subsample and replicate cross-validation, over the leave-one-out `cv`,
  as "much faster and may give more stable results". `cv_method`
  defaults to `randomcv`.

- **The grid.** Both ends of the smoothing grid are arguments, and a
  minimum on either end is reported as such, not returned as an optimum.
  See the section on curve shapes below.

- **Parsing.** The priming block is read by keyword, so an absent
  optional `moredetail` flag does not shift the values.

- **Smoothing check.** The smoothing that `treePL` reports in its log is
  compared with the value written to the configuration, and the run
  stops if they differ.

## Shape of the cross-validation curve

The selected smoothing value is the one with the lowest chi-square.
`cv_shape` reports what that value means, and a message or warning
states it when the function runs:

- `"interior"`: the curve turns upward inside the grid; the selected
  value is an optimum.

- `"edge"`: the minimum is the first or last value of the grid and the
  curve is still changing there. The optimum most likely lies beyond the
  grid: extend `cvstart` or `cvstop` and repeat before interpreting any
  age (Maurin, 2020). Reported as a warning.

- `"plateau"`: the three values at the end of the grid where the curve
  flattens differ by less than 1% in chi-square and lie within 1% of the
  minimum. Values on the plateau fit the data about equally well, so the
  selected value is a nominal choice among them; report node ages across
  the plateau with
  [`report_smoothing_sensitivity()`](https://beeamerino.github.io/PhyloCactus/reference/report_smoothing_sensitivity.md).
  Reported as a message.

## Failed subsamples

Under `randomcv`, `treePL` writes in the log of the cross-validation
(`cv_<label>.log`) the running chi-square sum after each random
subsample. A subsample whose contribution is 1e+30 or more has failed,
and the total at that smoothing value is then of the same order, so the
value drops out of the selection. The function counts the failed
subsamples at every value and warns when one or more failed at any value
of 0.1 or below, the part of the grid where a selection falls in
practice: the warning names each such value with how many of how many
subsamples failed there, says whether the selected value is among them,
and names the smoothing with the lowest chi-square among the values of
0.1 or below at which no subsample failed, which is also returned as
`cv_smoothing_without_failures`. With `cv_outlier_factor = NULL`, a
value with a failed subsample carries a total of the order of the
failure, so that smoothing is the selected value whenever the selection
falls at 0.1 or below, and the warning then says so. With the default,
the failed subsample is left out (see the section on anomalous
subsamples), the selection can fall on a value where it failed, and the
warning says that too. Values above 0.1 are left out: on the reference
dataset no run has selected a value above 0.1, and subsamples fail often
there. The running sum is printed with eight significant digits, so the
chi-square of the subsamples after a failure cannot be recovered from
the log, and the count of failures is a lower bound.

## Anomalous subsamples

A failed subsample is not the only one that can move the selection. In
the re-dating of the reference dataset under one fixed-age calibration
(Tutorial 2, Module 10), one or two of the ten subsamples of each
cross-validation fitted badly from 0.1 down to about 1e-07, or failed,
and normally below; summed over every subsample, the curves selected the
smallest values of the grid, and summed without those subsamples each
had its minimum at 1e-05. The default factor of 5 was set after those
curves had been read. With `cv_outlier_factor` set, a subsample that
fails at any value of 0.1 or below, or adds more than
`cv_outlier_factor` times the median increment of the subsamples at such
a value, is left out, and the smoothing is selected on the curve of the
remaining subsamples over the values of 0.1 or below. The values are
compared on one set of subsamples: because the subsamples after a
failure are lost at that value, either those values or those subsamples
are left out, whichever keeps more subsample-by-value cells. The curve
is classified as the plain sum is (see above). When no subsample is
anomalous the selection is that of the plain sum. When fewer than five
subsamples or five values would remain, the function warns and selects
on the plain sum. A warning names the subsamples left out and why, and
gives the selection over every subsample, returned as
`cv_smoothing_all_subsamples`; both curves are written to
`cv_<label>_selection.csv`, and the account to
`cv_<label>_selection.txt`.

## References

Maurin, K. J. L. (2020). An empirical guide for producing a dated
phylogeny with treePL in a maximum likelihood framework.
*arXiv*:2008.07054.
[doi:10.48550/arXiv.2008.07054](https://doi.org/10.48550/arXiv.2008.07054)

Sanderson, M. J. (2002). Estimating absolute rates of molecular
evolution and divergence times: a penalized likelihood approach.
*Molecular Biology and Evolution*, 19(1), 101-109.
[doi:10.1093/oxfordjournals.molbev.a003974](https://doi.org/10.1093/oxfordjournals.molbev.a003974)

Smith, S. A., & O'Meara, B. C. (2012). treePL: divergence time
estimation using penalized likelihood for large phylogenies.
*Bioinformatics*, 28(20), 2689-2690.
[doi:10.1093/bioinformatics/bts492](https://doi.org/10.1093/bioinformatics/bts492)
