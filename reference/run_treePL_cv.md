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
  dataset reaches a plateau below 1e-06, which only a grid extending
  well below that range can show. See the section on curve shapes below.

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

- treepl_bin:

  Name or path of the `treePL` binary.

- work_dir:

  Directory to run in. Defaults to the current one. `treePL` writes
  beside its configuration, so every intermediate lands here.

- quiet:

  Suppress the per-stage progress messages.

## Value

Invisibly, a list with `smoothing`, `cv_table`, `cv_at_edge`,
`cv_shape`, `prime_lines`, `primes`, `tree_file` (the chronogram
written) and the paths of the three configurations.

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
