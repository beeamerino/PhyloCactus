# Prime, Cross-Validate and Date a Tree with treePL

Runs the three stages of the `treePL` protocol of Maurin (2020)
(priming, cross-validation and dating) from R, and returns what each
stage chose. It replaces the shell wrapper this package shipped until
2026-09-02.

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
  cvstop = 1e-08,
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
  pairs with; the shell script's 100 was chosen for a rule that was
  stable under it. Raising `n_prime` while keeping `"lowest"` is a
  change of analysis, not a refinement of one, and it should be recorded
  as such.

- prime_rule:

  `"lowest"` (Maurin) or `"modal"` (the shell script's behaviour). See
  the note on `n_prime` above before changing either.

- cv_method:

  `"randomcv"` (Maurin's recommendation) or `"cv"` (leave-one-out).

- cvstart, cvstop:

  Ends of the smoothing grid. The defaults span 1e+03 down to 1e-08, low
  enough to contain the 1e-06 to 1e-08 range Maurin reports for a tree
  with rescaled branch lengths, which the shell script's floor of 1e-04
  was not.

- treepl_bin:

  Name or path of the `treePL` binary.

- work_dir:

  Directory to run in. Defaults to the current one. `treePL` writes
  beside its configuration, so every intermediate lands here.

- quiet:

  Suppress the per-stage progress messages.

## Value

Invisibly, a list with `smoothing`, `cv_table`, `prime_lines`, `primes`,
`tree_file` (the chronogram written) and the paths of the three
configurations.

## Why this exists in R

The shell script this replaces (github.com/tongjial/treepl_wrapper)
carries no licence, so it could not be redistributed inside a GPL-3
package. It also wrote `smoothing = ` into the final configuration, a
keyword `treePL` does not recognise: the line was discarded without a
message and every chronogram produced by this project before 2026-09-02
was dated at the built-in default of 10, with the cross-validation that
precedes it having no effect on any result. Writing the configuration
here puts that keyword under `.verify_treepl_smoothing()`, which
compares what was asked against what `treePL`'s own log reports it used.

## Where this follows the protocol and the wrapper did not

- **Priming selection.** Maurin instructs repeating the priming analysis
  and taking the lines with the *lowest* `opt` and `optad`. The shell
  script took the *most frequent* combination across repeats.
  `prime_rule` defaults to Maurin's rule; `"modal"` reproduces the old
  behaviour.

- **Cross-validation method.** Maurin recommends `randomcv`, random
  subsample and replicate cross-validation, over the leave-one-out `cv`,
  as "much faster and may give more stable results". The shell script
  hard-coded `cv`. `cv_method` defaults to `randomcv`.

- **The grid.** The shell script hard-coded a grid from 1e-04 to 1e+04.
  Both ends are arguments here, and a minimum landing on either end is
  reported rather than returned as if it were an optimum. See
  `.select_cv_smoothing()`.

- **Parsing.** The priming block is read by keyword rather than by
  taking the last character of each line and reassembling positionally,
  which misaligns whenever an optional `moredetail` flag is absent.

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
