# Automate treePL Divergence Time Estimation Pipeline Across Bootstrap Cohorts

Automates cross-validation parameter optimization, rate smoothing
selection, and chronogram estimation across temporal bootstrap
replicates using `treePL` (Sanderson, 2002; Smith & O'Meara, 2012).
Propagating temporal uncertainty across branch-length resampled
bootstrap trees yields empirical confidence intervals for node age
estimates. Every treePL output (maximum-likelihood chronogram and each
bootstrap chronogram) is validated after execution: the resulting tree
must exist, be non-empty, parse as a valid Newick topology, and be
ultrametric. A run that fails silently (e.g., because the underlying
`treePL` binary did not converge) is therefore reported as an explicit
error rather than propagated downstream as a corrupted chronogram. Note
that the optimal rate-smoothing parameter is cross-validated once on the
maximum-likelihood tree and reused, unmodified, across all bootstrap
replicates; this is a standard computational shortcut for treePL-based
dating pipelines (per-replicate cross-validation is prohibitively
expensive at typical bootstrap replicate counts), following the
empirical protocol of Maurin (2020).

## Usage

``` r
automate_treePL(
  cfg_file,
  ml_tree_file,
  bs_trees_file,
  results_dir,
  treePL_out,
  num_bs = NULL,
  numsites = NULL,
  outgroup = NULL,
  seed = NULL,
  rescale_factor = 100,
  n_prime = 10L,
  prime_rule = c("lowest", "modal"),
  cv_method = c("randomcv", "cv"),
  cvstart = 1000,
  cvstop = 1e-08,
  wrapper_sh = NULL
)
```

## Arguments

- cfg_file:

  Character. Path to primary `treePL` configuration file specifying
  calibration bounds and parameters.

- ml_tree_file:

  Character. Path to input maximum-likelihood reference tree file.

- bs_trees_file:

  Character. Path to input temporal bootstrap trees file.

- results_dir:

  Character. Directory path to save intermediate optimization results.

- treePL_out:

  Character. Destination directory path for final output chronograms.

- num_bs:

  Integer or NULL. Maximum number of temporal bootstrap trees to
  evaluate. If `NULL`, processes all available trees.

- numsites:

  Integer or NULL. Alignment length, in sites, of the supermatrix
  actually analysed. If `NULL`, the `numsites` line already present in
  `cfg_file` is used. This must be the length of the matrix `RAxML-NG`
  analysed (typically the `*.raxml.reduced.phy` produced by
  [`preprocess_partitions()`](https://beeamerino.github.io/PhyloCactus/reference/preprocess_partitions.md)),
  not a rounded figure: `treePL` uses it to convert branch lengths into
  expected substitution counts, so an incorrect value biases the rate
  smoothing. Declare the true alignment length here; the division by
  `rescale_factor` is applied internally.

- outgroup:

  Character vector of terminals used to root the maximum-likelihood tree
  and every bootstrap replicate before penalized-likelihood dating, or
  `NULL`. Defaults to `NULL`, in which case the set is derived from the
  maximum-likelihood tree with
  [`resolve_rooting_outgroup()`](https://beeamerino.github.io/PhyloCactus/reference/resolve_rooting_outgroup.md).
  Adapt `outgroup` or specify `resolve_rooting_outgroup(pattern = ...)`
  to the outgroup lineage sampled in your dataset (for example
  `"^(Talinum|Talinella)_"` when Talinaceae roots the tree). No fixed
  default is offered because a clade-level rooting set depends on what
  the supermatrix sampled and cannot be a package constant.

  **Why a clade rather than a terminal.** Rooting is imposed after the
  search: `RAxML-NG` returns an unrooted topology and the root is placed
  on a chosen edge. Naming a single terminal of a sampled outgroup clade
  places the root *inside* that clade, leaving it paraphyletic in the
  final tree and collapsing its crown node onto the root. Any
  calibration addressed by the MRCA of that clade then lands on the root
  instead of on the node it was written for. Supplying the whole clade
  places the root on its stem edge, which is the intended edge. Rooting
  is performed by
  [`root_on_clade()`](https://beeamerino.github.io/PhyloCactus/reference/root_on_clade.md),
  which also handles the basal polytomy `RAxML-NG` writes.

  Bootstrap replicates are handled with the intersection of this set and
  each replicate's tip labels, so a replicate missing some terminals is
  still rooted; only a replicate missing all of them is an error.

  **Topological assumption.** Placing the root on the outgroup lineage
  (Talinaceae or Portulacaceae) establishes the basal split for dating.
  In the reference dataset, Talinaceae roots the tree, placing the root
  on the stem of the ACP clade (Anacampserotaceae, Cactaceae,
  Portulacaceae) and allowing maximum-likelihood inference to test
  alternative topological resolutions among the three core families
  (Ramirez-Barahona et al., 2020; Zuntini et al., 2024). State this
  assumption in Methods, and treat the root age as conditional on it.

- seed:

  Integer or NULL. Seed for every stochastic step of the run: R's choice
  of which bootstrap replicates to date, and `treePL`'s own `seed`
  keyword, which is written into the maximum-likelihood configuration
  and, offset by the replicate index, into each replicate's. `treePL`
  seeds itself from the clock when the keyword is absent, which leaves
  both the cross-validation and the simulated annealing unreproducible.
  Defaults to `NULL`, which draws one and reports it; record the
  reported value, as it is required to reproduce the run.

- rescale_factor:

  Numeric. Multiplier applied to every branch length before dating, and
  the divisor applied to `numsites` in the configurations written for
  `treePL`. Defaults to `100`.

  **Why the two are one argument.** `treePL` reads a branch as
  `edge.length * numsites` expected substitutions. Branch lengths are
  rescaled because a substantial fraction of a low-divergence plastid
  supermatrix falls below the internal minimum `treePL` imposes on a
  branch, and those branches would otherwise be clamped to a common
  value, erasing the rate signal across them. Rescaling without dividing
  `numsites` by the same factor leaves `treePL` reading a matrix it
  believes to be `rescale_factor` times more informative than it is. The
  likelihood term grows with the substitution counts while the roughness
  penalty does not, so the effective smoothing becomes weaker than the
  nominal value by that factor, rates vary almost freely between
  branches, and node ages stop being determined by the data and start
  being determined by the edges of the region the calibrations leave
  feasible. The symptom is a chronogram whose calibrated nodes sit
  exactly on their bounds. Coupling the two here makes the

- n_prime:

  Integer. Priming repeats on the maximum-likelihood tree. See
  [`run_treePL_cv()`](https://beeamerino.github.io/PhyloCactus/reference/run_treePL_cv.md).

- prime_rule:

  `"lowest"` (Maurin 2020) or `"modal"` (the retired shell script's
  rule).

- cv_method:

  `"randomcv"`, recommended by Maurin (2020) as faster and more stable,
  or `"cv"` for leave-one-out, which is what the retired shell script
  used.

- cvstart, cvstop:

  Ends of the cross-validation smoothing grid. The defaults reach lower
  than the shell script's fixed floor of 1e-04, which the August 2026
  run of this project hit without the chi-square curve ever turning.

- wrapper_sh:

  Retired on 2026-09-02 and ignored, with a warning. Priming,
  cross-validation and dating of the maximum-likelihood tree now run
  through
  [`run_treePL_cv()`](https://beeamerino.github.io/PhyloCactus/reference/run_treePL_cv.md),
  which implements the protocol of Maurin (2020) in R. The shell script
  this argument used to point at carried no licence and wrote
  `smoothing = `, a keyword `treePL` discards without a message, so
  every chronogram produced before that date was dated at the built-in
  default of 10.

## Value

Invisible NULL upon completion.

## References

Sanderson, M. J. (2002). Estimating absolute rates of molecular
evolution and divergence times: a penalized likelihood approach.
*Molecular Biology and Evolution*, 19(1), 101-109.
[doi:10.1093/oxfordjournals.molbev.a003974](https://doi.org/10.1093/oxfordjournals.molbev.a003974)

Smith, S. A., & O’Meara, B. C. (2012). treePL: divergence time
estimation using penalized likelihood for large phylogenies.
*Bioinformatics*, 28(20), 2689-2690.
[doi:10.1093/bioinformatics/bts492](https://doi.org/10.1093/bioinformatics/bts492)

Maurin, K. J. (2020). An empirical guide for producing a dated phylogeny
with treePL in a maximum likelihood framework. *arXiv preprint
arXiv:2008.07054*.
[doi:10.48550/arXiv.2008.07054](https://doi.org/10.48550/arXiv.2008.07054)

## Examples

``` r
if (FALSE) { # \dontrun{
automate_treePL(
  cfg_file = "calibrations.cfg",
  ml_tree_file = "bestTree.tree",
  bs_trees_file = "temporal_bootstraps.tree",
  results_dir = "auto_results",
  treePL_out = "8_Dating",
  num_bs = 100
)
} # }
```
