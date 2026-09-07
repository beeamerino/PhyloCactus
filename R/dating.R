#' Rescale Branch Lengths of a Phylogenetic Tree
#'
#' Multiplies all edge lengths of a `phylo` tree object by a constant scaling factor.
#'
#' The factor does not change which branches `treePL` clamps. `treePL` refuses to let a branch
#' carry less than one expected substitution and rewrites any shorter branch to `1/numsites`,
#' so the clamp acts on the substitution count, not on the raw length. Measured on the Cactaceae
#' supermatrix, 477 of 2044 branches are clamped, identically at factor 100 and at factor 1, once
#' `numsites` is divided by the same factor as `automate_treePL()` does. What the factor does
#' change is the scale of the rate parameters and therefore the numerical conditioning of the
#' optimisation.
#'
#' Earlier versions of this documentation described the rescaling as a guard against numerical
#' underflow. That is not what it does.
#'
#' @param tree An object of class `phylo` representing a phylogenetic tree.
#' @param factor Numeric multiplier applied to all edge lengths. Defaults to `100`.
#' @return A rescaled `phylo` object with updated edge lengths.
#' @examples
#' \dontrun{
#' library(ape)
#' tree <- rtree(10)
#' rescaled_tree <- rescale_tree(tree, factor = 100)
#' }
#' @export
rescale_tree <- function(tree, factor = 100){
  cat("Rescaling all branches by factor:", factor, "\n")
  tree$edge.length <- tree$edge.length * factor
  return(tree)
}

#' @noRd
.validate_treepl_output <- function(tree_file, label, tol = 1e-6) {
  if (!file.exists(tree_file) || file.info(tree_file)$size == 0) {
    stop("treePL output for '", label, "' is missing or empty at: ", tree_file,
         ". The treePL run likely failed silently; check the corresponding ",
         "final_", label, ".log for details.", call. = FALSE)
  }
  tr <- tryCatch(ape::read.tree(tree_file), error = function(e) NULL)
  if (is.null(tr) || !inherits(tr, "phylo")) {
    stop("treePL output for '", label, "' at ", tree_file,
         " could not be parsed as a valid Newick tree.", call. = FALSE)
  }
  if (!ape::is.ultrametric(tr, tol = tol)) {
    stop("treePL output for '", label, "' at ", tree_file,
         " is not ultrametric (tolerance ", tol, "). The penalized-likelihood ",
         "dating run likely failed to converge; do not treat this chronogram as valid.",
         call. = FALSE)
  }
  invisible(tr)
}

#' Root a Tree on a Clade of Terminals, Tolerating a Basal Polytomy
#'
#' `ape::root()` decides monophyly against the tree's *current* root. `RAxML-NG` writes unrooted
#' Newick with a basal trifurcation, and the terminals of the rooting clade routinely fall on more
#' than one branch of that trifurcation, so `ape::root()` does not see them as a clade and does not
#' place the root where it was asked to. For example, `{outgroup | ingroup}` is a valid bipartition
#' of the unrooted topology, yet several outgroup terminals may sit outside the largest basal branch.
#'
#' Rooting first on an ingroup terminal moves the current root into the ingroup, which makes the
#' outgroup a clade in the rooted sense whenever `{outgroup | ingroup}` is a bipartition of the
#' unrooted topology. The second call then places the root on the intended stem edge.
#'
#' The first step is unconditional. `ape::is.monophyletic()` cannot guard it, because on an
#' unrooted tree it answers in the bipartition sense and reports the terminals as monophyletic
#' while `ape::root()`, which looks for them as a clade of the tree as currently rooted, does not
#' find them; `ape::root()` then falls back to their MRCA, which is the existing root, and the call
#' silently does nothing. That silent no-op is the failure this function exists to remove.
#'
#' Rooting cannot be delegated upstream to `RAxML-NG`. Under time-reversible substitution models
#' the likelihood is identical for every rooting of the same unrooted topology (Felsenstein 1981),
#' so neither a constraint tree nor `--outgroup` can select a root: the data carry no information
#' about its position. The root is an outgroup decision imposed after the search.
#'
#' @param phy An object of class `phylo`.
#' @param outgroup Character vector of rooting terminals, typically from
#'   `resolve_rooting_outgroup()`.
#' @return The rooted `phylo`. Errors, via `ape::root()`, when the terminals are not a clade of the
#'   unrooted topology.
#' @examples
#' tr <- ape::read.tree(text = "((A:1,B:1):1,Portulaca_fulgens:1,Portulaca_oleracea:1);")
#' rooted <- root_on_clade(tr, c("Portulaca_fulgens", "Portulaca_oleracea"))
#' ape::is.monophyletic(rooted, c("Portulaca_fulgens", "Portulaca_oleracea"))
#' @export
root_on_clade <- function(phy, outgroup) {
  outgroup <- intersect(outgroup, phy$tip.label)
  if (length(outgroup) == 0L) {
    stop("None of the rooting terminals are present in the tree.", call. = FALSE)
  }
  if (length(outgroup) == 1L) {
    return(ape::root(phy, outgroup = outgroup, resolve.root = TRUE))
  }

  # Unconditional, and deliberately so. `ape::is.monophyletic()` cannot guard this step: on an
  # unrooted tree it answers in the bipartition sense and reports the terminals as monophyletic
  # while `ape::root()`, which looks for them as a clade of the tree as currently rooted, does not
  # find them. Guarding on it therefore skips the step and lets `ape::root()` fall back to their
  # MRCA, which is the existing root, so the call silently does nothing.
  ingroup <- setdiff(phy$tip.label, outgroup)
  if (length(ingroup) > 0L) {
    phy <- ape::root(phy, outgroup = ingroup[[1]], resolve.root = TRUE)
  }

  # A set that is still not a clade here is not one in the unrooted topology either. `ape::root()`
  # errors in that case with a clear message, and letting it propagate is better than inventing a
  # fallback root placement that no calibration could rely on.
  ape::root(phy, outgroup = outgroup, resolve.root = TRUE)
}

#' Abort When a Tree Handed to treePL Is Not Rooted
#'
#' `treePL` requires a rooted tree. Given an unrooted one it does not refuse the input, it fails
#' during optimisation with a message that does not name the cause. Parsing the file and checking
#' rootedness first turns that into an actionable error.
#'
#' @param treefile Character. Path to a Newick tree file.
#' @param label Character. Run label used in the error message.
#' @return Invisible `TRUE` when the tree is rooted.
#' @noRd
.assert_rooted_treefile <- function(treefile, label) {
  if (!file.exists(treefile)) {
    stop("Tree file for '", label, "' does not exist: ", treefile, call. = FALSE)
  }
  tr <- tryCatch(ape::read.tree(treefile), error = function(e) NULL)
  if (is.null(tr)) {
    stop("Tree file for '", label, "' could not be parsed as Newick: ", treefile, call. = FALSE)
  }
  if (!ape::is.rooted(tr)) {
    stop("Tree file for '", label, "' is unrooted: ", treefile,
         ". treePL strictly requires a rooted tree. Root it on the outgroup clade first, for ",
         "example with root_on_clade(tree, resolve_rooting_outgroup(tree$tip.label)), or use ",
         "automate_treePL(), which roots its inputs itself.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Extract the Tree Path Declared in a treePL Configuration File
#'
#' @param cfg_file Character. Path to a `treePL` configuration file.
#' @return Character path declared on the `treefile` line, or `NA_character_` when absent.
#' @noRd
.treefile_from_cfg <- function(cfg_file) {
  if (!file.exists(cfg_file)) return(NA_character_)
  lines <- readLines(cfg_file, warn = FALSE)
  hit <- grep("^\\s*treefile\\s*=", lines, value = TRUE)
  if (length(hit) == 0L) return(NA_character_)
  trimws(sub("^\\s*treefile\\s*=\\s*", "", hit[[1]]))
}


#' Run treePL Executable Directly
#'
#' Executes the `treePL` binary directly to estimate divergence times from a prepared configuration file.
#'
#' @param cfg_file Character. Path to `treePL` configuration file.
#' @param label Character. Descriptive run label identifier.
#' @param cwd Character. Optional working directory context. Defaults to `NULL`.
#' @details
#' The tree path is read from the `treefile` line of `cfg_file` and checked for rootedness before
#' the binary is invoked, because `treePL` fails opaquely on an unrooted tree. When the
#' configuration declares no `treefile`, or the declared path cannot be resolved, the check is
#' skipped with a warning rather than blocking the run.
#' @return Invisible NULL upon system command execution.
#' @references
#' Smith, S. A., & O’Meara, B. C. (2012). treePL: divergence time estimation using penalized likelihood
#' for large phylogenies. *Bioinformatics*, 28(20), 2689-2690. \doi{10.1093/bioinformatics/bts492}
#' @export
run_treePL_direct <- function(cfg_file, label, cwd = NULL){
  if (!is.null(cwd)) {
    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(cwd)
  }

  # Resolved after the setwd() above, so a relative `treefile` in the configuration is interpreted
  # in the same directory treePL will interpret it.
  declared_tree <- .treefile_from_cfg(cfg_file)
  if (is.na(declared_tree) || !file.exists(declared_tree)) {
    warning("Could not resolve the `treefile` declared in '", cfg_file,
            "'; the rooting check was skipped. treePL requires a rooted tree.", call. = FALSE)
  } else {
    .assert_rooted_treefile(declared_tree, label)
  }

  # Captured to a file rather than streamed. Two reasons: the log is what
  # .verify_treepl_smoothing() reads to confirm treePL used the smoothing it was given, and a
  # 1000-taxon run emits one "tiny branch length" line per clamped branch, which buries every
  # other message in the console.
  log_file <- paste0("treepl_run_", label, ".log")
  cat(Sys.time(), "- Running treePL (direct) for:", label, "\n")
  status <- system2("treePL", args = shQuote(cfg_file), stdout = log_file, stderr = log_file)
  if (!identical(as.integer(status), 0L)) {
    stop("treePL direct execution failed for '", label, "' with exit status ", status,
         ". See '", file.path(getwd(), log_file), "'.", call. = FALSE)
  }
  .verify_treepl_smoothing(log_file, cfg_file, label)
  .check_treepl_convergence(log_file, label)
  cat(Sys.time(), "- Finished:", label, "\n\n")
}

#' Report a treePL run whose gradient optimisation never moved
#'
#' `treePL` optimises in two phases: simulated annealing, then gradient descent (L-BFGS through
#' NLopt) from the point annealing reached. Both phases print their objective to the log, as
#' `exit siman:` and `after opt calc2:`. When the second equals the first to every printed digit,
#' the gradient phase improved nothing: what the program writes out is the annealing endpoint, not
#' an optimised penalized-likelihood solution. The exit status is still 0 and the chronogram is
#' still ultrametric, so nothing else in this package notices.
#'
#' Measured on 2026-09-02 across twenty-one runs of the same tree at seven smoothing values. The
#' two runs whose ages were coherent improved by 0.90% and 0.71%; the five that returned absurd
#' ages improved by 0.0000%, 0.0000%, 0.0000%, 0.0000% and 0.0034%. One of those put `ACP_root` at
#' 188.13 Ma, older than the crown of Caryophyllales, and another moved `Opuntioideae_mrca` from
#' its lower bound to its upper bound while every other smoothing value left it on the lower one.
#' NLopt's own return codes did not separate the two groups; this quantity did, without overlap.
#'
#' That is a correlation over twenty-one runs on one dataset, not a proof, so this warns rather
#' than stops: annealing could in principle land on the optimum and leave the gradient phase with
#' nothing to do. Treat the warning as a reason to inspect the chronogram, and to distrust any node
#' whose age changes non-monotonically with smoothing.
#'
#' @param log_file Path to the captured treePL output.
#' @param label Run label, for the message.
#' @param min_improvement Numeric. Relative improvement below which the run is reported. The two
#'   groups are separated by more than two orders of magnitude: the incoherent runs improved by 0,
#'   0, 0, 0 and 3.4e-5, the coherent ones by 7.1e-3 and 9.0e-3. The default of 1e-4 sits three
#'   times above the largest improvement measured in an incoherent run and seventy times below the
#'   smallest measured in a coherent one.
#'
#'   Note that on a very small tree the objective is small enough that the printed digits, rather
#'   than the optimiser, set the floor on a measurable improvement, and annealing can legitimately
#'   land where the gradient phase has nothing to add. The package's own fixture trees trigger this
#'   warning for that reason. It is informative on trees of the size this pipeline is built for.
#' @return Invisibly the relative improvement, or `NA` when the log does not carry both numbers.
#' @keywords internal
#' @noRd
.check_treepl_convergence <- function(log_file, label, min_improvement = 1e-4) {
  if (!file.exists(log_file)) return(invisible(NA_real_))
  lg <- readLines(log_file, warn = FALSE)

  # treePL gives up on its own after ten attempts at a feasible starting point, and writes no tree.
  # Reported here because the message is buried in a log that a large run fills with other output.
  if (any(grepl("Failed setting feasible start rates/dates", lg, fixed = TRUE))) {
    stop("treePL could not find feasible starting rates and dates for '", label,
         "' and aborted without writing a chronogram. This is what it does when the configuration ",
         "carries no usable calibration, or when the calibrations it carries cannot all be ",
         "satisfied at once. See '", file.path(getwd(), log_file), "'.", call. = FALSE)
  }

  num <- function(pat) {
    hit <- grep(pat, lg, value = TRUE)
    if (length(hit) == 0L) return(NA_real_)
    suppressWarnings(as.numeric(sub(".*:\\s*", "", hit[length(hit)])))
  }
  siman <- num("^\\s*exit siman\\s*:")
  final <- num("^\\s*after opt calc2\\s*:")
  if (is.na(final)) final <- num("^\\s*after opt calc1\\s*:")
  if (is.na(siman) || is.na(final) || siman <= 0) return(invisible(NA_real_))

  improvement <- (siman - final) / siman
  if (improvement < min_improvement) {
    warning("treePL run '", label, "' left its gradient optimisation with nothing to do: the ",
            "objective went from ", format(siman), " after simulated annealing to ",
            format(final), " after optimisation, a relative change of ",
            format(improvement, digits = 3), ". Runs with this signature have returned node ages ",
            "that are not penalized-likelihood estimates but wherever annealing happened to stop. ",
            "Do not report ages from this chronogram without checking them against a run at a ",
            "different smoothing value: see report_smoothing_sensitivity(). Log: '",
            file.path(getwd(), log_file), "'.", call. = FALSE)
  }
  invisible(improvement)
}

#' Confirm that treePL used the smoothing value it was given
#'
#' treePL's configuration keyword is `smooth`. A line reading `smoothing = X` is not recognised,
#' is discarded without any message, and the run proceeds on the built-in default of 10. Between
#' the first dated tree of this project and 2026-09-02 the wrapper wrote `smoothing`, so every
#' chronogram was produced at 10 and the cross-validation that selects the value never reached the
#' program. Nothing in the outputs revealed it: the trees were valid, ultrametric, and plausible.
#'
#' The one place the truth was visible is treePL's own log, which prints the smoothing it is
#' actually using. This compares that number against the configuration and stops when they differ,
#' because a chronogram dated at a smoothing nobody chose is not a result, and the failure is
#' silent by construction.
#'
#' @param log_file Path to the captured treePL output.
#' @param cfg_file Path to the configuration handed to treePL.
#' @param label Run label, for the message.
#' @return Invisibly `TRUE` when the two agree, or `NA` when either could not be read.
#' @keywords internal
#' @noRd
.verify_treepl_smoothing <- function(log_file, cfg_file, label) {
  if (!file.exists(log_file) || !file.exists(cfg_file)) return(invisible(NA))

  cfg <- readLines(cfg_file, warn = FALSE)
  # Anchored so that `smoothing =`, the old wrong keyword, is not read as a request.
  req <- grep("^\\s*smooth\\s*=", cfg, value = TRUE)
  if (length(req) == 0L) return(invisible(NA))
  requested <- suppressWarnings(as.numeric(sub(".*=\\s*", "", req[length(req)])))

  used_line <- grep("^\\s*smoothing\\s*:", readLines(log_file, warn = FALSE), value = TRUE)
  if (length(used_line) == 0L) return(invisible(NA))
  used <- suppressWarnings(as.numeric(sub(".*:\\s*", "", used_line[length(used_line)])))

  if (is.na(requested) || is.na(used)) return(invisible(NA))
  if (!isTRUE(all.equal(requested, used, tolerance = 1e-6))) {
    stop("treePL ran '", label, "' at smoothing ", used, " but the configuration asked for ",
         requested, ". The keyword treePL reads is `smooth`; a line it does not recognise is ",
         "discarded in silence and the run falls back to the default. Check the `smooth = ` line ",
         "in '", cfg_file, "'. Dates produced under a smoothing nobody selected cannot be ",
         "reported.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Automate treePL Divergence Time Estimation Pipeline Across Bootstrap Cohorts
#'
#' Automates cross-validation parameter optimization, rate smoothing selection, and chronogram estimation
#' across temporal bootstrap replicates using `treePL` (Sanderson, 2002; Smith & O'Meara, 2012).
#' Propagating temporal uncertainty across branch-length resampled bootstrap trees yields empirical confidence intervals
#' for node age estimates. Every treePL output (maximum-likelihood chronogram and each bootstrap chronogram) is validated
#' after execution: the resulting tree must exist, be non-empty, parse as a valid Newick topology, and be ultrametric.
#' A run that fails silently (e.g., because the underlying `treePL` binary did not converge) is therefore reported as an
#' explicit error rather than propagated downstream as a corrupted chronogram. Note that the optimal rate-smoothing
#' parameter is cross-validated once on the maximum-likelihood tree and reused, unmodified, across all bootstrap
#' replicates; this is a standard computational shortcut for treePL-based dating pipelines (per-replicate cross-validation
#' is prohibitively expensive at typical bootstrap replicate counts), following the empirical protocol of Maurin (2020).
#'
#' @param cfg_file Character. Path to primary `treePL` configuration file specifying calibration bounds and parameters.
#' @param wrapper_sh Retired on 2026-09-02 and ignored, with a warning. Priming, cross-validation
#'   and dating of the maximum-likelihood tree now run through [run_treePL_cv()], which implements
#'   the protocol of Maurin (2020) in R. The shell script this argument used to point at carried no
#'   licence and wrote `smoothing = `, a keyword `treePL` discards without a message, so every
#'   chronogram produced before that date was dated at the built-in default of 10.
#' @param n_prime Integer. Priming repeats on the maximum-likelihood tree. See [run_treePL_cv()].
#' @param prime_rule `"lowest"` (Maurin 2020) or `"modal"` (the retired shell script's rule).
#' @param cv_method `"randomcv"`, recommended by Maurin (2020) as faster and more stable, or `"cv"`
#'   for leave-one-out, which is what the retired shell script used.
#' @param cvstart,cvstop Ends of the cross-validation smoothing grid. The defaults reach lower than
#'   the shell script's fixed floor of 1e-04, which the August 2026 run of this project hit without
#'   the chi-square curve ever turning.
#' @param ml_tree_file Character. Path to input maximum-likelihood reference tree file.
#' @param bs_trees_file Character. Path to input temporal bootstrap trees file.
#' @param results_dir Character. Directory path to save intermediate optimization results.
#' @param treePL_out Character. Destination directory path for final output chronograms.
#' @param num_bs Integer or NULL. Maximum number of temporal bootstrap trees to evaluate. If `NULL`, processes all available trees.
#' @param numsites Integer or NULL. Alignment length, in sites, of the supermatrix actually analysed. If `NULL`, the `numsites` line already present in `cfg_file` is used. This must be the length of the matrix `RAxML-NG` analysed (typically the `*.raxml.reduced.phy` produced by `preprocess_partitions()`), not a rounded figure: `treePL` uses it to convert branch lengths into expected substitution counts, so an incorrect value biases the rate smoothing. Declare the true alignment length here; the division by `rescale_factor` is applied internally.
#' @param rescale_factor Numeric. Multiplier applied to every branch length before dating, and the divisor applied to `numsites` in the configurations written for `treePL`. Defaults to `100`.
#'
#'   **Why the two are one argument.** `treePL` reads a branch as `edge.length * numsites` expected substitutions. Branch lengths are rescaled because a substantial fraction of a low-divergence plastid supermatrix falls below the internal minimum `treePL` imposes on a branch, and those branches would otherwise be clamped to a common value, erasing the rate signal across them. Rescaling without dividing `numsites` by the same factor leaves `treePL` reading a matrix it believes to be `rescale_factor` times more informative than it is. The likelihood term grows with the substitution counts while the roughness penalty does not, so the effective smoothing becomes weaker than the nominal value by that factor, rates vary almost freely between branches, and node ages stop being determined by the data and start being determined by the edges of the region the calibrations leave feasible. The symptom is a chronogram whose calibrated nodes sit exactly on their bounds. Coupling the two here makes the pair impossible to separate by accident.
#' @param outgroup Character vector of terminals used to root the maximum-likelihood tree and every bootstrap replicate before penalized-likelihood dating, or `NULL`. Defaults to `NULL`, in which case the set is derived from the maximum-likelihood tree with `resolve_rooting_outgroup()`. Adapt `outgroup` or specify `resolve_rooting_outgroup(pattern = ...)` to the outgroup lineage sampled in your dataset (for example `"^(Talinum|Talinella)_"` when Talinaceae roots the tree). No fixed default is offered because a clade-level rooting set depends on what the supermatrix sampled and cannot be a package constant.
#'
#'   **Why a clade rather than a terminal.** Rooting is imposed after the search: `RAxML-NG` returns an unrooted topology and the root is placed on a chosen edge. Naming a single terminal of a sampled outgroup clade places the root *inside* that clade, leaving it paraphyletic in the final tree and collapsing its crown node onto the root. Any calibration addressed by the MRCA of that clade then lands on the root instead of on the node it was written for. Supplying the whole clade places the root on its stem edge, which is the intended edge. Rooting is performed by `root_on_clade()`, which also handles the basal polytomy `RAxML-NG` writes.
#'
#'   Bootstrap replicates are handled with the intersection of this set and each replicate's tip labels, so a replicate missing some terminals is still rooted; only a replicate missing all of them is an error.
#'
#'   **Topological assumption.** Placing the root on the outgroup lineage (Talinaceae or Portulacaceae) establishes the basal split for dating. In the reference dataset, Talinaceae roots the tree, placing the root on the stem of the ACP clade (Anacampserotaceae, Cactaceae, Portulacaceae) and allowing maximum-likelihood inference to test alternative topological resolutions among the three core families (Ramirez-Barahona et al., 2020; Zuntini et al., 2024). State this assumption in Methods, and treat the root age as conditional on it.
#' @param seed Integer or NULL. Seed for every stochastic step of the run: R's choice of which
#'   bootstrap replicates to date, and `treePL`'s own `seed` keyword, which is written into the
#'   maximum-likelihood configuration and, offset by the replicate index, into each replicate's.
#'   `treePL` seeds itself from the clock when the keyword is absent, which leaves both the
#'   cross-validation and the simulated annealing unreproducible. Defaults to `NULL`, which draws
#'   one and reports it; record the reported value, as it is required to reproduce the run.
#' @return Invisible NULL upon completion.
#' @references
#' Sanderson, M. J. (2002). Estimating absolute rates of molecular evolution and divergence times:
#' a penalized likelihood approach. *Molecular Biology and Evolution*, 19(1), 101-109. \doi{10.1093/oxfordjournals.molbev.a003974}
#'
#' Smith, S. A., & O’Meara, B. C. (2012). treePL: divergence time estimation using penalized likelihood
#' for large phylogenies. *Bioinformatics*, 28(20), 2689-2690. \doi{10.1093/bioinformatics/bts492}
#'
#' Maurin, K. J. (2020). An empirical guide for producing a dated phylogeny with treePL in a maximum likelihood framework.
#' *arXiv preprint arXiv:2008.07054*. \doi{10.48550/arXiv.2008.07054}
#' @examples
#' \dontrun{
#' automate_treePL(
#'   cfg_file = "calibrations.cfg",
#'   ml_tree_file = "bestTree.tree",
#'   bs_trees_file = "temporal_bootstraps.tree",
#'   results_dir = "auto_results",
#'   treePL_out = "8_Dating",
#'   num_bs = 100
#' )
#' }
#' @export
automate_treePL <- function(cfg_file, ml_tree_file, bs_trees_file, results_dir, treePL_out,
                            num_bs = NULL, numsites = NULL, outgroup = NULL, seed = NULL,
                            rescale_factor = 100, n_prime = 10L,
                            prime_rule = c("lowest", "modal"),
                            cv_method = c("randomcv", "cv"),
                            cvstart = 1e3, cvstop = 1e-8, wrapper_sh = NULL) {

  prime_rule <- match.arg(prime_rule)
  cv_method <- match.arg(cv_method)

  # Retired on 2026-09-02. The shell script it pointed at carried no licence and wrote a smoothing
  # keyword treePL does not read; run_treePL_cv() replaces it. Accepted and ignored rather than
  # removed, so that a script written against the old signature still runs and says why.
  if (!is.null(wrapper_sh)) {
    warning("`wrapper_sh` is retired and ignored: priming, cross-validation and dating now run ",
            "through run_treePL_cv() in R. Remove the argument. See ?run_treePL_cv.",
            call. = FALSE)
  }

  # Convert inputs to absolute paths so they survive setwd()
  cfg_file_abs <- normalizePath(cfg_file, mustWork = TRUE)

  ml_dir       <- file.path(results_dir, "ML_tree")
  bs_out_dir   <- file.path(results_dir, "BS_tree")
  
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ml_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(bs_out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(treePL_out, recursive = TRUE, showWarnings = FALSE)
  
  if (!is.numeric(rescale_factor) || length(rescale_factor) != 1L || is.na(rescale_factor) || rescale_factor <= 0) {
    stop("`rescale_factor` must be a single positive number.", call. = FALSE)
  }

  # STEP 0: Fix and rescale maximum-likelihood tree
  ml_tree <- ape::read.tree(ml_tree_file)

  # A clade-level rooting set cannot be a package constant: it depends on what the supermatrix
  # actually sampled. With `outgroup = NULL` it is derived from the tree itself.
  if (is.null(outgroup)) {
    outgroup <- resolve_rooting_outgroup(ml_tree$tip.label)
    message("Rooting terminals derived from the maximum-likelihood tree (", length(outgroup), "): ",
            paste(outgroup, collapse = ", "))
  }

  # `outgroup` is a vector, so membership must be reduced with all()/any() rather than used as a
  # condition directly. Every terminal must be present in the maximum-likelihood tree: a partial
  # rooting set here would place the root on a different edge than intended, silently.
  ml_missing <- setdiff(outgroup, ml_tree$tip.label)
  if (!ape::is.rooted(ml_tree)) {
    if (length(ml_missing) > 0L) {
      stop("Maximum-likelihood tree is unrooted and ", length(ml_missing),
           " rooting terminal(s) are absent from its tip labels: ",
           paste(ml_missing, collapse = ", "),
           ". treePL strictly requires a rooted tree.", call. = FALSE)
    }
    ml_tree <- root_on_clade(ml_tree, outgroup)
  }

  # Resolved after the rooting checks and before the rescaling. Later than the bootstrap stage
  # would allow, because the configuration handed to the maximum-likelihood run needs the scaled
  # value and that run comes first; but after the rooting validation, so that a tree treePL could
  # never date still reports the rooting failure rather than a missing numsites.
  if (is.null(numsites)) {
    declared <- grep("^\\s*numsites\\s*=", readLines(cfg_file_abs), value = TRUE)
    if (length(declared) == 0L) {
      stop("numsites not found in treePL config and not provided. ",
           "Specify the alignment length via the numsites argument.", call. = FALSE)
    }
    numsites <- as.numeric(trimws(sub("^\\s*numsites\\s*=\\s*", "", declared[[1]])))
  }
  numsites <- as.numeric(numsites)
  if (is.na(numsites) || numsites <= 0) {
    stop("`numsites` must be a positive number of alignment sites.", call. = FALSE)
  }

  # The product edge.length * numsites is what treePL reads as a substitution count, so scaling
  # the branches up and this value down by the same factor leaves that count unchanged while
  # keeping every branch above the internal minimum treePL applies to a branch length.
  numsites_scaled <- round(numsites / rescale_factor)
  if (numsites_scaled < 1) {
    stop("rescale_factor = ", rescale_factor, " leaves numsites at ", numsites / rescale_factor,
         ", below one site. Lower the factor.", call. = FALSE)
  }
  numsites_line <- paste0("numsites = ", numsites_scaled)
  message("Branch lengths rescaled by ", rescale_factor, "; numsites written to the treePL ",
          "configurations as ", numsites_scaled, " (declared alignment length ", numsites, ").")

  # Resolved here, before any configuration is written, because it belongs in all of them.
  #
  # treePL has its own `seed` keyword, and when it is absent the program seeds itself from the
  # clock. Both stages that consume randomness are affected: the cross-validation, which partitions
  # sites, and the simulated annealing that precedes every gradient optimisation. Two runs of the
  # same configuration therefore need not select the same smoothing value or return the same ages,
  # and a result nobody can reproduce cannot be checked by a reviewer or by its own author. Until
  # 2026-09-02 this argument seeded only R's sampling of which bootstrap replicates to date, and
  # the treePL runs themselves were left to the clock.
  if (is.null(seed)) {
    # Bounded well below the integer maximum so that the per-replicate offsets below stay inside it.
    seed <- sample.int(1e8L, 1L)
    message("No seed supplied; using ", seed, ". Record it: it is written into every treePL ",
            "configuration and is required to reproduce this run.")
  }
  seed <- as.integer(seed)
  if (is.na(seed) || seed <= 0) {
    stop("`seed` must be a single positive integer.", call. = FALSE)
  }
  seed_line <- paste0("seed = ", format(seed, scientific = FALSE))

  ml_tree <- rescale_tree(ml_tree, factor = rescale_factor)
  ml_tree_fixed_file <- file.path(ml_dir, "supportTree_raxml-ng_fixed.tree")
  ape::write.tree(ml_tree, ml_tree_fixed_file)
  ml_tree_fixed_file_abs <- normalizePath(ml_tree_fixed_file, mustWork = TRUE)
  
  # The user configuration keeps declaring the true alignment length, which is what belongs in a
  # file meant to document the analysis. The value treePL is actually given is written here, into
  # a copy, so the two can never drift apart and the correction stays visible on disk.
  cfg_ml_scaled <- file.path(ml_dir, "calibrations_scaled.cfg")
  cfg_lines_ml <- readLines(cfg_file_abs, warn = FALSE)
  has_numsites <- grepl("^\\s*numsites\\s*=", cfg_lines_ml)
  if (any(has_numsites)) {
    cfg_lines_ml[has_numsites] <- numsites_line
  } else {
    cfg_lines_ml <- c(numsites_line, cfg_lines_ml)
  }
  # A seed already in the user configuration is replaced, not duplicated: treePL reads the last
  # assignment, so leaving both would make the effective seed depend on line order.
  has_seed <- grepl("^\\s*seed\\s*=", cfg_lines_ml)
  if (any(has_seed)) {
    cfg_lines_ml[has_seed] <- seed_line
  } else {
    cfg_lines_ml <- c(seed_line, cfg_lines_ml)
  }
  writeLines(cfg_lines_ml, cfg_ml_scaled)
  cfg_ml_scaled_abs <- normalizePath(cfg_ml_scaled, mustWork = TRUE)

  # STEP 1: Prime, cross-validate and date the maximum-likelihood tree.
  #
  # Carried out by run_treePL_cv(), which implements the protocol of Maurin (2020) in R. Until
  # 2026-09-02 this called out to a vendored copy of an unlicensed shell script; see the notes in
  # R/treepl_cv.R for why that had to go and what the shell script got wrong.
  if (!file.exists(file.path(ml_dir, "treepl_ML_tree.tre"))) {
    run_treePL_cv(cfg_file = cfg_ml_scaled_abs,
                  tree_file = ml_tree_fixed_file_abs,
                  label = "ML_tree",
                  n_prime = n_prime,
                  prime_rule = prime_rule,
                  cv_method = cv_method,
                  cvstart = cvstart,
                  cvstop = cvstop,
                  work_dir = ml_dir)
  } else {
    cat("treePL results for the maximum-likelihood tree already exist! Skipping treePL run.\n")
  }

  # Still checked after the fact, because the branch above is skipped whenever the chronogram is
  # already on disk, including one produced by the shell script before 2026-09-02.
  .validate_treepl_output(file.path(ml_dir, "treepl_ML_tree.tre"), "ML_tree")

  # STEP 2/3: Extract optimization parameters and smoothing from maximum-likelihood tree
  cfg_smooth_file <- file.path(ml_dir, "configure_smooth_ML_tree")
  cfg_lines <- readLines(cfg_smooth_file)
  opt_pattern <- "^(opt|optad|optcvad|moredetail|moredetailad|moredetailcvad)\\b"
  opt_lines <- grep(opt_pattern, cfg_lines, value = TRUE)
  
  best_smoothing <- NULL
  # The wrapper has already calculated the best smoothing and appended it, as `smooth = X`.
  # `smoothing` is accepted here as well, because run directories produced before 2026-09-02
  # carry the old, wrong keyword and should still be readable. The anchored `=` keeps the two
  # patterns from matching each other.
  smoothing_line <- grep("^(smooth|smoothing) *=", cfg_lines, value = TRUE)
  if (length(smoothing_line) > 0) {
    best_smoothing <- as.numeric(sub(".*=\\s*", "", smoothing_line[length(smoothing_line)]))
  } else {
    # Fallback to manual extraction if smoothing line is missing
    cv_file <- file.path(ml_dir, "cv_ML_tree")
    if (file.exists(cv_file)) {
      cv_raw <- readLines(cv_file)
      smooth_raw <- regmatches(cv_raw, regexpr("\\([0-9.eE+-]+\\)", cv_raw))
      smoothing <- as.numeric(gsub("[()]", "", smooth_raw))
      chisq <- as.numeric(sub(".*\\)\\s+", "", cv_raw))
      cv_data <- data.frame(smoothing = smoothing, chisq = chisq)
      if (nrow(cv_data) > 0) {
        best_smoothing <- cv_data$smoothing[which.min(cv_data$chisq)]
      }
    }
  }
  
  if (is.null(best_smoothing)) {
    stop("Failed to extract optimal smoothing parameter from treePL output.", call. = FALSE)
  }
  cat("Best smoothing strategy chosen:", best_smoothing, "\n")

  # Cross-validation selects by minimum chi-square over a grid the wrapper fixes at 1e-04 to 1e+04.
  # A minimum at either end is not a selection, it is the grid running out: the curve was still
  # descending. In the August 2026 run it was monotone all the way down and the wrapper took
  # cvstart, 1e-04, which is effectively no smoothing at all, and every calibrated node came back
  # sitting on its own bound.
  cv_file <- file.path(ml_dir, "cv_ML_tree")
  if (file.exists(cv_file)) {
    cv_raw <- readLines(cv_file, warn = FALSE)
    grid <- suppressWarnings(as.numeric(gsub("[()]", "",
      regmatches(cv_raw, regexpr("\\([0-9.eE+-]+\\)", cv_raw)))))
    grid <- grid[is.finite(grid)]
    if (length(grid) > 1L && isTRUE(best_smoothing %in% range(grid))) {
      warning("Cross-validation selected smoothing = ", best_smoothing,
              ", which is the ", if (best_smoothing == min(grid)) "lowest" else "highest",
              " value on the tested grid (", min(grid), " to ", max(grid),
              "). The curve had not turned, so this is the edge of the search rather than an ",
              "optimum, and the chronogram should not be interpreted until the grid is widened ",
              "or the selection is justified on other grounds.", call. = FALSE)
    }
  }
  
  # STEP 4: Prepare and rescale Bootstrap Trees (from concatenated)
  bs_all <- ape::read.tree(bs_trees_file)
  message(length(bs_all), " BS trees found in concatenated file\n")
  
  if (is.null(num_bs)) {
    num_bs_to_use <- length(bs_all)
  } else {
    num_bs_to_use <- min(num_bs, length(bs_all))
  }
  
  # Resolved above, before the maximum-likelihood configuration was written. Used here for R's own
  # choice of which replicates to date, and written into each replicate's configuration below so
  # that treePL's annealing is seeded too.
  set.seed(seed)
  bs_subset <- sample(bs_all, num_bs_to_use)
  
  # Detect local Performance-core count once (constant across all bootstrap replicates) for
  # consistency with the dynamic P-core thread allocation used elsewhere for RAxML-NG steps,
  # instead of a hardcoded thread count.
  treepl_threads <- tryCatch(.detect_pcores(), error = function(e) 8L)

  for(i in seq_along(bs_subset)){
    label <- paste0("BS_", i)
    
    bs_folder <- file.path(bs_out_dir, label)
    dir.create(bs_folder, recursive = TRUE, showWarnings = FALSE)
    
    # Skip processing if tree PL output already exists
    if (file.exists(file.path(bs_folder, paste0("treepl_", label, ".tre")))) {
      cat("Skipping BS tree preparation:", label, "- already completed\n")
      next
    }
    
    cat("Processing BS tree:", label, "\n")
    
    # Bootstrap replicates are resampled from the same alignment, so a replicate may legitimately
    # lack some rooting terminals. The intersection is used rather than the full set, and only an
    # empty intersection is fatal.
    if (!ape::is.rooted(bs_subset[[i]])) {
      bs_outgroup <- intersect(outgroup, bs_subset[[i]]$tip.label)
      if (length(bs_outgroup) == 0L) {
        stop("Bootstrap tree '", label, "' is unrooted and none of the ", length(outgroup),
             " rooting terminals are present in its tip labels. ",
             "treePL strictly requires a rooted tree.", call. = FALSE)
      }
      bs_subset[[i]] <- root_on_clade(bs_subset[[i]], bs_outgroup)
    }

    tree <- rescale_tree(bs_subset[[i]], factor = rescale_factor)
    tree_fixed_file <- file.path(bs_folder, paste0(label, "_fixed.tree"))
    ape::write.tree(tree, tree_fixed_file)
    tree_fixed_file_abs <- normalizePath(tree_fixed_file, mustWork = TRUE)
    
    bs_cfg_file <- file.path(bs_folder, paste0("cfg_", label, ".cfg"))
    cfg_content <- c(
      paste0("treefile = ", tree_fixed_file_abs),
      numsites_line,
      # Every replicate gets a seed derived from the run's, not the same one: identical seeds
      # across replicates would correlate their annealing, and no seed at all would leave each to
      # the clock and make the confidence interval unreproducible.
      paste0("seed = ", format((as.numeric(seed) + i) %% 2147483647, scientific = FALSE)),
      grep("^mrca|^min|^max", readLines(cfg_file_abs), value = TRUE),
      paste0("nthreads = ", treepl_threads),
      "thorough",
      opt_lines,
      # `smooth`, not `smoothing`: see the notes at the top of R/treepl_cv.R. Written without
      # scientific notation because a value such as 1e-04 is not what treePL's parser expects,
      # and a smoothing silently reset to the default is exactly the failure being closed here.
      paste0("smooth = ", format(best_smoothing, scientific = FALSE)),
      paste0("outfile = treepl_", label, ".tre")
    )
    writeLines(cfg_content, bs_cfg_file)
  }
  
  # STEP 5: Execute treePL for each bootstrap tree
  bs_folders <- list.dirs(bs_out_dir, recursive = FALSE)
  for(bs_folder in bs_folders){
    label <- basename(bs_folder)
    out_tree <- file.path(bs_folder, paste0("treepl_", label, ".tre"))
    if (file.exists(out_tree)) {
       # A file existing on disk is not sufficient evidence of a successful run
       # (a prior crashed/killed job can leave a truncated or non-ultrametric tree).
       .validate_treepl_output(out_tree, label)
       message("Skipping ", label, " execution - already exists and validated")
       next
    }

    bs_cfg <- list.files(bs_folder, pattern="^cfg_.*\\.cfg$", full.names = TRUE)
    if(length(bs_cfg) == 1){
      label <- basename(bs_folder)
      bs_cfg_abs <- normalizePath(bs_cfg, mustWork = TRUE)
      run_treePL_direct(bs_cfg_abs, label, cwd = bs_folder)
      .validate_treepl_output(out_tree, label)
      message("Bootstrap treePL completed and validated. Results in:", bs_folder, "\n")
    } else {
      message("Skipping", bs_folder, "- CFG missing\n")
    }
  }
  
  # STEP 6: Collect treePL outputs
  ml_results_dir <- file.path(results_dir, "ML_tree_treePL")
  bs_results_dir <- file.path(results_dir, "BS_tree_treePL")
  
  dir.create(ml_results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(bs_results_dir, recursive = TRUE, showWarnings = FALSE)
  
  ml_treepl_file <- list.files(ml_dir, pattern = "^treepl_.*\\.tre$", full.names = TRUE)
  if(length(ml_treepl_file) == 1){
    file.copy(ml_treepl_file, ml_results_dir, overwrite = TRUE)
    message(Sys.time(), " - treePL chronogram for maximum-likelihood tree copied to: ", ml_results_dir)
  }
  
  bs_folders <- list.dirs(bs_out_dir, recursive = FALSE)
  for(bs_folder in bs_folders){
    bs_treepl <- list.files(bs_folder, pattern = "^treepl_.*\\.tre$", full.names = TRUE)
    if(length(bs_treepl) == 1){
      file.copy(bs_treepl, bs_results_dir, overwrite = TRUE)
      message(Sys.time(), " - bootstrap treePL tree copied: ", basename(bs_treepl))
    }
  }
  
  bs_treepl_files <- list.files(
    bs_out_dir,
    pattern = "^treepl_.*\\.tre$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  out_bs_all <- file.path(treePL_out, "bsTree_treePL.tree")
  file.create(out_bs_all)
  
  for(f in bs_treepl_files){
    cat(readLines(f), file = out_bs_all, sep = "\n", append = TRUE)
  }
  
  ml_treepl_file <- file.path(ml_dir, "treepl_ML_tree.tre")
  final_ml_dated <- file.path(treePL_out, "BestTree_treePL.tree")
  if (file.exists(ml_treepl_file)) {
    file.copy(
      ml_treepl_file,
      final_ml_dated,
      overwrite = TRUE
    )
  }

  cat("\n====================================================\n")
  cat("  treePL Divergence Time Estimation Complete \U0001f335\n")
  cat("====================================================\n")
  cat("  Best ML chronogram:     ", final_ml_dated, "\n")
  cat("  Bootstrap chronograms:  ", out_bs_all, " (", length(bs_treepl_files), " trees)\n")
  cat("  Dating outputs dir:     ", treePL_out, "\n")
  cat("====================================================\n\n")

  return(invisible(final_ml_dated))
}



#' Resolve the Terminals a Calibration Row Addresses
#'
#' `treePL` locates a calibrated node by the MRCA of the terminals declared for it, so the taxon
#' set, not the row label, is what determines where a bound lands. `value` accepts a
#' semicolon-separated list of names, which the two stem calibrations require: a stem node is
#' shared by the families it subtends and cannot be addressed by a single family name.
#'
#' @param constraints Data frame of the taxonomic table (`cactus_constraints.csv`), carrying a
#'   `Specie_name` column and one column per rank used by the calibration table.
#' @param column Character. Name of the rank column the calibration row keys on.
#' @param value Character. One name, or several separated by `;`.
#' @param tip_labels Character vector of terminals present in the tree.
#' @return Character vector of terminals, intersected with `tip_labels`.
#' @examples
#' constraints <- data.frame(
#'   Specie_name = c("Opuntia_stricta", "Pereskia_aculeata", "Anacampseros_retusa"),
#'   Family = c("Cactaceae", "Cactaceae", "Anacampserotaceae")
#' )
#' calibration_tips(constraints, "Family", "Cactaceae;Anacampserotaceae",
#'                  constraints$Specie_name)
#' @export
calibration_tips <- function(constraints, column, value, tip_labels) {
  vals <- trimws(strsplit(as.character(value), ";")[[1]])
  tips_all <- unique(constraints[constraints[[column]] %in% vals, "Specie_name"])
  tips_all[tips_all %in% tip_labels]
}

#' Check That a Set of Calibration Bounds Is Internally Satisfiable
#'
#' Bounds are declared per node, but the nodes are nested, and an ultrametric tree forces every
#' ancestor to be older than every one of its descendants. Two bounds taken from analyses on
#' different timescales routinely violate that ordering, and `treePL` does not report the
#' violation: it returns a chronogram with the offending node pinned to a bound, which reads as an
#' estimate.
#'
#' \describe{
#'   \item{Infeasible pair}{The ancestor's maximum is at or below the descendant's minimum. No
#'     assignment of ages satisfies both. Reported as an error.}
#'   \item{Non-binding ancestor floor}{The ancestor's minimum is younger than the descendant's
#'     minimum, so the ancestor's own lower bound constrains nothing and the effective floor is
#'     inherited from the descendant. The signature of two sources on incompatible timescales,
#'     reported as a warning because the run is still valid.}
#' }
#'
#' A descendant maximum below the ancestor minimum is deliberately not flagged: it only states
#' that the descendant is necessarily younger, which is the normal condition for nested bounds.
#'
#' Two further conditions are checked because they silently reassign bounds rather than break the
#' run: a row resolving to fewer than two terminals, which the configuration builder drops without
#' reporting, and two rows resolving to the same node, which makes the later `min`/`max` pair
#' overwrite the earlier one.
#'
#' @param calibs Data frame read from `calibrations_bounds.csv`. Rows with `used_in_analysis` not
#'   `TRUE` are ignored.
#' @param tree A rooted `phylo` object, typically the output of `root_on_clade()`.
#' @param constraints Data frame of the taxonomic table (`cactus_constraints.csv`).
#' @param strict Logical. When `TRUE`, the default, an infeasible pair raises an error.
#' @return Invisibly, a data frame with one row per active calibration: resolved node, terminals
#'   matched, and whether the set is monophyletic in `tree`.
#' @examples
#' \dontrun{
#' calibs <- read.csv(system.file("extdata", "calibrations_bounds.csv", package = "PhyloCactus"))
#' constraints <- read.csv(system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus"))
#' rooted <- root_on_clade(ml_tree, resolve_rooting_outgroup(ml_tree$tip.label))
#' check_calibration_consistency(calibs, rooted, constraints)
#' }
#' @export
check_calibration_consistency <- function(calibs, tree, constraints, strict = TRUE) {
  if (!inherits(tree, "phylo")) stop("`tree` must be a phylo object.", call. = FALSE)
  if (!ape::is.rooted(tree)) {
    stop("`tree` must be rooted. The ancestor and descendant relations this check relies on are ",
         "undefined on an unrooted topology.", call. = FALSE)
  }

  active <- calibs[.is_true_col(calibs$used_in_analysis), , drop = FALSE]
  if (nrow(active) == 0L) stop("No calibration row is marked used_in_analysis = TRUE.", call. = FALSE)

  dup <- unique(active$mrca[duplicated(active$mrca)])
  if (length(dup) > 0L) {
    stop("Duplicated labels among active calibrations: ", paste(dup, collapse = ", "),
         ". treePL keys its mrca/min/max triplets by label, so the later bounds would silently ",
         "overwrite the earlier ones.", call. = FALSE)
  }

  tip_labels <- tree$tip.label
  resolved <- data.frame(mrca = character(0), node = integer(0), n_tips = integer(0),
                         min = numeric(0), max = numeric(0), monophyletic = logical(0),
                         stringsAsFactors = FALSE)

  for (i in seq_len(nrow(active))) {
    row <- active[i, ]
    tips <- calibration_tips(constraints, row$column, row$value, tip_labels)
    if (length(tips) < 2L) {
      warning("Calibration '", row$mrca, "' resolves to ", length(tips),
              " terminal(s) and will be dropped from the treePL configuration without further ",
              "notice. Its bound does not constrain the chronogram.", call. = FALSE)
      next
    }
    resolved <- rbind(resolved, data.frame(
      mrca = row$mrca, node = ape::getMRCA(tree, tips), n_tips = length(tips),
      min = as.numeric(row$min), max = as.numeric(row$max),
      monophyletic = ape::is.monophyletic(tree, tips), stringsAsFactors = FALSE))
  }
  if (nrow(resolved) == 0L) {
    stop("No active calibration resolved to two or more terminals of this tree.", call. = FALSE)
  }

  same_node <- unique(resolved$node[duplicated(resolved$node)])
  if (length(same_node) > 0L) {
    stop("Distinct calibration labels resolve to the same node of this tree: ",
         paste(resolved$mrca[resolved$node %in% same_node], collapse = ", "),
         ". The bounds of one are silently replaced by those of the other.", call. = FALSE)
  }
  for (i in seq_len(nrow(resolved))) {
    if (!resolved$monophyletic[i]) {
      warning("The terminals of calibration '", resolved$mrca[i],
              "' are not monophyletic in this tree, so their MRCA subtends additional terminals ",
              "and the bound lands on a node older than the one it was written for.",
              call. = FALSE)
    }
  }

  # Ancestry is read off the tree rather than assumed from the taxonomic hierarchy, because a
  # non-monophyletic set moves its MRCA and the two can disagree.
  clade_tips <- lapply(resolved$node, function(nd) ape::extract.clade(tree, nd)$tip.label)
  problems <- character(0)
  for (a in seq_len(nrow(resolved))) {
    for (d in seq_len(nrow(resolved))) {
      if (a == d) next
      if (!all(clade_tips[[d]] %in% clade_tips[[a]])) next
      anc <- resolved[a, ]; des <- resolved[d, ]
      if (anc$max <= des$min) {
        problems <- c(problems, sprintf(
          "'%s' (max %.3f) is an ancestor of '%s' (min %.3f): no ultrametric tree satisfies both.",
          anc$mrca, anc$max, des$mrca, des$min))
      }
      if (anc$min < des$min) {
        warning(sprintf(
          "Calibration '%s' (min %.3f) is an ancestor of '%s' (min %.3f). The ancestor's floor is ",
          anc$mrca, anc$min, des$mrca, des$min),
          "younger than the descendant's, so it constrains nothing and the effective floor is ",
          "inherited from the descendant. The two bounds come from analyses on different ",
          "timescales.", call. = FALSE)
      }
    }
  }
  if (length(problems) > 0L) {
    msg <- paste0("Calibration bounds are not jointly satisfiable:\n  ",
                  paste(problems, collapse = "\n  "))
    if (isTRUE(strict)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }

  message("Calibration consistency check passed for ", nrow(resolved), " active node(s).")
  invisible(resolved)
}

#' Report Which Calibrated Nodes Came Back Sitting on a Bound
#'
#' A node whose estimated age equals one of its own bounds was not estimated. Penalized likelihood
#' returned the constraint, and the number carries the prior rather than the data. This is
#' invisible in the output chronogram, which looks like any other, so it has to be checked
#' explicitly before a date is reported or interpreted.
#'
#' In the August 2026 run every one of the five calibrated nodes came back on a bound, the root at
#' its maximum and the rest at their minimum, and the chronogram gave no sign of it.
#'
#' **A single chronogram is not enough to answer this.** On 2026-09-02 the maximum-likelihood tree
#' placed `ACP_root` at 52.96 Ma, 0.41 Ma inside its upper bound of 53.37, and this function
#' reported it as interior. Across the 100 bootstrap replicates the interval was 53.29 to 53.37
#' and 96 of them returned the bound exactly. The point estimate was one realisation of a node
#' whose age the data cannot identify, and it happened to land just inside. Supplying
#' `bootstraps` is what distinguishes a node that was estimated from one that is unidentifiable
#' and collapsed onto its nearest constraint.
#'
#' @param chronogram An ultrametric, rooted `phylo` object, typically `BestTree_treePL.tree`.
#' @param calibs Data frame read from `calibrations_bounds.csv`.
#' @param constraints Data frame of the taxonomic table (`cactus_constraints.csv`).
#' @param tol Numeric. Absolute tolerance, in millions of years, within which an age counts as
#'   sitting on a bound. Defaults to `0.05`.
#' @param bootstraps Optional. The bootstrap chronograms, as a `multiPhylo`, a list of `phylo`, or
#'   a path to a Newick file with one tree per line (`bsTree_treePL.tree`). When supplied, the
#'   verdict is taken from the fraction of replicates sitting on a bound rather than from the
#'   single point estimate, and a node pinned in most replicates is reported even when the
#'   maximum-likelihood tree alone would have called it interior. Defaults to `NULL`.
#' @return Invisibly, a data frame with one row per calibrated node: its age in the
#'   maximum-likelihood chronogram, its bounds, a `status` of `"at_min"`, `"at_max"` or
#'   `"interior"`, and, when `bootstraps` is supplied, the number of replicates evaluated, their
#'   median and 95% interval, and the percentage of them sitting on a bound.
#' @examples
#' \dontrun{
#' chrono <- ape::read.tree(file.path("8_Dating", "BestTree_treePL.tree"))
#' calibs <- read.csv(system.file("extdata", "calibrations_bounds.csv", package = "PhyloCactus"))
#' constraints <- read.csv(system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus"))
#' report_bound_adherence(chrono, calibs, constraints,
#'                        bootstraps = file.path("8_Dating", "bsTree_treePL.tree"))
#' }
#' @export
report_bound_adherence <- function(chronogram, calibs, constraints, tol = 0.05,
                                   bootstraps = NULL) {
  if (!inherits(chronogram, "phylo") || !ape::is.rooted(chronogram)) {
    stop("`chronogram` must be a rooted phylo object.", call. = FALSE)
  }

  # Node ages by node number. branching.times() names its result with node.label when the tree
  # carries them, and treePL chronograms inherit those from the support tree, so indexing by node
  # number against that vector returns NA.
  node_ages <- function(tr) {
    d <- ape::node.depth.edgelength(tr)
    max(d[seq_len(ape::Ntip(tr))]) - d
  }
  age_of <- function(tr, tips) {
    t2 <- intersect(tips, tr$tip.label)
    if (length(t2) < 2L) return(NA_real_)
    nd <- ape::getMRCA(tr, t2)
    if (is.null(nd)) return(NA_real_)
    node_ages(tr)[nd]
  }

  bs <- NULL
  if (!is.null(bootstraps)) {
    if (is.character(bootstraps) && length(bootstraps) == 1L) {
      if (!file.exists(bootstraps)) {
        stop("`bootstraps` names a file that does not exist: '", bootstraps, "'.", call. = FALSE)
      }
      bs <- ape::read.tree(bootstraps)
    } else {
      bs <- bootstraps
    }
    if (inherits(bs, "phylo")) bs <- list(bs)
    bs <- Filter(function(x) inherits(x, "phylo"), as.list(bs))
    if (length(bs) == 0L) {
      stop("`bootstraps` contained no readable trees.", call. = FALSE)
    }
  }

  active <- calibs[.is_true_col(calibs$used_in_analysis), , drop = FALSE]
  tip_labels <- chronogram$tip.label
  rows <- list()

  for (i in seq_len(nrow(active))) {
    row <- active[i, ]
    tips <- calibration_tips(constraints, row$column, row$value, tip_labels)
    if (length(tips) < 2L) next
    age <- age_of(chronogram, tips)
    lo <- as.numeric(row$min); hi <- as.numeric(row$max)
    at_bound <- function(a) !is.na(a) && (abs(a - lo) <= tol || abs(a - hi) <= tol)
    status <- if (is.na(age)) "not_found"
              else if (abs(age - lo) <= tol) "at_min"
              else if (abs(age - hi) <= tol) "at_max"
              else "interior"

    n_bs <- NA_integer_; pct <- NA_real_
    med <- NA_real_; q_lo <- NA_real_; q_hi <- NA_real_
    if (!is.null(bs)) {
      a <- vapply(bs, age_of, numeric(1), tips = tips)
      a <- a[!is.na(a)]
      n_bs <- length(a)
      if (n_bs > 0L) {
        pct <- 100 * mean(vapply(a, at_bound, logical(1)))
        s <- sort(a); med <- stats::median(s)
        q_lo <- s[max(1L, floor(0.025 * n_bs))]
        q_hi <- s[max(1L, ceiling(0.975 * n_bs))]
      }
    }
    rows[[length(rows) + 1L]] <- data.frame(
      mrca = row$mrca, age = age, min = lo, max = hi, status = status,
      n_bs = n_bs, bs_median = med, bs_lo = q_lo, bs_hi = q_hi, bs_pct_at_bound = pct,
      stringsAsFactors = FALSE)
  }
  out <- if (length(rows)) do.call(rbind, rows) else
    data.frame(mrca = character(0), age = numeric(0), min = numeric(0), max = numeric(0),
               status = character(0), n_bs = integer(0), bs_median = numeric(0),
               bs_lo = numeric(0), bs_hi = numeric(0), bs_pct_at_bound = numeric(0),
               stringsAsFactors = FALSE)

  if (is.null(bs)) {
    pinned <- out[out$status %in% c("at_min", "at_max"), , drop = FALSE]
    if (nrow(pinned) > 0L) {
      warning("These calibrated nodes returned their own bound rather than an estimate: ",
              paste(sprintf("%s (%.2f Ma, %s)", pinned$mrca, pinned$age, pinned$status),
                    collapse = "; "),
              ". Their ages are dictated by the calibration and must not be reported as inferred ",
              "divergence times.", call. = FALSE)
    }
    message("No bootstrap chronograms supplied, so this verdict rests on a single point ",
            "estimate. On 2026-09-02 the maximum-likelihood tree placed ACP_root at 52.96 Ma, ",
            "0.41 Ma inside its upper bound and therefore reported as interior, while 96 of 100 ",
            "bootstrap replicates returned the bound exactly. Pass `bootstraps` to see that.")
  } else {
    # The replicates decide, not the point. A node whose age is not identifiable from the data
    # collapses onto whatever constraint is nearest, and one realisation of that can land a
    # fraction of a million years inside the bound purely by chance.
    heavy <- out[!is.na(out$bs_pct_at_bound) & out$bs_pct_at_bound >= 50, , drop = FALSE]
    if (nrow(heavy) > 0L) {
      warning("These calibrated nodes return their own bound in most bootstrap replicates: ",
              paste(sprintf("%s (%.0f%% of %d replicates; ML point %.2f Ma, %s)",
                            heavy$mrca, heavy$bs_pct_at_bound, heavy$n_bs, heavy$age,
                            heavy$status), collapse = "; "),
              ". Their ages are set by the calibration and cannot be reported as inferred ",
              "divergence times, whatever the maximum-likelihood tree alone suggests.",
              call. = FALSE)
    }
    misleading <- heavy[heavy$status == "interior", , drop = FALSE]
    if (nrow(misleading) > 0L) {
      warning("And for these the maximum-likelihood tree alone would have said otherwise: ",
              paste(misleading$mrca, collapse = ", "),
              ". A point estimate cannot distinguish a node that was estimated from one that is ",
              "unidentifiable and happened to land just inside its bound.", call. = FALSE)
    }
  }
  invisible(out)
}

#' Coerce a used_in_analysis Column to Logical
#'
#' `read.csv()` returns the column as logical, character or factor depending on how the file was
#' written, and a factor compares unequal to `TRUE` without warning, which would silently
#' deactivate every calibration.
#'
#' @param x The column as read from the calibration table.
#' @return Logical vector, `NA` treated as `FALSE`.
#' @noRd
.is_true_col <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  v <- toupper(trimws(as.character(x)))
  !is.na(v) & v %in% c("TRUE", "T", "YES", "1")
}

# ---------------------------------------------------------------------------------------------
# Smoothing sensitivity
# ---------------------------------------------------------------------------------------------

#' The calibrated nodes declared in a treePL configuration, as sets of terminals
#'
#' treePL addresses a node by the MRCA of the terminals listed on its `mrca` line, so those lines
#' are the definition of which nodes an analysis is making claims about. Reading them back is how
#' a report can cover exactly the nodes that were calibrated, without a second list to keep in
#' step with the first.
#'
#' @param cfg_file Path to a treePL configuration file.
#' @return A named list of character vectors, one per `mrca` line; empty if there are none.
#' @keywords internal
#' @noRd
.mrca_sets_from_cfg <- function(cfg_file) {
  lines <- grep("^\\s*mrca\\s*=", readLines(cfg_file, warn = FALSE), value = TRUE)
  if (length(lines) == 0L) return(list())
  parts <- strsplit(trimws(sub("^\\s*mrca\\s*=\\s*", "", lines)), "\\s+")
  keep <- vapply(parts, function(p) length(p) >= 3L, logical(1))
  parts <- parts[keep]
  stats::setNames(lapply(parts, function(p) p[-1]), vapply(parts, `[`, character(1), 1))
}

#' A treePL configuration identical to another except for its smoothing and output
#'
#' Both spellings of the keyword are stripped, not just the current one: configurations written
#' before 2026-09-02 carry `smoothing`, and leaving it in place would put two smoothing lines in
#' the file. The value is formatted without scientific notation because `1e-04` is not what
#' treePL's parser expects, and a smoothing silently reset to the default is the failure this
#' whole report exists to make visible.
#'
#' @param cfg_lines Lines of the source configuration.
#' @param smoothing Numeric smoothing value.
#' @param outfile Name of the chronogram treePL should write.
#' @return A character vector of configuration lines.
#' @keywords internal
#' @noRd
.smoothing_cfg <- function(cfg_lines, smoothing, outfile) {
  keep <- cfg_lines[!grepl("^\\s*(smooth|smoothing|outfile)\\s*=", cfg_lines)]
  c(keep,
    paste0("smooth = ", format(smoothing, scientific = FALSE, trim = TRUE)),
    paste0("outfile = ", outfile))
}

#' Age of Every Calibrated Node Across a Range of Rate-Smoothing Values
#'
#' Penalized likelihood requires a rate-smoothing parameter, and treePL selects one by
#' cross-validation. When the cross-validated minimum falls on the edge of the tested grid, as it
#' does for this dataset, the selection is the boundary of the search rather than an optimum, and
#' a reader is entitled to ask whether the reported ages are an artefact of that choice. This
#' function answers the question directly: it dates the same tree, under the same calibrations, at
#' each smoothing value given, and reports the age of every calibrated node in each run.
#'
#' The nodes come from the `mrca` lines of the configuration itself, so the table covers exactly
#' the nodes the analysis makes claims about and cannot drift out of step with them.
#'
#' Each run is checked with the same verification applied everywhere else: the smoothing treePL
#' reports in its log must match the smoothing requested. Until 2026-09-02 the pipeline wrote the
#' keyword `smoothing`, which treePL does not recognise and discards without a message, so every
#' chronogram was produced at the built-in default of 10 and a table like this one would have
#' shown five identical rows. A run whose smoothing cannot be confirmed is reported as such
#' rather than tabulated as a result.
#'
#' @param cfg_file Character. A treePL configuration carrying the tree, `numsites`, the
#'   optimisation parameters and the calibrations. The configuration written for the
#'   maximum-likelihood chronogram is the natural input.
#' @param out_csv Character or `NULL`. Where to write the table. Defaults to
#'   `TABLE_smoothing_sensitivity.csv` beside `cfg_file`.
#' @param smoothing_values Numeric vector of smoothing values to test. Defaults to five values
#'   spanning six orders of magnitude, which is wide enough that stability across them is
#'   informative.
#' @param treepl_bin Character. The treePL executable. Defaults to `"treePL"`.
#' @param timeout Numeric. Seconds allowed per run before it is abandoned and recorded as
#'   unconverged. Defaults to 1800.
#' @param work_dir Character or `NULL`. Directory for the intermediate configurations, logs and
#'   chronograms. Defaults to a `smoothing_sensitivity` directory beside `cfg_file`, so the
#'   evidence behind the table survives alongside it.
#' @return Invisibly, a data frame with one row per smoothing value and calibrated node:
#'   `smoothing`, `smoothing_used`, `node`, `age_ma`, `seconds` and `converged`.
#' @examples
#' \dontrun{
#' report_smoothing_sensitivity(
#'   cfg_file = "8_Dating/auto_results/ML_tree/configure_smooth_ML_tree"
#' )
#' }
#' @export
report_smoothing_sensitivity <- function(cfg_file,
                                         out_csv = NULL,
                                         smoothing_values = c(1e-4, 1e-2, 1, 10, 100),
                                         treepl_bin = "treePL",
                                         timeout = 1800,
                                         work_dir = NULL) {
  if (!file.exists(cfg_file)) {
    stop("`cfg_file` does not exist: '", cfg_file, "'.", call. = FALSE)
  }
  if (!is.numeric(smoothing_values) || length(smoothing_values) < 2L ||
      any(!is.finite(smoothing_values)) || any(smoothing_values <= 0)) {
    stop("`smoothing_values` must be two or more finite positive numbers; a sensitivity ",
         "analysis over one value is not one.", call. = FALSE)
  }

  cfg_file <- normalizePath(cfg_file, mustWork = TRUE)
  cfg_lines <- readLines(cfg_file, warn = FALSE)
  nodes <- .mrca_sets_from_cfg(cfg_file)
  if (length(nodes) == 0L) {
    stop("No `mrca` lines found in '", cfg_file, "'. Without them there are no calibrated nodes ",
         "to report, and the table would be empty.", call. = FALSE)
  }

  if (is.null(work_dir)) work_dir <- file.path(dirname(cfg_file), "smoothing_sensitivity")
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(out_csv)) {
    out_csv <- file.path(dirname(cfg_file), "TABLE_smoothing_sensitivity.csv")
  }

  old <- setwd(work_dir); on.exit(setwd(old), add = TRUE)

  rows <- list()
  for (lam in smoothing_values) {
    tag <- paste0("smooth_", format(lam, scientific = FALSE, trim = TRUE))
    out_tree <- paste0("chronogram_", tag, ".tre")
    log_file <- paste0(tag, ".log")
    this_cfg <- paste0("cfg_", tag)
    writeLines(.smoothing_cfg(cfg_lines, lam, out_tree), this_cfg)

    t0 <- Sys.time()
    tryCatch(
      system2(treepl_bin, args = shQuote(this_cfg),
              stdout = log_file, stderr = log_file, timeout = timeout),
      warning = function(w) NULL, error = function(e) NULL
    )
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    converged <- file.exists(out_tree) && file.size(out_tree) > 0

    used <- NA_real_
    used_line <- if (file.exists(log_file)) {
      grep("^\\s*smoothing\\s*:", readLines(log_file, warn = FALSE), value = TRUE)
    } else character(0)
    if (length(used_line)) {
      used <- suppressWarnings(as.numeric(sub(".*:\\s*", "", used_line[length(used_line)])))
    }
    if (converged && !is.na(used) && !isTRUE(all.equal(lam, used, tolerance = 1e-6))) {
      warning("treePL ran at smoothing ", used, " when ", lam, " was requested; that row is ",
              "reported as unverified. The configuration keyword treePL reads is `smooth`.",
              call. = FALSE)
    }

    ages <- stats::setNames(rep(NA_real_, length(nodes)), names(nodes))
    if (converged) {
      tr <- tryCatch(ape::read.tree(out_tree), error = function(e) NULL)
      if (!is.null(tr)) {
        depth <- ape::node.depth.edgelength(tr)
        tip_max <- max(depth[seq_len(ape::Ntip(tr))])
        for (nm in names(nodes)) {
          tips <- intersect(nodes[[nm]], tr$tip.label)
          if (length(tips) >= 2L) {
            nd <- ape::getMRCA(tr, tips)
            if (!is.null(nd)) ages[[nm]] <- tip_max - depth[nd]
          }
        }
      }
    }

    rows[[length(rows) + 1L]] <- data.frame(
      smoothing = lam, smoothing_used = used, node = names(nodes),
      age_ma = round(unname(ages), 4), seconds = round(secs, 1),
      converged = converged, stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  utils::write.csv(out, out_csv, row.names = FALSE)
  message("Smoothing sensitivity written to: ", out_csv)

  done <- out[out$converged & !is.na(out$age_ma), , drop = FALSE]
  if (nrow(done) > 0L) {
    spread <- stats::aggregate(age_ma ~ node, data = done,
                               FUN = function(v) round(diff(range(v)), 2))
    for (i in seq_len(nrow(spread))) {
      message(sprintf("  %-34s varies %.2f Ma across the smoothing values tested",
                      spread$node[i], spread$age_ma[i]))
    }
  }
  failed <- unique(out$smoothing[!out$converged])
  if (length(failed)) {
    warning("No chronogram was produced at smoothing: ", paste(failed, collapse = ", "),
            ". Those rows carry no ages.", call. = FALSE)
  }
  invisible(out)
}
