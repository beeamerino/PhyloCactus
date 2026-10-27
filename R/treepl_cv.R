# -------------------------------------------------------------------------------------------------
# treePL priming, cross-validation and dating, in R
# -------------------------------------------------------------------------------------------------
# These three stages are driven from R, not by a shell wrapper: the wrapper once used here
# (github.com/tongjial/treepl_wrapper) carries no licence, so it cannot be redistributed in a
# GPL-3 package, and it wrote `smoothing = ` where treePL reads `smooth = `, which silently dated
# every chronogram at the built-in default of 10.
#
# The procedure itself is not the script's: it is the published protocol of Maurin (2020), which is
# CC BY 4.0. What follows implements that protocol directly, which removes the licence problem,
# makes the smoothing keyword impossible to get wrong, and puts three of the protocol's own
# recommendations into effect that the shell script did not follow.
# -------------------------------------------------------------------------------------------------

#' Parse the optimization parameters treePL suggests after a priming run
#'
#' `treePL` ends a priming run with the line `PLACE THE LINES BELOW IN THE CONFIGURATION FILE`,
#' followed by an `opt`, `optad` and `optcvad` assignment and, for each, an optional `moredetail`
#' flag. The assignments are read by name, so an absent optional flag does not shift the values.
#'
#' @param lines Character vector: the captured output of one priming run.
#' @return A one-row data frame with integer `opt`, `optad`, `optcvad` and logical `moredetail`,
#'   `moredetailad`, `moredetailcvad`, or `NULL` when the block is absent.
#' @references
#' Maurin, K. J. L. (2020). An empirical guide for producing a dated phylogeny with treePL in a
#' maximum likelihood framework. *arXiv*:2008.07054. \doi{10.48550/arXiv.2008.07054}
#' @keywords internal
#' @noRd
.parse_treepl_prime <- function(lines) {
  start <- grep("PLACE THE LINES BELOW", lines, fixed = TRUE)
  if (length(start) == 0L) return(NULL)
  block <- lines[seq(start[length(start)] + 1L, length(lines))]

  num <- function(key) {
    hit <- grep(paste0("^\\s*", key, "\\s*="), block, value = TRUE)
    if (length(hit) == 0L) return(NA_integer_)
    suppressWarnings(as.integer(trimws(sub(".*=\\s*", "", hit[1]))))
  }
  # Anchored and terminated so that `optad` is not read as `opt`, and `moredetail` is not matched
  # by `moredetailad`.
  flag <- function(key) any(grepl(paste0("^\\s*", key, "\\s*$"), block))

  out <- data.frame(
    opt            = num("opt"),
    optad          = num("optad"),
    optcvad        = num("optcvad"),
    moredetail     = flag("moredetail"),
    moredetailad   = flag("moredetailad"),
    moredetailcvad = flag("moredetailcvad"),
    stringsAsFactors = FALSE
  )
  if (all(is.na(c(out$opt, out$optad, out$optcvad)))) return(NULL)
  out
}

#' Turn parsed priming results into the configuration lines treePL reads
#'
#' @param primes A data frame of parsed priming runs, one row each.
#' @param rule `"lowest"` follows Maurin (2020), who instructs the user to repeat priming and take
#'   the lines with the lowest `opt` and `optad`. `"modal"` takes the most frequent combination.
#' @return Character vector of configuration lines.
#' @keywords internal
#' @noRd
.prime_cfg_lines <- function(primes, rule = c("lowest", "modal")) {
  rule <- match.arg(rule)
  primes <- primes[stats::complete.cases(primes[, c("opt", "optad", "optcvad")]), , drop = FALSE]
  if (nrow(primes) == 0L) return(character(0))

  if (rule == "lowest") {
    chosen <- primes[order(primes$opt, primes$optad, primes$optcvad), ][1, ]
  } else {
    key <- apply(primes[, 1:6], 1, paste, collapse = "|")
    chosen <- primes[key == names(sort(table(key), decreasing = TRUE))[1], ][1, ]
  }

  c(paste0("opt = ", chosen$opt),
    if (isTRUE(chosen$moredetail)) "moredetail",
    paste0("optad = ", chosen$optad),
    if (isTRUE(chosen$moredetailad)) "moredetailad",
    paste0("optcvad = ", chosen$optcvad),
    if (isTRUE(chosen$moredetailcvad)) "moredetailcvad")
}

#' Read a treePL cross-validation output file, pick the smoothing value and classify the curve
#'
#' The file holds one `chisq: (smoothing) value` line per smoothing value evaluated. The selected
#' value is the one with the lowest chi-square, which is the criterion of Maurin (2020). What the
#' selected value means depends on the shape of the curve, which is classified as one of three:
#'
#' * `"interior"`: the curve turns upward inside the grid, and the selected value is an optimum.
#' * `"edge"`: the minimum is the first or last value of the grid and the curve is still changing
#'   there. The optimum most likely lies beyond the grid; Maurin (2020) instructs extending it.
#' * `"plateau"`: the three values at one end of the grid, that end being the one where the minimum
#'   lies or the low-smoothing end, differ by less than `plateau_tol` of their chi-square and lie
#'   within `plateau_tol` of the minimum. Smoothing values on the plateau fit the data about equally
#'   well, so the selected value is a nominal choice among them and not an optimum.
#'
#' @param cv_file Path to the `cvoutfile` treePL wrote.
#' @param plateau_tol Relative tolerance for the plateau criterion. Defaults to `0.01`.
#' @return A list with `smoothing`, the full `table` in ascending order of smoothing, `at_edge`
#'   and `shape`.
#' @keywords internal
#' @noRd
.select_cv_smoothing <- function(cv_file, plateau_tol = 0.01) {
  if (!file.exists(cv_file) || file.info(cv_file)$size == 0) {
    stop("The cross-validation produced no output at '", cv_file,
         "'. treePL writes this file only when the cross-validation completes.", call. = FALSE)
  }
  lines <- grep("chisq", readLines(cv_file, warn = FALSE), value = TRUE)
  sm <- suppressWarnings(as.numeric(gsub(".*\\(([^)]*)\\).*", "\\1", lines)))
  ch <- suppressWarnings(as.numeric(sub(".*\\)\\s*", "", lines)))
  keep <- is.finite(sm) & is.finite(ch)
  if (!any(keep)) {
    stop("No usable `chisq: (smoothing) value` lines in '", cv_file, "'.", call. = FALSE)
  }
  sm <- sm[keep]; ch <- ch[keep]
  tab <- data.frame(smoothing = sm, chisq = ch)[order(sm), ]
  best <- sm[which.min(ch)]
  best_chisq <- min(ch)
  # Position in the sorted grid, not all.equal(): all.equal() switches to an absolute tolerance of
  # 1.5e-8 when the values compared are that small, so on a grid reaching 1e-14 every smoothing value
  # below about 1e-8 would compare equal to the floor and be read as the edge of the grid.
  best_pos <- which.min(tab$chisq)
  at_low <- best_pos == 1L
  at_high <- best_pos == nrow(tab)
  at_edge <- at_low || at_high

  # Three values at one end of the grid that differ by less than `plateau_tol` and sit within
  # `plateau_tol` of the minimum. The second condition keeps a flat tail lying above an interior
  # optimum from being read as a plateau.
  flat_at <- function(end) {
    if (nrow(tab) < 4L) return(FALSE)
    x <- if (end == "low") utils::head(tab, 3L) else utils::tail(tab, 3L)
    span <- max(x$chisq) - min(x$chisq)
    is.finite(span) &&
      span / abs(min(x$chisq)) < plateau_tol &&
      (min(x$chisq) - best_chisq) / abs(best_chisq) < plateau_tol
  }

  shape <- if (at_high) {
    if (flat_at("high")) "plateau" else "edge"
  } else if (flat_at("low")) {
    "plateau"
  } else if (at_low) {
    "edge"
  } else {
    "interior"
  }

  if (shape == "plateau") {
    message("The cross-validation curve reaches a plateau: the three smoothing values at the ",
            if (at_high) "upper" else "lower", " end of the grid differ by less than ",
            100 * plateau_tol, "% in chi-square and lie within ", 100 * plateau_tol,
            "% of the minimum. The selected value, ", format(best, scientific = FALSE),
            ", is the lowest chi-square among values that fit the data about equally well, not an ",
            "optimum. Report node ages across the plateau: see report_smoothing_sensitivity().")
  } else if (shape == "edge") {
    warning("The cross-validation minimum falls on the edge of its own grid: smoothing ",
            format(best, scientific = FALSE), " with the grid running from ",
            format(min(sm), scientific = FALSE), " to ", format(max(sm), scientific = FALSE),
            ", and the chi-square curve was still changing there. The optimum most likely lies ",
            "beyond the grid; Maurin (2020) instructs extending the grid past it and repeating. ",
            "Widen `cvstart`/`cvstop` and rerun before reporting ages dated at this value.",
            call. = FALSE)
  }
  list(smoothing = best, table = tab, at_edge = at_edge, shape = shape)
}

#' Count the failed random subsamples of a treePL cross-validation
#'
#' Under `randomcv`, treePL prints in its log, for each smoothing value, the running chi-square sum
#' after each random subsample (`(chisq): v`) and then the total. The contribution of a subsample is
#' the difference between consecutive sums, and a subsample has failed when it adds `fail_at` or
#' more: in cross-validations of the reference Cactaceae dataset under five seeds, a failed
#' subsample added 1e+36 to 1e+41, against
#' totals of 30 to 1200 where none failed. The sum is printed with eight significant digits, so a
#' failure hides the contributions of the subsamples after it, and a later failure smaller than the
#' first by more than eight orders of magnitude is not seen; the count is a lower bound. A
#' `smoothing:` block without subsamples (treePL ends the log with the fit at the selected value)
#' is not read.
#'
#' @param cv_log Path to the log of the cross-validation stage.
#' @param fail_at Contribution at or above which a subsample is counted as failed.
#' @return A data frame with `smoothing`, `subsamples` and `failed`, in ascending order of
#'   smoothing, or `NULL` when the log is absent or carries no per-subsample line.
#' @keywords internal
#' @noRd
.read_cv_subsamples <- function(cv_log, fail_at = 1e30) {
  if (!file.exists(cv_log)) return(NULL)
  x <- readLines(cv_log, warn = FALSE)
  start <- grep("^\\s*smoothing\\s*:", x)
  rows <- lapply(seq_along(start), function(i) {
    end <- if (i < length(start)) start[i + 1L] - 1L else length(x)
    b <- x[start[i]:end]
    run <- suppressWarnings(as.numeric(sub("^\\s*\\(chisq\\)\\s*:\\s*", "",
                                           grep("^\\s*\\(chisq\\)\\s*:", b, value = TRUE))))
    run <- run[is.finite(run)]
    if (!length(run)) return(NULL)
    data.frame(smoothing = suppressWarnings(as.numeric(sub("^\\s*smoothing\\s*:\\s*", "", b[1L]))),
               subsamples = length(run), failed = sum(diff(c(0, run)) >= fail_at))
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(NULL)
  out <- out[is.finite(out$smoothing), , drop = FALSE]
  out <- out[order(out$smoothing), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Warn on the failed subsamples of a cross-validation at the values it can select
#'
#' A smoothing value at which a subsample failed carries a total chi-square of the order of the
#' failure, so it drops out of the selection, and the value selected is then the minimum over the
#' values that survived. On the reference Cactaceae dataset, seeds 20260919 and 8675309 lost one
#' subsample at 1e-07 to 1e-04 and at 1e-05 to 1e-03 and selected 1e-02 and 1e-08, with no message.
#' Values above `max_smoothing` are left out: every one of five seeds failed in four to seven
#' subsamples at 10 and above, and none selected a value above 0.1.
#'
#' @param failures The table of `.read_cv_subsamples()`, or `NULL`.
#' @param cv_table The cross-validation table of `.select_cv_smoothing()`.
#' @param selected The selected smoothing value.
#' @param label Run label, for the message.
#' @param max_smoothing Largest smoothing value checked.
#' @return A list with `failures` (as given) and `smoothing_without_failures`, the value with the
#'   lowest chi-square among those of `max_smoothing` or below at which no subsample failed (`NA`
#'   when there is none or the log could not be read).
#' @keywords internal
#' @noRd
.check_cv_subsamples <- function(failures, cv_table, selected, label, max_smoothing = 0.1) {
  if (is.null(failures) || nrow(failures) == 0L) {
    message("The cross-validation log of '", label, "' carries no per-subsample chi-square, so ",
            "failed subsamples were not checked.")
    return(list(failures = failures, smoothing_without_failures = NA_real_))
  }
  same <- function(a, b) match(signif(a, 6), signif(b, 6))
  low <- failures[failures$smoothing <= max_smoothing * (1 + 1e-9), , drop = FALSE]
  ok <- cv_table[!is.na(same(cv_table$smoothing, low$smoothing[low$failed == 0L])), , drop = FALSE]
  best_ok <- if (nrow(ok)) ok$smoothing[which.min(ok$chisq)] else NA_real_
  bad <- low[low$failed > 0L, , drop = FALSE]
  if (nrow(bad)) {
    bad <- bad[order(-bad$smoothing), , drop = FALSE]
    fmt <- function(v) vapply(v, function(s) format(s, digits = 6), character(1))
    each <- paste0(fmt(bad$smoothing), " (", bad$failed, " of ", bad$subsamples,
                   c(" subsamples", rep("", nrow(bad) - 1L)), ")")
    listed <- if (length(each) == 1L) each else
      paste0(paste(each[-length(each)], collapse = ", "), " and ", each[length(each)])
    among <- !is.na(same(selected, bad$smoothing))
    warning("The cross-validation of '", label, "' failed in part: at smoothing ", listed,
            ", a random subsample returned a chi-square of ", format(1e30), " or more. A value with a ",
            "failed subsample carries a total of that order, so the selection falls on the values ",
            "without failures whenever there are any. The selected value, ", fmt(selected), ", ",
            if (among) "is among them. " else "is not among them. ",
            if (is.na(best_ok)) {
              paste0("A subsample failed at every value of ", fmt(max_smoothing), " or below. ")
            } else {
              paste0("The smoothing with the lowest chi-square among the values of ", fmt(max_smoothing),
                     " or below at which no subsample failed is ", fmt(best_ok), " (chi-square ",
                     format(ok$chisq[which.min(ok$chisq)], digits = 8), "). ")
            },
            "Repeating the cross-validation with another `seed` shows whether the selection depends on ",
            "the subsamples drawn; the failures per value are in `cv_failures` of the result.",
            call. = FALSE)
  }
  list(failures = failures, smoothing_without_failures = best_ok)
}

#' Prime, Cross-Validate and Date a Tree with treePL
#'
#' Runs the three stages of the `treePL` protocol of Maurin (2020) (priming, cross-validation and
#' dating) from R, and returns what each stage chose.
#'
#' @section Choices relative to the protocol:
#' * **Priming selection.** Maurin instructs repeating the priming analysis and taking the lines
#'   with the *lowest* `opt` and `optad`. `prime_rule` defaults to this rule; `"modal"` takes the
#'   most frequent combination across repeats instead.
#' * **Cross-validation method.** Maurin recommends `randomcv`, random subsample and replicate
#'   cross-validation, over the leave-one-out `cv`, as "much faster and may give more stable
#'   results". `cv_method` defaults to `randomcv`.
#' * **The grid.** Both ends of the smoothing grid are arguments, and a minimum on either end is
#'   reported as such, not returned as an optimum. See the section on curve shapes below.
#' * **Parsing.** The priming block is read by keyword, so an absent optional `moredetail` flag
#'   does not shift the values.
#' * **Smoothing check.** The smoothing that `treePL` reports in its log is compared with the value
#'   written to the configuration, and the run stops if they differ.
#'
#' @param cfg_file Path to the configuration carrying `numsites`, the calibrations and any
#'   `nthreads`/`seed` lines. The `treefile` line is written by this function.
#' @param tree_file Path to the rooted tree to date.
#' @param label Run label. Every file this function writes carries it.
#' @param n_prime Integer. Priming repeats. The default is 10.
#'
#'   **`n_prime` and `prime_rule` interact, and not symmetrically.** Under `"modal"` the selection
#'   is stable in `n_prime`: more repeats sharpen an estimate of the most frequent combination.
#'   Under `"lowest"` it is not, because the minimum of a sample can only fall as the sample grows,
#'   so a hundred repeats will select lower parameters than ten of the same runs would. Maurin's
#'   instruction is to repeat the priming analysis a few times and take the lowest, which is what
#'   the default pairs with. Raising `n_prime` while
#'   keeping `"lowest"` is a change of analysis, not a refinement of one, and it should be recorded
#'   as such.
#' @param prime_rule `"lowest"` (Maurin) or `"modal"` (most frequent combination). See the note
#'   on `n_prime` above before changing either.
#' @param cv_method `"randomcv"` (Maurin's recommendation) or `"cv"` (leave-one-out).
#' @param cvstart,cvstop Ends of the smoothing grid, evaluated at one value per order of magnitude.
#'   The defaults span 1e+03 down to 1e-14. Maurin (2020) reports optima between 1e-06 and 1e-08
#'   for trees whose branch lengths are rescaled as this package rescales them; the reference
#'   Cactaceae dataset has an interior minimum at 1e-05. Only a grid extending well beyond that
#'   range tells an interior minimum from an edge or a plateau. See the section on curve shapes below.
#' @param cv_nthreads Integer. Number of threads for the cross-validation stage. Defaults to `1L`.
#'   `treePL` evaluates cross-validation replicates in an OpenMP loop, inside which simulated
#'   annealing calls the non-reentrant C `rand()`; the chi-square sum is also accumulated through an
#'   OpenMP reduction, whose order may vary between runs. With more than one thread, two runs with
#'   the same `seed` return different cross-validation curves and can select different smoothing
#'   values. With one thread the curve is reproducible for a given `seed`, at a cost of about 4% in
#'   wall-clock time on the reference dataset.
#' @param treepl_bin Name or path of the `treePL` binary.
#' @param work_dir Directory to run in. Defaults to the current one. `treePL` writes beside its
#'   configuration, so every intermediate lands here.
#' @param quiet Suppress the per-stage progress messages.
#' @return Invisibly, a list with `smoothing`, `cv_table`, `cv_at_edge`, `cv_shape`,
#'   `cv_failures` (per smoothing value, the subsamples evaluated and how many failed; `NULL` when
#'   the log carries no per-subsample line), `cv_smoothing_without_failures` (see the section on
#'   failed subsamples below), `prime_lines`, `primes`, `tree_file` (the chronogram written) and the
#'   paths of the three configurations.
#' @section Shape of the cross-validation curve:
#' The selected smoothing value is the one with the lowest chi-square. `cv_shape` reports what that
#' value means, and a message or warning states it when the function runs:
#' * `"interior"`: the curve turns upward inside the grid; the selected value is an optimum.
#' * `"edge"`: the minimum is the first or last value of the grid and the curve is still changing
#'   there. The optimum most likely lies beyond the grid: extend `cvstart` or `cvstop` and repeat
#'   before interpreting any age (Maurin, 2020). Reported as a warning.
#' * `"plateau"`: the three values at the end of the grid where the curve flattens differ by less
#'   than 1% in chi-square and lie within 1% of the minimum. Values on the plateau fit the data about
#'   equally well, so the selected value is a nominal choice among them; report node ages across the
#'   plateau with [report_smoothing_sensitivity()]. Reported as a message.
#' @section Failed subsamples:
#' Under `randomcv`, `treePL` writes in the log of the cross-validation (`cv_<label>.log`) the
#' running chi-square sum after each random subsample. A subsample whose contribution is 1e+30 or
#' more has failed, and the total at that smoothing value is then of the same order, so the value
#' drops out of the selection. The function counts the failed subsamples at every value and warns
#' when one or more failed at any value of 0.1 or below, the part of the grid where a selection
#' falls in practice: the warning names each such value with how many of how many subsamples failed
#' there, says whether the selected value is among them, and names the smoothing with the lowest
#' chi-square among the values of 0.1 or below at which no subsample failed, which is also returned
#' as `cv_smoothing_without_failures`. Values above 0.1 are left out because they fail under most
#' seeds on the reference dataset and are never selected there. On that dataset, two of five seeds
#' lost one subsample at values next to 1e-05 and selected 1e-02 and 1e-08, while the reference seed
#' lost none at 0.1 or below and selected 1e-05. The running sum is printed with eight significant
#' digits, so the chi-square of the subsamples after a failure cannot be recovered from the log, and
#' the count of failures is a lower bound.
#' @references
#' Maurin, K. J. L. (2020). An empirical guide for producing a dated phylogeny with treePL in a
#' maximum likelihood framework. *arXiv*:2008.07054. \doi{10.48550/arXiv.2008.07054}
#'
#' Sanderson, M. J. (2002). Estimating absolute rates of molecular evolution and divergence times:
#' a penalized likelihood approach. *Molecular Biology and Evolution*, 19(1), 101-109.
#' \doi{10.1093/oxfordjournals.molbev.a003974}
#'
#' Smith, S. A., & O'Meara, B. C. (2012). treePL: divergence time estimation using penalized
#' likelihood for large phylogenies. *Bioinformatics*, 28(20), 2689-2690.
#' \doi{10.1093/bioinformatics/bts492}
#' @export
run_treePL_cv <- function(cfg_file, tree_file, label,
                          n_prime = 10L,
                          prime_rule = c("lowest", "modal"),
                          cv_method = c("randomcv", "cv"),
                          cvstart = 1e3, cvstop = 1e-14,
                          cv_nthreads = 1L,
                          treepl_bin = "treePL",
                          work_dir = NULL,
                          quiet = FALSE) {
  prime_rule <- match.arg(prime_rule)
  cv_method <- match.arg(cv_method)
  stopifnot(is.character(label), length(label) == 1L, nzchar(label))
  cv_nthreads <- as.integer(cv_nthreads)
  if (is.na(cv_nthreads) || cv_nthreads < 1L) {
    stop("`cv_nthreads` must be a positive integer.", call. = FALSE)
  }
  if (!file.exists(cfg_file)) stop("Configuration not found: ", cfg_file, call. = FALSE)
  if (!file.exists(tree_file)) stop("Tree not found: ", tree_file, call. = FALSE)
  cfg_abs <- normalizePath(cfg_file, mustWork = TRUE)
  tree_abs <- normalizePath(tree_file, mustWork = TRUE)
  .assert_rooted_treefile(tree_abs, label)

  if (!is.null(work_dir)) {
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    old <- setwd(work_dir); on.exit(setwd(old), add = TRUE)
  }
  say <- function(...) if (!quiet) message(...)

  user_lines <- readLines(cfg_abs, warn = FALSE)
  # Anything that belongs to a later stage is dropped: this function writes each stage's own
  # directives, and a stale one inherited from the user configuration would run the wrong analysis.
  base_lines <- user_lines[!grepl(
    "^\\s*(treefile|outfile|cvoutfile|smooth|smoothing|prime|cv|randomcv|cvstart|cvstop)\\s*(=|$)",
    user_lines)]

  # ---- Stage 1: priming -------------------------------------------------------------------------
  prime_cfg <- paste0("configure_prime_", label)
  writeLines(c(paste0("treefile = ", tree_abs), base_lines, "thorough", "prime"), prime_cfg)

  say("Priming ", n_prime, " times (", label, ")")
  primes <- vector("list", n_prime)
  for (i in seq_len(n_prime)) {
    out <- suppressWarnings(tryCatch(
      system2(treepl_bin, args = shQuote(prime_cfg),
              stdout = TRUE, stderr = TRUE),
      error = function(e) character(0)
    ))
    primes[[i]] <- .parse_treepl_prime(out)
  }
  primes <- do.call(rbind, Filter(Negate(is.null), primes))
  if (is.null(primes) || nrow(primes) == 0L) {
    stop("No priming run of '", label, "' produced the block treePL prints after ",
         "`PLACE THE LINES BELOW IN THE CONFIGURATION FILE`. Check that '", treepl_bin,
         "' is on the path and that the configuration is valid.", call. = FALSE)
  }
  utils::write.csv(primes, paste0("prime_", label, ".csv"), row.names = FALSE)
  prime_lines <- .prime_cfg_lines(primes, rule = prime_rule)
  say("  ", nrow(primes), "/", n_prime, " runs parsed; selected (", prime_rule, "): ",
      paste(prime_lines, collapse = "; "))

  # ---- Stage 2: cross-validation ----------------------------------------------------------------
  if (cv_nthreads > 1L) {
    warning("`cv_nthreads` is set to ", cv_nthreads, ". treePL evaluates cross-validation replicates ",
            "in an OpenMP loop in which simulated annealing calls the non-reentrant C `rand()`, so two ",
            "runs with the same `seed` can return different curves and select different smoothing ",
            "values. Use `cv_nthreads = 1L` for a curve that is reproducible for a given `seed`.",
            call. = FALSE)
  }

  # Strip any inherited nthreads directive for the cross-validation configuration to ensure cv_nthreads governs
  cv_base_lines <- base_lines[!grepl("^\\s*nthreads\\s*(=|$)", base_lines)]
  cv_out <- paste0("cv_", label)
  cv_cfg <- paste0("configure_cv_", label)
  writeLines(c(paste0("treefile = ", tree_abs), cv_base_lines,
               paste0("nthreads = ", cv_nthreads), "thorough", prime_lines,
               cv_method, paste0("cvoutfile = ", cv_out),
               paste0("cvstart = ", format(cvstart, scientific = FALSE)),
               paste0("cvstop = ", format(cvstop, scientific = FALSE))), cv_cfg)

  say("Cross-validating with `", cv_method, "` over ", format(cvstop, scientific = FALSE),
      " to ", format(cvstart, scientific = FALSE), " (this is the slow stage, cv_nthreads = ", cv_nthreads, ")")
  cv_log <- paste0("cv_", label, ".log")
  status <- system2(treepl_bin, args = shQuote(cv_cfg), stdout = cv_log, stderr = cv_log)
  if (!identical(as.integer(status), 0L)) {
    stop("The treePL cross-validation for '", label, "' exited with status ", status,
         ". See '", file.path(getwd(), cv_log), "'.", call. = FALSE)
  }
  sel <- .select_cv_smoothing(cv_out)
  say("  Smoothing selected: ", format(sel$smoothing, scientific = FALSE),
      " (lowest chi-square of ", nrow(sel$table), " values; curve shape: ", sel$shape, ")")
  failed <- .check_cv_subsamples(.read_cv_subsamples(cv_log), sel$table, sel$smoothing, label)

  # ---- Stage 3: dating --------------------------------------------------------------------------
  # `smooth`, never `smoothing`. run_treePL_direct() confirms against treePL's own log that the
  # value asked for is the value used, and reports a run whose optimisation never moved.
  out_tree <- paste0("treepl_", label, ".tre")
  smooth_cfg <- paste0("configure_smooth_", label)
  writeLines(c(paste0("treefile = ", tree_abs), base_lines, "thorough", prime_lines,
               paste0("smooth = ", format(sel$smoothing, scientific = FALSE)),
               paste0("outfile = ", out_tree)), smooth_cfg)

  say("Dating at smooth = ", format(sel$smoothing, scientific = FALSE))
  run_treePL_direct(smooth_cfg, label, treepl_bin = treepl_bin)
  .validate_treepl_output(out_tree, label)

  invisible(list(smoothing = sel$smoothing, cv_table = sel$table, cv_at_edge = sel$at_edge,
                 cv_shape = sel$shape, cv_failures = failed$failures,
                 cv_smoothing_without_failures = failed$smoothing_without_failures,
                 prime_lines = prime_lines, primes = primes,
                 tree_file = normalizePath(out_tree, mustWork = TRUE),
                 prime_cfg = prime_cfg, cv_cfg = cv_cfg, smooth_cfg = smooth_cfg))
}
