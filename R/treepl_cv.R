# -------------------------------------------------------------------------------------------------
# treePL priming, cross-validation and dating, in R
# -------------------------------------------------------------------------------------------------
# Until 2026-09-02 these three stages were driven by a vendored copy of a shell script from
# github.com/tongjial/treepl_wrapper. That script carries no licence, so no permission to
# redistribute it exists, and shipping it inside a GPL-3 package was not tenable. It also wrote
# `smoothing = ` where treePL reads `smooth = `, which silently dated every chronogram this project
# produced at the built-in default of 10.
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
#'   Cactaceae dataset reaches a plateau below 1e-06, which only a grid extending well below that
#'   range can show. See the section on curve shapes below.
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
#'   `prime_lines`, `primes`, `tree_file` (the chronogram written) and the paths of the three
#'   configurations.
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

  # ---- Stage 3: dating --------------------------------------------------------------------------
  # `smooth`, never `smoothing`. run_treePL_direct() confirms against treePL's own log that the
  # value asked for is the value used, and reports a run whose optimisation never moved.
  out_tree <- paste0("treepl_", label, ".tre")
  smooth_cfg <- paste0("configure_smooth_", label)
  writeLines(c(paste0("treefile = ", tree_abs), base_lines, "thorough", prime_lines,
               paste0("smooth = ", format(sel$smoothing, scientific = FALSE)),
               paste0("outfile = ", out_tree)), smooth_cfg)

  say("Dating at smooth = ", format(sel$smoothing, scientific = FALSE))
  run_treePL_direct(smooth_cfg, label)
  .validate_treepl_output(out_tree, label)

  invisible(list(smoothing = sel$smoothing, cv_table = sel$table, cv_at_edge = sel$at_edge,
                 cv_shape = sel$shape,
                 prime_lines = prime_lines, primes = primes,
                 tree_file = normalizePath(out_tree, mustWork = TRUE),
                 prime_cfg = prime_cfg, cv_cfg = cv_cfg, smooth_cfg = smooth_cfg))
}
