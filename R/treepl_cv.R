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

#' Parse the optimisation parameters treePL suggests after a priming run
#'
#' `treePL` ends a priming run with the line `PLACE THE LINES BELOW IN THE CONFIGURATION FILE`,
#' followed by an `opt`, `optad` and `optcvad` assignment and, for each, an optional `moredetail`
#' flag. The shell wrapper this replaces read those lines by taking the *last character* of each
#' and reassembling them positionally, which breaks whenever an optional flag is absent and the
#' lines shift. This reads the assignments by name instead.
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
#'   the lines with the lowest `opt` and `optad`. `"modal"` reproduces the behaviour of the shell
#'   wrapper this replaces, which took the most frequent combination; it is kept so that a run made
#'   before 2026-09-02 can be reproduced.
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

#' Read a treePL cross-validation output file and pick the smoothing value
#'
#' The file holds one `chisq: (smoothing) value` line per smoothing value evaluated. The selected
#' value is the one with the lowest chi-square, which is the criterion Maurin (2020) states.
#'
#' When that value is the smallest or largest evaluated, the curve never turned: the minimum lies
#' at the edge of the grid, not inside it, and the true optimum is most likely beyond it. Maurin
#' gives this exact instruction (lower `cvstop` when the lowest chi-square falls on it) and this
#' is what happened to this project on 2026-09-02, when a grid running from 1e+04 down to 1e-04
#' returned 1e-04, its own floor, with the curve still descending. Maurin's own example settled at
#' 1e-06 to 1e-08 on a tree whose branch lengths had been rescaled the way this pipeline rescales
#' them, so a floor of 1e-04 was never low enough.
#'
#' @param cv_file Path to the `cvoutfile` treePL wrote.
#' @return A list with `smoothing`, the full `table`, and `at_edge`.
#' @keywords internal
#' @noRd
.select_cv_smoothing <- function(cv_file) {
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
  at_edge <- isTRUE(all.equal(best, min(sm))) || isTRUE(all.equal(best, max(sm)))

  # A minimum on the edge can mean two different things, and treating them alike is misleading.
  # Either the curve was still falling steeply when the grid ran out, which is the case Maurin
  # describes and which calls for extending it; or it has flattened onto a plateau and the edge is
  # simply where a curve that had already converged happened to stop. The second is not a problem
  # to fix. Measured here as the relative change across the last three grid points at that end: on
  # 2026-09-02 the run went 301.77, 301.08, 301.06, 300.61 over four orders of magnitude, a change
  # of 0.4%, while the ages moved by at most 1.25 Ma across those same four orders.
  plateau <- NA
  if (at_edge && nrow(tab) >= 4L) {
    edge <- if (isTRUE(all.equal(best, min(sm)))) utils::head(tab, 3L) else utils::tail(tab, 3L)
    span <- max(edge$chisq) - min(edge$chisq)
    plateau <- is.finite(span) && span / abs(min(edge$chisq)) < 0.01
  }

  if (isTRUE(at_edge) && isTRUE(plateau)) {
    message("The cross-validation minimum sits at the edge of the grid (smoothing ",
            format(best, scientific = FALSE), "), but the chi-square curve is flat there: it ",
            "changes by less than 1% across the last three values tested. That is a plateau ",
            "reached, not a search cut short, so extending the grid would move the selected value ",
            "without moving the result. Confirm by dating at two smoothing values inside the flat ",
            "region and comparing the ages: see report_smoothing_sensitivity().")
  } else if (at_edge) {
    warning("The cross-validation minimum falls on the edge of its own grid: smoothing ",
            format(best, scientific = FALSE), " with the grid running from ",
            format(min(sm), scientific = FALSE), " to ", format(max(sm), scientific = FALSE),
            ", and the chi-square curve was still changing there. This is where the search ",
            "stopped and not where the optimum lies; Maurin (2020) instructs extending the grid ",
            "past it and repeating. Widen `cvstart`/`cvstop` and rerun before reporting ages ",
            "dated at this value.", call. = FALSE)
  }
  list(smoothing = best, table = tab, at_edge = at_edge, at_plateau = plateau)
}

#' Prime, Cross-Validate and Date a Tree with treePL
#'
#' Runs the three stages of the `treePL` protocol of Maurin (2020) (priming, cross-validation and
#' dating) from R, and returns what each stage chose. It replaces the shell wrapper this package
#' shipped until 2026-09-02.
#'
#' @section Why this exists in R:
#' The shell script this replaces (github.com/tongjial/treepl_wrapper) carries no licence, so it
#' could not be redistributed inside a GPL-3 package. It also wrote `smoothing = ` into the final
#' configuration, a keyword `treePL` does not recognise: the line was discarded without a message
#' and every chronogram produced by this project before 2026-09-02 was dated at the built-in
#' default of 10, with the cross-validation that precedes it having no effect on any result.
#' Writing the configuration here puts that keyword under `.verify_treepl_smoothing()`, which
#' compares what was asked against what `treePL`'s own log reports it used.
#'
#' @section Where this follows the protocol and the wrapper did not:
#' * **Priming selection.** Maurin instructs repeating the priming analysis and taking the lines
#'   with the *lowest* `opt` and `optad`. The shell script took the *most frequent* combination
#'   across repeats. `prime_rule` defaults to Maurin's rule; `"modal"` reproduces the old behaviour.
#' * **Cross-validation method.** Maurin recommends `randomcv`, random subsample and replicate
#'   cross-validation, over the leave-one-out `cv`, as "much faster and may give more stable
#'   results". The shell script hard-coded `cv`. `cv_method` defaults to `randomcv`.
#' * **The grid.** The shell script hard-coded a grid from 1e-04 to 1e+04. Both ends are arguments
#'   here, and a minimum landing on either end is reported rather than returned as if it were an
#'   optimum. See `.select_cv_smoothing()`.
#' * **Parsing.** The priming block is read by keyword rather than by taking the last character of
#'   each line and reassembling positionally, which misaligns whenever an optional `moredetail`
#'   flag is absent.
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
#'   the default pairs with;
#'   the shell script's 100 was chosen for a rule that was stable under it. Raising `n_prime` while
#'   keeping `"lowest"` is a change of analysis, not a refinement of one, and it should be recorded
#'   as such.
#' @param prime_rule `"lowest"` (Maurin) or `"modal"` (the shell script's behaviour). See the note
#'   on `n_prime` above before changing either.
#' @param cv_method `"randomcv"` (Maurin's recommendation) or `"cv"` (leave-one-out).
#' @param cvstart,cvstop Ends of the smoothing grid. The defaults span 1e+03 down to 1e-08, low
#'   enough to contain the 1e-06 to 1e-08 range Maurin reports for a tree with rescaled branch
#'   lengths, which the shell script's floor of 1e-04 was not.
#' @param treepl_bin Name or path of the `treePL` binary.
#' @param work_dir Directory to run in. Defaults to the current one. `treePL` writes beside its
#'   configuration, so every intermediate lands here.
#' @param quiet Suppress the per-stage progress messages.
#' @return Invisibly, a list with `smoothing`, `cv_table`, `prime_lines`, `primes`, `tree_file`
#'   (the chronogram written) and the paths of the three configurations.
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
                          cvstart = 1e3, cvstop = 1e-8,
                          treepl_bin = "treePL",
                          work_dir = NULL,
                          quiet = FALSE) {
  prime_rule <- match.arg(prime_rule)
  cv_method <- match.arg(cv_method)
  stopifnot(is.character(label), length(label) == 1L, nzchar(label))
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
  cv_out <- paste0("cv_", label)
  cv_cfg <- paste0("configure_cv_", label)
  writeLines(c(paste0("treefile = ", tree_abs), base_lines, "thorough", prime_lines,
               cv_method, paste0("cvoutfile = ", cv_out),
               paste0("cvstart = ", format(cvstart, scientific = FALSE)),
               paste0("cvstop = ", format(cvstop, scientific = FALSE))), cv_cfg)

  say("Cross-validating with `", cv_method, "` over ", format(cvstop, scientific = FALSE),
      " to ", format(cvstart, scientific = FALSE), " (this is the slow stage)")
  cv_log <- paste0("cv_", label, ".log")
  status <- system2(treepl_bin, args = shQuote(cv_cfg), stdout = cv_log, stderr = cv_log)
  if (!identical(as.integer(status), 0L)) {
    stop("The treePL cross-validation for '", label, "' exited with status ", status,
         ". See '", file.path(getwd(), cv_log), "'.", call. = FALSE)
  }
  sel <- .select_cv_smoothing(cv_out)
  say("  Smoothing selected: ", format(sel$smoothing, scientific = FALSE),
      " (lowest chi-square of ", nrow(sel$table), " values)")

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
                 prime_lines = prime_lines, primes = primes,
                 tree_file = normalizePath(out_tree, mustWork = TRUE),
                 prime_cfg = prime_cfg, cv_cfg = cv_cfg, smooth_cfg = smooth_cfg))
}
