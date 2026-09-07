#' @noRd
.is_absolute_path <- function(p) grepl("^(/|~|[A-Za-z]:[/\\\\])", p)

#' Expand a Path to an Absolute Location Without Requiring It to Exist
#'
#' Downstream stages need a usable path for files a previous stage has not written yet, so
#' `normalizePath(mustWork = TRUE)` cannot be used here.
#' @param p Character. Path to expand.
#' @return Absolute path as a character scalar, or `NA_character_` for an empty input.
#' @noRd
.abs_path <- function(p) {
  if (length(p) != 1L || is.na(p) || !nzchar(p)) return(NA_character_)
  p <- path.expand(p)
  if (!.is_absolute_path(p)) p <- file.path(getwd(), p)
  suppressWarnings(normalizePath(p, winslash = "/", mustWork = FALSE))
}

#' Resolve Every File Path of a Phylogenetic Run from Its Naming Convention
#'
#' Reconstructs the absolute path of every file a `PhyloCactus` inference run produces, given only
#' the output directory, the run prefix and the supermatrix location. Because each path is derived
#' rather than carried in memory, any stage of the pipeline can be executed in a fresh R session
#' without first re-running the stages before it.
#'
#' @details
#' The pipeline is a sequence of stages whose outputs feed the next, and a linear script holds those
#' outputs in variables such as `analysed_phy` or `ml_results`. Restarting R and resuming at a later
#' stage clears those variables, and the stage fails on a missing object rather than on a missing
#' file. Hard-coding the filenames instead removes the dependency on session state but reintroduces
#' the problem this naming convention exists to solve: the same literal name is retyped at a dozen
#' call sites, drifts out of step with the run prefix, and continues to be read after it has stopped
#' describing its contents.
#'
#' Deriving the paths from the convention keeps a single definition of every filename while leaving
#' each stage independently runnable. Two entries are resolved by inspection rather than by
#' convention alone, because they record facts about the run rather than choices about naming:
#'
#' \itemize{
#'   \item `analysed_phy` is the reduced PHYLIP when `RAxML-NG` collapsed identical terminals during
#'     validation, and the supermatrix itself when it did not. `was_reduced` reports which occurred.
#'   \item `best_tree` and `ml_trees` resolve to the `ml_search/` subdirectory when the maximum
#'     likelihood search was submitted to a cluster, and to the run root when it was executed locally.
#' }
#'
#' Passing `require` converts a missing input into an error naming the stage that produces it,
#' instead of letting the failure surface several calls later as an unreadable file.
#'
#' @param output_dir Character. Run directory holding the inference outputs. Defaults to `"7_Phylogenetics"`.
#' @param prefix Character. Run prefix shared by every file of the run. Defaults to `"cactus"`.
#' @param supermatrix Character. Path to the concatenated supermatrix exported by
#'   [run_concatenation_pipeline()]. Defaults to the Module 6 location.
#' @param require Character vector. Names of entries that must already exist on disk. Any that do
#'   not trigger an error naming the module responsible for writing them. Defaults to `NULL`.
#' @return An object of class `cactus_run_paths`: a named list of absolute paths, plus `was_reduced`,
#'   `ml_on_cluster`, and an `exists` list of logicals reporting which files are present.
#' @seealso [preprocess_partitions()], [calculate_ml_tree()], [collect_bootstraps()].
#' @examples
#' \dontrun{
#' # Start of any stage, in a fresh R session
#' paths <- resolve_run_paths(output_dir = "7_Phylogenetics", prefix = "cactus")
#' paths
#'
#' # Refuse to start Module 9 unless Module 8 actually finished
#' paths <- resolve_run_paths(require = c("best_tree", "constraint_tree"))
#' }
#' @export
resolve_run_paths <- function(output_dir = "7_Phylogenetics",
                              prefix = "cactus",
                              supermatrix = file.path("6_Concatenated", "concatenated_alignments",
                                                      "ALIGNMENT_supermatrix.phy"),
                              require = NULL) {

  stopifnot(length(output_dir) == 1L, length(prefix) == 1L, nzchar(prefix))

  root <- .abs_path(output_dir)
  sm   <- .abs_path(supermatrix)

  reduced_phy <- file.path(root, paste0(prefix, ".raxml.reduced.phy"))
  was_reduced <- file.exists(reduced_phy)

  # The cluster script writes into ml_search/; a local run writes into the run root. Whichever
  # exists is the one the later stages must read.
  ml_dir_cluster  <- file.path(root, "ml_search")
  best_tree_local <- file.path(root, paste0(prefix, "_search.raxml.bestTree"))
  best_tree_hpc   <- file.path(ml_dir_cluster, paste0(prefix, "_search.raxml.bestTree"))
  ml_on_cluster   <- file.exists(best_tree_hpc) && !file.exists(best_tree_local)
  ml_root         <- if (ml_on_cluster) ml_dir_cluster else root

  out <- list(
    root        = root,
    prefix      = prefix,
    supermatrix = sm,

    # Module 7
    validated_part = file.path(root, paste0(prefix, "_partitions_validated.txt")),
    analysed_phy   = if (was_reduced) reduced_phy else sm,
    best_models    = file.path(root, paste0(prefix, "_modeltest.part.aicc")),

    # Module 8
    constraint_tree = file.path(root, paste0(prefix, "_constraints.tree")),
    best_tree       = file.path(ml_root, paste0(prefix, "_search.raxml.bestTree")),
    ml_trees        = file.path(ml_root, paste0(prefix, "_search.raxml.mlTrees")),
    ml_script       = file.path(ml_dir_cluster, "run_ml_search.sh"),

    # Module 9
    all_bootstraps   = file.path(root, paste0(prefix, "_ALL_bootstraps.tree")),
    support_tree     = file.path(root, paste0(prefix, "_support.raxml.support")),
    support_tree_tbe = file.path(root, paste0(prefix, "_support_tbe.raxml.support")),

    # Module 10
    # The prefix here must match what calculate_temporal_bootstraps() is actually given, not the
    # directory name. Reading "<dir>/<dir>.raxml.bootstraps" points at a file that never exists.
    temporal_bs_dir    = file.path(root, paste0(prefix, "_temporal_bs")),
    temporal_bs_prefix = paste0(prefix, "_temporal"),
    temporal_bs        = file.path(root, paste0(prefix, "_temporal_bs"),
                                   paste0(prefix, "_temporal.raxml.bootstraps")),

    was_reduced   = was_reduced,
    ml_on_cluster = ml_on_cluster
  )

  path_fields <- setdiff(names(out),
                         c("prefix", "temporal_bs_prefix", "was_reduced", "ml_on_cluster"))
  out$exists <- stats::setNames(
    lapply(path_fields, function(k) file.exists(out[[k]])),
    path_fields
  )

  if (!is.null(require)) {
    produced_by <- c(
      supermatrix      = "Module 6 (run_concatenation_pipeline)",
      validated_part   = "Module 7 (preprocess_partitions)",
      analysed_phy     = "Module 7 (preprocess_partitions)",
      best_models      = "Module 7 (run_modeltest_ng)",
      constraint_tree  = "Module 8 (build_constraint_scaffold)",
      best_tree        = "Module 8 (calculate_ml_tree or generate_ml_search_script)",
      ml_trees         = "Module 8 (calculate_ml_tree or generate_ml_search_script)",
      all_bootstraps   = "Module 9 (collect_bootstraps)",
      support_tree     = "Module 9 (map_branch_supports)",
      support_tree_tbe = "Module 9 (map_branch_supports)",
      temporal_bs      = "Module 10 (calculate_temporal_bootstraps)"
    )
    unknown <- setdiff(require, path_fields)
    if (length(unknown) > 0) {
      stop("Unknown path name(s) in `require`: ", paste(unknown, collapse = ", "),
           "\nAvailable: ", paste(path_fields, collapse = ", "), call. = FALSE)
    }
    missing <- require[!vapply(require, function(k) isTRUE(out$exists[[k]]), logical(1))]
    if (length(missing) > 0) {
      detail <- vapply(missing, function(k) {
        sprintf("  %-16s %s\n%18s expected at: %s", k,
                if (k %in% names(produced_by)) paste0("written by ", produced_by[[k]]) else "",
                "", out[[k]])
      }, character(1))
      stop("This stage needs files that are not on disk yet:\n",
           paste(detail, collapse = "\n"), call. = FALSE)
    }
  }

  class(out) <- "cactus_run_paths"
  out
}

#' @param x A `cactus_run_paths` object.
#' @param ... Ignored.
#' @rdname resolve_run_paths
#' @export
print.cactus_run_paths <- function(x, ...) {
  cat("PhyloCactus run paths\n")
  cat("  root  : ", x$root, "\n", sep = "")
  cat("  prefix: ", x$prefix, "\n", sep = "")
  cat("  matrix analysed: ",
      if (isTRUE(x$was_reduced)) "reduced (identical terminals collapsed)" else "full supermatrix",
      "\n", sep = "")
  cat("  ML search      : ",
      if (isTRUE(x$ml_on_cluster)) "cluster (ml_search/)" else "local (run root)",
      "\n\n", sep = "")
  for (k in names(x$exists)) {
    cat(sprintf("  [%s] %-16s %s\n",
                if (isTRUE(x$exists[[k]])) "x" else " ", k, x[[k]]))
  }
  invisible(x)
}

#' Resolve the Set of Rooting Terminals from a Tip Label Vector
#'
#' Selects every terminal belonging to the outgroup lineage that carries the root, returning them
#' as a vector suitable for `ape::root()` and for the `outgroup` argument of the inference and
#' dating functions.
#'
#' @details
#' Rooting is a property of the rooted tree, not of the search. `RAxML-NG` returns an unrooted
#' topology, and the root is imposed afterwards by choosing the edge it sits on. Choosing a single
#' terminal of a sampled outgroup clade places that root *inside* the clade, which renders the
#' clade paraphyletic in the final tree and collapses its crown node onto the root. Any calibration
#' addressed by the MRCA of that clade then silently lands on the root instead. Passing the whole
#' clade places the root on its stem edge and avoids both consequences.
#'
#' The trailing underscore in the default pattern is required rather than stylistic. Terminals use
#' the underscore as binomial separator, so a bare `"^Portulaca"` also matches *Portulacaria*
#' (Didiereaceae), which is neither Portulacaceae nor part of the intended rooting sample. The same
#' reasoning is applied in `run_marker_screening()` and `run_concatenation_pipeline()`.
#'
#' @param tip_labels Character vector of terminal labels, typically `phylo$tip.label` or the row
#'   names of the supermatrix.
#' @param pattern Character. Regular expression matched against `tip_labels`. Defaults to
#'   `"^Portulaca_"`, the Portulacaceae sample of the reference Cactaceae dataset.
#' @return Character vector of matching terminals, sorted. Errors when no terminal matches, since a
#'   silently empty rooting set would leave every downstream tree unrooted.
#' @examples
#' resolve_rooting_outgroup(c("Portulaca_oleracea", "Portulacaria_afra", "Opuntia_ficus-indica"))
#' @export
resolve_rooting_outgroup <- function(tip_labels, pattern = "^Portulaca_") {
  if (!is.character(tip_labels) || length(tip_labels) == 0L) {
    stop("`tip_labels` must be a non-empty character vector.", call. = FALSE)
  }
  tips <- sort(unique(grep(pattern, tip_labels, value = TRUE)))
  if (length(tips) == 0L) {
    stop("No rooting terminal matched '", pattern, "' among the ", length(tip_labels),
         " labels supplied. Set `pattern` to the outgroup lineage sampled in this dataset.",
         call. = FALSE)
  }
  tips
}
