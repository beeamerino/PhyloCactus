# -------------------------------------------------------------
# Separating imposed topology from estimated topology
# -------------------------------------------------------------
# A topological constraint is reproduced by every bootstrap replicate that was run under it, so the
# support value returned for a constrained bipartition is the constraint restated, not a measurement.
# The distinction matters twice in this package: once when tabulating support for a manuscript, and
# again when placing a calibration, because a bound addressed to a node whose existence was imposed
# dates a node the analysis assumed rather than one it recovered.
#
# The functions here recompute bipartitions from the Newick files rather than trusting node labels or
# clade names, because a calibration is resolved by the MRCA of its terminals and a correct label
# guarantees nothing about which node that is.

# Unrooted bipartitions of `tree`, keyed as a canonical string over `tip_universe`.
#
# A bipartition splits the terminals into two sets, and rooting decides only which of the two is
# written as the clade. Both sides are therefore reduced to one key by always keeping the side that
# excludes `ref_tip`. Comparing rooted clade sets instead is the error that made two earlier audits
# report a difference where the topologies agreed.
.bipartition_keys <- function(tree, tip_universe, ref_tip) {
  parts <- ape::prop.part(tree)
  labels <- attr(parts, "labels")
  vapply(parts, function(idx) {
    side <- labels[idx]
    side <- side[side %in% tip_universe]
    if (ref_tip %in% side) side <- setdiff(tip_universe, side)
    if (length(side) < 2L || length(side) > length(tip_universe) - 2L) return(NA_character_)
    paste(sort(side), collapse = "\r")
  }, character(1))
}

.bipartition_key_of <- function(tips, tip_universe, ref_tip) {
  side <- intersect(tips, tip_universe)
  if (ref_tip %in% side) side <- setdiff(tip_universe, side)
  if (length(side) < 2L || length(side) > length(tip_universe) - 2L) return(NA_character_)
  paste(sort(side), collapse = "\r")
}

#' Classify Every Internal Node as Constrained or Estimated
#'
#' Compares the bipartitions of a reference topology against those of the topological constraint it
#' was inferred under, and reports for each internal node whether the split was imposed or estimated.
#' Support values attached to imposed splits are returned unchanged but must not be reported as
#' support: when the bootstrap replicates were themselves run with `--tree-constraint`, every
#' replicate reproduces those splits and the value is 1 by construction.
#'
#' @param tree An `ape` `phylo` object, or a path to a Newick file. Pass the annotated support tree
#'   (`.raxml.support`) to carry support values through into the result.
#' @param constraint_tree An `ape` `phylo` object, or a path to a Newick file, holding the topology
#'   passed to `RAxML-NG --tree-constraint`. Terminals absent from `tree` (collapsed as identical
#'   during parsing) are dropped before comparison.
#' @return A `data.frame` with one row per internal node: `node` (index in the `ape` numbering),
#'   `n_tips` (size of the smaller side of the bipartition), `constrained` (logical) and `support`
#'   (numeric, `NA` when the tree carries no node labels).
#' @examples
#' \dontrun{
#' nodes <- classify_constrained_nodes(
#'   tree            = "7_Phylogenetics/cactus_support_fbp.raxml.support",
#'   constraint_tree = "7_Phylogenetics/cactus_constraints.tree"
#' )
#' table(nodes$constrained)
#' summary(nodes$support[!nodes$constrained])
#' }
#' @export
classify_constrained_nodes <- function(tree, constraint_tree) {
  if (is.character(tree)) tree <- ape::read.tree(tree)
  if (is.character(constraint_tree)) constraint_tree <- ape::read.tree(constraint_tree)

  tip_universe <- tree$tip.label
  ref_tip <- sort(tip_universe)[[1L]]

  constrained_keys <- stats::na.omit(.bipartition_keys(constraint_tree, tip_universe, ref_tip))
  parts <- ape::prop.part(tree)
  if (length(parts) != tree$Nnode) {
    stop("prop.part() returned ", length(parts), " clades for ", tree$Nnode,
         " internal nodes; node indexing cannot be assumed.", call. = FALSE)
  }
  keys <- .bipartition_keys(tree, tip_universe, ref_tip)

  support <- rep(NA_real_, tree$Nnode)
  if (!is.null(tree$node.label)) {
    support <- suppressWarnings(as.numeric(tree$node.label))
  }

  n_tips <- vapply(parts, function(idx) {
    k <- length(idx)
    min(k, length(tip_universe) - k)
  }, numeric(1))

  data.frame(
    node        = seq_len(tree$Nnode) + ape::Ntip(tree),
    n_tips      = n_tips,
    constrained = !is.na(keys) & keys %in% constrained_keys,
    support     = support,
    stringsAsFactors = FALSE
  )
}

#' Audit the Nodes a Calibration Table Addresses
#'
#' For every row of a calibration table, resolves the node its terminals address in the reference
#' topology and reports whether that node was imposed by the topological constraint and what support
#' it carries. This is the check that has to pass before a chronogram is interpreted: a bound placed
#' on a constrained node dates a node the analysis assumed, and a bound placed on a weakly supported
#' node propagates that uncertainty into every age downstream of it without recording that it did.
#'
#' @param calibrations A `data.frame` in the format of `inst/extdata/calibrations_bounds.csv`, with
#'   columns `mrca`, `column`, `value`, `min`, `max` and `used_in_analysis`.
#' @param constraints A `data.frame` in the format of `inst/extdata/cactus_constraints.csv`, used to
#'   translate a taxon-set name into terminals.
#' @param tree An `ape` `phylo` object, or a path to a Newick file, holding the reference topology.
#' @param constraint_tree An `ape` `phylo` object, or a path to a Newick file, holding the topology
#'   passed to `RAxML-NG --tree-constraint`.
#' @param support_trees Named list of annotated support trees (`phylo` objects or paths), one per
#'   metric, for example `list(fbp = "cactus_support_fbp.raxml.support", tbe = "...")`. One column is
#'   added per element. Defaults to an empty list.
#' @param active_only Logical. Restrict to rows with `used_in_analysis == TRUE`. Defaults to `TRUE`.
#' @return A `data.frame` with one row per calibration: `mrca`, `n_tips`, `node`, `constrained`,
#'   `is_root`, one support column per entry of `support_trees`, and `min` and `max`.
#' @examples
#' \dontrun{
#' audit_calibration_support(
#'   calibrations    = read.csv(system.file("extdata", "calibrations_bounds.csv",
#'                                          package = "PhyloCactus")),
#'   constraints     = read.csv(system.file("extdata", "cactus_constraints.csv",
#'                                          package = "PhyloCactus")),
#'   tree            = "7_Phylogenetics/ml_search/cactus_search.raxml.bestTree",
#'   constraint_tree = "7_Phylogenetics/cactus_constraints.tree",
#'   support_trees   = list(fbp = "7_Phylogenetics/cactus_support_fbp.raxml.support",
#'                          tbe = "7_Phylogenetics/cactus_support_tbe.raxml.support")
#' )
#' }
#' @export
audit_calibration_support <- function(calibrations, constraints, tree, constraint_tree,
                                      support_trees = list(), active_only = TRUE) {
  if (is.character(tree)) tree <- ape::read.tree(tree)
  if (is.character(constraint_tree)) constraint_tree <- ape::read.tree(constraint_tree)

  tip_universe <- tree$tip.label
  ref_tip <- sort(tip_universe)[[1L]]
  constrained_keys <- stats::na.omit(.bipartition_keys(constraint_tree, tip_universe, ref_tip))
  root_node <- ape::Ntip(tree) + 1L

  rows <- if (isTRUE(active_only)) {
    calibrations[calibrations$used_in_analysis == TRUE, , drop = FALSE]
  } else {
    calibrations
  }

  support_lookup <- lapply(support_trees, function(st) {
    if (is.character(st)) st <- ape::read.tree(st)
    keys <- .bipartition_keys(st, tip_universe, ref_tip)
    vals <- suppressWarnings(as.numeric(st$node.label))
    stats::setNames(vals, keys)
  })

  out <- lapply(seq_len(nrow(rows)), function(i) {
    row <- rows[i, ]
    tips <- calibration_tips(constraints, row$column, row$value, tip_universe)
    if (length(tips) < 2L) {
      base <- data.frame(mrca = row$mrca, n_tips = length(tips), node = NA_integer_,
                         constrained = NA, is_root = NA, stringsAsFactors = FALSE)
    } else {
      key <- .bipartition_key_of(tips, tip_universe, ref_tip)
      node <- ape::getMRCA(tree, tips)
      base <- data.frame(
        mrca        = row$mrca,
        n_tips      = length(tips),
        node        = node,
        constrained = !is.na(key) && key %in% constrained_keys,
        is_root     = identical(as.integer(node), as.integer(root_node)),
        stringsAsFactors = FALSE
      )
    }
    for (m in names(support_lookup)) {
      key <- if (length(tips) >= 2L) .bipartition_key_of(tips, tip_universe, ref_tip) else NA_character_
      base[[m]] <- if (is.na(key)) NA_real_ else unname(support_lookup[[m]][key])
    }
    base$min <- row$min
    base$max <- row$max
    base
  })

  do.call(rbind, out)
}
