#' Collapse Weakly Supported Nodes into Soft Polytomies
#'
#' Internal helper used by `integrate_publication_tree()`. Normalizes node support labels to a
#' 0-1 scale (support values stored as 0-100 percentages, e.g. RAxML-NG's default TBE/FBP output,
#' are detected and divided by 100; values already on a 0-1 scale are left unchanged), then collapses
#' every internal node whose support falls below `collapse_cutoff` into a soft polytomy by zeroing its
#' incident edge length and applying `ape::di2multi()`. Extracted as a standalone, independently
#' testable package-internal function (previously inline code inside `integrate_publication_tree()`).
#'
#' @param ml_tree A `phylo` object with numeric (or percentage) support values in `node.label`.
#' @param collapse_cutoff Numeric. Support threshold (0-1 scale) below which a node is collapsed. Defaults to `0.70`.
#' @return A `phylo` object with weakly supported nodes collapsed into soft polytomies.
#' @noRd
.collapse_weak_support_nodes <- function(ml_tree, collapse_cutoff = 0.70) {
  # Set supports numeric
  supp_vals <- suppressWarnings(as.numeric(ml_tree$node.label))
  if (any(!is.na(supp_vals)) && max(supp_vals, na.rm = TRUE) > 1.0) {
    supp_vals <- supp_vals / 100
  }

  # Nodes with a missing or non-numeric support label must not be silently treated as
  # well-supported: `which(NA < collapse_cutoff)` would otherwise exclude them from
  # collapsing by default. Flag them explicitly instead.
  if (any(is.na(supp_vals))) {
    warning(sum(is.na(supp_vals)), " internal node(s) have a missing or non-numeric support ",
            "label and cannot be evaluated against collapse_cutoff = ", collapse_cutoff,
            "; these nodes will NOT be collapsed regardless of their true support.", call. = FALSE)
  }

  # Identify nodes below cutoff support (FBP)
  nodes_to_collapse <- which(supp_vals < collapse_cutoff) + ape::Ntip(ml_tree)

  # Dynamic editorial collapse (collapses weak nodes)
  collapse_tree <- ml_tree
  for (nd in nodes_to_collapse) {
    edge_idx <- which(collapse_tree$edge[, 2] == nd)
    if (length(edge_idx) == 1 && !is.null(collapse_tree$edge.length)) {
      collapse_tree$edge.length[edge_idx] <- 0
    }
  }

  collapsed_ml <- ape::di2multi(collapse_tree, tol = 1e-10)
  collapsed_ml <- ape::collapse.singles(collapsed_ml)

  # ape::is.rooted() infers rootedness from the root's out-degree: a root left with
  # exactly 2 children reads as rooted, but collapsing a node that sits directly under
  # the root (as above) can leave the root with 3+ children. ape then reinterprets that
  # as an ordinary unrooted trichotomy rather than a soft polytomy, so ape::is.binary()
  # silently reports the tree as still fully resolved. An explicit (zero-length)
  # root.edge keeps the tree's rootedness (and therefore the polytomy) unambiguous.
  if (!ape::is.rooted(collapsed_ml)) {
    collapsed_ml$root.edge <- 0
  }

  collapsed_ml
}

#' Render Final Publication Figures and Registry
#'
#' Maps statistical support values (e.g., Felsenstein Bootstrap Proportions, FBP) onto nodes of the
#' final chronogram and maximum-likelihood phylogeny. Nodes failing to meet the minimum support threshold
#' (`collapse_cutoff`) are systematically collapsed into soft polytomies (analytical uncertainty)
#' to produce conservative, publication-ready figures.
#'
#' @param ml_support_tree_path Character. Best ML support tree path.
#' @param summary_chronogram_path Character. Chronogram path with HPD annotations.
#' @param constraints_path Character. Taxonomy constraints CSV path.
#' @param out_dir Character. Publication figures directory.
#' @param collapse_cutoff Numeric. Felsenstein Bootstrap Proportion (FBP) threshold below which weakly supported nodes are collapsed into soft polytomies (0.0 to 1.0; default 0.70; collapsing 700 of 986 free nodes in the empirical dataset).
#' @return A data frame listing exported figure files and threshold parameters.
#' @export
integrate_publication_tree <- function(ml_support_tree_path, summary_chronogram_path, constraints_path = NULL, out_dir, collapse_cutoff = 0.70) {

  ml_tree <- ape::read.tree(ml_support_tree_path)
  chrono_tree <- ape::read.tree(summary_chronogram_path)

  collapsed_ml <- .collapse_weak_support_nodes(ml_tree, collapse_cutoff = collapse_cutoff)

  # Plot using ggtree
  p <- ggtree::ggtree(collapsed_ml, layout = "rectangular", size = 0.35) +
    ggtree::geom_tiplab(size = 1.0, colour = "black") +
    ggtree::theme_tree()
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(out_dir, "Figure_1_ML_Collapsed.pdf")
  
  ggplot2::ggsave(
    filename = pdf_file,
    plot = p,
    width = 8.5,
    height = 11,
    device = grDevices::pdf
  )
  
  # Chronogram layout plots
  p_chrono <- ggtree::ggtree(chrono_tree, layout = "rectangular", size = 0.35) +
    ggtree::theme_tree2() +
    ggplot2::labs(x = "Time before present (Ma)")
  
  pdf_chrono <- file.path(out_dir, "Figure_2_Chronogram_Dating.pdf")
  ggplot2::ggsave(
    filename = pdf_chrono,
    plot = p_chrono,
    width = 8.5,
    height = 11,
    device = grDevices::pdf
  )
  
  return(data.frame(
    figure = c("Figure_1_ML", "Figure_2_Chrono"),
    path = c(pdf_file, pdf_chrono),
    cutoff_used = collapse_cutoff,
    stringsAsFactors = FALSE
  ))
}
