#' Render Final Publication Figures and Registry
#'
#' Maps statistical support values (e.g., Transfer Bootstrap Expectation, TBE) onto nodes of the
#' final chronogram and maximum-likelihood phylogeny. Nodes failing to meet the minimum support threshold
#' (`collapse_cutoff`) are systematically collapsed into soft polytomies (analytical uncertainty)
#' to produce conservative, publication-ready figures.
#'
#' @param ml_support_tree_path Character. Best ML support tree path.
#' @param summary_chronogram_path Character. Chronogram path with HPD annotations.
#' @param constraints_path Character. Taxonomy constraints CSV path.
#' @param out_dir Character. Publication figures directory.
#' @param collapse_cutoff Numeric. TBE threshold to collapse weak nodes into soft polytomies (0.0 to 1.0; default 0.70).
#' @return A data frame listing exported figure files and threshold parameters.
#' @export
integrate_publication_tree <- function(ml_support_tree_path, summary_chronogram_path, constraints_path, out_dir, collapse_cutoff = 0.70) {
  
  ml_tree <- ape::read.tree(ml_support_tree_path)
  chrono_tree <- ape::read.tree(summary_chronogram_path)
  
  # Set supports numeric
  supp_vals <- suppressWarnings(as.numeric(ml_tree$node.label))
  if (all(is.na(supp_vals)) && !is.null(ml_tree$node.label)) {
    # If supports are % based (0 to 100) translate
    supp_vals <- as.numeric(ml_tree$node.label) / 100
  }
  
  # Identify nodes below cutoff TBE
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
