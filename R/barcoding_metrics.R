# Per-class metrics of the molecular diagnostic branch (11_barcoding/11_metrics/).
# Phase 6D of PhyloCactus 0.5.0, decisions M1 to M4 of BMM (26-09); section 5 of the validation plan.
# It reads the tables of steps 7 and 9 and the library summary, and runs no classifier.

#' The predictions of one classifier at its operating point of step 9
#'
#' The nearest neighbour through the identification path is cut at `t*` of each locus and scheme
#' (6A, D1); IdTaxa at the row `idtaxa_rule` of its operating table, 60 by the decision of 26-09
#' (6B proposal, section 12). The operating point is read, never recomputed.
#' @noRd
.bc_metrics_operating <- function(pred, operating, method, idtaxa_rule) {
  key_p <- paste(pred$locus, pred$scheme)
  if (method == "idtaxa") operating <- operating[operating$rule == idtaxa_rule, , drop = FALSE]
  t_op <- operating$threshold[match(key_p, paste(operating$locus, operating$scheme))]
  if (anyNA(t_op)) {
    miss <- unique(key_p[is.na(t_op)])
    stop("No operating threshold in step 9 for ", paste(miss, collapse = ", "), " (", method,
         "). Run sweep_barcoding_threshold() first.", call. = FALSE)
  }
  out <- lapply(split(seq_len(nrow(pred)), key_p), function(i) {
    p <- pred[i, , drop = FALSE]
    p <- if (method == "idtaxa") .bc_apply_confidence(p, t_op[i[1]]) else .bc_apply_remoteness(p, t_op[i[1]])
    p$operating_threshold <- t_op[i[1]]
    p
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' One row per class of one locus and scheme: queries, correct, assigned, recall and precision
#'
#' The class is the species in scheme E and the genus in scheme G. An assignment is a state 1 at the
#' species (E) or a state 2 at the genus (G): a state 1 of scheme G names an absent species and is a
#' misassignment, never an assignment to its genus. Precision is NA for a class never assigned.
#' @noRd
.bc_metrics_classes <- function(p, scheme) {
  truth <- if (scheme == "species") p$true_species else p$true_genus
  named <- if (scheme == "species") ifelse(p$state == 1L, p$predicted_species, NA_character_) else
    ifelse(p$state == 2L, p$predicted_genus, NA_character_)
  hit <- !is.na(named) & named == truth
  classes <- sort(unique(truth), method = "radix")
  queries <- vapply(classes, function(k) sum(truth == k), integer(1))
  correct <- vapply(classes, function(k) sum(truth == k & hit), integer(1))
  assigned <- vapply(classes, function(k) sum(!is.na(named) & named == k), integer(1))
  data.frame(class = classes, queries = unname(queries), correct = unname(correct),
             assigned = unname(assigned), recall = unname(correct / queries),
             precision = ifelse(assigned > 0L, unname(correct / pmax(assigned, 1L)), NA_real_),
             stringsAsFactors = FALSE)
}

#' Summarise the Per-Class Metrics of the Molecular Diagnostic Branch
#'
#' Step 11 of the branch (Phase 6D). For each classifier at its operating point of step 9, per locus
#' and per scheme, it writes the metrics of section 5 of the validation plan: per-class recall and
#' precision as distributions (median, interquartile range, the classes with recall 0 and the
#' classes never assigned), the unassigned rate and the two kinds of misassignment as separate
#' figures, and the five context columns of section 3 of the plan on every row. It runs no
#' classifier and computes no operating point: both come from steps 7 and 9.
#'
#' @details
#' The class is the species in scheme E and the genus in scheme G. In scheme E an assignment is a
#' state 1; a state 2 with the right genus is an honest genus, reported as its own rate, not as a
#' correct species. In scheme G an assignment is a state 2; a state 1 names a species that is absent
#' from the training set and is a misassignment even when its genus is right. The unassigned, the
#' wrong species and the wrong genus are never summed into one error.
#'
#' Precision is NA, not 0, for a class that is never assigned, and those classes are counted. The
#' interquartile range uses `stats::quantile()` type 7.
#'
#' With `figures = TRUE` it draws, for each locus, the curve of step 9 as the share of correct
#' assignments at the rank of the scheme against the unassigned share, one panel per scheme, both
#' classifiers, their operating points marked.
#'
#' @param classifier_dir Character. Directory of step 7, with the prediction tables.
#' @param threshold_dir Character. Directory of step 9, with the operating and curve tables.
#' @param library_dir Character. Directory of step 4, with `TABLE_barcoding_library_summary.csv`.
#' @param output_dir Character. Where the tables and figures are written.
#' @param methods Character vector, any of `"nn_add"` (the nearest neighbour through the
#'   identification path, at `t*`) and `"idtaxa"`.
#' @param idtaxa_rule Character. The row of the IdTaxa operating table used as its operating
#'   point. Defaults to `"default_60"`, as decided on 2026-09-26.
#' @param figures Logical. Draw the curve figures. Defaults to `TRUE`.
#' @return Invisibly, a list with the class and summary tables of each method.
#' @seealso [sweep_barcoding_threshold()], [classify_barcoding_folds()].
#' @examples
#' \dontrun{
#' summarise_barcoding_metrics()
#' }
#' @export
summarise_barcoding_metrics <- function(classifier_dir = file.path("11_barcoding", "7_classifier"),
                                        threshold_dir = file.path("11_barcoding", "9_threshold"),
                                        library_dir = file.path("11_barcoding", "4_library"),
                                        output_dir = file.path("11_barcoding", "11_metrics"),
                                        methods = c("nn_add", "idtaxa"),
                                        idtaxa_rule = "default_60",
                                        figures = TRUE) {
  methods <- match.arg(methods, c("nn_add", "idtaxa"), several.ok = TRUE)
  .bc_assert_output_dir(output_dir)
  f_lib <- file.path(library_dir, "TABLE_barcoding_library_summary.csv")
  if (!file.exists(f_lib)) {
    stop("No TABLE_barcoding_library_summary.csv in ", library_dir, ". Run finalize_barcoding_library() first.",
         call. = FALSE)
  }
  ctx_cols <- c("total_species", "species_with_replicate", "total_accessions", "total_genera",
                "genera_with_2plus_species")
  lib <- utils::read.csv(f_lib, stringsAsFactors = FALSE)[, c("locus", ctx_cols)]

  op_file <- c(nn_add = "TABLE_barcoding_threshold_operating.csv",
               idtaxa = "TABLE_barcoding_threshold_operating_idtaxa.csv")
  inputs <- list()
  for (m in methods) {
    f_op <- file.path(threshold_dir, op_file[[m]])
    if (!file.exists(f_op)) {
      stop("No ", op_file[[m]], " in ", threshold_dir, ". Run sweep_barcoding_threshold(",
           if (m == "idtaxa") "method = \"idtaxa\"" else "alignment = \"add\"", ") first.", call. = FALSE)
    }
    f_pred <- file.path(classifier_dir, sprintf("TABLE_barcoding_predictions_%s_%s.csv", c("species", "genus"), m))
    if (!all(file.exists(f_pred))) {
      stop("No ", basename(f_pred[!file.exists(f_pred)][1]), " in ", classifier_dir,
           ". Run classify_barcoding_folds() first.", call. = FALSE)
    }
    pred <- do.call(rbind, lapply(f_pred, utils::read.csv, stringsAsFactors = FALSE))
    inputs[[m]] <- .bc_metrics_operating(pred, utils::read.csv(f_op, stringsAsFactors = FALSE), m, idtaxa_rule)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  q7 <- function(x, p) if (length(x)) unname(stats::quantile(x, p, type = 7)) else NA_real_
  out <- list()
  for (m in methods) {
    p_all <- inputs[[m]]
    class_rows <- list(); summ_rows <- list()
    for (k in split(seq_len(nrow(p_all)), paste(p_all$locus, p_all$scheme, sep = "\r"))) {
      p <- p_all[k, , drop = FALSE]
      l <- p$locus[1]; sc <- p$scheme[1]
      cl <- .bc_metrics_classes(p, sc)
      class_rows[[length(class_rows) + 1L]] <- cbind(data.frame(method = m, locus = l, scheme = sc,
                                                                stringsAsFactors = FALSE), cl)
      cat5 <- factor(.bc_outcome_category(p, sc), levels = .bc_outcome_levels())
      n5 <- as.list(as.integer(table(cat5))); names(n5) <- .bc_outcome_levels()
      n <- nrow(p)
      zero <- cl$class[cl$recall == 0]
      pr <- cl$precision[!is.na(cl$precision)]
      op_row <- data.frame(method = m, locus = l, scheme = sc, threshold = p$operating_threshold[1],
                           rule = if (m == "idtaxa") idtaxa_rule else "quantile", stringsAsFactors = FALSE)
      ctx <- lib[match(l, lib$locus), ctx_cols, drop = FALSE]
      if (nrow(ctx) == 0L || anyNA(ctx)) {
        stop("Locus ", l, " is not in the library summary of ", library_dir, ".", call. = FALSE)
      }
      summ_rows[[length(summ_rows) + 1L]] <- cbind(
        op_row, ctx,
        data.frame(queries = n, classes = nrow(cl)), as.data.frame(n5),
        data.frame(correct_rate = (if (sc == "species") n5$correct_species else n5$correct_genus) / n,
                   unassigned_rate = n5$unassigned / n,
                   wrong_species_rate = n5$wrong_species / n,
                   wrong_genus_rate = n5$wrong_genus / n,
                   honest_genus_rate = if (sc == "species") n5$correct_genus / n else NA_real_,
                   recall_median = stats::median(cl$recall), recall_q1 = q7(cl$recall, 0.25),
                   recall_q3 = q7(cl$recall, 0.75),
                   classes_recall_zero = length(zero),
                   classes_recall_zero_list = if (length(zero)) paste(zero, collapse = ";") else NA_character_,
                   precision_median = if (length(pr)) stats::median(pr) else NA_real_,
                   precision_q1 = q7(pr, 0.25), precision_q3 = q7(pr, 0.75),
                   classes_never_predicted = sum(is.na(cl$precision)),
                   stringsAsFactors = FALSE))
    }
    class_tab <- do.call(rbind, class_rows); rownames(class_tab) <- NULL
    summ_tab <- do.call(rbind, summ_rows); rownames(summ_tab) <- NULL
    o <- order(summ_tab$locus, match(summ_tab$scheme, c("species", "genus")), method = "radix")
    summ_tab <- summ_tab[o, , drop = FALSE]; rownames(summ_tab) <- NULL
    o <- order(class_tab$locus, match(class_tab$scheme, c("species", "genus")), class_tab$class, method = "radix")
    class_tab <- class_tab[o, , drop = FALSE]; rownames(class_tab) <- NULL
    utils::write.csv(class_tab, file.path(output_dir, sprintf("TABLE_barcoding_metrics_class_%s.csv", m)), row.names = FALSE)
    utils::write.csv(summ_tab, file.path(output_dir, sprintf("TABLE_barcoding_metrics_summary_%s.csv", m)), row.names = FALSE)
    out[[m]] <- list(class = class_tab, summary = summ_tab)
  }

  if (isTRUE(figures)) .bc_metrics_figures(threshold_dir, output_dir, out, methods)
  .bc_metrics_banner(out, output_dir)
  invisible(out)
}

#' Curve figures of step 11: correct share against unassigned share, both classifiers per locus
#' @noRd
.bc_metrics_figures <- function(threshold_dir, output_dir, out, methods) {
  curve_file <- c(nn_add = "TABLE_barcoding_threshold_curve.csv", idtaxa = "TABLE_barcoding_threshold_curve_idtaxa.csv")
  label <- c(nn_add = "Nearest neighbour (t*)", idtaxa = "IdTaxa (60)")
  curves <- do.call(rbind, lapply(methods, function(m) {
    f <- file.path(threshold_dir, curve_file[[m]])
    if (!file.exists(f)) return(NULL)
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    d$correct <- ifelse(d$scheme == "species", d$correct_species, d$correct_genus)
    data.frame(classifier = label[[m]], locus = d$locus, scheme = d$scheme, threshold = d$threshold,
               unassigned_share = d$unassigned / d$queries, correct_share = d$correct / d$queries,
               stringsAsFactors = FALSE)
  }))
  if (is.null(curves) || nrow(curves) == 0L) return(invisible(NULL))
  points <- do.call(rbind, lapply(methods, function(m) {
    s <- out[[m]]$summary
    data.frame(classifier = label[[m]], locus = s$locus, scheme = s$scheme, unassigned_share = s$unassigned_rate,
               correct_share = s$correct_rate, stringsAsFactors = FALSE)
  }))
  scheme_lab <- c(species = "Scheme E (species)", genus = "Scheme G (genus)")
  for (l in sort(unique(curves$locus), method = "radix")) {
    cv <- curves[curves$locus == l, , drop = FALSE]
    cv <- cv[order(cv$classifier, cv$scheme, cv$unassigned_share, cv$threshold), , drop = FALSE]
    pt <- points[points$locus == l, , drop = FALSE]
    cv$scheme <- factor(scheme_lab[cv$scheme], levels = scheme_lab)
    pt$scheme <- factor(scheme_lab[pt$scheme], levels = scheme_lab)
    g <- ggplot2::ggplot(cv, ggplot2::aes(x = .data$unassigned_share, y = .data$correct_share,
                                          colour = .data$classifier)) +
      ggplot2::geom_path(linewidth = 0.5) +
      ggplot2::geom_point(data = pt, size = 2.5, shape = 21, fill = "white", stroke = 1) +
      ggplot2::facet_wrap(~ scheme) +
      ggplot2::scale_x_continuous(limits = c(0, 1)) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::scale_colour_manual(values = c("Nearest neighbour (t*)" = "#1b6ca8", "IdTaxa (60)" = "#c0392b")) +
      ggplot2::labs(title = l, x = "Unassigned share", y = "Correct at the rank of the scheme",
                    colour = NULL, caption = "Circles: operating points. Correct: species in E, genus in G.") +
      ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "bottom")
    .bc_gap_save(g, file.path(output_dir, paste0("FIG_barcoding_metrics_curve_", l)), width = 180, height = 100)
  }
  invisible(NULL)
}

#' Closing banner of step 11
#' @noRd
.bc_metrics_banner <- function(out, output_dir) {
  cat("\n====================================================\n")
  cat("  Barcoding Metrics Complete \U0001f335\n")
  cat("====================================================\n")
  for (m in names(out)) {
    s <- out[[m]]$summary
    cat("  ", m, ": ", length(unique(s$locus)), " loci, ", nrow(s), " locus-scheme rows\n", sep = "")
  }
  cat("  Output directory:      ", output_dir, "\n")
  cat("====================================================\n\n")
  invisible(TRUE)
}
