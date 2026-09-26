# ------------------------------------------------------------------------------
# Molecular diagnostic branch (11_barcoding/9_threshold/), threshold of remoteness. PhyloCactus
# 0.5.0, Phase 6A.
#
# The nearest neighbour has no rule of remoteness: in CN2, 288 of 312 alien queries came out in state
# 1 with the name of a cactus species. The rule added here sends a query whose distance to its
# nearest neighbour exceeds a threshold to state 3. The curve reports, for every threshold, what the
# legitimate queries become and how many alien queries are rejected. The operating threshold of a
# locus is the quantile 0.99, type 1, of the distances of the legitimate queries of scheme G, a cactus
# whose species is missing from the library and whose genus is not, which is the query the rule must
# not reject. It is fixed without the alien queries, which then report how many of them it rejects
# (decisions D1 and D2 of BMM, 2026-09-26).
# ------------------------------------------------------------------------------

#' The rule of remoteness applied to a table of predictions
#'
#' A query whose `nn_distance` exceeds `threshold` becomes state 3 with the reason `remoteness`,
#' with no species and no genus. A query already in state 3 keeps its own reason, a query with no
#' distance is left as it is, and the distance itself is never rewritten.
#' @noRd
.bc_apply_remoteness <- function(pred, threshold) {
  far <- !is.na(pred$nn_distance) & pred$nn_distance > threshold & pred$state != 3L
  if (!any(far)) return(pred)
  pred$state[far] <- 3L
  pred$predicted_species[far] <- NA_character_
  pred$predicted_genus[far] <- NA_character_
  pred$reason[far] <- "remoteness"
  pred
}

#' The confidence rule of IdTaxa applied to a table of predictions
#'
#' IdTaxa ran at threshold 0 and wrote the confidence of each rank (Phase 6B, K3). At a threshold
#' `threshold`: species confidence at or above it is state 1; genus confidence alone at or above it is
#' state 2; neither is state 3 with the reason `low_confidence`. A query already in state 3 for another
#' reason (`no_match`) keeps it. The confidences are never rewritten, and `candidates` is left as
#' step 7 wrote it. Applied at the threshold of step 7 it gives back the table of step 7.
#' @noRd
.bc_apply_confidence <- function(pred, threshold) {
  keep3 <- !is.na(pred$reason) & pred$reason != "low_confidence" & pred$state == 3L
  sp <- !is.na(pred$species_confidence) & pred$species_confidence >= threshold & !keep3
  ge <- !sp & !is.na(pred$genus_confidence) & pred$genus_confidence >= threshold & !keep3
  no <- !sp & !ge & !keep3
  pred$state[sp] <- 1L
  pred$predicted_species[sp] <- pred$species_idtaxa[sp]
  pred$predicted_genus[sp] <- pred$genus_idtaxa[sp]
  pred$reason[sp] <- NA_character_
  pred$state[ge] <- 2L
  pred$predicted_species[ge] <- NA_character_
  pred$predicted_genus[ge] <- pred$genus_idtaxa[ge]
  pred$reason[ge] <- NA_character_
  pred$state[no] <- 3L
  pred$predicted_species[no] <- NA_character_
  pred$predicted_genus[no] <- NA_character_
  pred$reason[no] <- "low_confidence"
  pred
}

#' Curve of IdTaxa for one locus and one scheme: one row per confidence threshold
#'
#' The grid is 0, 60, the operating threshold and every confidence observed in the locus, genus and
#' species, legitimate or comparable outgroup. An outgroup query counts as rejected when it ends in
#' state 3.
#' @noRd
.bc_confidence_curve <- function(pred, scheme, outgroup = NULL, extra = NULL) {
  aj <- .bc_comparable_aliens(outgroup)
  n_aj <- if (is.null(aj)) 0L else nrow(aj)
  observed <- c(pred$genus_confidence, pred$species_confidence, extra,
                if (n_aj > 0L) c(aj$genus_confidence, aj$species_confidence))
  grid <- sort(unique(c(0, 60, observed[!is.na(observed)])))
  rows <- lapply(grid, function(t) .bc_confidence_row(pred, scheme, aj, t))
  do.call(rbind, rows)
}

.bc_confidence_row <- function(pred, scheme, aj, t) {
  n_aj <- if (is.null(aj)) 0L else nrow(aj)
  cat5 <- factor(.bc_outcome_category(.bc_apply_confidence(pred, t), scheme), levels = .bc_outcome_levels())
  counts <- as.list(as.integer(table(cat5)))
  names(counts) <- .bc_outcome_levels()
  n_rejected <- if (n_aj > 0L) sum(.bc_apply_confidence(aj, t)$state == 3L) else NA_integer_
  cbind(data.frame(locus = pred$locus[1], scheme = scheme, threshold = t, queries = nrow(pred),
                   stringsAsFactors = FALSE),
        as.data.frame(counts),
        data.frame(outgroup_comparable = n_aj, outgroup_rejected = as.integer(n_rejected)))
}

#' Five disjoint outcomes of a legitimate query (validation plan, sec. 5)
#'
#' A wrong species and an honest genus are not the same error and are never summed. In scheme G the
#' true species is not in the training set, so every state 1 names a false species.
#' @noRd
.bc_outcome_category <- function(pred, scheme) {
  out <- rep("unassigned", nrow(pred))
  e1 <- pred$state == 1L
  e2 <- pred$state == 2L
  sp_ok <- !is.na(pred$predicted_species) & pred$predicted_species == pred$true_species
  ge_ok <- !is.na(pred$predicted_genus) & pred$predicted_genus == pred$true_genus
  if (scheme == "species") {
    out[e1 & sp_ok] <- "correct_species"
    out[e1 & !sp_ok] <- "wrong_species"
  } else {
    out[e1] <- "wrong_species"
  }
  out[e2 & ge_ok] <- "correct_genus"
  out[e2 & !ge_ok] <- "wrong_genus"
  out
}

.bc_outcome_levels <- function() {
  c("correct_species", "wrong_species", "correct_genus", "wrong_genus", "unassigned")
}

#' Exact binomial interval at 95 %, as stats::binom.test() computes it; NA for an empty denominator
#' @noRd
.bc_clopper_pearson <- function(x, n, level = 0.95) {
  if (is.na(n) || n == 0L) return(c(NA_real_, NA_real_))
  a <- (1 - level) / 2
  lo <- if (x == 0) 0 else stats::qbeta(a, x, n - x + 1)
  hi <- if (x == n) 1 else stats::qbeta(1 - a, x + 1, n - x)
  c(lo, hi)
}

#' Alien queries of a locus that can be compared at all: the ones that match it in some direction
#' @noRd
.bc_comparable_aliens <- function(outgroup) {
  if (is.null(outgroup) || nrow(outgroup) == 0L) return(outgroup)
  outgroup[outgroup$orientation != "no_match", , drop = FALSE]
}

#' Curve of one locus and one scheme: one row per threshold
#'
#' The grid is every distance observed in the locus, legitimate or comparable alien, plus zero, so
#' the curve is exact and no grid is chosen. `extra` adds the distances of the other scheme of the
#' locus, so both schemes are read on the same grid. An alien query counts as rejected when it ends
#' in state 3, whatever the reason. With no comparable alien query the rejection is left as NA.
#' @noRd
.bc_threshold_curve <- function(pred, scheme, outgroup = NULL, extra = NULL) {
  aj <- .bc_comparable_aliens(outgroup)
  n_aj <- if (is.null(aj)) 0L else nrow(aj)
  observed <- c(pred$nn_distance, extra, if (n_aj > 0L) aj$nn_distance)
  threshold_grid <- sort(unique(c(0, observed[!is.na(observed)])))
  rows <- lapply(threshold_grid, function(t) {
    cat5 <- factor(.bc_outcome_category(.bc_apply_remoteness(pred, t), scheme),
                   levels = .bc_outcome_levels())
    counts <- as.list(as.integer(table(cat5)))
    names(counts) <- .bc_outcome_levels()
    n_rejected <- if (n_aj > 0L) sum(.bc_apply_remoteness(aj, t)$state == 3L) else NA_integer_
    cbind(data.frame(locus = pred$locus[1], scheme = scheme, threshold = t, queries = nrow(pred),
                     stringsAsFactors = FALSE),
          as.data.frame(counts),
          data.frame(outgroup_comparable = n_aj, outgroup_rejected = as.integer(n_rejected)))
  })
  do.call(rbind, rows)
}

#' Operating threshold of a locus: a declared quantile of the legitimate queries of scheme G
#'
#' It takes the predictions of scheme G and nothing else. The alien queries are not an argument,
#' which is what keeps them out of the choice.
#' @noRd
.bc_operating_point <- function(pred_genus, q = 0.99, type = 1L) {
  if (!all(pred_genus$scheme == "genus")) {
    stop("The operating threshold is derived from scheme genus (G) only (decision D1).", call. = FALSE)
  }
  d <- pred_genus$nn_distance[!is.na(pred_genus$nn_distance)]
  if (length(d) == 0L) return(NA_real_)
  unname(stats::quantile(d, q, type = type))
}

#' Sweep the Threshold of Remoteness of the Molecular Diagnostic Branch
#'
#' Step 9 of the branch. Without a rule of remoteness the nearest neighbour names a cactus species
#' for any query: in the negative control CN2, 288 of 312 outgroup queries came out in state 1. This
#' step applies the rule (a query whose distance to its nearest neighbour exceeds a threshold is
#' state 3, reason `remoteness`) over every threshold observed, and fixes an operating threshold per
#' locus.
#'
#' @details
#' For every locus and scheme the curve reports the legitimate queries in five disjoint outcomes
#' (correct species, wrong species, correct genus, wrong genus, not assigned; a wrong species and an
#' honest genus are never summed) and, when CN2 has comparable outgroup queries for the locus, how
#' many of them the rule rejects. Outgroup queries that match the locus in neither direction are not
#' comparable and stay out.
#'
#' The operating threshold of a locus is the quantile `q` (type `quantile_type` of
#' [stats::quantile()]) of the distances of the legitimate queries of scheme G: a cactus whose
#' species is missing from the library while its genus is present, which is the query the rule must
#' not reject. It is fixed without the outgroup queries, which then report the rejection with its
#' exact binomial interval; a locus with no comparable outgroup query reports it as not measured.
#'
#' @param classifier_dir Character. Output directory of [classify_barcoding_folds()].
#' @param controls_dir Character. Output directory of [run_barcoding_controls()]. When it holds no
#'   `TABLE_barcoding_cn2_queries.csv` the rejection is reported as not measured everywhere.
#' @param output_dir Character. Destination; refused inside a directory of the phylogeny.
#' @param alignment Character. Which predictions of step 7 are swept: `"add"`, the default, reads
#'   the tables of `alignment = "add"`, measured the way a user query is measured; `"library"` reads
#'   those of the joint alignment.
#' @param method Character. `"nn"`, the nearest neighbour and its rule of remoteness, or `"idtaxa"`,
#'   IdTaxa and its confidence threshold (Phase 6B, decision K7 of BMM). With `"idtaxa"` the tables
#'   read are `TABLE_barcoding_predictions_<scheme>_idtaxa.csv` and
#'   `TABLE_barcoding_cn2_queries_idtaxa.csv`, the curve runs over the confidence threshold, the
#'   operating threshold of a locus is the quantile `q` of `genus_confidence` of its scheme G queries
#'   (the rule of D1 turned around, since confidence rises where distance falls), and it is reported
#'   next to 60, DECIPHER's default, in the column `rule`. The tables written carry the suffix
#'   `_idtaxa`.
#' @param q Numeric. Quantile of scheme G that sets the operating threshold. Defaults to `0.99` of
#'   the distances for `"nn"` and `0.01` of the genus confidences for `"idtaxa"`.
#' @param quantile_type Integer. Type of [stats::quantile()]. Defaults to `1L`, whose result is always
#'   a distance observed in the data.
#' @param figures Logical. Write one figure per locus, PDF and PNG at 300 dpi, 180 mm wide.
#' @return Invisibly, a list with `curve` and `operating`. Writes
#'   `TABLE_barcoding_threshold_curve.csv`, `TABLE_barcoding_threshold_operating.csv` and, with
#'   `figures = TRUE`, `FIG_barcoding_threshold_<locus>.pdf` and `.png`.
#' @examples
#' \dontrun{
#' sweep_barcoding_threshold(
#'   classifier_dir = "11_barcoding/7_classifier",
#'   controls_dir = "11_barcoding/8_controls",
#'   output_dir = "11_barcoding/9_threshold"
#' )
#' }
#' @export
sweep_barcoding_threshold <- function(classifier_dir = file.path("11_barcoding", "7_classifier"),
                                      controls_dir = file.path("11_barcoding", "8_controls"),
                                      output_dir = file.path("11_barcoding", "9_threshold"),
                                      alignment = c("add", "library"),
                                      q = NULL,
                                      quantile_type = 1L,
                                      figures = TRUE,
                                      method = c("nn", "idtaxa")) {
  alignment <- match.arg(alignment)
  method <- match.arg(method)
  .bc_assert_output_dir(output_dir)
  if (method == "idtaxa") {
    return(.bc_sweep_idtaxa(classifier_dir, controls_dir, output_dir,
                            q = if (is.null(q)) 0.01 else q, quantile_type = quantile_type,
                            figures = figures))
  }
  if (is.null(q)) q <- 0.99
  suffix <- if (alignment == "add") "nn_add" else "nn"
  files <- file.path(classifier_dir, paste0("TABLE_barcoding_predictions_", c("species", "genus"),
                                               "_", suffix, ".csv"))
  missing <- files[!file.exists(files)]
  if (length(missing) > 0L) {
    stop("No ", basename(missing[1]), " in ", classifier_dir, ". Run classify_barcoding_folds(method = \"nn\"",
         if (alignment == "add") ", alignment = \"add\"", ") first.", call. = FALSE)
  }
  pred <- list(species = utils::read.csv(files[1], stringsAsFactors = FALSE),
               genus = utils::read.csv(files[2], stringsAsFactors = FALSE))
  f_cn2 <- file.path(controls_dir, "TABLE_barcoding_cn2_queries.csv")
  outgroup <- if (file.exists(f_cn2)) utils::read.csv(f_cn2, stringsAsFactors = FALSE) else NULL
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  loci <- sort(unique(c(pred$species$locus, pred$genus$locus)), method = "radix")
  curve_list <- list()
  operating_rows <- list()
  for (l in loci) {
    aj_l <- if (is.null(outgroup)) NULL else outgroup[outgroup$locus == l, , drop = FALSE]
    aj_c <- .bc_comparable_aliens(aj_l)
    n_aj <- if (is.null(aj_c)) 0L else nrow(aj_c)
    t_op <- .bc_operating_point(pred$genus[pred$genus$locus == l, , drop = FALSE], q = q, type = quantile_type)
    all_distances <- c(pred$species$nn_distance[pred$species$locus == l],
               pred$genus$nn_distance[pred$genus$locus == l])
    for (sc in c("species", "genus")) {
      p <- pred[[sc]][pred[[sc]]$locus == l, , drop = FALSE]
      if (nrow(p) == 0L) next
      curve_list[[length(curve_list) + 1L]] <- .bc_threshold_curve(p, sc, outgroup = aj_l, extra = all_distances)
      cat5 <- factor(.bc_outcome_category(if (is.na(t_op)) p else .bc_apply_remoteness(p, t_op), sc),
                     levels = .bc_outcome_levels())
      counts <- as.list(as.integer(table(cat5)))
      names(counts) <- .bc_outcome_levels()
      n_rejected <- if (n_aj > 0L && !is.na(t_op)) sum(.bc_apply_remoteness(aj_c, t_op)$state == 3L) else NA_integer_
      ci <- if (is.na(n_rejected)) c(NA_real_, NA_real_) else .bc_clopper_pearson(n_rejected, n_aj)
      operating_rows[[length(operating_rows) + 1L]] <- cbind(
        data.frame(locus = l, scheme = sc, threshold = t_op, q = q, quantile_type = as.integer(quantile_type),
                   queries = nrow(p), stringsAsFactors = FALSE),
        as.data.frame(counts),
        data.frame(outgroup_comparable = n_aj, outgroup_rejected = as.integer(n_rejected),
                   rejection = if (is.na(n_rejected)) NA_real_ else n_rejected / n_aj,
                   rejection_ci_low = ci[1], rejection_ci_high = ci[2],
                   rejection_status = if (is.na(n_rejected)) "not measured" else "measured",
                   stringsAsFactors = FALSE))
    }
  }
  curve_tab <- do.call(rbind, curve_list)
  oper_tab <- do.call(rbind, operating_rows)
  rownames(curve_tab) <- NULL
  rownames(oper_tab) <- NULL
  utils::write.csv(curve_tab, file.path(output_dir, "TABLE_barcoding_threshold_curve.csv"), row.names = FALSE)
  utils::write.csv(oper_tab, file.path(output_dir, "TABLE_barcoding_threshold_operating.csv"), row.names = FALSE)

  if (isTRUE(figures)) {
    for (l in unique(curve_tab$locus)) {
      g <- .bc_threshold_plot(curve_tab[curve_tab$locus == l, , drop = FALSE], oper_tab[oper_tab$locus == l, , drop = FALSE], l)
      .bc_gap_save(g, file.path(output_dir, paste0("FIG_barcoding_threshold_", l)), width = 180, height = 110)
    }
  }
  .bc_threshold_banner(oper_tab, output_dir)
  invisible(list(curve = curve_tab, operating = oper_tab))
}

#' Step 9 with IdTaxa (Phase 6B, K7)
#' @noRd
.bc_sweep_idtaxa <- function(classifier_dir, controls_dir, output_dir, q, quantile_type, figures) {
  files <- file.path(classifier_dir, paste0("TABLE_barcoding_predictions_", c("species", "genus"), "_idtaxa.csv"))
  missing <- files[!file.exists(files)]
  if (length(missing) > 0L) {
    stop("No ", basename(missing[1]), " in ", classifier_dir,
         ". Run classify_barcoding_folds(method = \"idtaxa\") first.", call. = FALSE)
  }
  pred <- list(species = utils::read.csv(files[1], stringsAsFactors = FALSE),
               genus = utils::read.csv(files[2], stringsAsFactors = FALSE))
  need <- c("genus_idtaxa", "species_idtaxa", "genus_confidence", "species_confidence")
  if (!all(need %in% names(pred$genus))) {
    stop("The IdTaxa tables of ", classifier_dir, " have no confidence columns: they were written ",
         "before Phase 6B. Run classify_barcoding_folds(method = \"idtaxa\") again.", call. = FALSE)
  }
  f_cn2 <- file.path(controls_dir, "TABLE_barcoding_cn2_queries_idtaxa.csv")
  outgroup <- if (file.exists(f_cn2)) utils::read.csv(f_cn2, stringsAsFactors = FALSE) else NULL
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  loci <- sort(unique(c(pred$species$locus, pred$genus$locus)), method = "radix")
  curve_list <- list()
  operating_rows <- list()
  for (l in loci) {
    aj_l <- if (is.null(outgroup)) NULL else outgroup[outgroup$locus == l, , drop = FALSE]
    aj_c <- .bc_comparable_aliens(aj_l)
    n_aj <- if (is.null(aj_c)) 0L else nrow(aj_c)
    g_l <- pred$genus[pred$genus$locus == l, , drop = FALSE]
    gc <- g_l$genus_confidence[!is.na(g_l$genus_confidence)]
    t_q <- if (length(gc) == 0L) NA_real_ else unname(stats::quantile(gc, q, type = quantile_type))
    extra <- c(pred$species$genus_confidence[pred$species$locus == l],
               pred$species$species_confidence[pred$species$locus == l],
               g_l$genus_confidence, g_l$species_confidence, t_q)
    for (sc in c("species", "genus")) {
      p <- pred[[sc]][pred[[sc]]$locus == l, , drop = FALSE]
      if (nrow(p) == 0L) next
      curve_list[[length(curve_list) + 1L]] <- .bc_confidence_curve(p, sc, outgroup = aj_l, extra = extra)
      for (rule in c("quantile", "default_60")) {
        t_op <- if (rule == "quantile") t_q else 60
        r <- if (is.na(t_op)) NULL else .bc_confidence_row(p, sc, aj_c, t_op)
        n_rej <- if (is.null(r) || n_aj == 0L) NA_integer_ else r$outgroup_rejected
        ci <- if (is.na(n_rej)) c(NA_real_, NA_real_) else .bc_clopper_pearson(n_rej, n_aj)
        counts <- if (is.null(r)) stats::setNames(as.list(rep(NA_integer_, 5)), .bc_outcome_levels()) else
          as.list(r[, .bc_outcome_levels()])
        operating_rows[[length(operating_rows) + 1L]] <- cbind(
          data.frame(locus = l, scheme = sc, rule = rule, threshold = t_op,
                     q = if (rule == "quantile") q else NA_real_,
                     quantile_type = if (rule == "quantile") as.integer(quantile_type) else NA_integer_,
                     queries = nrow(p), stringsAsFactors = FALSE),
          as.data.frame(counts),
          data.frame(outgroup_comparable = n_aj, outgroup_rejected = as.integer(n_rej),
                     rejection = if (is.na(n_rej)) NA_real_ else n_rej / n_aj,
                     rejection_ci_low = ci[1], rejection_ci_high = ci[2],
                     rejection_status = if (is.na(n_rej)) "not measured" else "measured",
                     stringsAsFactors = FALSE))
      }
    }
  }
  curve_tab <- do.call(rbind, curve_list)
  oper_tab <- do.call(rbind, operating_rows)
  rownames(curve_tab) <- NULL
  rownames(oper_tab) <- NULL
  utils::write.csv(curve_tab, file.path(output_dir, "TABLE_barcoding_threshold_curve_idtaxa.csv"), row.names = FALSE)
  utils::write.csv(oper_tab, file.path(output_dir, "TABLE_barcoding_threshold_operating_idtaxa.csv"), row.names = FALSE)
  if (isTRUE(figures)) {
    for (l in unique(curve_tab$locus)) {
      o_l <- oper_tab[oper_tab$locus == l & oper_tab$rule == "quantile", , drop = FALSE]
      g <- .bc_threshold_plot(curve_tab[curve_tab$locus == l, , drop = FALSE], o_l, l, method = "idtaxa")
      .bc_gap_save(g, file.path(output_dir, paste0("FIG_barcoding_threshold_idtaxa_", l)), width = 180, height = 110)
    }
  }
  .bc_threshold_banner(oper_tab[oper_tab$rule == "quantile", , drop = FALSE], output_dir)
  invisible(list(curve = curve_tab, operating = oper_tab))
}

#' Figure of one locus: what the legitimate queries become and how many alien queries are rejected
#'
#' Editorial rule of the branch (acta of the figures of Phase 4): English, the figures in the
#' subtitle and read from the table, nothing written inside the plotting area.
#' @noRd
.bc_threshold_plot <- function(curve_tab, oper_tab, locus, method = "nn") {
  idt <- identical(method, "idtaxa")
  x_scale <- if (idt) 1 else 100
  col <- .bc_gap_colors()
  label <- c(species = "Scheme E (species in library)", genus = "Scheme G (species left out)")
  long_tab <- do.call(rbind, lapply(split(curve_tab, curve_tab$scheme), function(d) {
    series <- data.frame(
      scheme = label[[d$scheme[1]]], threshold = d$threshold * x_scale,
      Unassigned = d$unassigned / d$queries,
      `Wrong species named` = d$wrong_species / d$queries,
      `Outgroup rejected` = d$outgroup_rejected / d$outgroup_comparable,
      check.names = FALSE, stringsAsFactors = FALSE)
    stats::reshape(series, direction = "long", varying = c("Unassigned", "Wrong species named", "Outgroup rejected"),
                   v.names = "value", timevar = "series",
                   times = c("Unassigned", "Wrong species named", "Outgroup rejected"), idvar = c("scheme", "threshold"))
  }))
  long_tab <- long_tab[!is.na(long_tab$value), , drop = FALSE]
  t_op <- oper_tab$threshold[1]
  o_g <- oper_tab[oper_tab$scheme == "genus", , drop = FALSE]
  measured <- nrow(o_g) == 1L && o_g$rejection_status == "measured"
  subtitle_text <- paste0(
    if (idt) paste0("Confidence threshold ", format(round(t_op, 1), nsmall = 1), " (quantile ", oper_tab$q[1],
                    " of scheme G genus confidence). ")
    else paste0("Threshold ", format(round(t_op * 100, 2), nsmall = 2), " % (quantile ", oper_tab$q[1], " of scheme G). "),
    if (nrow(o_g) == 1L) paste0("Scheme G unassigned: ", o_g$unassigned, " of ", o_g$queries, ". ") else "",
    if (measured) paste0("Outgroup rejected: ", o_g$outgroup_rejected, " of ", o_g$outgroup_comparable, ".")
    else "Outgroup rejection not measured.")
  ggplot2::ggplot(long_tab, ggplot2::aes(x = .data$threshold, y = .data$value, colour = .data$series)) +
    ggplot2::geom_step(linewidth = 0.5, direction = "hv") +
    ggplot2::geom_vline(xintercept = t_op * x_scale, linetype = "dashed", colour = col[["neutral"]], linewidth = 0.4) +
    ggplot2::facet_wrap(~ .data$scheme, nrow = 1) +
    ggplot2::scale_colour_manual(values = c(Unassigned = col[["intra"]], `Wrong species named` = col[["inter"]],
                                            `Outgroup rejected` = col[["highlight"]])) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::labs(title = locus, subtitle = subtitle_text,
                  x = if (idt) "IdTaxa confidence threshold" else "Distance to the nearest neighbour (% divergence)",
                  y = "Proportion of queries",
                  colour = NULL, caption = "Dashed: operating threshold.") +
    .bc_gap_theme() +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 7, colour = "grey30", hjust = 0))
}

#' Closing banner of step 9
#' @noRd
.bc_threshold_banner <- function(oper_tab, output_dir) {
  g <- oper_tab[oper_tab$scheme == "genus", , drop = FALSE]
  cat("\n====================================================\n")
  cat("  Barcoding Remoteness Threshold Complete \U0001f335\n")
  cat("====================================================\n")
  cat("  Loci:                  ", length(unique(oper_tab$locus)), "\n")
  cat("  Outgroup rejection measured in ", sum(g$rejection_status == "measured"), " loci\n", sep = "")
  cat("  Output directory:      ", output_dir, "\n")
  cat("====================================================\n\n")
  invisible(TRUE)
}
