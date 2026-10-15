# ------------------------------------------------------------------------------
# Molecular diagnostic branch (11_barcoding/), figures of an identification drawn from the tables of
# its report. PhyloCactus 0.5.0.
#
# The report of an identification (report_barcoding_identification()) writes its answer per locus and
# the sampling of the named genera as CSV files. These functions draw the same figures from those
# tables, so that a figure can be drawn again without the library, as Tutorial 6 does with the tables
# distributed in extdata.
#
# Governing documents: 01_minutes_and_phases/04_barcoding_0.5.0/ of the audit repository, proposal of
# Tutorial 6 of 2026-10-14 (decisions T6-5 and F2 to F4 of BMM).
# ------------------------------------------------------------------------------

#' A table given as a data frame or as the path of a CSV file, with the columns a figure needs
#' @noRd
.bc_plot_table <- function(x, need, what) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop("Not found: ", x, ".", call. = FALSE)
    x <- utils::read.csv(x, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(x)) stop("The ", what, " must be a data frame or the path of a CSV file.", call. = FALSE)
  miss <- setdiff(need, names(x))
  if (length(miss)) stop("The ", what, " lacks the column", if (length(miss) > 1L) "s", " ", paste(miss, collapse = ", "), ".",
                         call. = FALSE)
  x
}

#' The rows of loci found in the query (the others have no region to show)
#' @noRd
.bc_plot_found <- function(t) {
  if (!"reason" %in% names(t)) return(t)
  t[!t$reason %in% c("no_overlap", "assembly_failed", "single_species_library"), , drop = FALSE]
}

#' Draw the answer of an identification
#'
#' Draws, from the tables the identification report writes, the answer per locus of an
#' identification: for each locus found, the IdTaxa confidence at genus rank (striped bar) and at
#' species rank (solid bar), coloured by state (1 species, 2 genus, 3 no answer), with the threshold
#' as a dashed line. When the sampling table is given, the sampling of the named genera in the library
#' (sequences per species and locus) is drawn below. The figure is the one of section 2 and section 6
#' of the report ([report_barcoding_identification()]), drawn without the library.
#'
#' @param answer Data frame or path. `REPORT_<run_name>_answer.csv`, or a table with its columns
#'   `query`, `locus`, `compartment`, `state`, `predicted_species`, `predicted_genus`,
#'   `genus_confidence` and `species_confidence` (`reason`, if present, leaves out the loci not
#'   found).
#' @param sampling Data frame, path or `NULL`. `REPORT_<run_name>_sampling.csv`, or a table with its
#'   columns `genus`, `species`, `locus` and `sequences`; restricted to the genera named in `answer`.
#' @param threshold Numeric. The IdTaxa threshold drawn as a dashed line. Defaults to 60.
#' @param max_species Integer. Above this number of species in the sampling, the species not named
#'   are gathered into one row per genus, whose cells count the species with sequences at each locus.
#'   Defaults to 40.
#' @return A `ggplot` object, or a `patchwork` object with the sampling below when `sampling` is
#'   given. Nothing is written.
#' @seealso [plot_barcoding_identification_summary()], [report_barcoding_identification()].
#' @examples
#' \dontrun{
#' plot_barcoding_identification("11_barcoding/10_identify/REPORT_MN517611.1_answer.csv",
#'                               sampling = "11_barcoding/10_identify/REPORT_MN517611.1_sampling.csv")
#' }
#' @export
plot_barcoding_identification <- function(answer, sampling = NULL, threshold = 60, max_species = 40L) {
  a <- .bc_plot_table(answer, c("query", "locus", "compartment", "state", "predicted_species", "predicted_genus",
                                "genus_confidence", "species_confidence"), "answer table")
  a$state <- suppressWarnings(as.integer(a$state))
  f <- .bc_plot_found(a)
  p <- .bc_report_answer_plot(if (nrow(f)) f else a, threshold)
  if (is.null(sampling)) return(p)
  s <- .bc_plot_table(sampling, c("genus", "species", "locus", "sequences"), "sampling table")
  named_genera <- unique(stats::na.omit(f$predicted_genus[f$state %in% 1:2]))
  named_species <- unique(stats::na.omit(f$predicted_species[f$state %in% 1L]))
  if (length(named_genera) && any(s$genus %in% named_genera)) s <- s[s$genus %in% named_genera, , drop = FALSE]
  s <- .bc_plot_gather_sampling(s, named_species, max_species)
  n_rows <- length(unique(s$species))
  patchwork::wrap_plots(p, .bc_report_sampling_plot(s, named_species), ncol = 1, heights = c(1, max(0.5, n_rows / 20)))
}

#' Sampling of many species: the species not named gathered into one row per genus
#' @noRd
.bc_plot_gather_sampling <- function(s, named_species, max_species) {
  s$sequences <- suppressWarnings(as.integer(s$sequences))
  if (length(unique(s$species)) <= max_species) return(s)
  keep <- s[s$species %in% named_species, , drop = FALSE]
  rest <- s[!s$species %in% named_species, , drop = FALSE]
  if (!nrow(rest)) return(keep)
  g <- do.call(rbind, lapply(split(rest, list(rest$genus, rest$locus), drop = TRUE), function(x) {
    n <- length(unique(x$species))
    data.frame(genus = x$genus[1], species = paste0(x$genus[1], ": ", n, " other species"), locus = x$locus[1],
               sequences = length(unique(x$species[x$sequences > 0])), stringsAsFactors = FALSE)
  }))
  rbind(keep[, names(g)], g)
}

#' The answer of one identification in words: one species, one genus, a conflict or no answer
#' @noRd
.bc_plot_answer_words <- function(f) {
  sp <- unique(stats::na.omit(f$predicted_species[f$state %in% 1L]))
  ge <- unique(stats::na.omit(f$predicted_genus[f$state %in% 1:2]))
  if (length(sp) == 1L) return(paste0(gsub("_", " ", sp), " (species)"))
  if (length(sp) > 1L || length(ge) > 1L) return("conflict")
  if (length(ge) == 1L) return(paste0(ge, " (genus)"))
  "no answer"
}

#' Draw the answers of several identifications side by side
#'
#' One row per identification and one column per locus found in any of them, each cell coloured by
#' the state of that locus (1 species, 2 genus, 3 no answer); a cross marks a locus named at genus rank
#' only by rule (reason `genus_only_locus`, argument `genus_only_loci` of
#' [identify_barcoding_query()]). The answer of each identification (one species, one genus, a
#' conflict, or no answer) is written next to its name. When a locus was found in more than one
#' region of a query (contigs of reads), the cell shows its most informative row.
#'
#' @param answers Named list. Each element is the answer table of one identification
#'   (`REPORT_<run_name>_answer.csv` or its path, as in [plot_barcoding_identification()]); the names
#'   label the rows.
#' @return A `ggplot` object; its data hold one row per identification and locus found, with the
#'   columns `example`, `locus`, `state`, `genus_only` and `answer`. Nothing is written.
#' @seealso [plot_barcoding_identification()].
#' @examples
#' \dontrun{
#' plot_barcoding_identification_summary(list(
#'   "MN517611.1" = "11_barcoding/10_identify/REPORT_MN517611.1_answer.csv",
#'   "LT745702.1" = "11_barcoding/10_identify/REPORT_LT745702.1_answer.csv"))
#' }
#' @export
plot_barcoding_identification_summary <- function(answers) {
  if (!is.list(answers) || is.data.frame(answers)) stop("answers must be a named list of answer tables.", call. = FALSE)
  if (is.null(names(answers)) || any(!nzchar(names(answers)))) {
    stop("answers must be a named list: the names label the rows.", call. = FALSE)
  }
  rows <- lapply(names(answers), function(n) {
    a <- .bc_plot_table(answers[[n]], c("locus", "state", "predicted_species", "predicted_genus", "species_confidence"),
                        paste0("answer table '", n, "'"))
    a$state <- suppressWarnings(as.integer(a$state))
    if (!"reason" %in% names(a)) a$reason <- NA_character_
    f <- .bc_plot_found(a)
    if (!nrow(f)) return(NULL)
    conf <- ifelse(is.na(f$species_confidence), -1, f$species_confidence)
    f <- f[order(f$locus, f$state, -conf, method = "radix"), , drop = FALSE]
    f <- f[!duplicated(f$locus), , drop = FALSE]
    data.frame(example = n, locus = f$locus, state = as.character(f$state),
               genus_only = f$reason %in% "genus_only_locus", answer = .bc_plot_answer_words(f), stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows)
  if (is.null(d)) stop("No locus was found in any of the identifications.", call. = FALSE)
  d$example <- factor(d$example, levels = rev(names(answers)))
  d$state <- factor(d$state, levels = c("1", "2", "3"))
  lab <- unique(d[, c("example", "answer")])
  ylab <- stats::setNames(paste0(lab$example, ": ", lab$answer), as.character(lab$example))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$locus, y = .data$example)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$state), colour = "white") +
    ggplot2::geom_point(data = d[d$genus_only, , drop = FALSE], shape = 4, size = 3) +
    ggplot2::scale_fill_manual(values = .bc_state_colours(), name = "State", drop = FALSE,
                               labels = c("1" = "1 species", "2" = "2 genus", "3" = "3 no answer")) +
    ggplot2::scale_y_discrete(labels = ylab) +
    ggplot2::labs(x = NULL, y = NULL, caption = "Cross: locus named at genus rank only (genus_only_locus).") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
}
