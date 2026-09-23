# ------------------------------------------------------------------------------
# Molecular diagnostic branch (11_barcoding/), barcode gap. PhyloCactus 0.5.0, Phase 4.
#
# Rule of the validation plan (sec. 6.1), fixed before any data was seen: the candidate threshold of
# a locus is the 95th percentile of its intraspecific distribution, published next to the 5th
# percentile of the interspecific nearest-neighbour distribution. When the first is above the second
# there is no gap in that locus and it is declared; no other percentile is looked for. Both distance
# models, `raw` and K80, are computed and reported; neither is chosen afterwards for separating
# better. Pairs without a value are counted and left out of the percentiles, never imputed.
# ------------------------------------------------------------------------------

#' Distance matrix of one locus
#' @param dna `DNAbin` of the locus.
#' @param model Character. Model of [ape::dist.dna()].
#' @noRd
.bc_distance_matrix <- function(dna, model) {
  as.matrix(ape::dist.dna(dna, model = model, pairwise.deletion = TRUE))
}

#' Comparable positions of every pair of the locus
#'
#' Counts, for every pair, the positions where both sequences carry A, C, G or T. A distance resting
#' on too few of them says nothing, and the fourth amendment of the validation plan (E13, 2026-09-22)
#' takes those pairs out of the distributions.
#' @param dna `DNAbin` of the locus.
#' @noRd
.bc_comparable_sites <- function(dna) {
  x <- tolower(as.character(as.matrix(dna)))
  m <- matrix(x %in% c("a", "c", "g", "t"), nrow = nrow(x), ncol = ncol(x))
  out <- m %*% t(m)
  dimnames(out) <- list(rownames(x), rownames(x))
  out
}

#' Intraspecific and nearest-neighbour distributions of one locus, and the gap decision
#'
#' @param dmat Distance matrix with `sid` in its dimnames.
#' @param species Character vector of species named by `sid`.
#' @param locus,model Character. Copied into every row.
#' @return List with `intra` (one row per pair with a value), `inter` (one row per sequence, with
#'   `NA` when no neighbour of another species has a value) and `resumen`.
#' @noRd
.compute_barcode_gap <- function(dmat, species, locus, model, comparable = NULL, min_comparable = 100L) {
  sids <- rownames(dmat)
  species <- species[sids]
  genus <- .bc_genus(unname(species))
  names(genus) <- names(species)

  # E13: a pair resting on fewer than `min_comparable` positions has no value
  short <- matrix(FALSE, nrow(dmat), ncol(dmat), dimnames = dimnames(dmat))
  if (!is.null(comparable)) {
    cmp <- comparable[sids, sids, drop = FALSE]
    short <- !is.na(cmp) & cmp < min_comparable
    dmat[short] <- NA_real_
  }
  same_sp <- outer(species, species, "==")
  upper <- upper.tri(short)
  n_short_intra <- sum(short & same_sp & upper)
  n_short_inter <- sum(short & !same_sp & upper)
  empty_intra <- data.frame(locus = character(0), modelo = character(0), especie = character(0),
                            sid_a = character(0), sid_b = character(0), distancia = numeric(0),
                            stringsAsFactors = FALSE)

  n_intra_na <- 0L
  intra_rows <- list()
  for (sp in sort(unique(unname(species)), method = "radix")) {
    ids <- sort(names(species)[species == sp], method = "radix")
    if (length(ids) < 2L) next
    cb <- utils::combn(ids, 2L)
    dd <- dmat[cbind(cb[1, ], cb[2, ])]
    n_intra_na <- n_intra_na + sum(is.na(dd))
    keep <- !is.na(dd)
    if (any(keep)) {
      intra_rows[[sp]] <- data.frame(locus = locus, modelo = model, especie = sp,
                                     sid_a = cb[1, keep], sid_b = cb[2, keep],
                                     distancia = unname(dd[keep]), stringsAsFactors = FALSE)
    }
  }
  intra <- if (length(intra_rows) > 0) do.call(rbind, intra_rows) else empty_intra
  rownames(intra) <- NULL

  inter <- do.call(rbind, lapply(sids, function(s) {
    other <- names(species)[species != species[[s]]]
    if (length(other) == 0) return(NULL)
    d <- dmat[s, other]
    row <- data.frame(locus = locus, modelo = model, sid = s, especie = unname(species[[s]]),
                      sid_vecino = NA_character_, especie_vecino = NA_character_,
                      distancia = NA_real_, stringsAsFactors = FALSE)
    if (all(is.na(d))) return(row)
    nb <- sort(other[!is.na(d) & d == min(d, na.rm = TRUE)], method = "radix")[1]
    row$sid_vecino <- nb
    row$especie_vecino <- unname(species[[nb]])
    row$distancia <- unname(d[nb])
    row
  }))
  if (is.null(inter)) {
    inter <- data.frame(locus = character(0), modelo = character(0), sid = character(0),
                        especie = character(0), sid_vecino = character(0),
                        especie_vecino = character(0), distancia = numeric(0), stringsAsFactors = FALSE)
  }
  rownames(inter) <- NULL

  q <- function(x, p) if (length(x) == 0) NA_real_ else stats::quantile(x, p, type = 7, names = FALSE)
  med <- function(x) if (length(x) == 0) NA_real_ else stats::median(x)
  iv <- intra$distancia
  nv <- inter$distancia[!is.na(inter$distancia)]
  p95 <- q(iv, 0.95)
  p5 <- q(nv, 0.05)

  # E14: three descriptive columns, because a 5th percentile of zero hides how much the two
  # distributions overlap
  cero <- !is.na(inter$distancia) & inter$distancia == 0
  congenere <- cero & !is.na(inter$especie_vecino) &
    .bc_genus(inter$especie) == .bc_genus(inter$especie_vecino)

  resumen <- data.frame(
    locus = locus, modelo = model,
    pares_intra = nrow(intra), mediana_intra = med(iv), p95_intra = p95,
    consultas_inter = nrow(inter), mediana_inter = med(nv), p5_inter = p5,
    pares_intra_sin_valor = as.integer(n_intra_na),
    consultas_inter_sin_valor = as.integer(sum(is.na(inter$distancia))),
    minimo_comparables = as.integer(min_comparable),
    pares_intra_solapamiento_corto = as.integer(n_short_intra),
    pares_inter_solapamiento_corto = as.integer(n_short_inter),
    consultas_inter_cero = as.integer(sum(cero)),
    prop_inter_cero = if (length(nv) == 0) NA_real_ else sum(cero) / length(nv),
    inter_cero_congenere = as.integer(sum(congenere)),
    inter_cero_otro_genero = as.integer(sum(cero & !congenere)),
    prop_inter_bajo_umbral = if (length(nv) == 0 || is.na(p95)) NA_real_ else
      mean(inter$distancia < p95, na.rm = TRUE),
    umbral_candidato = p95,
    hay_gap = isTRUE(p95 < p5),
    stringsAsFactors = FALSE
  )
  list(intra = intra, inter = inter, resumen = resumen)
}

#' Editorial palette of the branch figures
#'
#' The two series of the barcode gap, in the Dark2 hues the package already uses for its figures
#' (`R/validation.R`). Checked for colour vision deficiency: the adjacent pair keeps a separation of
#' 11.6 in deuteranopia, above the threshold of 8, and both hold contrast over a light surface.
#' @noRd
.bc_gap_colors <- function() {
  c(intra = "#1b9e77", inter = "#d95f02", resaltado = "#7570b3", neutro = "grey35")
}

#' Editorial theme of the branch figures
#'
#' The style the package uses elsewhere: `theme_minimal`, bold titles, black axis text, no minor
#' grid and the legend at the bottom. Sizes are set for a double column figure (180 mm), where the
#' base size lands between 8 and 10 points once printed.
#' @noRd
.bc_gap_theme <- function(base_size = 9) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1, hjust = 0,
                                         margin = ggplot2::margin(t = 0, r = 0, b = 2, l = 14)),
      plot.subtitle = ggplot2::element_text(size = base_size - 1, colour = "grey25",
                                            margin = ggplot2::margin(t = 0, r = 0, b = 4, l = 14)),
      plot.title.position = "plot",
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1, colour = "black"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey88"),
      legend.position = "bottom",
      legend.key.size = grid::unit(3.5, "mm"),
      legend.text = ggplot2::element_text(size = base_size - 1),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 2)
    )
}

#' Largest intraspecific distance of every species against its nearest neighbour of another species
#'
#' The per-species view of the barcode gap: a species without gap is one whose nearest neighbour of
#' another species is no farther than its own most distant pair. Species without replica have no
#' intraspecific distance and do not appear.
#' @noRd
.bc_gap_species_summary <- function(intra, inter) {
  if (is.null(intra) || nrow(intra) == 0) {
    return(data.frame(locus = character(0), modelo = character(0), especie = character(0),
                      intra_max = numeric(0), inter_min = numeric(0), pares_intra = integer(0),
                      sin_gap = logical(0), stringsAsFactors = FALSE))
  }
  key <- paste(intra$modelo, intra$especie, sep = "\r")
  out <- do.call(rbind, lapply(sort(unique(key), method = "radix"), function(k) {
    d <- intra[key == k, , drop = FALSE]
    m <- d$modelo[1]
    sp <- d$especie[1]
    vecinos <- inter$distancia[inter$modelo == m & inter$especie == sp]
    vecinos <- vecinos[!is.na(vecinos)]
    intra_max <- max(d$distancia)
    inter_min <- if (length(vecinos) == 0) NA_real_ else min(vecinos)
    data.frame(locus = d$locus[1], modelo = m, especie = sp, intra_max = intra_max,
               inter_min = inter_min, pares_intra = nrow(d),
               sin_gap = isTRUE(inter_min <= intra_max), stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' A distance as a short percentage, for the labels of the figures
#' @noRd
.bc_gap_pct <- function(x) {
  if (length(x) == 0 || is.na(x)) return("NA")
  paste0(signif(x * 100, 3), "%")
}

#' Histogram of the two distributions of one model, as proportions of each distribution
#'
#' Editorial rule of the review of 2026-09-22: every visible string is English and short, and
#' nothing is written inside the plotting area. The two percentiles are read in the subtitle; the
#' dashed and dotted lines only mark where they fall.
#' @noRd
.bc_gap_distribution_panel <- function(intra, inter, r, locus, bins = 40L) {
  col <- .bc_gap_colors()
  etiquetas <- c(intra = "Intraspecific", inter = "Interspecific, nearest neighbour")
  x_intra <- intra$distancia * 100
  x_inter <- inter$distancia[!is.na(inter$distancia)] * 100
  rango <- range(c(x_intra, x_inter, 0), na.rm = TRUE)
  if (!is.finite(diff(rango)) || diff(rango) == 0) rango <- c(0, max(1, rango[2]))
  brk <- seq(rango[1], rango[2], length.out = bins + 1L)
  ancho <- diff(brk)[1]
  centros <- (brk[-1] + brk[-length(brk)]) / 2
  cuenta <- function(x, etiqueta) {
    n <- tabulate(cut(x, breaks = brk, include.lowest = TRUE, labels = FALSE), nbins = bins)
    data.frame(centro = centros, proporcion = if (sum(n) > 0) n / sum(n) else as.numeric(n),
               distribucion = etiqueta, stringsAsFactors = FALSE)
  }
  df <- rbind(cuenta(x_intra, etiquetas[["intra"]]), cuenta(x_inter, etiquetas[["inter"]]))
  df$distribucion <- factor(df$distribucion, levels = unname(etiquetas))
  y_max <- max(df$proporcion, na.rm = TRUE)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$centro, y = .data$proporcion, fill = .data$distribucion)) +
    ggplot2::geom_col(width = ancho, alpha = 0.6, position = "identity") +
    ggplot2::geom_vline(xintercept = r$p95_intra * 100, colour = col[["intra"]], linewidth = 0.5,
                        linetype = "dashed") +
    ggplot2::geom_vline(xintercept = r$p5_inter * 100, colour = col[["inter"]], linewidth = 0.5,
                        linetype = "dotted") +
    ggplot2::scale_fill_manual(values = stats::setNames(c(col[["intra"]], col[["inter"]]), unname(etiquetas)),
                               drop = FALSE) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = c(0, 0),
                                limits = c(0, y_max * 1.05)) +
    ggplot2::labs(
      title = paste0(locus, " (", r$modelo, ")"),
      subtitle = paste0("P95 intra ", .bc_gap_pct(r$p95_intra), " | P5 inter ", .bc_gap_pct(r$p5_inter),
                        " | ", if (isTRUE(r$hay_gap)) "gap" else "no gap"),
      x = "Distance (% divergence)", y = "Share of distribution", fill = NULL
    ) +
    .bc_gap_theme()
}

#' Per-species panel: most distant conspecific pair against the nearest neighbour of another species
#'
#' The 1:1 line is the gap: a species below it has a neighbour of another species closer than its
#' own most distant pair. The reading of the line goes in the caption, not inside the panel.
#' @noRd
.bc_gap_species_panel <- function(sp, r, locus) {
  col <- .bc_gap_colors()
  d <- sp[sp$modelo == r$modelo & !is.na(sp$inter_min), , drop = FALSE]
  techo <- max(c(d$intra_max, d$inter_min, 0), na.rm = TRUE) * 100
  if (!is.finite(techo) || techo == 0) techo <- 1
  # Room below zero: most species have their neighbour of another species at distance 0 and their
  # point would be cut in half against the axis
  lim <- c(-0.03 * techo, techo * 1.05)
  sin_gap <- if (nrow(d) == 0) NA_real_ else mean(d$sin_gap)

  ggplot2::ggplot(d, ggplot2::aes(x = .data$intra_max * 100, y = .data$inter_min * 100)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_point(shape = 21, size = 1.6, stroke = 0.3, colour = "white",
                        fill = col[["neutro"]], alpha = 0.75) +
    ggplot2::coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
    ggplot2::labs(
      title = paste0("Per species (", r$modelo, ")"),
      subtitle = if (is.na(sin_gap)) "No replicated species" else
        paste0(nrow(d), " replicated species | ", round(100 * sin_gap), "% without gap"),
      x = "Widest conspecific pair (%)", y = "Nearest other species (%)"
    ) +
    .bc_gap_theme()
}

#' Figure of one locus: the two distributions and the per-species view, one column per model
#' @noRd
.bc_gap_plot <- function(intra, inter, resumen, locus) {
  sp <- .bc_gap_species_summary(intra, inter)
  paneles <- list()
  for (m in unique(resumen$modelo)) {
    r <- resumen[resumen$modelo == m, , drop = FALSE][1, ]
    paneles[[length(paneles) + 1L]] <- .bc_gap_distribution_panel(
      intra[intra$modelo == m, , drop = FALSE], inter[inter$modelo == m, , drop = FALSE], r, locus)
  }
  for (m in unique(resumen$modelo)) {
    r <- resumen[resumen$modelo == m, , drop = FALSE][1, ]
    paneles[[length(paneles) + 1L]] <- .bc_gap_species_panel(sp, r, locus)
  }
  minimo <- if (!is.null(resumen$minimo_comparables)) resumen$minimo_comparables[1] else 100L
  patchwork::wrap_plots(paneles, ncol = length(unique(resumen$modelo))) +
    patchwork::plot_annotation(
      tag_levels = "A",
      caption = paste0("Dashed: P95 intra. Dotted: P5 inter. Below the 1:1 line: no gap. Pairs under ",
                       minimo, " shared sites excluded."),
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(size = 7, colour = "grey30", hjust = 0))
    ) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

#' Overview figure: the two percentiles of every locus, and how much the distributions overlap
#' @noRd
.bc_gap_overview_plot <- function(summary_tab) {
  col <- .bc_gap_colors()
  d <- summary_tab
  orden <- d$locus[d$modelo == d$modelo[1]][order(d$p95_intra[d$modelo == d$modelo[1]])]
  d$locus <- factor(d$locus, levels = orden)
  # What the figure says is counted here, never asserted: a locus-model pair has a gap when its
  # P95 intraspecific falls below its P5 interspecific (column hay_gap of the summary, when present)
  con_gap <- if (!is.null(d$hay_gap)) !is.na(d$hay_gap) & as.logical(d$hay_gap) else d$p95_intra < d$p5_inter
  loci <- unique(as.character(d$locus))
  loci_con_gap <- unique(as.character(d$locus[con_gap]))
  titulo <- if (length(loci_con_gap) == 0) "No locus separates species by distance" else
    paste0(length(loci_con_gap), " of ", length(loci), " loci separate species by distance")
  subtitulo <- paste0("P95 intra above P5 inter in ", sum(!con_gap), " of ", length(con_gap),
                      " locus-model pairs")
  largo <- rbind(
    data.frame(locus = d$locus, modelo = d$modelo, valor = d$p95_intra * 100,
               serie = "P95 intraspecific", stringsAsFactors = FALSE),
    data.frame(locus = d$locus, modelo = d$modelo, valor = d$p5_inter * 100,
               serie = "P5 interspecific", stringsAsFactors = FALSE)
  )
  p1 <- ggplot2::ggplot(largo, ggplot2::aes(x = .data$valor, y = .data$locus)) +
    ggplot2::geom_line(ggplot2::aes(group = interaction(.data$locus, .data$modelo)),
                       colour = "grey70", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$serie), size = 1.8) +
    ggplot2::facet_wrap(~ .data$modelo, nrow = 1) +
    ggplot2::scale_colour_manual(values = stats::setNames(c(col[["intra"]], col[["inter"]]),
                                                          c("P95 intraspecific", "P5 interspecific"))) +
    ggplot2::labs(title = titulo, subtitle = subtitulo,
                  x = "Distance (% divergence)", y = NULL, colour = NULL) +
    .bc_gap_theme()

  z <- d[d$modelo == d$modelo[1], , drop = FALSE]
  p2 <- ggplot2::ggplot(z, ggplot2::aes(x = .data$prop_inter_cero * 100, y = .data$locus)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$prop_inter_cero * 100,
                                       y = .data$locus, yend = .data$locus),
                          colour = "grey80", linewidth = 0.5) +
    ggplot2::geom_point(colour = col[["inter"]], size = 1.8) +
    ggplot2::scale_x_continuous(limits = c(0, 100)) +
    ggplot2::labs(title = "Identical to another species",
                  subtitle = "Queries with a neighbour at distance 0",
                  x = "% of queries", y = NULL) +
    .bc_gap_theme()

  patchwork::wrap_plots(p1, p2, widths = c(2, 1)) +
    patchwork::plot_annotation(tag_levels = "A",
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(size = 7, colour = "grey30", hjust = 0)))
}

#' Barcode Gap of the Molecular Diagnostic Branch
#'
#' Step 6 of the branch. For every locus of the library, computes the intraspecific distances by
#' pairs and, for every sequence, the distance to its nearest neighbour of another species, with the
#' `raw` and K80 models. Derives the candidate threshold with the rule of the validation plan: the
#' 95th percentile of the intraspecific distribution, published next to the 5th percentile of the
#' interspecific one. When the first is above the second the locus is declared to have no gap.
#'
#' @param library_dir Character. Output directory of [finalize_barcoding_library()].
#' @param output_dir Character. Destination; refused inside a directory of the phylogeny.
#' @param models Character vector. Distance models of [ape::dist.dna()]; both are reported.
#' @param min_comparable Integer. Minimum number of positions where both sequences carry A, C, G or
#'   T for a pair to have a value (E13 of the validation plan, 2026-09-22). Pairs below it are counted
#'   apart and left out of every distribution and percentile; they are never imputed.
#' @param paralog_loci Character vector. Loci flagged in `posible_paralogo`.
#' @param figures Logical. Write one figure per locus, with one panel per model.
#' @return Invisibly, a list with `summary` and `distances` (by locus). Writes
#'   `TABLE_barcoding_gap_summary.csv`, `TABLE_barcoding_distances_<locus>.csv` and
#'   `FIG_barcoding_gap_<locus>.png`.
#' @examples
#' \dontrun{
#' analyze_barcode_gap(
#'   library_dir = "11_barcoding/4_library",
#'   output_dir = "11_barcoding/6_gap"
#' )
#' }
#' @export
analyze_barcode_gap <- function(library_dir = file.path("11_barcoding", "4_library"),
                                output_dir = file.path("11_barcoding", "6_gap"),
                                models = c("raw", "K80"),
                                min_comparable = 100L,
                                paralog_loci = "pepC_like",
                                figures = TRUE) {
  .bc_assert_output_dir(output_dir)
  files <- list.files(library_dir, pattern = "^LIB_.*\\.fasta$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No LIB_<locus>.fasta in ", library_dir, ". Run finalize_barcoding_library() first.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  summaries <- list()
  distances <- list()
  for (f in sort(files, method = "radix")) {
    locus <- sub("^LIB_(.*)\\.fasta$", "\\1", basename(f))
    dna <- ape::read.dna(f, format = "fasta")
    h <- .bc_parse_header(labels(dna))
    rownames(dna) <- h$sid
    species <- stats::setNames(h$species, h$sid)

    comparable <- .bc_comparable_sites(dna)
    per_model <- lapply(models, function(m) {
      .compute_barcode_gap(.bc_distance_matrix(dna, m), species, locus, m,
                           comparable = comparable, min_comparable = min_comparable)
    })
    intra <- do.call(rbind, lapply(per_model, function(x) x$intra))
    inter <- do.call(rbind, lapply(per_model, function(x) x$inter))
    resumen <- do.call(rbind, lapply(per_model, function(x) x$resumen))
    summaries[[locus]] <- resumen

    tab <- rbind(
      data.frame(modelo = intra$modelo, tipo = "intra", sid = intra$sid_a, especie = intra$especie,
                 sid_comparado = intra$sid_b, especie_comparada = intra$especie,
                 distancia = intra$distancia, stringsAsFactors = FALSE),
      data.frame(modelo = inter$modelo, tipo = "inter_vecino", sid = inter$sid, especie = inter$especie,
                 sid_comparado = inter$sid_vecino, especie_comparada = inter$especie_vecino,
                 distancia = inter$distancia, stringsAsFactors = FALSE)
    )
    rownames(tab) <- NULL
    distances[[locus]] <- tab
    utils::write.csv(tab, file.path(output_dir, paste0("TABLE_barcoding_distances_", locus, ".csv")),
                     row.names = FALSE)

    if (isTRUE(figures)) {
      fig <- .bc_gap_plot(intra, inter, resumen, locus)
      .bc_gap_save(fig, file.path(output_dir, paste0("FIG_barcoding_gap_", locus)),
                   width = 180, height = 155)
    }
    for (i in seq_len(nrow(resumen))) {
      message(sprintf("Locus '%s', model %s: %d intraspecific pairs, %d queries; P95 intra = %s, P5 inter = %s (%s); %s%% of queries have a neighbour of another species at distance 0; pairs left out for short overlap: %d intra, %d inter.",
                      locus, resumen$modelo[i], resumen$pares_intra[i], resumen$consultas_inter[i],
                      signif(resumen$p95_intra[i], 3), signif(resumen$p5_inter[i], 3),
                      if (isTRUE(resumen$hay_gap[i])) "gap" else "no gap",
                      signif(100 * resumen$prop_inter_cero[i], 3),
                      resumen$pares_intra_solapamiento_corto[i], resumen$pares_inter_solapamiento_corto[i]))
    }
  }

  summary_tab <- do.call(rbind, summaries)
  rownames(summary_tab) <- NULL
  summary_tab <- .bc_flag_provisional(.bc_flag_paralog(summary_tab, paralog_loci))
  utils::write.csv(summary_tab, file.path(output_dir, "TABLE_barcoding_gap_summary.csv"), row.names = FALSE)

  if (isTRUE(figures) && nrow(summary_tab) > 0) {
    .bc_gap_save(.bc_gap_overview_plot(summary_tab), file.path(output_dir, "FIG_barcoding_gap_overview"),
                 width = 180, height = 110)
  }

  invisible(list(summary = summary_tab, distances = distances))
}

#' Write a figure as vector (PDF, for the manuscript) and raster (PNG, to read), at column width
#' @param plot A `ggplot` or `patchwork` object.
#' @param path_base Character. Path without extension.
#' @param width,height Numeric. Millimetres; 180 mm is the usual double column.
#' @noRd
.bc_gap_save <- function(plot, path_base, width = 180, height = 155) {
  ggplot2::ggsave(paste0(path_base, ".pdf"), plot, width = width, height = height, units = "mm")
  ggplot2::ggsave(paste0(path_base, ".png"), plot, width = width, height = height, units = "mm", dpi = 300)
  invisible(c(paste0(path_base, ".pdf"), paste0(path_base, ".png")))
}
