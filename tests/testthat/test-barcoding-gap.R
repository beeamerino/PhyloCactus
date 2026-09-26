# Tests of the barcode gap of the molecular diagnostic branch (11_barcoding/6_gap/).
# Written before the functions they test. Phase 4 of PhyloCactus 0.5.0.
#
# Rule of the validation plan (sec. 6.1), fixed before any data was seen: the candidate threshold of
# a locus is the 95th percentile of its intraspecific distribution, published next to the 5th
# percentile of the interspecific nearest-neighbour distribution. When the first is above the second
# there is no gap in that locus and it is declared; no other percentile is looked for.

# Distance matrix written by hand: two species with replica and one species with a single sequence.
.gap_matrix <- function(values) {
  sids <- c("A1.1", "A2.1", "B1.1", "B2.1", "C1.1")
  m <- matrix(NA_real_, 5, 5, dimnames = list(sids, sids))
  m[lower.tri(m)] <- values
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  diag(m) <- 0
  m
}

.gap_species <- function() {
  stats::setNames(c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_stricta", "Opuntia_stricta",
                    "Cereus_jamacaru"),
                  c("A1.1", "A2.1", "B1.1", "B2.1", "C1.1"))
}

test_that("the two distributions are built from the matrix, and a species without replica still gives a query", {
  # Order of the lower triangle: (A2,A1) (B1,A1) (B2,A1) (C1,A1) (B1,A2) (B2,A2) (C1,A2) (B2,B1) (C1,B1) (C1,B2)
  m <- .gap_matrix(c(0.01, 0.20, 0.22, 0.30, 0.21, 0.23, 0.31, 0.02, 0.25, 0.26))

  g <- .compute_barcode_gap(m, .gap_species(), locus = "matK", model = "raw")

  expect_setequal(paste(g$intra$sid_a, g$intra$sid_b), c("A1.1 A2.1", "B1.1 B2.1"))
  expect_equal(g$intra$distance[paste(g$intra$sid_a, g$intra$sid_b) == "A1.1 A2.1"], 0.01)
  expect_true(all(g$intra$locus == "matK" & g$intra$model == "raw"))
  # Nearest neighbour of another species, per sequence
  nn <- stats::setNames(g$inter$distance, g$inter$sid)
  expect_equal(unname(nn[c("A1.1", "A2.1", "B1.1", "B2.1", "C1.1")]), c(0.20, 0.21, 0.20, 0.22, 0.25))
  expect_equal(g$inter$nn_sid[g$inter$sid == "A1.1"], "B1.1")
  expect_equal(g$inter$nn_species[g$inter$sid == "C1.1"], "Opuntia_stricta")
  # The species with a single sequence gives no intraspecific pair, but is a query and can be a neighbour
  expect_false("Cereus_jamacaru" %in% g$intra$species)
  expect_true("C1.1" %in% g$inter$sid)
  expect_equal(g$summary_df$intra_pairs, 2L)
  expect_equal(g$summary_df$inter_queries, 5L)
})

test_that("the candidate threshold is the 95th intraspecific percentile, with a gap and without it", {
  sp <- .gap_species()
  con_gap <- .compute_barcode_gap(
    .gap_matrix(c(0.01, 0.20, 0.22, 0.30, 0.21, 0.23, 0.31, 0.02, 0.25, 0.26)), sp, "matK", "raw")
  intra <- con_gap$intra$distance
  inter <- con_gap$inter$distance
  expect_equal(con_gap$summary_df$candidate_threshold, stats::quantile(intra, 0.95, type = 7, names = FALSE))
  expect_equal(con_gap$summary_df$p5_inter, stats::quantile(inter, 0.05, type = 7, names = FALSE))
  expect_true(con_gap$summary_df$has_gap)

  # Overlap: the intraspecific distances are larger than the closest interspecific ones
  no_gap <- .compute_barcode_gap(
    .gap_matrix(c(0.30, 0.05, 0.22, 0.28, 0.06, 0.23, 0.29, 0.31, 0.25, 0.26)), sp, "matK", "raw")
  expect_false(no_gap$summary_df$has_gap)
  expect_equal(no_gap$summary_df$candidate_threshold,
               stats::quantile(no_gap$intra$distance, 0.95, type = 7, names = FALSE))
  expect_gt(no_gap$summary_df$candidate_threshold, no_gap$summary_df$p5_inter)
})

test_that("pairs without a value are declared, excluded from the percentiles and never imputed", {
  sp <- .gap_species()
  v <- c(NA, 0.20, 0.22, 0.30, 0.21, 0.23, 0.31, 0.02, 0.25, 0.26)  # the pair of Opuntia_robusta has no value
  g <- .compute_barcode_gap(.gap_matrix(v), sp, "matK", "K80")

  expect_equal(g$summary_df$intra_pairs, 1L)
  expect_equal(g$summary_df$intra_pairs_no_value, 1L)
  expect_false(any(is.na(g$intra$distance)))
  expect_false("A1.1 A2.1" %in% paste(g$intra$sid_a, g$intra$sid_b))
  # The percentile uses the values that exist, and the maximum is not imputed anywhere
  expect_equal(g$summary_df$p95_intra, 0.02)
  expect_equal(g$summary_df$candidate_threshold, 0.02)
  expect_equal(g$summary_df$inter_queries_no_value, 0L)
})

test_that("both models are reported for every locus, and neither is chosen for separating better", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir)
  set.seed(4L)
  base <- paste(sample(c("A", "C", "G", "T"), 200, replace = TRUE), collapse = "")
  vary <- function(s, k) { for (p in sample(seq_len(200), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1); s }
  seqs <- c(`Opuntia_robusta|A1.1` = vary(base, 2), `Opuntia_robusta|A2.1` = vary(base, 2),
            `Opuntia_stricta|B1.1` = vary(base, 30), `Opuntia_stricta|B2.1` = vary(base, 32),
            `Cereus_jamacaru|C1.1` = vary(base, 60))
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(seqs), file.path(lib_dir, "LIB_matK.fasta"))

  out <- suppressMessages(analyze_barcode_gap(library_dir = lib_dir, output_dir = file.path(tmp, "6_gap")))

  s <- utils::read.csv(file.path(tmp, "6_gap", "TABLE_barcoding_gap_summary.csv"), stringsAsFactors = FALSE)
  expect_setequal(s$model, c("raw", "K80"))
  expect_equal(nrow(s), 2L)
  expect_true(all(s$locus == "matK"))
  expect_true(all(s$candidate_threshold == s$p95_intra))
  expect_equal(nrow(out$summary), 2L)
})

test_that("step 6 writes the distances and the figure, flags the loci and names the step that is missing", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir)
  set.seed(5L)
  base <- paste(sample(c("A", "C", "G", "T"), 200, replace = TRUE), collapse = "")
  vary <- function(s, k) { for (p in sample(seq_len(200), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1); s }
  for (l in c("matK", "pepC_like")) {
    seqs <- c(`Opuntia_robusta|A1.1` = vary(base, 2), `Opuntia_robusta|A2.1` = vary(base, 3),
              `Opuntia_stricta|B1.1` = vary(base, 30), `Opuntia_stricta|B2.1` = vary(base, 33))
    names(seqs) <- paste0(sub("\\|.*$", "", names(seqs)), "|", l, "_", sub("^.*\\|", "", names(seqs)))
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(seqs), file.path(lib_dir, paste0("LIB_", l, ".fasta")))
  }

  suppressMessages(analyze_barcode_gap(library_dir = lib_dir, output_dir = file.path(tmp, "6_gap"),
                                       paralog_loci = "pepC_like"))

  out_dir <- file.path(tmp, "6_gap")
  d <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_distances_matK.csv"), stringsAsFactors = FALSE)
  expect_setequal(unique(d$type), c("intra", "inter_nn"))
  expect_setequal(unique(d$model), c("raw", "K80"))
  expect_true(all(file.exists(file.path(out_dir, paste0("FIG_barcoding_gap_", c("matK", "pepC_like"), ".png")))))

  s <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_gap_summary.csv"), stringsAsFactors = FALSE)
  expect_true(all(s$possible_paralog[s$locus == "pepC_like"]))
  expect_false(any(s$possible_paralog[s$locus == "matK"]))
  expect_false(any(s$provisional_locus))

  expect_error(analyze_barcode_gap(library_dir = file.path(tmp, "none"), output_dir = file.path(tmp, "6_gap")),
               "finalize_barcoding_library")
  expect_error(analyze_barcode_gap(library_dir = lib_dir, output_dir = file.path("4_Cleaned", "gap")),
               "phylogeny")
})

test_that("the gap does not depend on the random state", {
  sp <- .gap_species()
  m <- .gap_matrix(c(0.01, 0.20, 0.22, 0.30, 0.21, 0.23, 0.31, 0.02, 0.25, 0.26))
  set.seed(1L)
  a <- .compute_barcode_gap(m, sp, "matK", "raw")
  set.seed(77L)
  b <- .compute_barcode_gap(m, sp, "matK", "raw")
  expect_identical(a, b)
})

# Fourth amendment of the validation plan (E13 and E14, 2026-09-22): pairs resting on fewer than
# `min_comparable` comparable positions are treated as having no value, and the summary carries three
# descriptive columns next to the gap.

.gap_comparable <- function(values) {
  sids <- c("A1.1", "A2.1", "B1.1", "B2.1", "C1.1")
  m <- matrix(NA_integer_, 5, 5, dimnames = list(sids, sids))
  m[lower.tri(m)] <- as.integer(values)
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  diag(m) <- 500L
  m
}

test_that("a pair with too few comparable positions has no value, and is counted apart", {
  sp <- .gap_species()
  d <- .gap_matrix(c(0.01, 0.20, 0.22, 0.30, 0.21, 0.23, 0.31, 0.02, 0.25, 0.26))
  # The pair of Opuntia_robusta and the pair (A1.1, B1.1) rest on 3 and 5 positions
  cmp <- .gap_comparable(c(3, 5, 500, 500, 500, 500, 500, 500, 500, 500))

  g <- .compute_barcode_gap(d, sp, "matK", "raw", comparable = cmp, min_comparable = 100L)

  expect_equal(g$summary_df$intra_pairs, 1L)
  expect_equal(g$summary_df$intra_pairs_short_overlap, 1L)
  expect_false("A1.1 A2.1" %in% paste(g$intra$sid_a, g$intra$sid_b))
  # A1.1 no longer has B1.1 as its neighbour: that pair has no value
  expect_equal(g$inter$nn_sid[g$inter$sid == "A1.1"], "B2.1")
  expect_equal(g$inter$distance[g$inter$sid == "A1.1"], 0.22)
  expect_equal(g$summary_df$inter_pairs_short_overlap, 1L)
  expect_equal(g$summary_df$min_comparable, 100L)
})

test_that("without a comparable matrix nothing is excluded, and the counters are zero", {
  sp <- .gap_species()
  g <- .compute_barcode_gap(.gap_matrix(c(0.01, 0.20, 0.22, 0.30, 0.21, 0.23, 0.31, 0.02, 0.25, 0.26)),
                            sp, "matK", "raw")
  expect_equal(g$summary_df$intra_pairs, 2L)
  expect_equal(g$summary_df$intra_pairs_short_overlap, 0L)
  expect_equal(g$summary_df$inter_pairs_short_overlap, 0L)
})

test_that("the three descriptive columns measure the overlap that a percentile of zero hides", {
  sp <- .gap_species()
  # A1.1 and B1.1 are identical (0), so are A2.1 and C1.1; the rest are far apart
  d <- .gap_matrix(c(0.01, 0.00, 0.22, 0.30, 0.21, 0.23, 0.00, 0.02, 0.25, 0.26))

  g <- .compute_barcode_gap(d, sp, "matK", "raw")

  # Queries whose nearest neighbour of another species is at distance 0: A1.1, B1.1, A2.1, C1.1
  expect_equal(g$summary_df$inter_queries_zero, 4L)
  expect_equal(g$summary_df$prop_inter_zero, 4 / 5)
  # Opuntia_robusta with Opuntia_stricta is within the genus; with Cereus_jamacaru it is not
  expect_equal(g$summary_df$inter_zero_congeneric, 2L)
  expect_equal(g$summary_df$inter_zero_other_genus, 2L)
  # Queries whose neighbour is closer than the candidate threshold
  expect_equal(g$summary_df$prop_inter_below_threshold,
               mean(g$inter$distance < g$summary_df$candidate_threshold, na.rm = TRUE))
})

test_that("step 6 declares the minimum of comparable positions it used", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir)
  set.seed(6L)
  base <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  short <- paste0(substr(base, 1, 40), paste(rep("-", 260), collapse = ""))
  seqs <- c(`Opuntia_robusta|A1.1` = base, `Opuntia_robusta|A2.1` = base,
            `Opuntia_stricta|B1.1` = short,
            `Cereus_jamacaru|C1.1` = paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = ""))
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(seqs), file.path(lib_dir, "LIB_matK.fasta"))

  suppressMessages(analyze_barcode_gap(library_dir = lib_dir, output_dir = file.path(tmp, "6_gap"),
                                       models = "raw", min_comparable = 100L, figures = FALSE))

  s <- utils::read.csv(file.path(tmp, "6_gap", "TABLE_barcoding_gap_summary.csv"), stringsAsFactors = FALSE)
  expect_equal(s$min_comparable, 100L)
  # B1.1 shares only 40 positions with the rest: every pair of B1.1 has no value
  expect_true(s$inter_pairs_short_overlap > 0L)
  expect_true(all(c("prop_inter_zero", "inter_zero_congeneric", "inter_zero_other_genus",
                    "prop_inter_below_threshold") %in% names(s)))
})

# Figures (BMM, 2026-09-22): the editorial design of the package and Q1 standards. Vector output for
# the manuscript and a raster copy for reading, the two distributions and the per-species view in the
# same figure, and one overview figure across loci.

test_that("the per-species view pairs each species' largest intraspecific distance with its nearest neighbour", {
  intra <- data.frame(
    locus = "matK", model = "raw",
    species = c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_stricta"),
    sid_a = c("A1.1", "A1.1", "B1.1"), sid_b = c("A2.1", "A3.1", "B2.1"),
    distance = c(0.01, 0.03, 0.05), stringsAsFactors = FALSE
  )
  inter <- data.frame(
    locus = "matK", model = "raw",
    sid = c("A1.1", "A2.1", "A3.1", "B1.1", "B2.1", "C1.1"),
    species = c(rep("Opuntia_robusta", 3), "Opuntia_stricta", "Opuntia_stricta", "Cereus_jamacaru"),
    nn_sid = NA_character_, nn_species = NA_character_,
    distance = c(0.08, 0.09, 0.07, 0.02, 0.04, 0.10), stringsAsFactors = FALSE
  )

  sp <- .bc_gap_species_summary(intra, inter)

  r <- sp[sp$species == "Opuntia_robusta", ]
  expect_equal(c(r$intra_max, r$inter_min), c(0.03, 0.07))
  expect_false(r$no_gap)
  s <- sp[sp$species == "Opuntia_stricta", ]
  expect_equal(c(s$intra_max, s$inter_min), c(0.05, 0.02))
  expect_true(s$no_gap)
  # A species without replica has no intraspecific distance and does not appear
  expect_false("Cereus_jamacaru" %in% sp$species)
  expect_true(all(sp$model == "raw"))
})

test_that("step 6 writes a vector figure and a raster copy per locus, plus one overview", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir)
  set.seed(7L)
  base <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  vary <- function(s, k) { for (p in sample(seq_len(300), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1); s }
  for (l in c("matK", "rbcL")) {
    seqs <- c(`Opuntia_robusta|A1.1` = vary(base, 2), `Opuntia_robusta|A2.1` = vary(base, 3),
              `Opuntia_stricta|B1.1` = vary(base, 25), `Opuntia_stricta|B2.1` = vary(base, 28),
              `Cereus_jamacaru|C1.1` = vary(base, 60))
    names(seqs) <- sub("\\|", paste0("|", l, "_"), names(seqs))
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(seqs), file.path(lib_dir, paste0("LIB_", l, ".fasta")))
  }
  out_dir <- file.path(tmp, "6_gap")

  suppressMessages(analyze_barcode_gap(library_dir = lib_dir, output_dir = out_dir))

  for (l in c("matK", "rbcL")) {
    expect_true(file.exists(file.path(out_dir, paste0("FIG_barcoding_gap_", l, ".pdf"))), info = l)
    expect_true(file.exists(file.path(out_dir, paste0("FIG_barcoding_gap_", l, ".png"))), info = l)
  }
  expect_true(file.exists(file.path(out_dir, "FIG_barcoding_gap_overview.pdf")))
  expect_true(file.exists(file.path(out_dir, "FIG_barcoding_gap_overview.png")))

  # figures = FALSE writes none of them
  out2 <- file.path(tmp, "6_gap_sin_figuras")
  suppressMessages(analyze_barcode_gap(library_dir = lib_dir, output_dir = out2, figures = FALSE))
  expect_length(list.files(out2, pattern = "^FIG_"), 0L)
})

test_that("the figures use the palette of the package and both models, with the two views", {
  intra <- data.frame(locus = "matK", model = rep(c("raw", "K80"), each = 2),
                      species = "Opuntia_robusta", sid_a = "A1.1", sid_b = c("A2.1", "A3.1"),
                      distance = c(0.01, 0.02, 0.011, 0.021), stringsAsFactors = FALSE)
  inter <- data.frame(locus = "matK", model = rep(c("raw", "K80"), each = 3),
                      sid = c("A1.1", "A2.1", "B1.1"), species = c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_stricta"),
                      nn_sid = "B1.1", nn_species = "Opuntia_stricta",
                      distance = c(0.05, 0.06, 0.05, 0.051, 0.061, 0.051), stringsAsFactors = FALSE)
  summary_df <- data.frame(locus = "matK", model = c("raw", "K80"), p95_intra = c(0.02, 0.021),
                        p5_inter = c(0.05, 0.051), candidate_threshold = c(0.02, 0.021),
                        has_gap = TRUE, intra_pairs = 2L, inter_queries = 3L,
                        prop_inter_zero = 0, stringsAsFactors = FALSE)

  p <- .bc_gap_plot(intra, inter, summary_df, "matK")

  expect_s3_class(p, "patchwork")
  expect_equal(length(p$patches$plots) + 1L, 4L)
  expect_true(all(c(.bc_gap_colors()[["intra"]], .bc_gap_colors()[["inter"]]) %in%
                    c("#1b9e77", "#d95f02")))
})

# Editorial review of 2026-09-22 (BMM): the figures were written in Spanish, with long descriptions
# and labels drawn inside the panels that collided with the bars. The rule fixed now: every visible
# string is English and short, and nothing is written inside the plotting area; the two percentiles
# travel in the subtitle, not in annotations next to their own vertical lines.

.gap_fig_strings <- function(q) {
  label_set <- unlist(q$labels, use.names = FALSE)
  annotations <- unlist(lapply(q$layers, function(l) {
    if (!is.null(l$aes_params$label)) as.character(l$aes_params$label) else NULL
  }), use.names = FALSE)
  out <- c(label_set, annotations)
  out[!is.na(out) & nzchar(out)]
}

.gap_fig_has_text_layer <- function(q) {
  any(vapply(q$layers, function(l) inherits(l$geom, "GeomText") || inherits(l$geom, "GeomLabel"),
             logical(1)))
}

test_that("every visible string of the figures is English and short, and no text is drawn inside the panels", {
  intra <- data.frame(locus = "matK", model = rep(c("raw", "K80"), each = 2),
                      species = "Opuntia_robusta", sid_a = "A1.1", sid_b = c("A2.1", "A3.1"),
                      distance = c(0.01, 0.02, 0.011, 0.021), stringsAsFactors = FALSE)
  inter <- data.frame(locus = "matK", model = rep(c("raw", "K80"), each = 3),
                      sid = c("A1.1", "A2.1", "B1.1"),
                      species = c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_stricta"),
                      nn_sid = "B1.1", nn_species = "Opuntia_stricta",
                      distance = c(0.05, 0.06, 0.05, 0.051, 0.061, 0.051), stringsAsFactors = FALSE)
  summary_df <- data.frame(locus = "matK", model = c("raw", "K80"), p95_intra = c(0.02, 0.021),
                        p5_inter = c(0.05, 0.051), candidate_threshold = c(0.02, 0.021),
                        has_gap = TRUE, intra_pairs = 2L, inter_queries = 3L,
                        prop_inter_zero = 0, min_comparable = 100L, stringsAsFactors = FALSE)
  vision <- data.frame(locus = c("matK", "rbcL"), model = "raw", p95_intra = c(0.02, 0.03),
                       p5_inter = c(0, 0), prop_inter_zero = c(0.2, 0.5), stringsAsFactors = FALSE)

  r <- summary_df[summary_df$model == "raw", , drop = FALSE][1, ]
  d <- .bc_gap_distribution_panel(intra[intra$model == "raw", ], inter[inter$model == "raw", ], r, "matK")
  s <- .bc_gap_species_panel(.bc_gap_species_summary(intra, inter), r, "matK")
  p <- .bc_gap_plot(intra, inter, summary_df, "matK")
  o <- .bc_gap_overview_plot(vision)

  legend_labels <- levels(d$data$distribution)
  series <- unique(as.character(o$patches$plots[[1]]$data$series))
  pie <- p$patches$annotation$caption
  visibles <- c(.gap_fig_strings(d), .gap_fig_strings(s), legend_labels, series,
                .gap_fig_strings(o$patches$plots[[1]]), .gap_fig_strings(o), pie)

  # No Spanish and no accented characters anywhere the reader can see
  spanish_pattern <- paste0("(?i)\\b(distancia|distancias|especie|especies|consulta|consultas|pares|",
                       "modelo|proporcion|intraespecific[ao]s?|interespecific[ao]s?|umbral|vecino|",
                       "ningun|identica|porcentaje|replica|conspecifico|proximo|diagonal|bajo|",
                       "tiene|queda|siempre|secuencias|percentil)\\b")
  bad <- visibles[grepl(spanish_pattern, visibles, perl = TRUE)]
  expect_equal(bad, character(0))
  expect_identical(visibles, iconv(visibles, "UTF-8", "ASCII", sub = "?"))

  # Short: titles, subtitles, axes and legend keys of one panel; the caption is the only long string
  expect_lte(max(nchar(c(.gap_fig_strings(d), .gap_fig_strings(s), legend_labels, series))), 60L)
  expect_lte(max(nchar(c(.gap_fig_strings(o$patches$plots[[1]]), .gap_fig_strings(o)))), 60L)
  expect_lte(nchar(pie), 125L)

  # Nothing written inside the plotting area: that is where the labels collided
  expect_false(.gap_fig_has_text_layer(d))
  expect_false(.gap_fig_has_text_layer(s))
  expect_false(.gap_fig_has_text_layer(o$patches$plots[[1]]))
  expect_false(.gap_fig_has_text_layer(o))

  # The two percentiles are still published, now in the subtitle of the panel
  expect_match(d$labels$subtitle, "P95")
  expect_match(d$labels$subtitle, "P5")
  expect_match(d$labels$title, "matK")
  expect_match(d$labels$title, "raw")
})

# Visual review of the regenerated figures (BMM, 2026-09-22): three defects left. In the overview the
# tag of the panel is printed on top of the title; the title and the subtitle assert a result instead
# of reading it from the table, so they would keep saying the same with other data; and in the
# per-species panel the species whose nearest neighbour is at distance 0 are cut by the axis.

test_that("the overview reads its title from the table instead of asserting a result", {
  no_gap <- data.frame(locus = c("matK", "rbcL"), model = "raw", p95_intra = c(0.02, 0.03),
                        p5_inter = c(0, 0), prop_inter_zero = c(0.2, 0.5), has_gap = FALSE,
                        stringsAsFactors = FALSE)
  con_gap <- data.frame(locus = c("matK", "rbcL"), model = "raw", p95_intra = c(0.005, 0.03),
                        p5_inter = c(0.01, 0), prop_inter_zero = c(0.2, 0.5), has_gap = c(TRUE, FALSE),
                        stringsAsFactors = FALSE)

  a <- .bc_gap_overview_plot(no_gap)$patches$plots[[1]]$labels
  b <- .bc_gap_overview_plot(con_gap)$patches$plots[[1]]$labels

  expect_match(a$title, "^No locus")
  expect_match(b$title, "1 of 2")
  expect_false(identical(a$title, b$title))
  expect_match(a$subtitle, "2 of 2")
  expect_match(b$subtitle, "1 of 2")
  expect_lte(max(nchar(c(a$title, b$title, a$subtitle, b$subtitle))), 60L)

  # A table without the has_gap column is read from the two percentiles, not refused
  sin_col <- no_gap[, setdiff(names(no_gap), "has_gap")]
  expect_match(.bc_gap_overview_plot(sin_col)$patches$plots[[1]]$labels$title, "^No locus")
})

test_that("the title is indented so the tag of the panel does not fall on it", {
  th <- .bc_gap_theme()
  expect_gte(as.numeric(th$plot.title$margin)[4], 10)
  expect_gte(as.numeric(th$plot.subtitle$margin)[4], 10)
})

test_that("a species whose nearest neighbour is at distance 0 is drawn whole, not cut by the axis", {
  intra <- data.frame(locus = "matK", model = "raw", species = c("Opuntia_robusta", "Opuntia_stricta"),
                      sid_a = c("A1.1", "B1.1"), sid_b = c("A2.1", "B2.1"),
                      distance = c(0.01, 0.02), stringsAsFactors = FALSE)
  inter <- data.frame(locus = "matK", model = "raw", sid = c("A1.1", "B1.1"),
                      species = c("Opuntia_robusta", "Opuntia_stricta"),
                      nn_sid = c("B1.1", "A1.1"), nn_species = c("Opuntia_stricta", "Opuntia_robusta"),
                      distance = c(0, 0.05), stringsAsFactors = FALSE)
  r <- data.frame(locus = "matK", model = "raw", p95_intra = 0.02, p5_inter = 0,
                  has_gap = FALSE, intra_pairs = 2L, inter_queries = 2L, prop_inter_zero = 0.5,
                  stringsAsFactors = FALSE)

  s <- .bc_gap_species_panel(.bc_gap_species_summary(intra, inter), r, "matK")

  expect_lt(s$coordinates$limits$y[1], 0)
  expect_lt(s$coordinates$limits$x[1], 0)
  expect_gt(s$coordinates$limits$y[2], 0)
})

test_that("the title of the narrow panel of the overview fits its third of the page", {
  # In the figure of 2026-09-22 the title of the right panel was cut against the edge: it holds a
  # third of the 180 mm, so it takes a shorter title than the left one.
  d <- data.frame(locus = c("matK", "rbcL"), model = "raw", p95_intra = c(0.02, 0.03),
                  p5_inter = c(0, 0), prop_inter_zero = c(0.2, 0.5), has_gap = FALSE,
                  stringsAsFactors = FALSE)

  o <- .bc_gap_overview_plot(d)

  expect_lte(nchar(o$labels$title), 30L)
  expect_lte(nchar(o$labels$subtitle), 40L)
  expect_lte(nchar(o$patches$plots[[1]]$labels$title), 40L)
})
