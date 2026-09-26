# Tests of the remoteness threshold of the molecular diagnostic branch (11_barcoding/9_threshold/).
# Written before the code. Phase 6A of PhyloCactus 0.5.0, proposal of 2026-09-25, decisions D1 to D4.
#
# Without a rule of remoteness the nearest neighbour names a cactus species for anything: in CN2, 288
# of 312 alien queries came out in state 1. The rule sends a query whose distance to its nearest
# neighbour exceeds a threshold to state 3. The threshold of each locus is the quantile 0.99, type 1,
# of the distances of the legitimate queries of scheme G (a cactus whose species is missing from the
# library and whose genus is not), and it is fixed without looking at the alien queries, which only
# report afterwards how many it rejects.

.thr_pred <- function(scheme = "species") {
  # Legitimate queries of one locus by hand, with every category the plan separates
  data.frame(
    locus = "matK", scheme = scheme, fold = 1:6, stratum = NA_character_,
    sid = paste0("q", 1:6),
    true_species = c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_robusta",
                          "Opuntia_stricta", "Cereus_jamacaru", "Cereus_jamacaru"),
    true_genus = c("Opuntia", "Opuntia", "Opuntia", "Opuntia", "Cereus", "Cereus"),
    method = "nn", alignment = "add", orientation = "forward",
    state = c(1L, 1L, 2L, 3L, 2L, 3L),
    predicted_species = c("Opuntia_robusta", "Opuntia_stricta", NA, NA, NA, NA),
    predicted_genus = c("Opuntia", "Opuntia", "Opuntia", NA, "Opuntia", NA),
    candidates = NA_character_,
    nn_distance = c(0.01, 0.02, 0.03, 0.04, 0.05, NA),
    margin = 0, confidence = NA_real_,
    reason = c(NA, NA, NA, "tie_across_genera", NA, "no_comparable_positions"),
    stringsAsFactors = FALSE
  )
}

.thr_alien <- function() {
  # Alien queries as CN2 writes them, one of them with no homology to the locus
  data.frame(
    locus = "matK", sid = paste0("o", 1:5),
    query_species = c("Portulaca_amilis", "Portulaca_amilis", "Talinum_paniculatum",
                         "Anacampseros_filamentosa", "Talinopsis_frutescens"),
    orientation = c("forward", "forward", "reverse", "forward", "no_match"),
    state = c(1L, 1L, 1L, 3L, 3L),
    predicted_species = c("Opuntia_robusta", "Opuntia_stricta", "Cereus_jamacaru", NA, NA),
    predicted_genus = c("Opuntia", "Opuntia", "Cereus", NA, NA),
    candidates = NA_character_,
    nn_distance = c(0.015, 0.06, 0.08, 0.09, NA),
    margin = 0, confidence = NA_real_,
    reason = c(NA, NA, NA, "tie_across_genera", "no_match"),
    stringsAsFactors = FALSE
  )
}

.thr_dirs <- function(tmp, with_cn2 = TRUE) {
  cls <- file.path(tmp, "7_classifier")
  ctl <- file.path(tmp, "8_controls")
  dir.create(cls, recursive = TRUE, showWarnings = FALSE)
  dir.create(ctl, recursive = TRUE, showWarnings = FALSE)
  e <- .thr_pred("species")
  g <- .thr_pred("genus")
  g$nn_distance <- c(0.02, 0.03, 0.06, 0.07, 0.08, 0.10)
  g$reason[6] <- NA
  g$state[6] <- 2L
  g$predicted_genus[6] <- "Cereus"
  # A second locus with no alien query at all
  e2 <- e; e2$locus <- "rbcL"
  g2 <- g; g2$locus <- "rbcL"
  utils::write.csv(rbind(e, e2), file.path(cls, "TABLE_barcoding_predictions_species_nn_add.csv"), row.names = FALSE)
  utils::write.csv(rbind(g, g2), file.path(cls, "TABLE_barcoding_predictions_genus_nn_add.csv"), row.names = FALSE)
  if (with_cn2) {
    utils::write.csv(.thr_alien(), file.path(ctl, "TABLE_barcoding_cn2_queries.csv"), row.names = FALSE)
  }
  list(classifier_dir = cls, controls_dir = ctl, out = file.path(tmp, "9_threshold"))
}

test_that("the rule of remoteness sends only the distant queries to state 3 and keeps every earlier reason", {
  p <- .thr_pred()
  r <- .bc_apply_remoteness(p, 0.025)
  # At or below the threshold nothing changes
  expect_identical(r[1:2, ], p[1:2, ])
  # Beyond it: state 3, reason remoteness, no species and no genus
  expect_equal(r$state[3], 3L)
  expect_equal(r$reason[3], "remoteness")
  expect_true(is.na(r$predicted_genus[3]))
  expect_equal(r$state[5], 3L)
  expect_equal(r$reason[5], "remoteness")
  # A query already in state 3 keeps its own reason, and one with no distance is left as it was
  expect_equal(r$reason[4], "tie_across_genera")
  expect_equal(r$reason[6], "no_comparable_positions")
  # The distance is never altered: the rule reads it, it does not rewrite it
  expect_identical(r$nn_distance, p$nn_distance)
})

test_that("the five categories are disjoint and add up to the queries, in both schemes", {
  for (sc in c("species", "genus")) {
    p <- .thr_pred(sc)
    cat5 <- .bc_outcome_category(p, sc)
    expect_length(cat5, nrow(p))
    expect_true(all(cat5 %in% c("correct_species", "wrong_species", "correct_genus",
                                "wrong_genus", "unassigned")))
  }
  e <- .bc_outcome_category(.thr_pred("species"), "species")
  expect_equal(e, c("correct_species", "wrong_species", "correct_genus", "unassigned",
                    "wrong_genus", "unassigned"))
})

test_that("in scheme G no query is ever a correct species, even when it is forced to state 1", {
  p <- .thr_pred("genus")
  # The fixture names the true species in state 1; in scheme G that species is not in the training
  # set, so a species name there is always a false one
  g <- .bc_outcome_category(p, "genus")
  expect_false(any(g == "correct_species"))
  expect_equal(g[1:2], c("wrong_species", "wrong_species"))
})

test_that("above every distance the curve reproduces step 7 row by row, and below every distance nothing is assigned", {
  p <- .thr_pred()
  expect_identical(.bc_apply_remoteness(p, 1), p)
  below <- .bc_apply_remoteness(p, -1)
  expect_true(all(below$state == 3L))
  expect_true(all(is.na(below$predicted_species)))
})

test_that("the share of unassigned queries never falls as the threshold goes down", {
  cv <- .bc_threshold_curve(.thr_pred(), "species", outgroup = .thr_alien())
  cv <- cv[order(cv$threshold, decreasing = TRUE), ]
  expect_true(all(diff(cv$unassigned) >= 0))
  expect_true(all(diff(cv$outgroup_rejected) >= 0))
  # One row per distance observed, legitimate or alien, plus zero
  observed <- c(.thr_pred()$nn_distance, .thr_alien()$nn_distance[.thr_alien()$orientation != "no_match"])
  expect_setequal(cv$threshold, sort(unique(c(0, observed[!is.na(observed)]))))
  # The five categories add up to the queries at every threshold
  expect_true(all(cv$correct_species + cv$wrong_species + cv$correct_genus +
                    cv$wrong_genus + cv$unassigned == cv$queries))
})

test_that("alien queries with no homology stay out of the curve, and a locus with only those is not measured", {
  a <- .thr_alien()
  cv <- .bc_threshold_curve(.thr_pred(), "species", outgroup = a)
  expect_true(all(cv$outgroup_comparable == 4L))
  # An alien query already in state 3 counts as rejected at any threshold
  expect_true(all(cv$outgroup_rejected >= 1L))

  solo <- a[a$orientation == "no_match", , drop = FALSE]
  cv0 <- .bc_threshold_curve(.thr_pred(), "species", outgroup = solo)
  expect_true(all(cv0$outgroup_comparable == 0L))
  expect_true(all(is.na(cv0$outgroup_rejected)))
})

test_that("the operating threshold is fixed without the alien queries", {
  skip_if_not_installed("withr")
  tmp1 <- withr::local_tempdir()
  tmp2 <- withr::local_tempdir()
  d1 <- .thr_dirs(tmp1, with_cn2 = TRUE)
  d2 <- .thr_dirs(tmp2, with_cn2 = FALSE)
  suppressMessages(sweep_barcoding_threshold(classifier_dir = d1$classifier_dir, controls_dir = d1$controls_dir,
                                             output_dir = d1$out, figures = FALSE))
  suppressMessages(sweep_barcoding_threshold(classifier_dir = d2$classifier_dir, controls_dir = d2$controls_dir,
                                             output_dir = d2$out, figures = FALSE))
  o1 <- utils::read.csv(file.path(d1$out, "TABLE_barcoding_threshold_operating.csv"), stringsAsFactors = FALSE)
  o2 <- utils::read.csv(file.path(d2$out, "TABLE_barcoding_threshold_operating.csv"), stringsAsFactors = FALSE)
  expect_equal(o1$threshold, o2$threshold)
  # Without CN2 the rejection is not measured anywhere, and it says so
  expect_true(all(o2$rejection_status == "not measured"))
  # With CN2, matK is measured and rbcL, which has no alien query, is not
  expect_equal(unique(o1$rejection_status[o1$locus == "matK"]), "measured")
  expect_equal(unique(o1$rejection_status[o1$locus == "rbcL"]), "not measured")
  expect_true(all(is.na(o1$rejection[o1$locus == "rbcL"])))
})

test_that("the operating threshold is the declared quantile of scheme G, type 1, and the table says which", {
  g <- .thr_pred("genus")
  g$nn_distance <- c(0.02, 0.03, 0.06, 0.07, 0.08, 0.10)
  op <- .bc_operating_point(g, q = 0.5, type = 1L)
  # Type 1 returns a distance that exists in the data
  expect_equal(op, unname(stats::quantile(g$nn_distance, 0.5, type = 1)))
  expect_true(op %in% g$nn_distance)
  expect_error(.bc_operating_point(.thr_pred("species"), q = 0.99, type = 1L), "genus")

  skip_if_not_installed("withr")
  tmp <- withr::local_tempdir()
  d <- .thr_dirs(tmp)
  suppressMessages(sweep_barcoding_threshold(classifier_dir = d$classifier_dir, controls_dir = d$controls_dir,
                                             output_dir = d$out, q = 0.99, quantile_type = 1L, figures = FALSE))
  o <- utils::read.csv(file.path(d$out, "TABLE_barcoding_threshold_operating.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("locus", "scheme", "threshold", "q", "quantile_type", "queries", "correct_species",
                    "wrong_species", "correct_genus", "wrong_genus", "unassigned",
                    "outgroup_comparable", "outgroup_rejected", "rejection", "rejection_ci_low",
                    "rejection_ci_high", "rejection_status") %in% names(o)))
  expect_true(all(o$q == 0.99))
  expect_true(all(o$quantile_type == 1L))
  # One row per locus and scheme, and the threshold is the same in both schemes of a locus
  expect_equal(nrow(o), 4L)
  expect_equal(length(unique(o$threshold[o$locus == "matK"])), 1L)
})

test_that("the interval of the rejection is the exact binomial one", {
  for (xn in list(c(0L, 7L), c(7L, 7L), c(12L, 13L))) {
    ci <- .bc_clopper_pearson(xn[1], xn[2])
    ref <- stats::binom.test(xn[1], xn[2])$conf.int
    expect_equal(unname(ci), as.numeric(ref), tolerance = 1e-12)
  }
  expect_true(all(is.na(.bc_clopper_pearson(0L, 0L))))
})

test_that("the sweep writes its tables and figures, and names the step that is missing", {
  skip_if_not_installed("withr")
  tmp <- withr::local_tempdir()
  d <- .thr_dirs(tmp)
  suppressMessages(sweep_barcoding_threshold(classifier_dir = d$classifier_dir, controls_dir = d$controls_dir,
                                             output_dir = d$out, figures = TRUE))
  expect_true(file.exists(file.path(d$out, "TABLE_barcoding_threshold_curve.csv")))
  expect_true(file.exists(file.path(d$out, "TABLE_barcoding_threshold_operating.csv")))
  for (l in c("matK", "rbcL")) {
    expect_true(file.exists(file.path(d$out, paste0("FIG_barcoding_threshold_", l, ".pdf"))))
    expect_true(file.exists(file.path(d$out, paste0("FIG_barcoding_threshold_", l, ".png"))))
  }
  cv <- utils::read.csv(file.path(d$out, "TABLE_barcoding_threshold_curve.csv"), stringsAsFactors = FALSE)
  expect_setequal(unique(paste(cv$locus, cv$scheme)),
                  c("matK species", "matK genus", "rbcL species", "rbcL genus"))

  expect_error(sweep_barcoding_threshold(classifier_dir = file.path(tmp, "none"),
                                         controls_dir = d$controls_dir, output_dir = d$out),
               "classify_barcoding_folds")
  expect_error(sweep_barcoding_threshold(classifier_dir = d$classifier_dir, controls_dir = d$controls_dir,
                                         output_dir = file.path("4_Cleaned", "thr")),
               "phylogeny")
})
