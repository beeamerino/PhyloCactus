# Tests of the per-class metrics of the molecular diagnostic branch (11_barcoding/11_metrics/).
# Written before the code. Phase 6D of PhyloCactus 0.5.0, decisions M1 to M4 of BMM (26-09): the
# metrics of section 5 of the validation plan, per locus and scheme, for the nearest neighbour
# through the identification path at t* and for IdTaxa at 60, side by side and never combined.
#
# The fixture is written by hand so that every figure can be checked by hand. None of these tests
# runs a classifier: the function reads the tables of steps 7 and 9 and the library summary.

.met_write <- function(d, path) utils::write.csv(d, path, row.names = FALSE)

.met_nn_row <- function(scheme, sid, true_species, state, species, genus, candidates, distance, reason = NA) {
  data.frame(locus = "matK", scheme = scheme, fold = sid, stratum = true_species, sid = paste0("s", sid),
             true_species = true_species, true_genus = sub("_.*", "", true_species),
             method = "nn", alignment = "add", orientation = "forward", state = state,
             predicted_species = species, predicted_genus = genus, candidates = candidates,
             nn_distance = distance, margin = 0.01, confidence = NA_real_, reason = reason,
             stringsAsFactors = FALSE)
}

.met_idt_row <- function(scheme, sid, true_species, genus_idtaxa, species_idtaxa, gc, sc) {
  data.frame(locus = "matK", scheme = scheme, fold = sid, stratum = true_species, sid = paste0("s", sid),
             true_species = true_species, true_genus = sub("_.*", "", true_species),
             method = "idtaxa", alignment = "library", state = 3L, predicted_species = NA_character_,
             predicted_genus = NA_character_, candidates = NA_character_, nn_distance = NA_real_,
             margin = NA_real_, confidence = NA_real_, reason = "low_confidence",
             genus_idtaxa = genus_idtaxa, species_idtaxa = species_idtaxa,
             genus_confidence = gc, species_confidence = sc, stringsAsFactors = FALSE)
}

.met_fixture <- function(tmp) {
  lib <- file.path(tmp, "4_library"); cls <- file.path(tmp, "7_classifier"); thr <- file.path(tmp, "9_threshold")
  for (d in c(lib, cls, thr)) dir.create(d, recursive = TRUE)
  .met_write(data.frame(locus = "matK", total_species = 40L, species_with_replicate = 3L, total_accessions = 90L,
                        total_genera = 5L, genera_with_2plus_species = 2L), file.path(lib, "TABLE_barcoding_library_summary.csv"))
  # Nearest neighbour, scheme E. t* = 0.1: s5 is farther and goes to state 3 at the operating point
  nn_e <- rbind(.met_nn_row("species", 1, "A_a", 1L, "A_a", "A", "A_a", 0.01),
                .met_nn_row("species", 2, "A_a", 1L, "A_b", "A", "A_b", 0.01),
                .met_nn_row("species", 3, "A_b", 1L, "A_b", "A", "A_b", 0.02),
                .met_nn_row("species", 4, "A_b", 2L, NA, "A", "A_a|A_b", 0.01),
                .met_nn_row("species", 5, "B_c", 1L, "B_c", "B", "B_c", 0.50),
                .met_nn_row("species", 6, "B_c", 3L, NA, NA, "A_a|B_c", 0.01, "tie_across_genera"))
  # Nearest neighbour, scheme G: a state 1 names an absent species
  nn_g <- rbind(.met_nn_row("genus", 1, "A_a", 2L, NA, "A", "A_b", 0.01),
                .met_nn_row("genus", 2, "A_b", 1L, "A_a", "A", "A_a", 0.01),
                .met_nn_row("genus", 3, "B_c", 2L, NA, "A", "A_a|A_b", 0.02),
                .met_nn_row("genus", 4, "B_d", 2L, NA, "B", "B_c", 0.03))
  .met_write(nn_e, file.path(cls, "TABLE_barcoding_predictions_species_nn_add.csv"))
  .met_write(nn_g, file.path(cls, "TABLE_barcoding_predictions_genus_nn_add.csv"))
  # IdTaxa, scheme E and G. At 60: i1 correct species, i2 wrong species, i3 correct genus,
  # i4 unassigned (it would be state 2 at the quantile threshold 10)
  id_e <- rbind(.met_idt_row("species", 1, "A_a", "A", "A_a", 95, 90),
                .met_idt_row("species", 2, "A_b", "A", "A_a", 99, 70),
                .met_idt_row("species", 3, "B_c", "B", "B_d", 80, 40),
                .met_idt_row("species", 4, "B_c", "B", "B_c", 30, 20))
  id_g <- rbind(.met_idt_row("genus", 1, "A_a", "A", "A_b", 90, 50),
                .met_idt_row("genus", 2, "B_c", "A", "A_a", 70, 20))
  .met_write(id_e, file.path(cls, "TABLE_barcoding_predictions_species_idtaxa.csv"))
  .met_write(id_g, file.path(cls, "TABLE_barcoding_predictions_genus_idtaxa.csv"))
  # Step 9: the operating points
  cats <- data.frame(correct_species = 0L, wrong_species = 0L, correct_genus = 0L, wrong_genus = 0L, unassigned = 0L)
  .met_write(cbind(data.frame(locus = "matK", scheme = c("species", "genus"), threshold = 0.1, q = 0.99,
                              quantile_type = 1L, queries = c(6L, 4L)), cats),
             file.path(thr, "TABLE_barcoding_threshold_operating.csv"))
  .met_write(cbind(data.frame(locus = "matK", scheme = rep(c("species", "genus"), each = 2),
                              rule = rep(c("quantile", "default_60"), 2), threshold = rep(c(10, 60), 2),
                              q = rep(c(0.01, NA), 2), quantile_type = rep(c(1L, NA), 2), queries = rep(c(4L, 2L), each = 2)), cats),
             file.path(thr, "TABLE_barcoding_threshold_operating_idtaxa.csv"))
  curve <- function(n) cbind(data.frame(locus = "matK", scheme = rep(c("species", "genus"), each = 2),
                                        threshold = rep(c(0, 1), 2), queries = rep(n, each = 2)),
                             data.frame(correct_species = c(1L, 0L, 0L, 0L), wrong_species = c(1L, 0L, 1L, 0L),
                                        correct_genus = c(1L, 1L, 1L, 1L), wrong_genus = 0L,
                                        unassigned = c(n[1] - 3L, n[1] - 1L, n[2] - 2L, n[2] - 1L)))
  .met_write(curve(c(6L, 4L)), file.path(thr, "TABLE_barcoding_threshold_curve.csv"))
  .met_write(curve(c(4L, 2L)), file.path(thr, "TABLE_barcoding_threshold_curve_idtaxa.csv"))
  list(library_dir = lib, classifier_dir = cls, threshold_dir = thr, output_dir = file.path(tmp, "11_metrics"))
}

.met_run <- function(f, ...) {
  suppressMessages(utils::capture.output(
    res <- summarise_barcoding_metrics(classifier_dir = f$classifier_dir, threshold_dir = f$threshold_dir,
                                       library_dir = f$library_dir, output_dir = f$output_dir, ...)))
  res
}

.met_read <- function(f, what, method) {
  utils::read.csv(file.path(f$output_dir, sprintf("TABLE_barcoding_metrics_%s_%s.csv", what, method)),
                  stringsAsFactors = FALSE)
}

test_that("per-class recall and precision of scheme E are those computed by hand", {
  tmp <- withr::local_tempdir(); f <- .met_fixture(tmp)
  .met_run(f, figures = FALSE)
  cl <- .met_read(f, "class", "nn_add")
  e <- cl[cl$scheme == "species", ]
  e <- e[order(e$class), ]
  expect_identical(e$class, c("A_a", "A_b", "B_c"))
  expect_identical(e$queries, c(2L, 2L, 2L))
  expect_identical(e$correct, c(1L, 1L, 0L))
  expect_equal(e$recall, c(0.5, 0.5, 0))
  # A_a named once and right; A_b named twice (s2, wrongly, and s3); B_c never named at t*
  expect_identical(e$assigned, c(1L, 2L, 0L))
  expect_equal(e$precision, c(1, 0.5, NA))
})

test_that("a class never predicted has NA precision, counted apart, never 0", {
  tmp <- withr::local_tempdir(); f <- .met_fixture(tmp)
  .met_run(f, figures = FALSE)
  s <- .met_read(f, "summary", "nn_add")
  e <- s[s$scheme == "species", ]
  expect_identical(e$classes_never_predicted, 1L)
  expect_equal(e$precision_median, 0.75)
  expect_identical(e$classes_recall_zero, 1L)
  expect_identical(e$classes_recall_zero_list, "B_c")
  expect_equal(e$recall_median, 0.5)
  expect_equal(c(e$recall_q1, e$recall_q3), unname(stats::quantile(c(0.5, 0.5, 0), c(0.25, 0.75))))
})

test_that("the unassigned and the two kinds of misassignment are separate columns, never summed", {
  tmp <- withr::local_tempdir(); f <- .met_fixture(tmp)
  .met_run(f, figures = FALSE)
  s <- .met_read(f, "summary", "nn_add")
  e <- s[s$scheme == "species", ]
  expect_identical(c(e$correct_species, e$wrong_species, e$correct_genus, e$wrong_genus, e$unassigned),
                   c(2L, 1L, 1L, 0L, 2L))
  expect_equal(c(e$unassigned_rate, e$wrong_species_rate, e$wrong_genus_rate), c(2, 1, 0) / 6)
  expect_equal(e$honest_genus_rate, 1 / 6)
  expect_false(any(grepl("error|misassign|accuracy_overall", names(s))))
})

test_that("in scheme G a state 1 is a misassignment and the class is the genus", {
  tmp <- withr::local_tempdir(); f <- .met_fixture(tmp)
  .met_run(f, figures = FALSE)
  cl <- .met_read(f, "class", "nn_add")
  g <- cl[cl$scheme == "genus", ]
  g <- g[order(g$class), ]
  expect_identical(g$class, c("A", "B"))
  expect_identical(g$correct, c(1L, 1L))
  expect_equal(g$recall, c(0.5, 0.5))
  expect_equal(g$precision, c(0.5, 1))
  s <- .met_read(f, "summary", "nn_add")
  gs <- s[s$scheme == "genus", ]
  expect_identical(c(gs$correct_species, gs$wrong_species, gs$correct_genus, gs$wrong_genus, gs$unassigned),
                   c(0L, 1L, 2L, 1L, 0L))
  expect_true(is.na(gs$honest_genus_rate))
})

test_that("the five context columns of the plan are on every summary row", {
  tmp <- withr::local_tempdir(); f <- .met_fixture(tmp)
  .met_run(f, figures = FALSE)
  ctx <- c("total_species", "species_with_replicate", "total_accessions", "total_genera", "genera_with_2plus_species")
  for (m in c("nn_add", "idtaxa")) {
    s <- .met_read(f, "summary", m)
    expect_true(all(ctx %in% names(s)))
    expect_identical(unname(unlist(s[1, ctx])), c(40L, 3L, 90L, 5L, 2L))
    expect_false(anyNA(s[, ctx]))
  }
})

test_that("the operating points are read from step 9: t* for the nearest neighbour, 60 for IdTaxa", {
  tmp <- withr::local_tempdir(); f <- .met_fixture(tmp)
  .met_run(f, figures = FALSE)
  s <- .met_read(f, "summary", "idtaxa")
  e <- s[s$scheme == "species", ]
  expect_equal(e$threshold, 60)
  expect_identical(e$rule, "default_60")
  # i4 has genus confidence 30: unassigned at 60, state 2 at the quantile row (10)
  expect_identical(c(e$correct_species, e$wrong_species, e$correct_genus, e$wrong_genus, e$unassigned),
                   c(1L, 1L, 1L, 0L, 1L))
  g <- s[s$scheme == "genus", ]
  expect_identical(c(g$correct_genus, g$wrong_genus), c(1L, 1L))
  n <- .met_read(f, "summary", "nn_add")
  expect_equal(n$threshold, c(0.1, 0.1))
})

test_that("both classifiers are written apart and never combined", {
  tmp <- withr::local_tempdir(); f <- .met_fixture(tmp)
  .met_run(f, figures = FALSE)
  files <- list.files(f$output_dir)
  expect_setequal(files, c("TABLE_barcoding_metrics_class_nn_add.csv", "TABLE_barcoding_metrics_summary_nn_add.csv",
                           "TABLE_barcoding_metrics_class_idtaxa.csv", "TABLE_barcoding_metrics_summary_idtaxa.csv"))
  for (m in c("nn_add", "idtaxa")) expect_true(all(.met_read(f, "summary", m)$method == m))
})

test_that("a missing step 9 table stops with a message that names the step", {
  tmp <- withr::local_tempdir(); f <- .met_fixture(tmp)
  file.remove(file.path(f$threshold_dir, "TABLE_barcoding_threshold_operating_idtaxa.csv"))
  expect_error(.met_run(f, figures = FALSE), "sweep_barcoding_threshold")
})

test_that("with figures, one curve figure per locus with both classifiers", {
  skip_if_not_installed("ggplot2")
  tmp <- withr::local_tempdir(); f <- .met_fixture(tmp)
  .met_run(f, figures = TRUE)
  expect_true(file.exists(file.path(f$output_dir, "FIG_barcoding_metrics_curve_matK.pdf")))
  expect_true(file.exists(file.path(f$output_dir, "FIG_barcoding_metrics_curve_matK.png")))
})
