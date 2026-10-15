# Decisions T6-5 of BMM (14-10, Phase 12, proposal of Tutorial 6): two exported functions draw, from the
# tables the identification report already writes (REPORT_<run>_answer.csv, REPORT_<run>_sampling.csv),
# the answer of one identification and the summary of several. Tutorial 6 draws its figures with them,
# so the vignette and the report show the same figure. Written before the code.

.plt_answer <- function(query = "q1") {
  data.frame(query = query,
             locus = c("ITS", "matK", "rbcL", "trnL-trnF"),
             compartment = c("nuclear", "plastid", "plastid", "plastid"),
             state = c(3L, 1L, 2L, 3L),
             predicted_species = c(NA, "Opuntia_alpha", NA, NA),
             predicted_genus = c(NA, "Opuntia", "Opuntia", NA),
             genus_confidence = c(NA, 95, 85, 30),
             species_confidence = c(NA, 80, 75, 20),
             species_threshold = c(60.9, 66.9, 60, 63.7),
             reason = c("no_overlap", NA, "genus_only_locus", "low_confidence"),
             stringsAsFactors = FALSE)
}

.plt_sampling <- function() {
  data.frame(genus = "Opuntia", species = rep(c("Opuntia_alpha", "Opuntia_beta"), each = 2),
             locus = rep(c("matK", "rbcL"), 2), sequences = c(3L, 1L, 0L, 2L), stringsAsFactors = FALSE)
}

test_that("plot_barcoding_identification() draws the answer per locus, and the sampling below it when given", {
  p <- plot_barcoding_identification(.plt_answer())
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
  q <- plot_barcoding_identification(.plt_answer(), sampling = .plt_sampling())
  expect_true(inherits(q, "patchwork"))
})

test_that("plot_barcoding_identification() draws the threshold given as a dashed line", {
  p <- plot_barcoding_identification(.plt_answer(), threshold = 60)
  g <- if (inherits(p, "patchwork")) p[[1]] else p
  hl <- Filter(function(l) inherits(l$geom, "GeomHline"), g$layers)
  expect_gte(length(hl), 1L)
  expect_true(any(vapply(hl, function(l) isTRUE(all.equal(l$data$yintercept %||% l$aes_params$yintercept, 60)), logical(1))))
})

test_that("plot_barcoding_identification() stops with the name of a missing column", {
  a <- .plt_answer()
  a$state <- NULL
  expect_error(plot_barcoding_identification(a), "state")
  s <- .plt_sampling()
  s$sequences <- NULL
  expect_error(plot_barcoding_identification(.plt_answer(), sampling = s), "sequences")
})

test_that("plot_barcoding_identification() reads the CSV files the report writes", {
  tmp <- withr::local_tempdir()
  fa <- file.path(tmp, "REPORT_run1_answer.csv"); fs <- file.path(tmp, "REPORT_run1_sampling.csv")
  utils::write.csv(.plt_answer(), fa, row.names = FALSE)
  utils::write.csv(.plt_sampling(), fs, row.names = FALSE)
  p <- plot_barcoding_identification(fa, sampling = fs)
  expect_true(inherits(p, "patchwork"))
})

test_that("the summary has one row per example and one column per locus found in any of them", {
  b <- .plt_answer("q2")
  b$state <- c(3L, 3L, 3L, 2L); b$predicted_species <- NA; b$predicted_genus <- c(NA, NA, NA, "Cereus")
  b$reason <- c("no_overlap", "low_confidence", "low_confidence", NA)
  p <- plot_barcoding_identification_summary(list(first = .plt_answer(), second = b))
  expect_true(inherits(p, "ggplot"))
  d <- p$data
  expect_true(all(c("example", "locus", "state", "genus_only", "answer") %in% names(d)))
  expect_setequal(unique(as.character(d$example)), c("first", "second"))
  # ITS was not found (no_overlap) in either example: no column for it
  expect_setequal(unique(as.character(d$locus)), c("matK", "rbcL", "trnL-trnF"))
  expect_equal(nrow(d), 6L)
})

test_that("the summary marks a row named at genus rank only and writes the answer of each example", {
  b <- .plt_answer("q2")
  b$state <- 3L; b$predicted_species <- NA; b$predicted_genus <- NA; b$reason <- "low_confidence"
  p <- plot_barcoding_identification_summary(list(first = .plt_answer(), second = b))
  d <- p$data
  expect_true(d$genus_only[d$example == "first" & d$locus == "rbcL"])
  expect_false(any(d$genus_only[d$locus != "rbcL"]))
  a1 <- unique(d$answer[d$example == "first"]); a2 <- unique(d$answer[d$example == "second"])
  expect_length(a1, 1L); expect_length(a2, 1L)
  expect_match(a1, "Opuntia alpha", fixed = TRUE)
  expect_match(a2, "no answer", fixed = TRUE)
})

test_that("the summary needs a named list of answer tables", {
  expect_error(plot_barcoding_identification_summary(list(.plt_answer())), "named")
  expect_error(plot_barcoding_identification_summary(.plt_answer()), "list")
})

test_that("neither function writes a file", {
  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)
  plot_barcoding_identification(.plt_answer(), sampling = .plt_sampling())
  plot_barcoding_identification_summary(list(first = .plt_answer()))
  expect_length(list.files(tmp, recursive = TRUE, all.files = TRUE, no.. = TRUE), 0L)
})
