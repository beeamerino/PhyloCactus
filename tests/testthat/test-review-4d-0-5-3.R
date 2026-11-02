# Closure of 0.5.2 into 0.5.3, phase 4d (review of the text of 0.5.3 by a session that did not write
# it, minutes of 27-10-2026 in the audit repository): the findings BMM accepted on the distributed
# registry and on the cross-validation warning (M5, m19). Written before the corrections.

.d4_extdata <- function(name) {
  f <- system.file("extdata", name, package = "PhyloCactus")
  skip_if(!nzchar(f), paste(name, "not installed"))
  f
}

test_that("the distributed registry carries no checklist content beyond the names and the match (M5)", {
  sp <- utils::read.csv(.d4_extdata("phylocactus_table_species.csv"), nrows = 1, check.names = FALSE)
  dropped <- c("checklist_fullName", "checklist_author", "checklist_rank", "checklist_taxon", "checklist_uuid")
  expect_length(intersect(dropped, names(sp)), 0L)
  expect_true(all(c("species", "accepted_in_checklist", "status_pipeline") %in% names(sp)))
})

test_that("the cross-validation warning says when the lowest value without failures is the selected one (m19)", {
  tmp <- withr::local_tempdir()
  grid <- 10^(1:-8)
  per <- matrix(300 / 10, nrow = length(grid), ncol = 10)
  per[grid %in% c(1e-5, 1e-4), 4] <- 1e36
  log <- unlist(lapply(seq_along(grid), function(i) {
    run <- cumsum(per[i, ])
    c(paste0("smoothing:", grid[i]), paste0("(chisq): ", signif(run, 8)),
      paste0("chisq: (", grid[i], ") ", signif(run[10], 8)))
  }))
  totals <- rowSums(per) + rev(seq_along(grid)) * 0.01
  writeLines(log, file.path(tmp, "cv.log"))
  writeLines(paste0("chisq: (", grid, ") ", signif(totals, 8)), file.path(tmp, "cv"))
  sel <- PhyloCactus:::.select_cv_smoothing(file.path(tmp, "cv"))
  w <- NULL
  withCallingHandlers(
    PhyloCactus:::.check_cv_subsamples(PhyloCactus:::.read_cv_subsamples(file.path(tmp, "cv.log")),
                                       sel$table, sel$smoothing, "run"),
    warning = function(cnd) { w <<- conditionMessage(cnd); invokeRestart("muffleWarning") })
  expect_false(is.null(w))
  expect_match(w, "), the selected value.", fixed = TRUE)
})

test_that("summarise_iucn_registry() leaves the checklist content out of the registry it writes (M5)", {
  sp <- data.frame(
    species = c("A_a", "A_b"), status_pipeline = c("sampled_and_dated", "unsampled_checklist_gap"),
    in_ml_tree = c(TRUE, FALSE), is_ingroup = c(TRUE, TRUE), in_checklist = c(TRUE, TRUE),
    accepted_in_checklist = c(TRUE, TRUE),
    checklist_pureName = c("A a", "A b"), checklist_fullName = c("A a L.", "A b Mill."),
    checklist_author = c("L.", "Mill."), checklist_rank = c("species", "species"),
    checklist_taxon = c("Taxon", "Taxon"), checklist_uuid = c("u1", "u2"),
    iucn_found = c(TRUE, TRUE), iucn_category = c("CR", "VU"),
    iucn_category_name = c("Critically Endangered", "Vulnerable"), iucn_is_endemic = c(TRUE, FALSE),
    iucn_endemic_locations = c("Mexico", NA), iucn_threat_names = c("Roads", NA),
    iucn_habitat_names = c("Desert", NA), stringsAsFactors = FALSE)
  out <- summarise_iucn_registry(sp, redlist_version = "2026-1", access_date = "2026-10-15")
  dropped <- c("checklist_fullName", "checklist_author", "checklist_rank", "checklist_taxon", "checklist_uuid")
  expect_length(intersect(dropped, names(out$registry)), 0L)
  expect_true(all(c("species", "in_checklist", "accepted_in_checklist", "checklist_pureName") %in% names(out$registry)))
})
