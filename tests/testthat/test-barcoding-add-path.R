# Tests of the alignment path of step 7 (11_barcoding/7_classifier/). Written before the code.
# Phase 6A of PhyloCactus 0.5.0, decision D5 of BMM and the probe of 2026-09-25.
#
# The distances of the legitimate queries of Phase 5A come from the matrix of the joint alignment of
# the library, in which the query took part. The distances of the alien queries of CN2, and of any
# query a user will bring, come from MAFFT --add --keeplength against the library, one call per
# query. The probe of 2026-09-25 measured both on 720 legitimate queries: 699 identical to rounding,
# 21 with a different distance, and 4 that changed their answer. In pepC_like the joint alignment
# left JN387243.1 at distance 0 from several species of Opuntia (state 2, correct genus); added on
# its own it sits at 0.0086 and the tie breaks towards a wrong species. The threshold of Phase 6A has
# to be derived from distances measured the way the identification will measure them, so step 7
# learns to classify the folds by that path. The tables of Phase 5A are not touched.

.add_mafft_ok <- function() {
  isTRUE(tryCatch(.bc_assert_mafft("mafft"), error = function(e) FALSE))
}

.add_fixture <- function(tmp, indels = FALSE) {
  # Six sequences per locus, two of a species, two of another of the same genus, two of another
  # genus. Without indels the alignment is the identity and both paths have to agree exactly.
  # Divergences of 3 to 7 %, not the 8 to 20 % of the fixture of step 7: the orientation of the add
  # path looks for shared 20-mers, and at 25 substitutions in 300 bases almost none survives, so a
  # legitimate query came back as no_match. On the real data none of the 720 legitimate
  # queries of the probe of 2026-09-25 did; the first version of this fixture was the outlier.
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(21L)
  base <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  vary <- function(s, k) { for (p in sample(seq_len(nchar(s)), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1); s }
  x <- c(`Opuntia_robusta|A1.1` = vary(base, 2), `Opuntia_robusta|A2.1` = vary(base, 3),
         `Opuntia_stricta|B1.1` = vary(base, 9), `Opuntia_stricta|B2.1` = vary(base, 10),
         `Cereus_jamacaru|C1.1` = vary(base, 18), `Cereus_horrida|C2.1` = vary(base, 20))
  if (indels) {
    # A deletion in one sequence of each species, so the joint alignment has gaps to place
    x[c(1, 3, 5)] <- vapply(x[c(1, 3, 5)], function(s) paste0(substr(s, 1, 120), substr(s, 131, 300)),
                            character(1))
    raw <- file.path(tmp, "raw.fasta")
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), raw)
    run_mafft(raw, file.path(lib_dir, "LIB_matK.fasta"), mafft_exec = "mafft", mafft_opts = "--auto")
  } else {
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(lib_dir, "LIB_matK.fasta"))
  }
  folds_dir <- file.path(tmp, "5_folds")
  suppressMessages(build_barcoding_folds(library_dir = lib_dir, output_dir = folds_dir))
  list(library_dir = lib_dir, folds_dir = folds_dir)
}

test_that("the training alignment of a fold holds only the training set, and no column of gaps only", {
  skip_if_not_installed("ape")
  m <- rbind(a = c("a", "c", "-", "g", "t"),
             b = c("a", "c", "-", "g", "-"),
             q = c("a", "c", "t", "g", "t"))
  dna <- ape::as.DNAbin(m)
  train_aln <- .bc_training_alignment(dna, c("a", "b"))
  expect_equal(rownames(train_aln), c("a", "b"))
  # The third column was open only for the query: with the query out it has no base left
  expect_equal(ncol(train_aln), 4L)
  expect_true(all(colSums(as.character(train_aln) != "-") > 0))
  expect_error(.bc_training_alignment(dna, c("a", "z")), "z")
})

test_that("without indels both paths give the same state, species and distance for every query", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  skip_if_not(.add_mafft_ok(), "MAFFT is not available")
  tmp <- withr::local_tempdir()
  f <- .add_fixture(tmp)
  out_dir <- file.path(tmp, "7_classifier")

  suppressMessages(classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                            output_dir = out_dir, method = "nn"))
  suppressMessages(classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                            output_dir = out_dir, method = "nn", alignment = "add"))

  for (sc in c("species", "genus")) {
    lib <- utils::read.csv(file.path(out_dir, paste0("TABLE_barcoding_predictions_", sc, "_nn.csv")),
                           stringsAsFactors = FALSE)
    add <- utils::read.csv(file.path(out_dir, paste0("TABLE_barcoding_predictions_", sc, "_nn_add.csv")),
                           stringsAsFactors = FALSE)
    expect_equal(add$sid, lib$sid)
    expect_equal(add$state, lib$state)
    expect_equal(add$predicted_species, lib$predicted_species)
    expect_equal(add$predicted_genus, lib$predicted_genus)
    expect_equal(add$nn_distance, lib$nn_distance, tolerance = 1e-9)
  }
})

test_that("the add path writes its own tables and leaves the tables of the library path untouched", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  skip_if_not(.add_mafft_ok(), "MAFFT is not available")
  tmp <- withr::local_tempdir()
  f <- .add_fixture(tmp, indels = TRUE)
  out_dir <- file.path(tmp, "7_classifier")

  suppressMessages(classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                            output_dir = out_dir, method = "nn"))
  previous <- list.files(out_dir, pattern = "_nn\\.csv$", full.names = TRUE)
  md5_before <- tools::md5sum(previous)

  suppressMessages(classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                            output_dir = out_dir, method = "nn", alignment = "add"))

  expect_identical(tools::md5sum(previous), md5_before)
  for (sc in c("species", "genus")) {
    add <- utils::read.csv(file.path(out_dir, paste0("TABLE_barcoding_predictions_", sc, "_nn_add.csv")),
                           stringsAsFactors = FALSE)
    folds <- utils::read.csv(file.path(f$folds_dir, paste0("TABLE_barcoding_folds_", sc, ".csv")),
                             stringsAsFactors = FALSE)
    # One row per query of the folds, and the path declared in every row
    expect_setequal(paste(add$locus, add$sid), paste(folds$locus, folds$sid))
    expect_true(all(c("alignment", "orientation") %in% names(add)))
    expect_true(all(add$alignment == "add"))
    expect_true(all(add$orientation %in% c("forward", "reverse", "no_match")))
  }
})

test_that("a legitimate query with no homology to its training set is state 3, as in CN2", {
  skip_if_not_installed("ape")
  skip_if_not(.add_mafft_ok(), "MAFFT is not available")
  set.seed(22L)
  base <- paste(sample(c("a", "c", "g", "t"), 300, replace = TRUE), collapse = "")
  vary <- function(s, k) { for (p in sample(seq_len(300), k)) substr(s, p, p) <- sample(c("a", "c", "g", "t"), 1); s }
  train_aln <- c(A1 = vary(base, 2), A2 = vary(base, 3), B1 = vary(base, 25), B2 = vary(base, 28))
  train_aln <- ape::as.DNAbin(do.call(rbind, strsplit(train_aln, "")))
  species <- c(A1 = "Opuntia_robusta", A2 = "Opuntia_robusta", B1 = "Opuntia_stricta", B2 = "Opuntia_stricta")
  alien <- paste(sample(c("a", "c", "g", "t"), 300, replace = TRUE), collapse = "")

  r <- .bc_classify_by_add(train_aln, species, alien, model = "raw", min_comparable = 100L)
  expect_equal(r$orientation, "no_match")
  expect_equal(r$state, 3L)
  expect_equal(r$reason, "no_match")
  expect_true(is.na(r$predicted_species))
  expect_true(is.na(r$nn_distance))
})

test_that("CN2 and the add path of step 7 are one path: the same query gets the same answer", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  skip_if_not(.add_mafft_ok(), "MAFFT is not available")
  tmp <- withr::local_tempdir()
  f <- .add_fixture(tmp, indels = TRUE)
  lib <- ape::read.dna(file.path(f$library_dir, "LIB_matK.fasta"), format = "fasta", as.matrix = TRUE)
  h <- .bc_parse_header(rownames(lib))
  rownames(lib) <- h$sid
  species <- stats::setNames(h$species, h$sid)

  set.seed(23L)
  q <- gsub("-", "", paste(as.character(lib["A2.1", ]), collapse = ""), fixed = TRUE)
  for (p in sample(seq_len(nchar(q)), 10)) substr(q, p, p) <- sample(c("a", "c", "g", "t"), 1)
  og <- ape::as.DNAbin(list(`Portulaca_amilis|X1.1` = strsplit(q, "")[[1]]))

  a <- .bc_classify_by_add(lib, species, q, model = "raw", min_comparable = 100L)
  b <- .bc_cn2_queries(og, lib, data.frame(sid = h$sid, species = h$species, stringsAsFactors = FALSE),
                       model = "raw", min_comparable = 100L, locus = "matK")
  expect_equal(a$state, b$state)
  expect_equal(a$predicted_species, b$predicted_species)
  expect_equal(a$nn_distance, b$nn_distance)
  expect_equal(a$orientation, b$orientation)
})

test_that("the add path needs MAFFT and says so before classifying anything, and refuses IdTaxa", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  f <- .add_fixture(tmp)
  out_dir <- file.path(tmp, "7_classifier")

  expect_error(classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                        output_dir = out_dir, method = "nn", alignment = "add",
                                        mafft_exec = "no_such_mafft_binary"),
               "MAFFT")
  expect_false(any(file.exists(file.path(out_dir, c("TABLE_barcoding_predictions_species_nn_add.csv",
                                                    "TABLE_barcoding_predictions_genus_nn_add.csv")))))
  # IdTaxa does not align, so the add path means nothing for it
  expect_error(classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                        output_dir = out_dir, method = "idtaxa", alignment = "add"),
               "nn")
})

test_that("step 7 writes its running time per locus and scheme, for both paths", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  f <- .add_fixture(tmp)
  out_dir <- file.path(tmp, "7_classifier")

  suppressMessages(classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                            output_dir = out_dir, method = "nn"))
  t <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_timing_nn.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("locus", "scheme", "method", "alignment", "folds", "queries",
                    "seconds") %in% names(t)))
  expect_equal(nrow(t), 2L)
  expect_setequal(t$scheme, c("species", "genus"))
  expect_true(all(t$seconds >= 0))
  expect_true(all(t$alignment == "library"))

  skip_if_not(.add_mafft_ok(), "MAFFT is not available")
  suppressMessages(classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                            output_dir = out_dir, method = "nn", alignment = "add"))
  ta <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_timing_nn_add.csv"), stringsAsFactors = FALSE)
  expect_true(all(ta$alignment == "add"))
  # The timing of one path does not overwrite the other
  expect_true(file.exists(file.path(out_dir, "TABLE_barcoding_timing_nn.csv")))
})
