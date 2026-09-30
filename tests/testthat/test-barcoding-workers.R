# Tests of the parallel run of step 7 on one machine (PA1 to PA3 of BMM, 29-09). Written before the
# code. With workers > 1, classify_barcoding_folds() starts that many R processes, one chunk each
# (the chunks of 6B), watches them, merges them with merge_barcoding_chunks(), writes a progress
# file, sends one email at the end, and stops the others if one fails.

.wk_fixture <- function(tmp) {
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir, showWarnings = FALSE)
  set.seed(9L)
  base <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  vary <- function(s, k) { for (p in sample(seq_len(300), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1); s }
  for (l in c("matK", "rbcL")) {
    x <- c(`Opuntia_robusta|A1.1` = vary(base, 2), `Opuntia_robusta|A2.1` = vary(base, 3),
           `Opuntia_stricta|B1.1` = vary(base, 25), `Opuntia_stricta|B2.1` = vary(base, 28),
           `Cereus_jamacaru|C1.1` = vary(base, 60), `Cereus_jamacaru|C3.1` = vary(base, 61),
           `Cereus_horrida|C2.1` = vary(base, 62))
    names(x) <- sub("\\|", paste0("|", l, "_"), names(x))
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(lib_dir, paste0("LIB_", l, ".fasta")))
  }
  folds_dir <- file.path(tmp, "5_folds")
  suppressMessages(build_barcoding_folds(library_dir = lib_dir, output_dir = folds_dir))
  list(library_dir = lib_dir, folds_dir = folds_dir)
}

.wk_skip <- function() {
  skip_on_cran()
  skip_if(!nzchar(Sys.which(file.path(R.home("bin"), "Rscript"))) && !file.exists(file.path(R.home("bin"), "Rscript")),
          "Rscript not found")
  # The workers load the installed package; skip when it lacks the chunks of 6B
  skip_if(!"n_chunks" %in% names(formals(get("classify_barcoding_folds", envir = asNamespace("PhyloCactus")))))
}

test_that("with workers, step 7 gives the same tables as one process, and says it merged them (PA1)", {
  .wk_skip()
  withr::local_options(PhyloCactus.poll_seconds = 1)
  tmp <- withr::local_tempdir()
  f <- .wk_fixture(tmp)
  one <- file.path(tmp, "7_one"); par <- file.path(tmp, "7_par")
  suppressMessages(utils::capture.output(classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                                                  output_dir = one, method = "nn")))
  suppressMessages(utils::capture.output(classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                                                  output_dir = par, method = "nn", workers = 2L)))
  for (sc in c("species", "genus")) {
    a <- utils::read.csv(file.path(one, paste0("TABLE_barcoding_predictions_", sc, "_nn.csv")))
    b <- utils::read.csv(file.path(par, paste0("TABLE_barcoding_predictions_", sc, "_nn.csv")))
    key <- c("locus", "scheme", "fold", "sid")
    a <- a[do.call(order, a[key]), ]; b <- b[do.call(order, b[key]), ]
    rownames(a) <- rownames(b) <- NULL
    expect_identical(a, b)
  }
  expect_true(file.exists(file.path(par, "chunks", "TABLE_barcoding_predictions_species_nn_chunk1of2.csv")))
})

test_that("with workers, step 7 writes a progress file that ends as finished (PA2)", {
  .wk_skip()
  withr::local_options(PhyloCactus.poll_seconds = 1)
  tmp <- withr::local_tempdir()
  f <- .wk_fixture(tmp)
  out <- file.path(tmp, "7_par")
  suppressMessages(utils::capture.output(classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                                                  output_dir = out, method = "nn", workers = 2L)))
  p <- file.path(out, "PROGRESS_nn.txt")
  expect_true(file.exists(p))
  txt <- paste(readLines(p), collapse = "\n")
  expect_match(txt, "finished", fixed = TRUE)
  expect_match(txt, "2 of 2 workers", fixed = TRUE)
})

test_that("a worker that fails stops the run with its error and the notification says FAILED (PA3)", {
  .wk_skip()
  withr::local_options(PhyloCactus.poll_seconds = 1)
  withr::local_envvar(PHYLOCACTUS_TEST_FAIL_CHUNK = "2", MY_EMAIL = "")
  tmp <- withr::local_tempdir()
  f <- .wk_fixture(tmp)
  out <- file.path(tmp, "7_par")
  expect_error(suppressMessages(utils::capture.output(
    classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir, output_dir = out,
                             method = "nn", workers = 2L, notify = TRUE))), "chunk 2")
  expect_match(paste(readLines(file.path(out, "PROGRESS_nn.txt")), collapse = "\n"), "failed", fixed = TRUE)
  expect_false(file.exists(file.path(out, "TABLE_barcoding_predictions_species_nn.csv")))
})

test_that("workers = 1 is the run of one process, and workers cannot be combined with chunk (PA1)", {
  expect_identical(formals(classify_barcoding_folds)$workers, 1L)
  expect_error(classify_barcoding_folds(workers = 2L, chunk = 1L, n_chunks = 2L), "workers")
})

test_that("Tutorial 5 runs nn_add with workers and declares n_workers (PA1)", {
  f <- system.file("scripts", "tutorial-5-cactus-phylogeny-barcoding.R", package = "PhyloCactus")
  skip_if(!nzchar(f), "Tutorial 5 not installed")
  ex <- as.list(parse(f))
  lhs <- vapply(Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")), ex), function(e) deparse(e[[2]]), "")
  expect_true("n_workers" %in% lhs)
  calls <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("classify_barcoding_folds")), ex)
  add <- Filter(function(e) identical(e$alignment, "add"), calls)
  expect_length(add, 1L)
  expect_identical(add[[1]]$workers, as.name("n_workers"))
})
