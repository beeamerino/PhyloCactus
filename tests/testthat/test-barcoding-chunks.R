# Tests of the chunked run of step 7 and of the SLURM job that runs it (11_barcoding/7_classifier/).
# Written before the code. Phase 6B of PhyloCactus 0.5.0, decisions K1, K4 and K5 of BMM (26-09).
#
# The full IdTaxa run takes about 34 hours of one core and runs on Leftraru as a job array: the
# folds are dealt to N chunks, each chunk is one task, and a last job merges the chunk tables. The
# merged tables must be the tables a single run would have written, byte for byte, and the merge
# must refuse anything else. None of this depends on the classifier, so it is tested here with the
# nearest neighbour, which runs under devtools::test(); the same checks with IdTaxa are in
# test-barcoding-idtaxa.R.

.chk_fixture <- function(tmp) {
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(9L)
  base <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  vary <- function(s, k) { for (p in sample(seq_len(300), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1); s }
  for (l in c("matK", "rbcL")) {
    x <- c(`Opuntia_robusta|A1.1` = vary(base, 2), `Opuntia_robusta|A2.1` = vary(base, 3),
           `Opuntia_robusta|A3.1` = vary(base, 4),
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

.chk_run <- function(f, out_dir, ...) {
  suppressMessages(utils::capture.output(
    classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                             output_dir = out_dir, method = "nn", ...)))
  invisible(out_dir)
}

.chk_md5 <- function(dir, suffix = "nn") {
  unname(tools::md5sum(file.path(dir, paste0("TABLE_barcoding_predictions_", c("species", "genus"),
                                             "_", suffix, ".csv"))))
}

# ---- K1: the seed of a query ------------------------------------------------------------------

# Reference of the contract, written independently of the package code: a polynomial hash of the
# four fields joined by a carriage return, base 31, modulo the Mersenne prime 2^31 - 1, combined
# with the seed of the run. Pinning the value in the test is what makes the seed of a query the
# same on the Mac, on Leftraru and in a later version of R.
.chk_reference_seed <- function(seed, locus, scheme, fold, sid) {
  m <- 2147483647
  key <- paste(locus, scheme, as.character(fold), sid, sep = "\r")
  h <- Reduce(function(acc, b) (acc * 31 + b) %% m, as.numeric(utf8ToInt(enc2utf8(key))), 0)
  s <- ((as.numeric(seed) %% m) * 1000003 + h) %% m
  as.integer(if (s == 0) 1 else s)
}

test_that("the seed of a query is a pure function of the seed, the locus, the scheme, the fold and the sid", {
  a <- .bc_query_seed(1L, "matK", "species", 3L, "matK_A1.1")

  expect_true(is.integer(a))
  expect_length(a, 1L)
  expect_true(a >= 1L)
  expect_identical(a, .chk_reference_seed(1L, "matK", "species", 3L, "matK_A1.1"))

  # Whatever the state of the session's generator, and without touching it
  set.seed(99L); stats::runif(3)
  before <- get(".Random.seed", envir = globalenv())
  expect_identical(.bc_query_seed(1L, "matK", "species", 3L, "matK_A1.1"), a)
  expect_true(identical(get(".Random.seed", envir = globalenv()), before))

  # Any field that changes gives another seed
  expect_false(identical(.bc_query_seed(2L, "matK", "species", 3L, "matK_A1.1"), a))
  expect_false(identical(.bc_query_seed(1L, "rbcL", "species", 3L, "matK_A1.1"), a))
  expect_false(identical(.bc_query_seed(1L, "matK", "genus", 3L, "matK_A1.1"), a))
  expect_false(identical(.bc_query_seed(1L, "matK", "species", 4L, "matK_A1.1"), a))
  expect_false(identical(.bc_query_seed(1L, "matK", "species", 3L, "matK_A2.1"), a))
})

# ---- K4: chunks ---------------------------------------------------------------------------------

test_that("the folds of a locus and scheme are dealt round-robin, each to one chunk, sizes within one", {
  for (n in 0:10) for (k in 1:4) {
    a <- .bc_chunk_assign(n, k)
    expect_length(a, n)
    expect_true(all(a %in% seq_len(k)))
    expect_identical(a, as.integer((seq_len(n) - 1L) %% k + 1L))
    if (n > 0L) {
      sizes <- tabulate(a, nbins = k)
      expect_lte(max(sizes) - min(sizes), 1L)
    }
  }
})

test_that("the chunk arguments are checked before anything is read", {
  tmp <- withr::local_tempdir()
  expect_error(classify_barcoding_folds(library_dir = tmp, folds_dir = tmp, output_dir = file.path(tmp, "o"),
                                        chunk = 1L), "n_chunks")
  expect_error(classify_barcoding_folds(library_dir = tmp, folds_dir = tmp, output_dir = file.path(tmp, "o"),
                                        chunk = 4L, n_chunks = 3L), "chunk")
  expect_error(classify_barcoding_folds(library_dir = tmp, folds_dir = tmp, output_dir = file.path(tmp, "o"),
                                        chunk = 1L, n_chunks = 0L), "n_chunks")
})

test_that("three chunks merged give the tables of one run without chunks, byte for byte", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  f <- .chk_fixture(tmp)
  one <- .chk_run(f, file.path(tmp, "one"))
  parts <- file.path(tmp, "parts")
  for (k in 1:3) .chk_run(f, parts, chunk = k, n_chunks = 3L)

  # Each chunk writes its own tables, its timing and its session, in chunks/, and nothing else
  for (k in 1:3) {
    tag <- paste0("_chunk", k, "of3")
    expect_true(file.exists(file.path(parts, "chunks", paste0("TABLE_barcoding_predictions_species_nn", tag, ".csv"))))
    expect_true(file.exists(file.path(parts, "chunks", paste0("TABLE_barcoding_predictions_genus_nn", tag, ".csv"))))
    expect_true(file.exists(file.path(parts, "chunks", paste0("TABLE_barcoding_timing_nn", tag, ".csv"))))
    expect_true(file.exists(file.path(parts, "chunks", paste0("SESSION_barcoding_nn", tag, ".txt"))))
  }
  expect_false(file.exists(file.path(parts, "TABLE_barcoding_predictions_species_nn.csv")))

  suppressMessages(utils::capture.output(
    merge_barcoding_chunks(library_dir = f$library_dir, folds_dir = f$folds_dir, output_dir = parts,
                           method = "nn", n_chunks = 3L)))

  expect_identical(.chk_md5(parts), .chk_md5(one))
  tim <- utils::read.csv(file.path(parts, "TABLE_barcoding_timing_nn.csv"), stringsAsFactors = FALSE)
  tim1 <- utils::read.csv(file.path(one, "TABLE_barcoding_timing_nn.csv"), stringsAsFactors = FALSE)
  expect_identical(tim[, c("locus", "scheme", "folds", "queries")], tim1[, c("locus", "scheme", "folds", "queries")])
})

test_that("the merge refuses a missing chunk, a fold twice, a fold missing and a fold not in step 5", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  f <- .chk_fixture(tmp)
  fresh <- function(name) {
    d <- file.path(tmp, name)
    for (k in 1:3) .chk_run(f, d, chunk = k, n_chunks = 3L)
    d
  }
  merge <- function(d) merge_barcoding_chunks(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                              output_dir = d, method = "nn", n_chunks = 3L)
  chunk_file <- function(d, k, sc = "species") {
    file.path(d, "chunks", paste0("TABLE_barcoding_predictions_", sc, "_nn_chunk", k, "of3.csv"))
  }

  d <- fresh("missing")
  file.remove(chunk_file(d, 2L, "genus"))
  expect_error(suppressMessages(merge(d)), "chunk 2")

  d <- fresh("twice")
  c1 <- utils::read.csv(chunk_file(d, 1L), stringsAsFactors = FALSE)
  c2 <- utils::read.csv(chunk_file(d, 2L), stringsAsFactors = FALSE)
  utils::write.csv(rbind(c2, c1[1, ]), chunk_file(d, 2L), row.names = FALSE)
  expect_error(suppressMessages(merge(d)), "more than once")

  d <- fresh("absent")
  c1 <- utils::read.csv(chunk_file(d, 1L), stringsAsFactors = FALSE)
  gone <- c1$fold[1]
  utils::write.csv(c1[c1$fold != gone | c1$locus != c1$locus[1], ], chunk_file(d, 1L), row.names = FALSE)
  expect_error(suppressMessages(merge(d)), "missing")

  d <- fresh("foreign")
  c1 <- utils::read.csv(chunk_file(d, 1L), stringsAsFactors = FALSE)
  c1$fold[1] <- 999L
  utils::write.csv(c1, chunk_file(d, 1L), row.names = FALSE)
  expect_error(suppressMessages(merge(d)), "999")
})

# ---- K5: the SLURM job --------------------------------------------------------------------------

test_that("the SLURM job runs one chunk per task and merges after the array, with nothing hard-coded", {
  tmp <- withr::local_tempdir()
  job_dir <- file.path(tmp, "job dir")

  files <- write_barcoding_slurm_job(
    job_dir = job_dir, n_chunks = 12L,
    library_dir = "/data/11_barcoding/4_library", folds_dir = "/data/11_barcoding/5_folds",
    output_dir = "/data/11_barcoding/7_classifier", method = "idtaxa",
    partition = "general", account = "acc123", time = "04:00:00", mem = "4G",
    r_setup = "module load R/4.4.1", lib_path = "/home/u/R/lib", email = "user@example.org")

  expect_setequal(basename(files), c("job_array.sh", "job_merge.sh", "submit.sh", "run_chunk.R", "run_merge.R"))
  expect_true(all(file.exists(files)))
  read <- function(x) paste(readLines(file.path(job_dir, x)), collapse = "\n")
  arr <- read("job_array.sh"); mrg <- read("job_merge.sh"); sub <- read("submit.sh")
  rc <- read("run_chunk.R"); rm_ <- read("run_merge.R")

  expect_match(arr, "#SBATCH --array=1-12", fixed = TRUE)
  expect_match(arr, "#SBATCH --partition=general", fixed = TRUE)
  expect_match(arr, "#SBATCH --account=acc123", fixed = TRUE)
  expect_match(arr, "#SBATCH --time=04:00:00", fixed = TRUE)
  expect_match(arr, "#SBATCH --mem=4G", fixed = TRUE)
  expect_match(arr, "#SBATCH --cpus-per-task=1", fixed = TRUE)
  expect_match(arr, "module load R/4.4.1", fixed = TRUE)
  expect_match(arr, "SLURM_ARRAY_TASK_ID", fixed = TRUE)
  expect_match(arr, "user@example.org", fixed = TRUE)
  # The path with a space reaches bash quoted
  expect_match(arr, shQuote(file.path(job_dir, "run_chunk.R")), fixed = TRUE)
  expect_match(sub, "--parsable", fixed = TRUE)
  expect_match(sub, "--dependency=afterok:", fixed = TRUE)
  expect_match(mrg, "run_merge.R", fixed = TRUE)

  # The R scripts are valid R and carry the arguments of the run
  expect_silent(parse(file.path(job_dir, "run_chunk.R")))
  expect_silent(parse(file.path(job_dir, "run_merge.R")))
  expect_match(rc, "classify_barcoding_folds", fixed = TRUE)
  expect_match(rc, "n_chunks = 12L", fixed = TRUE)
  expect_match(rc, "\"idtaxa\"", fixed = TRUE)
  expect_match(rc, "/data/11_barcoding/5_folds", fixed = TRUE)
  expect_match(rc, "/home/u/R/lib", fixed = TRUE)
  expect_match(rm_, "merge_barcoding_chunks", fixed = TRUE)

  # Nothing about a particular cluster is written unless it was passed
  all_text <- paste(arr, mrg, sub, rc, rm_)
  expect_false(grepl("leftraru|nlhpc", all_text, ignore.case = TRUE))

  # An optional value that is not passed leaves no line behind
  files2 <- write_barcoding_slurm_job(
    job_dir = file.path(tmp, "job2"), n_chunks = 2L, library_dir = "l", folds_dir = "f",
    output_dir = "o", method = "nn", partition = "general", time = "01:00:00", mem = "2G",
    r_setup = "module load R")
  arr2 <- paste(readLines(file.path(tmp, "job2", "job_array.sh")), collapse = "\n")
  expect_false(grepl("--account", arr2, fixed = TRUE))
  expect_false(grepl("--mail-user", arr2, fixed = TRUE))

  expect_error(write_barcoding_slurm_job(job_dir = file.path(tmp, "job3"), n_chunks = 0L,
                                         library_dir = "l", folds_dir = "f", output_dir = "o",
                                         partition = "general", time = "01:00:00", mem = "2G",
                                         r_setup = "module load R"),
               "n_chunks")
})
