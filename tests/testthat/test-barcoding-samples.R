# Tests of the sample sheet and of the matrix of accepted inputs (Phase 9B; decisions IR7 and IR8 of
# BMM, 28-09). Written before the code. The matrix is a table of the package; the sample sheet routes
# each row by it: an admitted input to its identification function, any other to its message. The
# tests that identify run IdTaxa and skip where it does not run (as in test-barcoding-identify.R).

.smp_matrix <- function() {
  f <- system.file("extdata", "barcoding_accepted_inputs.csv", package = "PhyloCactus")
  utils::read.csv(f, stringsAsFactors = FALSE)
}

test_that("the matrix of accepted inputs is a table of the package with a status, a route or a message per row", {
  m <- .smp_matrix()
  expect_true(all(c("input_type", "description", "format", "status", "function_name", "message", "evidence") %in% names(m)))
  expect_false(anyDuplicated(m$input_type) > 0)
  expect_true(all(m$status %in% c("admitted", "not_tested", "not_admitted", "rejected", "out_of_scope")))
  ad <- m[m$status == "admitted", ]
  expect_setequal(ad$input_type, c("sanger", "plastome", "genbank", "assembly", "reads_paired"))
  expect_true(all(ad$function_name %in% getNamespaceExports("PhyloCactus")))
  expect_true(all(nzchar(m$message[m$status != "admitted"])))
  expect_true(all(c("reads_single", "target_capture", "long_reads", "aligned_reads", "rad_gbs", "rnaseq",
                    "amplicon_mix", "ab1", "skim_nuclear") %in% m$input_type))
})

.smp_sheet <- function(tmp, rows) {
  f <- file.path(tmp, "samples.csv")
  utils::write.csv(rows, f, row.names = FALSE)
  f
}

.smp_row <- function(sample_id, input_type, file_1 = "x", file_2 = NA, declared_species = NA, voucher = "Coll 1") {
  data.frame(sample_id = sample_id, declared_species = declared_species, voucher = voucher, input_type = input_type,
             file_1 = file_1, file_2 = file_2, stringsAsFactors = FALSE)
}

.smp_run <- function(sheet, ...) {
  res <- NULL; msg <- character(0)
  txt <- suppressWarnings(utils::capture.output(withCallingHandlers(
    res <- identify_barcoding_samples(sheet, ...),
    message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") })))
  list(index = res, text = c(txt, msg))
}

test_that("every input that is not admitted is not run and gets the message of its row", {
  tmp <- withr::local_tempdir()
  m <- .smp_matrix()
  other <- m[m$status != "admitted", ]
  sheet <- .smp_sheet(tmp, do.call(rbind, lapply(seq_len(nrow(other)), function(i) .smp_row(paste0("s", i), other$input_type[i]))))
  x <- .smp_run(sheet, library_dir = file.path(tmp, "lib"), output_dir = file.path(tmp, "out"))
  idx <- x$index
  expect_equal(nrow(idx), nrow(other))
  expect_true(all(idx$status == "not_run"))
  expect_identical(idx$message[match(other$input_type, idx$input_type)], other$message)
  expect_true(file.exists(file.path(tmp, "out", "TABLE_barcoding_samples_index.csv")))
})

test_that("a sheet without a required column, or with an unknown input type, stops and says what is valid", {
  tmp <- withr::local_tempdir()
  bad <- .smp_row("s1", "sanger")[, c("sample_id", "input_type")]
  expect_error(identify_barcoding_samples(.smp_sheet(tmp, bad), library_dir = tmp, output_dir = file.path(tmp, "o1")), "file_1")
  unk <- .smp_row("s1", "nanopore_magic")
  expect_error(identify_barcoding_samples(.smp_sheet(tmp, unk), library_dir = tmp, output_dir = file.path(tmp, "o2")),
               "sanger")
  dup <- rbind(.smp_row("s1", "ab1"), .smp_row("s1", "ab1"))
  expect_error(identify_barcoding_samples(.smp_sheet(tmp, dup), library_dir = tmp, output_dir = file.path(tmp, "o3")),
               "sample_id")
})

# ---- Identification of admitted rows (IdTaxa) ------------------------------------------------------

.smp_idtaxa_runs <- function() {
  tryCatch({
    tr <- Biostrings::DNAStringSet(c(a = "ACGTACGTAACCGGTTAACC", b = "ACGTACGTACCCGGTTAACC",
                                     c = "TTTTGGGGAATTCCGGAATT", d = "TTTTGGGGAATTCCGGAATG"))
    lt <- DECIPHER::LearnTaxa(tr, c("Root;Opuntia;Opuntia_robusta", "Root;Opuntia;Opuntia_stricta",
                                    "Root;Cereus;Cereus_jamacaru", "Root;Cereus;Cereus_horrida"), verbose = FALSE)
    DECIPHER::IdTaxa(Biostrings::DNAStringSet("ACGTACGTAACCGGTTAACC"), lt, strand = "top", threshold = 60,
                     processors = 1L, verbose = FALSE)
    TRUE
  }, error = function(e) FALSE)
}

.smp_skip <- function() {
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("Biostrings")
  if (!suppressWarnings(suppressMessages(.smp_idtaxa_runs()))) {
    skip("DECIPHER::IdTaxa() does not run in this loading context: match() over XStringSet fails to dispatch inside IdTaxa. The same call works in a plain R session.")
  }
}

.smp_random <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
.smp_vary <- function(s, k) {
  for (p in sample(seq_len(nchar(s)), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1)
  s
}
.smp_library <- function(tmp) {
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir, recursive = TRUE)
  set.seed(41L)
  seqs <- list()
  for (l in c("matK", "rbcL")) {
    root <- .smp_random(300); x <- character(0)
    for (g in c("Opuntia", "Cereus")) {
      gb <- .smp_vary(root, 6)
      for (e in c("alpha", "beta")) {
        sb <- .smp_vary(gb, 3)
        for (r in 1:3) x[paste0(g, "_", e, "|", l, "_", substr(g, 1, 1), e, r, ".1")] <- .smp_vary(sb, 1)
      }
    }
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(lib_dir, paste0("LIB_", l, ".fasta")))
    seqs[[l]] <- x
  }
  list(dir = lib_dir, seqs = seqs)
}

test_that("a Sanger sample is identified as one sample, with its report, and compared with the declared species", {
  .smp_skip()
  tmp <- withr::local_tempdir()
  lib <- .smp_library(tmp)
  fa <- file.path(tmp, "spec1.fasta")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(a = unname(lib$seqs$matK[1]), b = unname(lib$seqs$rbcL[1]))), fa)
  fb <- file.path(tmp, "spec2.fasta")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(a = unname(lib$seqs$matK[1]))), fb)
  sheet <- .smp_sheet(tmp, rbind(.smp_row("spec1", "sanger", fa, declared_species = "Opuntia_alpha"),
                                 .smp_row("spec2", "sanger", fb, declared_species = "Cereus_beta", voucher = NA),
                                 .smp_row("spec3", "target_capture")))
  x <- .smp_run(sheet, library_dir = lib$dir, metrics_dir = NULL, output_dir = file.path(tmp, "out"))
  idx <- x$index
  expect_identical(idx$status, c("identified", "identified", "not_run"))
  expect_true(file.exists(file.path(tmp, "out", "spec1", "REPORT_spec1.html")))
  expect_equal(idx$loci_found[idx$sample_id == "spec1"], 2L)
  # The states are counted over the regions found, not over the loci absent from the sequences
  expect_equal(idx$state1 + idx$state2 + idx$state3, idx$loci_found)
  expect_false(idx$voucher_missing[idx$sample_id == "spec1"])
  expect_true(idx$voucher_missing[idx$sample_id == "spec2"])
  # The answer against the declared species
  expect_true(idx$declared_comparison[idx$sample_id == "spec1"] %in% c("species agrees", "genus agrees"))
  expect_true(idx$declared_comparison[idx$sample_id == "spec2"] %in% c("disagrees", "not assigned"))
  rec <- utils::read.csv(file.path(tmp, "out", "spec1", "RUN_spec1.csv"), stringsAsFactors = FALSE)
  expect_identical(rec$value[rec$name == "declared_species"], "Opuntia_alpha")
  expect_identical(rec$value[rec$name == "voucher"], "Coll 1")
  h <- paste(readLines(file.path(tmp, "out", "spec1", "REPORT_spec1.html"), warn = FALSE), collapse = "\n")
  top <- sub("id=\"section-2\".*", "", h)
  expect_true(grepl("Opuntia_alpha", top, fixed = TRUE))
  expect_true(grepl("Coll 1", top, fixed = TRUE))
  h2 <- paste(readLines(file.path(tmp, "out", "spec2", "REPORT_spec2.html"), warn = FALSE), collapse = "\n")
  expect_true(grepl("no voucher", h2, fixed = TRUE))
})

test_that("an error in one sample is recorded and the other samples are still identified", {
  .smp_skip()
  tmp <- withr::local_tempdir()
  lib <- .smp_library(tmp)
  fa <- file.path(tmp, "ok.fasta")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(a = unname(lib$seqs$matK[1]))), fa)
  sheet <- .smp_sheet(tmp, rbind(.smp_row("bad", "sanger", file.path(tmp, "missing.fasta")), .smp_row("ok", "sanger", fa)))
  x <- .smp_run(sheet, library_dir = lib$dir, metrics_dir = NULL, output_dir = file.path(tmp, "out"))
  expect_identical(x$index$status, c("error", "identified"))
  expect_true(nzchar(x$index$message[1]))
})
