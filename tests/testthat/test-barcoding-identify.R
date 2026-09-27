# Tests of step 10 of the molecular diagnostic branch: identification of a user query. Written
# before the code. Phase 6 bis of PhyloCactus 0.5.0, decisions I1 to I9 of BMM (27-09).
#
# IdTaxa does not run under devtools::test() (5A minutes, section 4). These tests carry the same
# probe as test-barcoding-idtaxa.R and run under R CMD check and with the installed package.

.idn_runs <- function() {
  tryCatch({
    tr <- Biostrings::DNAStringSet(c(a = "ACGTACGTAACCGGTTAACC", b = "ACGTACGTACCCGGTTAACC",
                                     c = "TTTTGGGGAATTCCGGAATT", d = "TTTTGGGGAATTCCGGAATG"))
    lt <- DECIPHER::LearnTaxa(tr, c("Root;Opuntia;Opuntia_robusta", "Root;Opuntia;Opuntia_stricta",
                                    "Root;Cereus;Cereus_jamacaru", "Root;Cereus;Cereus_horrida"),
                              verbose = FALSE)
    DECIPHER::IdTaxa(Biostrings::DNAStringSet("ACGTACGTAACCGGTTAACC"), lt, strand = "top",
                     threshold = 60, processors = 1L, verbose = FALSE)
    TRUE
  }, error = function(e) FALSE)
}

.idn_skip <- function() {
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  if (!suppressWarnings(suppressMessages(.idn_runs()))) {
    skip("DECIPHER::IdTaxa() does not run in this loading context: match() over XStringSet fails to dispatch inside IdTaxa. The same call works in a plain R session.")
  }
}

.idn_random <- function(n, gc = 0.5) {
  paste(sample(c("A", "C", "G", "T"), n, replace = TRUE,
               prob = c(1 - gc, gc, gc, 1 - gc) / 2), collapse = "")
}

.idn_vary <- function(s, k) {
  for (p in sample(seq_len(nchar(s)), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1)
  s
}

.idn_rc <- function(s) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(s)))

# Three loci of 300 bases, each from its own random root, so that no locus shares 20-mers with
# another; three genera, two species each, three accessions per species.
.idn_fixture <- function(tmp) {
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(41L)
  genera <- c("Opuntia", "Cereus", "Mammillaria")
  epithets <- c("alpha", "beta")
  seqs <- list()
  for (l in c("matK", "rbcL", "trnL-trnF")) {
    root <- .idn_random(300)
    x <- character(0)
    for (g in seq_along(genera)) {
      gb <- .idn_vary(root, 40)
      for (e in seq_along(epithets)) {
        sb <- .idn_vary(gb, 10)
        for (r in 1:3) {
          x[paste0(genera[g], "_", epithets[e], "|", l, "_", substr(genera[g], 1, 1), e, r, ".1")] <-
            .idn_vary(sb, 2)
        }
      }
    }
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(lib_dir, paste0("LIB_", l, ".fasta")))
    seqs[[l]] <- x
  }
  list(library_dir = lib_dir, seqs = seqs, out = file.path(tmp, "10_identify"))
}

.idn_run <- function(f, query, ...) {
  res <- NULL
  txt <- suppressWarnings(utils::capture.output(
    msg <- utils::capture.output(
      res <- identify_barcoding_query(query, library_dir = f$library_dir,
                                      metrics_dir = file.path(dirname(f$library_dir), "11_metrics"),
                                      output_dir = f$out, ...),
      type = "message")))
  list(tab = res, text = c(txt, msg))
}

.idn_columns <- c("query", "locus", "query_length", "path", "region_start", "region_end",
                  "region_length", "other_windows", "orientation", "state", "predicted_species",
                  "predicted_genus", "candidates", "genus_idtaxa", "species_idtaxa",
                  "genus_confidence", "species_confidence", "reason", "threshold",
                  "validation_ws_rate_species_present", "validation_ws_rate_species_absent")

# ---- I3 to I5: the Sanger path --------------------------------------------------------------------

test_that("a query of a library species gets one row per locus, with the three states only", {
  .idn_skip()
  f <- .idn_fixture(withr::local_tempdir())
  q <- c(q1 = unname(f$seqs$matK[1]))
  r <- .idn_run(f, q)$tab
  expect_true(all(.idn_columns %in% names(r)))
  expect_setequal(r$locus, c("matK", "rbcL", "trnL-trnF"))
  expect_true(all(r$state %in% 1:3))
  m <- r[r$locus == "matK", ]
  expect_identical(m$path, "whole")
  expect_identical(m$orientation, "forward")
  expect_true(m$state %in% 1:2)
  expect_identical(m$predicted_genus, "Opuntia")
  expect_true(all(r$reason[r$locus != "matK"] == "no_overlap"))
  expect_true(all(r$path[r$locus != "matK"] == "none"))
  expect_true(file.exists(file.path(f$out, "TABLE_barcoding_identify_query.csv")))
})

test_that("the reverse complement of a query gets the same answer", {
  .idn_skip()
  f <- .idn_fixture(withr::local_tempdir())
  s <- unname(f$seqs$rbcL[4])
  a <- .idn_run(f, c(q = s), locus = "rbcL")$tab
  b <- .idn_run(f, c(q = .idn_rc(s)), locus = "rbcL")$tab
  expect_identical(b$orientation, "reverse")
  for (col in c("state", "predicted_species", "predicted_genus", "genus_confidence",
                "species_confidence")) {
    expect_identical(b[[col]], a[[col]])
  }
})

test_that("a sequence foreign to every locus is state 3, no_overlap, with the wording of section 3.6", {
  .idn_skip()
  f <- .idn_fixture(withr::local_tempdir())
  set.seed(7L)
  out <- .idn_run(f, c(alien = .idn_random(500)))
  expect_true(all(out$tab$state == 3L))
  expect_true(all(out$tab$reason == "no_overlap"))
  expect_true(all(is.na(out$tab$predicted_genus)))
  expect_true(any(grepl("Not assignable: no overlap with the library", out$text, fixed = TRUE)))
})

test_that("a fragment shorter than min_overlap is state 3, short_overlap", {
  .idn_skip()
  f <- .idn_fixture(withr::local_tempdir())
  frag <- substr(unname(f$seqs$matK[7]), 101, 160)
  r <- .idn_run(f, c(frag = frag), locus = "matK")$tab
  expect_equal(nrow(r), 1L)
  expect_equal(r$state, 3L)
  expect_identical(r$reason, "short_overlap")
  expect_equal(r$region_length, 60L)
})

# ---- I3: the crop -----------------------------------------------------------------------------------

test_that("a locus inside 20 kb of flanks is cropped and answered as the bare locus", {
  .idn_skip()
  f <- .idn_fixture(withr::local_tempdir())
  s <- unname(f$seqs$matK[10])
  set.seed(8L)
  contig <- paste0(.idn_random(20000, 0.37), s, .idn_random(20000, 0.37))
  bare <- .idn_run(f, c(q = s), locus = "matK")$tab
  crop <- .idn_run(f, c(q = contig), locus = "matK")$tab
  expect_identical(crop$path, "crop")
  expect_equal(crop$query_length, 40300L)
  expect_gte(crop$region_start, 20001L)
  expect_lte(crop$region_end, 20300L)
  expect_equal(crop$other_windows, 0L)
  for (col in c("state", "predicted_species", "predicted_genus")) {
    expect_identical(crop[[col]], bare[[col]])
  }
})

test_that("a contig with two loci gets one row per locus and no combined row", {
  .idn_skip()
  f <- .idn_fixture(withr::local_tempdir())
  set.seed(9L)
  contig <- paste0(.idn_random(5000, 0.37), unname(f$seqs$matK[13]), .idn_random(5000, 0.37),
                   unname(f$seqs$rbcL[1]), .idn_random(5000, 0.37))
  r <- .idn_run(f, c(contig = contig))$tab
  expect_equal(nrow(r), 3L)
  expect_false(anyDuplicated(r$locus) > 0L)
  expect_identical(r$path[r$locus == "matK"], "crop")
  expect_identical(r$path[r$locus == "rbcL"], "crop")
  expect_identical(r$predicted_genus[r$locus == "matK"], "Mammillaria")
  expect_identical(r$predicted_genus[r$locus == "rbcL"], "Opuntia")
  expect_identical(r$reason[r$locus == "trnL-trnF"], "no_overlap")
})

test_that("the reverse complement of a contig gets the same answers", {
  .idn_skip()
  f <- .idn_fixture(withr::local_tempdir())
  set.seed(10L)
  contig <- paste0(.idn_random(5000, 0.37), unname(f$seqs$matK[13]), .idn_random(5000, 0.37),
                   unname(f$seqs$rbcL[1]), .idn_random(5000, 0.37))
  a <- .idn_run(f, c(contig = contig))$tab
  b <- .idn_run(f, c(contig = .idn_rc(contig)))$tab
  a <- a[order(a$locus), ]
  b <- b[order(b$locus), ]
  expect_identical(b$state, a$state)
  expect_identical(b$predicted_species, a$predicted_species)
  expect_identical(b$predicted_genus, a$predicted_genus)
  expect_identical(b$region_length, a$region_length)
})

# ---- I7: inputs ------------------------------------------------------------------------------------

test_that("a FASTQ file stops with the requirement of Phase 8", {
  f <- .idn_fixture(withr::local_tempdir())
  fq <- file.path(dirname(f$library_dir), "reads.fastq")
  writeLines(c("@read1", "ACGTACGTACGT", "+", "IIIIIIIIIIII"), fq)
  expect_error(identify_barcoding_query(fq, library_dir = f$library_dir, output_dir = f$out),
               "raw reads")
})

test_that("a declared locus gives one row per query, for that locus only", {
  .idn_skip()
  f <- .idn_fixture(withr::local_tempdir())
  q <- Biostrings::DNAStringSet(c(a = unname(f$seqs$rbcL[2]), b = unname(f$seqs$rbcL[5])))
  r <- .idn_run(f, q, locus = "rbcL")$tab
  expect_equal(nrow(r), 2L)
  expect_true(all(r$locus == "rbcL"))
  expect_identical(r$query, c("a", "b"))
})

# ---- I2: the model cache ---------------------------------------------------------------------------

test_that("a second call reuses the model and a changed library retrains it", {
  .idn_skip()
  f <- .idn_fixture(withr::local_tempdir())
  q <- c(q = unname(f$seqs$matK[1]))
  .idn_run(f, q, locus = "matK")
  model <- file.path(f$out, "models", "MODEL_idtaxa_matK.rds")
  expect_true(file.exists(model))
  md5_1 <- unname(tools::md5sum(model))
  time_1 <- file.mtime(model)
  Sys.sleep(1.1)
  .idn_run(f, q, locus = "matK")
  expect_identical(file.mtime(model), time_1)
  lib <- file.path(f$library_dir, "LIB_matK.fasta")
  x <- Biostrings::readDNAStringSet(lib)
  Biostrings::writeXStringSet(x[-length(x)], lib)
  .idn_run(f, q, locus = "matK")
  expect_false(identical(unname(tools::md5sum(model)), md5_1))
})

# ---- I8: equivalence with CN2 ----------------------------------------------------------------------

test_that("the identification reproduces the CN2 table with IdTaxa row by row", {
  .idn_skip()
  tmp <- withr::local_tempdir()
  f <- .idn_fixture(tmp)
  og <- file.path(tmp, "outgroup")
  dir.create(og)
  set.seed(11L)
  s <- .idn_vary(unname(f$seqs$matK[1]), 30)
  x <- stats::setNames(c(s, .idn_rc(s), .idn_random(300)),
                       paste0(c("Portulaca_amilis", "Portulaca_amilis", "Talinum_paniculatum"),
                              "|matK_h", 1:3, ".1"))
  og_fasta <- file.path(og, "matK.fasta")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), og_fasta)
  folds_dir <- file.path(tmp, "5_folds")
  suppressMessages(build_barcoding_folds(library_dir = f$library_dir, output_dir = folds_dir))
  suppressWarnings(suppressMessages(utils::capture.output(
    run_barcoding_controls(library_dir = f$library_dir, folds_dir = folds_dir,
                           output_dir = file.path(tmp, "8_controls"), outgroup_dir = og,
                           method = "idtaxa", controls = "CN2", loci = "matK"))))
  cn2 <- utils::read.csv(file.path(tmp, "8_controls", "TABLE_barcoding_cn2_queries_idtaxa.csv"),
                         stringsAsFactors = FALSE)
  r <- .idn_run(f, og_fasta, locus = "matK")$tab
  expect_identical(r$query, names(x))
  expect_identical(r$state, cn2$state)
  expect_identical(r$predicted_species, cn2$predicted_species)
  expect_identical(r$predicted_genus, cn2$predicted_genus)
  expect_identical(r$genus_idtaxa, cn2$genus_idtaxa)
  expect_identical(r$species_idtaxa, cn2$species_idtaxa)
  expect_equal(r$genus_confidence, cn2$genus_confidence, tolerance = 0)
  expect_equal(r$species_confidence, cn2$species_confidence, tolerance = 0)
  expect_identical(r$orientation[1:2], cn2$orientation[1:2])
  expect_identical(r$reason[3], "no_overlap")
  expect_identical(cn2$reason[3], "no_match")
})

# ---- Reproducibility -------------------------------------------------------------------------------

test_that("two runs are byte-identical and the random state of the session is left alone", {
  .idn_skip()
  f <- .idn_fixture(withr::local_tempdir())
  q <- c(a = unname(f$seqs$matK[4]), b = unname(f$seqs[["trnL-trnF"]][1]))
  set.seed(99L)
  before <- .Random.seed
  .idn_run(f, q, run_name = "one")
  expect_identical(.Random.seed, before)
  .idn_run(f, q, run_name = "two")
  expect_identical(unname(tools::md5sum(file.path(f$out, "TABLE_barcoding_identify_one.csv"))),
                   unname(tools::md5sum(file.path(f$out, "TABLE_barcoding_identify_two.csv"))))
})

# ---- I6: the context of 6D -------------------------------------------------------------------------

test_that("each row carries the 6D wrong-species rates of its locus, or NA with the reason", {
  .idn_skip()
  f <- .idn_fixture(withr::local_tempdir())
  md <- file.path(dirname(f$library_dir), "11_metrics")
  dir.create(md)
  utils::write.csv(data.frame(method = "idtaxa", locus = c("matK", "matK", "rbcL", "rbcL"),
                              scheme = c("species", "genus", "species", "genus"), threshold = 60,
                              wrong_species_rate = c(0.02, 0.04, 0.03, 0.02)),
                   file.path(md, "TABLE_barcoding_metrics_summary_idtaxa.csv"), row.names = FALSE)
  q <- c(q = unname(f$seqs$matK[1]))
  r <- .idn_run(f, q)$tab
  m <- r[r$locus == "matK", ]
  expect_equal(m$validation_ws_rate_species_present, 0.02)
  expect_equal(m$validation_ws_rate_species_absent, 0.04)
  expect_true(is.na(r$validation_ws_rate_species_present[r$locus == "trnL-trnF"]))
  out <- .idn_run(f, q, threshold = 50, run_name = "t50")
  expect_true(all(is.na(out$tab$validation_ws_rate_species_present)))
  expect_true(any(grepl("threshold", out$text)))
})
