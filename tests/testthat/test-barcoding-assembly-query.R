# Tests of the identification of a genome assembly (Phase 8, T3; decisions A1 to A5 and T3a of
# BMM, 28-09). Written before the code. The loci are located with minimap2 and cut with samtools
# faidx, then identified by step 10: the tests skip when either tool is missing from the PATH, and
# IdTaxa runs only under R CMD check and with the installed package (as in
# test-barcoding-identify.R).

.ias_runs <- function() {
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

.ias_skip_tools <- function() {
  if (!nzchar(Sys.which("minimap2"))) skip("minimap2 is not on the PATH.")
  if (!nzchar(Sys.which("samtools"))) skip("samtools is not on the PATH.")
}

.ias_skip <- function() {
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("Biostrings")
  .ias_skip_tools()
  if (!suppressWarnings(suppressMessages(.ias_runs()))) {
    skip("DECIPHER::IdTaxa() does not run in this loading context: match() over XStringSet fails to dispatch inside IdTaxa. The same call works in a plain R session.")
  }
}

.ias_random <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
.ias_vary <- function(s, k) {
  for (p in sample(seq_len(nchar(s)), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1)
  s
}
.ias_rc <- function(s) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(s)))

# The library of test-barcoding-identify.R: three loci of 300 bases, three genera, two species,
# three accessions each.
.ias_fixture <- function(tmp) {
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(41L)
  genera <- c("Opuntia", "Cereus", "Mammillaria")
  epithets <- c("alpha", "beta")
  seqs <- list()
  for (l in c("matK", "rbcL", "trnL-trnF")) {
    root <- .ias_random(300)
    x <- character(0)
    for (g in seq_along(genera)) {
      gb <- .ias_vary(root, 6)
      for (e in seq_along(epithets)) {
        sb <- .ias_vary(gb, 3)
        for (r in 1:3) {
          x[paste0(genera[g], "_", epithets[e], "|", l, "_", substr(genera[g], 1, 1), e, r, ".1")] <-
            .ias_vary(sb, 1)
        }
      }
    }
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(lib_dir, paste0("LIB_", l, ".fasta")))
    seqs[[l]] <- x
  }
  list(library_dir = lib_dir, seqs = seqs, out = file.path(tmp, "10_identify"), tmp = tmp)
}

# An assembly of Opuntia alpha: matK on scaffold_1 (forward), rbcL on scaffold_2 (reverse), a third
# scaffold without loci; trnL-trnF is absent. With copies = TRUE, rbcL is also on scaffold_4.
.ias_assembly <- function(f, copies = FALSE, gz = FALSE) {
  set.seed(42L)
  sc <- c(scaffold_1 = paste0(.ias_random(8000), unname(f$seqs$matK[1]), .ias_random(8000)),
          scaffold_2 = paste0(.ias_random(5000), .ias_rc(unname(f$seqs$rbcL[1])), .ias_random(6000)),
          scaffold_3 = .ias_random(12000))
  if (copies) sc["scaffold_4"] <- paste0(.ias_random(7000), unname(f$seqs$rbcL[1]), .ias_random(7000))
  path <- file.path(f$tmp, if (gz) "assembly.fna.gz" else "assembly.fna")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(sc), path, compress = gz)
  list(path = path, scaffolds = sc)
}

.ias_run <- function(f, assembly, ...) {
  res <- NULL
  msg <- character(0)
  txt <- suppressWarnings(utils::capture.output(withCallingHandlers(
    res <- identify_barcoding_assembly(assembly, library_dir = f$library_dir, metrics_dir = NULL,
                                       output_dir = f$out, ...),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    })))
  list(tab = res, text = c(txt, msg))
}

test_that("each locus is located on its scaffold, cut with its flanks and identified", {
  .ias_skip()
  f <- .ias_fixture(withr::local_tempdir())
  a <- .ias_assembly(f)
  r <- .ias_run(f, a$path)$tab
  expect_setequal(r$locus, c("matK", "rbcL", "trnL-trnF"))
  expect_true(all(c("scaffold", "scaffold_length", "hit_start", "hit_end", "mapq", "matching_bases",
                    "scaffolds_hit", "low_mapq") %in% names(r)))
  m <- r[r$locus == "matK", ]
  expect_identical(m$scaffold, "scaffold_1")
  expect_equal(m$scaffold_length, nchar(a$scaffolds[["scaffold_1"]]))
  expect_equal(m$scaffolds_hit, 1L)
  expect_false(m$low_mapq)
  expect_identical(m$predicted_genus, "Opuntia")
  expect_true(m$state %in% 1:2)
  # The region handed to step 10 is the hit and 500 bases on each side; step 10 cuts it to the core
  expect_equal(m$query_length, 300L + 2L * 500L)
  b <- r[r$locus == "rbcL", ]
  expect_identical(b$scaffold, "scaffold_2")
  expect_identical(b$predicted_genus, "Opuntia")
  expect_true(file.exists(file.path(f$out, "TABLE_barcoding_identify_assembly.csv")))
})

test_that("a locus absent from the assembly is state 3, no_overlap, with no scaffold", {
  .ias_skip()
  f <- .ias_fixture(withr::local_tempdir())
  r <- .ias_run(f, .ias_assembly(f)$path)$tab
  t <- r[r$locus == "trnL-trnF", ]
  expect_equal(nrow(t), 1L)
  expect_equal(t$state, 3L)
  expect_identical(t$reason, "no_overlap")
  expect_true(is.na(t$scaffold))
  expect_equal(t$scaffolds_hit, 0L)
})

test_that("a locus with two copies keeps the best hit and warns of the low mapping quality (T3a)", {
  .ias_skip()
  f <- .ias_fixture(withr::local_tempdir())
  x <- .ias_run(f, .ias_assembly(f, copies = TRUE)$path)
  b <- x$tab[x$tab$locus == "rbcL", ]
  expect_equal(b$scaffolds_hit, 2L)
  expect_true(b$mapq < 20L)
  expect_true(b$low_mapq)
  expect_true(any(grepl("rbcL.*mapping quality", x$text)))
  expect_false(x$tab$low_mapq[x$tab$locus == "matK"])
})

test_that("a gzipped assembly gives the same table as the plain file", {
  .ias_skip()
  f <- .ias_fixture(withr::local_tempdir())
  p <- .ias_run(f, .ias_assembly(f)$path)$tab
  g <- .ias_run(f, .ias_assembly(f, gz = TRUE)$path)$tab
  expect_identical(g, p)
})

test_that("the assembly folder is left as it was: no index is written next to it", {
  .ias_skip()
  f <- .ias_fixture(withr::local_tempdir())
  a <- .ias_assembly(f)
  before <- list.files(dirname(a$path))
  .ias_run(f, a$path)
  expect_setequal(setdiff(list.files(dirname(a$path)), basename(f$out)), before)
})

test_that("a missing external tool stops with its name before anything runs", {
  skip_if_not_installed("Biostrings")
  f <- .ias_fixture(withr::local_tempdir())
  a <- .ias_assembly(f)
  expect_error(identify_barcoding_assembly(a$path, library_dir = f$library_dir, metrics_dir = NULL,
                                           output_dir = f$out, minimap2 = "no_such_minimap2"),
               "no_such_minimap2")
  expect_false(dir.exists(f$out))
})

test_that("a missing assembly file stops with its path", {
  f <- .ias_fixture(withr::local_tempdir())
  expect_error(identify_barcoding_assembly(file.path(f$tmp, "none.fna"), library_dir = f$library_dir,
                                           metrics_dir = NULL, output_dir = f$out),
               "none.fna")
})
