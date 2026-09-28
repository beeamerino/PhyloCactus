# Tests of the extra records cropped from GenBank plastomes (Phase 8, F2; decisions PL1 to PL7 of
# BMM, 28-09). Written before the code. No network: the download is an argument of the function and
# the tests hand it GenBank flat files written here. No IdTaxa: the extraction does not classify.

.pex_random <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
.pex_vary <- function(s, k) {
  for (p in sample(seq_len(nchar(s)), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1)
  s
}
.pex_rc <- function(s) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(s)))

# Three loci of 300 bases; three genera, two species, three accessions each.
.pex_library <- function(tmp) {
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir, recursive = TRUE)
  set.seed(61L)
  seqs <- list()
  for (l in c("matK", "rbcL", "trnL-trnF")) {
    root <- .pex_random(300)
    x <- character(0)
    for (g in c("Opuntia", "Cereus", "Mammillaria")) {
      gb <- .pex_vary(root, 6)
      for (e in c("alpha", "beta")) {
        sb <- .pex_vary(gb, 3)
        for (r in 1:3) {
          x[paste0(g, "_", e, "|", l, "_", substr(g, 1, 1), substr(e, 1, 1), r, ".1")] <- .pex_vary(sb, 1)
        }
      }
    }
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(lib_dir, paste0("LIB_", l, ".fasta")))
    seqs[[l]] <- x
  }
  list(dir = lib_dir, seqs = seqs)
}

# One GenBank flat file
.pex_gb <- function(accession, organism, sequence, voucher = NA) {
  s <- tolower(sequence)
  starts <- seq(1, nchar(s), by = 60)
  origin <- vapply(starts, function(i) {
    chunk <- substr(s, i, min(i + 59, nchar(s)))
    blocks <- substring(chunk, seq(1, nchar(chunk), 10), pmin(seq(10, nchar(chunk) + 9, 10), nchar(chunk)))
    sprintf("%9d %s", i, paste(blocks, collapse = " "))
  }, "")
  c(sprintf("LOCUS       %s %d bp    DNA     circular PLN 01-JAN-2026", sub("\\..*$", "", accession), nchar(s)),
    sprintf("DEFINITION  %s chloroplast, complete genome.", organism),
    sprintf("VERSION     %s", accession),
    "FEATURES             Location/Qualifiers",
    sprintf("     source          1..%d", nchar(s)),
    sprintf("                     /organism=\"%s\"", organism),
    if (!is.na(voucher)) sprintf("                     /specimen_voucher=\"%s\"", voucher),
    "ORIGIN", origin, "//")
}

.pex_checklist <- c("Opuntia alpha", "Opuntia beta", "Cereus alpha", "Cereus beta", "Mammillaria alpha",
                    "Mammillaria beta", "Opuntia alpha subsp. gamma")

# A plastome of Opuntia alpha holding matK and rbcL (not trnL-trnF) between random flanks
.pex_plastome <- function(f, reverse = FALSE) {
  set.seed(62L)
  s <- paste0(.pex_random(20000), unname(f$seqs$matK[1]), .pex_random(20000), unname(f$seqs$rbcL[2]), .pex_random(20000))
  if (reverse) .pex_rc(s) else s
}

# The fetcher the function receives: GenBank text for a query or for accessions
.pex_fetch <- function(plastome_gb, library_gb = character(0)) {
  function(query = NULL, ids = NULL) {
    if (!is.null(ids)) return(paste(library_gb, collapse = "\n"))
    paste(plastome_gb, collapse = "\n")
  }
}

.pex_run <- function(f, fetch, tmp) {
  out <- file.path(tmp, "0_genomic_raw", "plastomes")
  suppressMessages(invisible(utils::capture.output(
    r <- extract_barcoding_plastome_loci(library_dir = f$dir, output_dir = out,
                                         checklist_path = .pex_checklist, fetch = fetch))))
  list(res = r, out = out)
}

test_that("the loci of a plastome are cropped to the extent of the library, one record per locus", {
  tmp <- withr::local_tempdir()
  f <- .pex_library(tmp)
  gb <- .pex_gb("NC_000001.1", "Opuntia alpha", .pex_plastome(f))
  x <- .pex_run(f, .pex_fetch(gb), tmp)
  rec <- utils::read.csv(file.path(x$out, "TABLE_genomic_extra_records.csv"), stringsAsFactors = FALSE)
  expect_setequal(rec$locus, c("matK", "rbcL"))
  expect_setequal(rec$sid, c("NC_000001.1__matK", "NC_000001.1__rbcL"))
  expect_true(all(rec$species == "Opuntia_alpha"))
  fa <- Biostrings::readDNAStringSet(file.path(x$out, "extra_records.fasta"))
  expect_setequal(names(fa), c("Opuntia_alpha|NC_000001.1__matK", "Opuntia_alpha|NC_000001.1__rbcL"))
  expect_identical(as.character(fa[["Opuntia_alpha|NC_000001.1__matK"]]), toupper(unname(f$seqs$matK[1])))
  expect_identical(as.character(fa[["Opuntia_alpha|NC_000001.1__rbcL"]]), toupper(unname(f$seqs$rbcL[2])))
})

test_that("an infraspecific name is cut to the binomial and an unmatched name is left out and listed", {
  tmp <- withr::local_tempdir()
  f <- .pex_library(tmp)
  gb <- c(.pex_gb("NC_000002.1", "Opuntia alpha subsp. gamma", .pex_plastome(f)),
          .pex_gb("NC_000003.1", "Opuntia inexistens", .pex_plastome(f)))
  x <- .pex_run(f, .pex_fetch(gb), tmp)
  rec <- utils::read.csv(file.path(x$out, "TABLE_genomic_extra_records.csv"), stringsAsFactors = FALSE)
  expect_true(all(rec$species == "Opuntia_alpha"))
  expect_false(any(grepl("NC_000003.1", rec$sid)))
  un <- utils::read.csv(file.path(x$out, "TABLE_unmatched_names.csv"), stringsAsFactors = FALSE)
  expect_true("NC_000003.1" %in% un$accession)
})

test_that("a plastome of a specimen already in a locus of the library does not add that locus", {
  tmp <- withr::local_tempdir()
  f <- .pex_library(tmp)
  gb <- .pex_gb("NC_000004.1", "Opuntia alpha", .pex_plastome(f), voucher = "Majure, L. 5610 (DES)")
  lib_gb <- .pex_gb("matK_Oa1.1", "Opuntia alpha", unname(f$seqs$matK[1]), voucher = "Majure 5610")
  x <- .pex_run(f, .pex_fetch(gb, lib_gb), tmp)
  rec <- utils::read.csv(file.path(x$out, "TABLE_genomic_extra_records.csv"), stringsAsFactors = FALSE)
  expect_identical(rec$locus, "rbcL")
  sk <- utils::read.csv(file.path(x$out, "TABLE_skipped_same_specimen.csv"), stringsAsFactors = FALSE)
  expect_identical(sk$locus, "matK")
  expect_identical(sk$library_sid, "matK_Oa1.1")
})

test_that("a plastome in the reverse strand gives the same records, with the strand recorded", {
  tmp <- withr::local_tempdir()
  f <- .pex_library(tmp)
  gb <- .pex_gb("NC_000005.1", "Opuntia alpha", .pex_plastome(f, reverse = TRUE))
  x <- .pex_run(f, .pex_fetch(gb), tmp)
  rec <- utils::read.csv(file.path(x$out, "TABLE_genomic_extra_records.csv"), stringsAsFactors = FALSE)
  expect_true(all(rec$strand == "reverse"))
  fa <- Biostrings::readDNAStringSet(file.path(x$out, "extra_records.fasta"))
  expect_identical(as.character(fa[["Opuntia_alpha|NC_000005.1__matK"]]), toupper(unname(f$seqs$matK[1])))
})

test_that("the provenance table holds accession, locus, coordinates, strand and md5", {
  tmp <- withr::local_tempdir()
  f <- .pex_library(tmp)
  gb <- .pex_gb("NC_000006.1", "Opuntia alpha", .pex_plastome(f))
  x <- .pex_run(f, .pex_fetch(gb), tmp)
  rec <- utils::read.csv(file.path(x$out, "TABLE_genomic_extra_records.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("sid", "accession", "species", "organism", "locus", "region_start", "region_end",
                    "strand", "region_length", "other_windows", "md5") %in% names(rec)))
  m <- rec[rec$locus == "matK", ]
  expect_equal(m$region_start, 20001L)
  expect_equal(m$region_end, 20300L)
  expect_true(file.exists(file.path(x$out, "MANIFEST.csv")))
  man <- utils::read.csv(file.path(x$out, "MANIFEST.csv"), stringsAsFactors = FALSE)
  expect_identical(man$md5[man$accession == "NC_000006.1"], m$md5)
})

test_that("step 1 takes the extra records only when asked, and marks their source", {
  registry <- data.frame(sid = c("A1.1", "B1.1"), species = c("Opuntia_alpha", "Cereus_beta"),
                         genus = c("Opuntia", "Cereus"), locus = c("matK", "matK"),
                         cluster_id = c(1L, 1L), genbank_name = c("Opuntia_alpha", "Cereus_beta"),
                         stringsAsFactors = FALSE)
  sequences <- c(A1.1 = "ACGT", B1.1 = "TTGG")
  same <- .bc_add_extra_records(registry, sequences, NULL)
  expect_identical(same$registry, registry)
  expect_identical(same$sequences, sequences)

  tmp <- withr::local_tempdir()
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c("Opuntia_beta|NC_1.1__matK" = "GGCC")),
                              file.path(tmp, "extra_records.fasta"))
  utils::write.csv(data.frame(sid = "NC_1.1__matK", accession = "NC_1.1", species = "Opuntia_beta",
                              organism = "Opuntia beta", locus = "matK"),
                   file.path(tmp, "TABLE_genomic_extra_records.csv"), row.names = FALSE)
  both <- .bc_add_extra_records(registry, sequences, tmp)
  expect_equal(nrow(both$registry), 3L)
  expect_identical(both$registry$source[both$registry$sid == "NC_1.1__matK"], "genbank_plastome")
  expect_true(all(both$registry$source[both$registry$sid != "NC_1.1__matK"] == "phylotaR"))
  expect_identical(unname(both$sequences["NC_1.1__matK"]), "GGCC")
  expect_false(anyDuplicated(both$registry$sid) > 0)
})
