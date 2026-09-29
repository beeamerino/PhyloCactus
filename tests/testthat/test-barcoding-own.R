# Tests of the own data of step 1 (decisions O1 to O3 of BMM, 29-09). Written before the code.
# Own data (Sanger, loci taken from own plastomes, anything of the laboratory) live in a folder of
# their own, `0_own/`, apart from the mined data: a sheet `SAMPLES.csv` (sample, species, voucher,
# data_type, published, notes) and one FASTA per library locus, named as the locus, with headers
# `>sample`. Step 1 reads it through `own_dir`; the records enter the registry with
# `source = "own"` and the sid `own:<sample>:<locus>`, and what cannot enter is listed.

.own_checklist <- function(dir) {
  f <- file.path(dir, "checklist.csv")
  utils::write.csv(data.frame(pureName = c("Eriosyce aurata", "Copiapoa cinerea", "Opuntia ficus-indica")),
                   f, row.names = FALSE)
  f
}

.own_dir <- function() {
  d <- tempfile("own_")
  dir.create(d)
  utils::write.csv(data.frame(
    sample = c("PG101", "PG102", "ST7", "BV3"),
    species = c("Eriosyce aurata", "Copiapoa cinerea", "Eriosyce informalis", "Opuntia ficus-indica"),
    voucher = c("PG 101 (SGO)", "PG 102 (SGO)", "ST 7 (CONC)", ""),
    data_type = c("plastome", "plastome", "sanger", "sanger"),
    published = c(FALSE, FALSE, FALSE, TRUE),
    notes = "", stringsAsFactors = FALSE), file.path(d, "SAMPLES.csv"), row.names = FALSE)
  writeLines(c(">PG101", "ACGTACGTAA", ">PG102", "ACGTACGTTT", ">ST7", "ACGTACGTCC", ">XX9", "ACGTACGTGG"),
             file.path(d, "matK.fasta"))
  writeLines(c(">PG101 rbcL from the supermatrix", "TTGGCCAA", ">BV3", "TTGGCCTT"), file.path(d, "rbcL.fasta"))
  writeLines(c(">PG101", "GGGGCCCC"), file.path(d, "accD.fasta"))
  d
}

test_that("own records are read from the sheet and the FASTA of each locus (O1)", {
  d <- .own_dir()
  own <- .bc_own_records(d, checklist_path = .own_checklist(d), loci = c("matK", "rbcL", "ITS"))
  r <- own$registry
  expect_true(all(c("sid", "species", "genus", "locus", "cluster_id", "genbank_name", "source") %in% names(r)))
  expect_setequal(r$sid, c("own:PG101:matK", "own:PG102:matK", "own:PG101:rbcL", "own:BV3:rbcL"))
  expect_true(all(r$source == "own"))
  expect_true(all(is.na(r$cluster_id)))
  expect_identical(r$species[r$sid == "own:PG101:matK"], "Eriosyce_aurata")
  expect_identical(r$genus[r$sid == "own:BV3:rbcL"], "Opuntia")
  expect_identical(r$genbank_name[r$sid == "own:PG102:matK"], "PG102")
  # The header is the sample; text after the first space is ignored
  expect_identical(unname(own$sequences[["own:PG101:rbcL"]]), "TTGGCCAA")
  expect_identical(names(own$sequences), r$sid)
})

test_that("what cannot enter is listed with its reason, not dropped silently (O1)", {
  d <- .own_dir()
  own <- .bc_own_records(d, checklist_path = .own_checklist(d), loci = c("matK", "rbcL", "ITS"))
  lo <- own$left_out
  expect_true(all(c("sample", "locus", "reason") %in% names(lo)))
  expect_identical(lo$reason[lo$sample == "ST7" & lo$locus == "matK"], "species_not_in_checklist")
  expect_identical(lo$reason[lo$sample == "XX9" & lo$locus == "matK"], "sample_not_in_sheet")
  expect_identical(lo$reason[lo$sample == "PG101" & lo$locus == "accD"], "locus_not_in_library")
  expect_identical(nrow(lo), 3L)
})

test_that("a malformed folder stops with a message (O1)", {
  d <- .own_dir()
  ck <- .own_checklist(d)
  expect_error(.bc_own_records(file.path(d, "none"), ck, "matK"), "SAMPLES.csv")
  s <- utils::read.csv(file.path(d, "SAMPLES.csv"), stringsAsFactors = FALSE)
  utils::write.csv(s[, c("sample", "voucher")], file.path(d, "SAMPLES.csv"), row.names = FALSE)
  expect_error(.bc_own_records(d, ck, "matK"), "species")
  utils::write.csv(rbind(s, s[1, ]), file.path(d, "SAMPLES.csv"), row.names = FALSE)
  expect_error(.bc_own_records(d, ck, "matK"), "PG101")
  utils::write.csv(s, file.path(d, "SAMPLES.csv"), row.names = FALSE)
  writeLines(c(">PG101", "ACGT", ">PG101", "ACGA"), file.path(d, "matK.fasta"))
  expect_error(.bc_own_records(d, ck, "matK"), "PG101")
})

test_that("the extras keep the source of own records already in the registry (O1)", {
  reg <- data.frame(sid = c("AB1.1", "own:PG101:matK"), species = "Eriosyce_aurata", genus = "Eriosyce",
                    locus = "matK", cluster_id = c(1L, NA), genbank_name = c("Eriosyce_aurata", "PG101"),
                    source = c("phylotaR", "own"), stringsAsFactors = FALSE)
  d <- tempfile("extra_"); dir.create(d)
  utils::write.csv(data.frame(sid = "NC1.1__matK", species = "Copiapoa_cinerea", locus = "matK",
                              organism = "Copiapoa cinerea"), file.path(d, "TABLE_genomic_extra_records.csv"), row.names = FALSE)
  writeLines(c(">Copiapoa_cinerea|NC1.1__matK", "ACGT"), file.path(d, "extra_records.fasta"))
  out <- .bc_add_extra_records(reg, c(AB1.1 = "ACGT", "own:PG101:matK" = "ACGA"), d)
  expect_identical(out$registry$source[out$registry$sid == "own:PG101:matK"], "own")
  expect_identical(out$registry$source[out$registry$sid == "AB1.1"], "phylotaR")
  expect_identical(out$registry$source[out$registry$sid == "NC1.1__matK"], "genbank_plastome")
})

test_that("step 1 takes own_dir, adds the own records before the cut of long records and lists what is left out (O1)", {
  expect_true("own_dir" %in% names(formals(assemble_barcoding_dataset)))
  expect_null(formals(assemble_barcoding_dataset)$own_dir)
  b <- paste(deparse(body(assemble_barcoding_dataset)), collapse = "\n")
  expect_match(b, ".bc_own_records(", fixed = TRUE)
  expect_match(b, "TABLE_barcoding_own_left_out.csv", fixed = TRUE)
  expect_lt(regexpr(".bc_own_records(", b, fixed = TRUE), regexpr(".bc_cut_long_records(", b, fixed = TRUE))
})

test_that("step 1 of Tutorial 5 shows own_dir (O1)", {
  f <- system.file("scripts", "tutorial-5-cactus-phylogeny-barcoding.R", package = "PhyloCactus")
  skip_if(!nzchar(f), "Tutorial 5 not installed")
  calls <- as.list(parse(f))
  step1 <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("assemble_barcoding_dataset")) &&
                    identical(e$wd_path, "0_phylotaR_raw_Ingroup"), calls)
  expect_length(step1, 1L)
  expect_true("own_dir" %in% names(step1[[1]]))
})
