# Decision P10a of BMM (13-10, Phase 12): the cut of records over 2000 bases (L2 of 28-09) applies to
# Module 1 of the phylogeny as it does to step 1 of the barcoding branch, so that both branches use
# the same records and a locus keeps the extent its shorter records cover. Folder 0 stays raw: the
# cut is made where the clusters of one marker are merged into its FASTA file. Written before the code.

.lr_random_dna <- function(n, seed) {
  withr::with_seed(seed, paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = ""))
}

.lr_write <- function(seqs, file) {
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(seqs), file)
  file
}

.lr_fixture <- function(dir) {
  core <- .lr_random_dna(800, 1)
  short <- vapply(1:10, function(i) substr(core, 1 + i * 3, 760 + i * 3), "")
  names(short) <- sprintf("Genus_species%02d", 1:10)
  long <- paste0(.lr_random_dna(700, 2), core, .lr_random_dna(700, 3))
  stray <- .lr_random_dna(2500, 4)
  c1 <- c(short[1:6], Genus_longcut = long)
  c2 <- c(short[7:10], Genus_stray = stray, Genus_stray = short[1])
  list(core = core, long = long,
       files = c(.lr_write(c1, file.path(dir, "CLUSTER_1.fasta")), .lr_write(c2, file.path(dir, "CLUSTER_2.fasta"))))
}

test_that("a record over 2000 bases is cut to the extent of the shorter records of its marker (P10a)", {
  skip_if_not_installed("Biostrings")
  dir <- withr::local_tempdir()
  f <- .lr_fixture(dir)
  out <- file.path(dir, "matK.fasta")
  tab <- PhyloCactus:::.cp_merge_marker_fasta(f$files, out)
  seqs <- as.character(Biostrings::readDNAStringSet(out))
  expect_true("Genus_longcut" %in% names(seqs))
  cut <- seqs[["Genus_longcut"]]
  expect_lte(nchar(cut), 2000L)
  expect_lt(abs(nchar(cut) - 800L), 60L)
  expect_true(grepl(substr(f$core, 100, 600), cut, fixed = TRUE))
  row <- tab[tab$header == "Genus_longcut", ]
  expect_equal(row$action, "cut")
  expect_equal(row$length, nchar(f$long))
  expect_equal(row$region_length, nchar(cut))
})

test_that("a long record sharing nothing with the shorter records is left out, and its species keeps its next record (P10a)", {
  skip_if_not_installed("Biostrings")
  dir <- withr::local_tempdir()
  f <- .lr_fixture(dir)
  out <- file.path(dir, "matK.fasta")
  tab <- PhyloCactus:::.cp_merge_marker_fasta(f$files, out)
  seqs <- as.character(Biostrings::readDNAStringSet(out))
  expect_equal(tab$action[tab$header == "Genus_stray"], "no_overlap")
  expect_equal(sum(names(seqs) == "Genus_stray"), 1L)
  expect_lte(nchar(seqs[["Genus_stray"]]), 2000L)
})

test_that("records of 2000 bases or fewer are untouched and one record per species is written, in file order (P10a)", {
  skip_if_not_installed("Biostrings")
  dir <- withr::local_tempdir()
  f <- .lr_fixture(dir)
  out <- file.path(dir, "matK.fasta")
  PhyloCactus:::.cp_merge_marker_fasta(f$files, out)
  seqs <- as.character(Biostrings::readDNAStringSet(out))
  ins <- unlist(lapply(f$files, function(x) as.character(Biostrings::readDNAStringSet(x))))
  expect_false(anyDuplicated(names(seqs)) > 0)
  for (nm in sprintf("Genus_species%02d", 1:10)) expect_identical(seqs[[nm]], ins[[nm]])
})

test_that("a marker without records over 2000 bases is written as before and gives an empty table (P10a)", {
  skip_if_not_installed("Biostrings")
  dir <- withr::local_tempdir()
  short <- c(Genus_a = .lr_random_dna(500, 5), Genus_b = .lr_random_dna(600, 6), Genus_a = .lr_random_dna(550, 7))
  file <- .lr_write(short, file.path(dir, "CLUSTER_3.fasta"))
  out <- file.path(dir, "rbcL.fasta")
  tab <- PhyloCactus:::.cp_merge_marker_fasta(file, out)
  seqs <- as.character(Biostrings::readDNAStringSet(out))
  expect_equal(nrow(tab), 0L)
  expect_identical(unname(seqs), unname(short[1:2]))
})

test_that("both phylogeny assemblies merge their marker FASTA files through the cut and write its table (P10a)", {
  for (fun in c("assemble_ingroup_phylotar", "assemble_outgroup_phylotar")) {
    body_txt <- paste(deparse(body(get(fun, envir = asNamespace("PhyloCactus")))), collapse = "\n")
    expect_match(body_txt, ".cp_merge_marker_fasta(", fixed = TRUE)
    expect_match(body_txt, "TABLE_LONG_RECORDS_CUT_", fixed = TRUE)
  }
})
