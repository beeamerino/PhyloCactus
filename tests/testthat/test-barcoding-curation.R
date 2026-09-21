# Tests of the molecular diagnostic branch (11_barcoding/), curation and screening steps.
# Written before the functions they test.

test_that("Genus_species|sid headers survive the alignment and curation of the phylogeny", {
  skip_if_not_installed("Biostrings")
  skip_if(!nzchar(Sys.which("mafft")), "MAFFT not installed")

  tmp <- withr::local_tempdir()
  input_dir <- file.path(tmp, "1_assembly")
  dir.create(input_dir)
  set.seed(11L)
  base <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  mutate <- function(s, pos) { substr(s, pos, pos) <- ifelse(substr(s, pos, pos) == "A", "C", "A"); s }
  seqs <- c(base, mutate(base, 50), mutate(base, 120), mutate(base, 200))
  names(seqs) <- c("Opuntia_ficus-indica|AB1.1", "Opuntia_ficus-indica|AB2.1",
                   "Opuntia_robusta|CD3.1", "Opuntia_robusta|CD4.1")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(seqs), file.path(input_dir, "matK.fasta"))

  suppressMessages(curate_barcoding_markers(input_dir, file.path(tmp, "2_curated"), loci = "matK"))

  # run_alignment_pipeline() writes the curated alignment as alignments/ALN_masked_final_<locus>.fasta
  aln_file <- file.path(tmp, "2_curated", "alignments", "ALN_masked_final_matK.fasta")
  expect_true(file.exists(aln_file))
  expect_setequal(names(Biostrings::readDNAStringSet(aln_file)), names(seqs))
})

test_that("curation only processes the loci it is given", {
  skip_if_not_installed("Biostrings")
  skip_if(!nzchar(Sys.which("mafft")), "MAFFT not installed")

  tmp <- withr::local_tempdir()
  input_dir <- file.path(tmp, "1_assembly")
  dir.create(input_dir)
  s <- paste(rep("ACGT", 60), collapse = "")
  for (l in c("matK", "NHX")) {
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(`Opuntia_robusta|A1.1` = s, `Opuntia_robusta|A2.1` = s,
                                                          `Opuntia_stricta|B1.1` = s, `Opuntia_stricta|B2.1` = s)),
                                file.path(input_dir, paste0(l, ".fasta")))
  }

  suppressMessages(curate_barcoding_markers(input_dir, file.path(tmp, "2_curated"), loci = "matK"))

  written <- basename(list.files(file.path(tmp, "2_curated"), pattern = "\\.fasta$", recursive = TRUE))
  expect_true(any(grepl("matK", written)))
  expect_false(any(grepl("NHX", written)))
})

test_that("the replication threshold is evaluated on the curated set, not on the assembled one", {
  assembled <- data.frame(
    sid     = c("a", "b", "c", "d"),
    species = c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_stricta", "Opuntia_stricta"),
    genus   = "Opuntia",
    locus   = "matK",
    stringsAsFactors = FALSE
  )
  # Curation removed one accession of each species: before, 2 species with replica; after, none.
  curated_sids <- c("a", "c")

  scr <- .bc_screen_threshold(assembled, curated_sids, min_species_with_replica = 2L)

  expect_equal(scr$especies_con_replica_ensamblado, 2L)
  expect_equal(scr$especies_con_replica_curado, 0L)
  expect_false(scr$pasa_umbral)
})

test_that("saturation is flagged and never excludes a locus", {
  scr <- .bc_screen_decision(
    data.frame(locus = c("matK", "ITS"), pasa_umbral = c(TRUE, TRUE), saturado = c(FALSE, TRUE),
               stringsAsFactors = FALSE)
  )
  expect_true(all(scr$entra))
  expect_equal(scr$marca_saturacion, c(FALSE, TRUE))
})

test_that("identical sequences are counted per locus and species on the curated alignment", {
  aln <- c(`Opuntia_robusta|a` = "ACGT-ACGT", `Opuntia_robusta|b` = "ACGT-ACGT", `Opuntia_robusta|c` = "ACGTTACGT",
           `Opuntia_stricta|d` = "ACGT-ACGA")

  tab <- .bc_identical_sequences(list(matK = aln))

  r <- tab[tab$species == "Opuntia_robusta", ]
  expect_equal(r$accesiones, 3L)
  expect_equal(r$secuencias_distintas, 2L)
  s <- tab[tab$species == "Opuntia_stricta", ]
  expect_equal(s$accesiones, 1L)
  expect_equal(s$secuencias_distintas, 1L)
})
