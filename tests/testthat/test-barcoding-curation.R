# Tests of the molecular diagnostic branch (11_barcoding/), curation and screening steps.
# Written before the functions they test.
# The screening reads its inputs from the directories of steps 1 and 2 (BMM, 2026-09-22).

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

test_that("screening reads the registry from the assembly directory and runs without the assembly object", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")

  tmp <- withr::local_tempdir()
  headers <- c("Opuntia_robusta|A1.1", "Opuntia_robusta|A2.1", "Opuntia_stricta|B1.1", "Opuntia_stricta|B2.1")
  sp <- sub("\\|.*$", "", headers)
  asm_dir <- file.path(tmp, "1_assembly")
  dir.create(asm_dir)
  utils::write.csv(data.frame(sid = sub("^.*\\|", "", headers), species = sp, genus = "Opuntia", locus = "matK",
                              cluster_id = 1L, nombre_genbank = sp, stringsAsFactors = FALSE),
                   file.path(asm_dir, "TABLE_barcoding_accession_registry.csv"), row.names = FALSE)

  set.seed(3L)
  base <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  mutate <- function(s, pos) { substr(s, pos, pos) <- ifelse(substr(s, pos, pos) == "A", "C", "A"); s }
  aln <- stats::setNames(c(base, mutate(base, 10), mutate(base, 100), mutate(mutate(base, 100), 200)), headers)
  cur_dir <- file.path(tmp, "2_curated")
  dir.create(file.path(cur_dir, "alignments"), recursive = TRUE)
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(aln), file.path(cur_dir, "alignments", "ALN_masked_final_matK.fasta"))

  suppressMessages(screen_barcoding_markers(assembly_dir = asm_dir, curated_dir = cur_dir, loci = "matK",
                                            output_dir = file.path(tmp, "3_screening"), min_species_with_replica = 2L))

  tab <- utils::read.csv(file.path(tmp, "3_screening", "TABLE_barcoding_marker_screening.csv"), stringsAsFactors = FALSE)
  expect_equal(tab$locus, "matK")
  expect_equal(tab$especies_con_replica, 2L)
  expect_equal(tab$entra, TRUE)
  # Paralogy surveillance (Phase 2, item 5): flagged in every screening table, never excluding
  expect_equal(tab$posible_paralogo, FALSE)

  suppressMessages(screen_barcoding_markers(assembly_dir = asm_dir, curated_dir = cur_dir, loci = "matK",
                                            output_dir = file.path(tmp, "3_screening_flag"),
                                            min_species_with_replica = 2L, paralog_loci = "matK"))
  for (f in c("TABLE_barcoding_marker_screening.csv", "TABLE_barcoding_identical_sequences.csv")) {
    t <- utils::read.csv(file.path(tmp, "3_screening_flag", f), stringsAsFactors = FALSE)
    expect_true(all(t$posible_paralogo), info = f)
  }
  t <- utils::read.csv(file.path(tmp, "3_screening_flag", "TABLE_barcoding_marker_screening.csv"), stringsAsFactors = FALSE)
  expect_equal(t$entra, TRUE)
})

test_that("screening says which earlier step is missing when the registry is absent", {
  tmp <- withr::local_tempdir()
  expect_error(
    screen_barcoding_markers(assembly_dir = file.path(tmp, "1_assembly"), curated_dir = file.path(tmp, "2_curated"),
                             loci = "matK", output_dir = file.path(tmp, "3_screening")),
    "assemble_barcoding_dataset"
  )
})
