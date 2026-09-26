# Tests of the molecular diagnostic branch (11_barcoding/), final library step.
# Written before the function they test (decisions of BMM, 2026-09-21, Phase 2 proposal sec. 8):
# sequences that match their locus in neither direction are excluded and declared; identical aligned
# sequences are collapsed within species and locus; replica means two or more distinct sequences;
# the threshold is evaluated again after the collapse.
# Like the phylogeny, the branch passes work between steps through directories, not R objects
# (BMM, 2026-09-22): the library reads the registry, the curated alignments and the screening
# decision from the folders of steps 1, 2 and 3, so it runs without re-running them.

# Minimal curated directory, laid out as run_alignment_pipeline() writes it.
.lib_fixture <- function(tmp, alns, strand_logs) {
  cur <- file.path(tmp, "2_curated")
  dir.create(file.path(cur, "alignments"), recursive = TRUE)
  dir.create(file.path(cur, "tables"), recursive = TRUE)
  for (l in names(alns)) {
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(alns[[l]]),
                                file.path(cur, "alignments", paste0("ALN_masked_final_", l, ".fasta")))
  }
  for (l in names(strand_logs)) {
    utils::write.csv(strand_logs[[l]], file.path(cur, "tables", paste0("LOG_STRAND_", l, ".csv")), row.names = FALSE)
  }
  cur
}

# Step outputs the final library reads from disk: the registry of step 1 and the screening table of step 3.
.lib_step_files <- function(tmp, registry, screening) {
  asm <- file.path(tmp, "1_assembly")
  scr <- file.path(tmp, "3_screening")
  dir.create(asm)
  dir.create(scr)
  utils::write.csv(registry, file.path(asm, "TABLE_barcoding_accession_registry.csv"), row.names = FALSE)
  utils::write.csv(screening, file.path(scr, "TABLE_barcoding_marker_screening.csv"), row.names = FALSE)
  list(assembly_dir = asm, screening_dir = scr)
}

.lib_registry <- function(headers, locus) {
  sp <- sub("\\|.*$", "", headers)
  data.frame(sid = sub("^.*\\|", "", headers), species = sp, genus = sub("_.*$", "", sp),
             locus = locus, cluster_id = 1L, genbank_name = sp, stringsAsFactors = FALSE)
}

.strand_log <- function(headers, unmatched = character(0)) {
  data.frame(Seq = headers, Length = 9L,
             ShareForward = ifelse(headers %in% unmatched, 0.1, 1),
             ShareReverse = 0,
             Action = ifelse(headers %in% unmatched, "no_match_either_direction", "kept"),
             stringsAsFactors = FALSE)
}

test_that("identical aligned sequences collapse within a species to the alphabetically smallest sid", {
  aln <- c(`Opuntia_robusta|B2.1` = "ACGT-ACGT", `Opuntia_robusta|A9.1` = "ACGT-ACGT",
           `Opuntia_robusta|C1.1` = "ACGTTACGT")

  out <- .bc_collapse_identical(aln, locus = "matK")

  expect_setequal(names(out$aln), c("Opuntia_robusta|A9.1", "Opuntia_robusta|C1.1"))
  m <- out$map
  expect_equal(nrow(m), 3L)
  expect_equal(m$representative[m$sid == "B2.1"], "A9.1")
  expect_equal(m$representative[m$sid == "A9.1"], "A9.1")
  expect_equal(m$representative[m$sid == "C1.1"], "C1.1")
  expect_equal(m$collapsed[m$sid == "B2.1"], TRUE)
  expect_equal(sum(m$collapsed), 1L)
  expect_true(all(m$locus == "matK"))
})

test_that("identical sequences of different species are not collapsed", {
  aln <- c(`Opuntia_robusta|A1.1` = "ACGT-ACGT", `Opuntia_stricta|B1.1` = "ACGT-ACGT")

  out <- .bc_collapse_identical(aln, locus = "matK")

  expect_setequal(names(out$aln), names(aln))
  expect_false(any(out$map$collapsed))
})

test_that("sequence identity ignores letter case", {
  aln <- c(`Opuntia_robusta|A1.1` = "acgt-acgt", `Opuntia_robusta|A2.1` = "ACGT-ACGT")

  out <- .bc_collapse_identical(aln, locus = "matK")

  expect_equal(names(out$aln), "Opuntia_robusta|A1.1")
})

test_that("the final library excludes non-homologous sequences, collapses identical ones and re-evaluates the threshold", {
  skip_if_not_installed("Biostrings")
  tmp <- withr::local_tempdir()

  # matK: two species with two distinct sequences each; one extra identical copy; one non-homologous sequence.
  matk <- c(`Opuntia_robusta|A1.1` = "ACGTACGTAC", `Opuntia_robusta|A2.1` = "ACGTACGTAA",
            `Opuntia_robusta|A3.1` = "ACGTACGTAC",
            `Opuntia_stricta|B1.1` = "TCGTACGTAC", `Opuntia_stricta|B2.1` = "TCGTACGTCC",
            `Opuntia_stricta|B3.1` = "GGGGCCCCTT")
  # ITS: two species whose replication is identical copies only; it passes before the collapse, not after.
  its <- c(`Opuntia_robusta|C1.1` = "ACGTACGTAC", `Opuntia_robusta|C2.1` = "ACGTACGTAC",
           `Opuntia_stricta|D1.1` = "TCGTACGTAC", `Opuntia_stricta|D2.1` = "TCGTACGTAC")
  # A non-homologous sequence the curation had already removed from the alignment.
  its_log_headers <- c(names(its), "Opuntia_stricta|D9.1")

  cur <- .lib_fixture(tmp, list(matK = matk, ITS = its),
                      list(matK = .strand_log(names(matk), "Opuntia_stricta|B3.1"),
                           ITS = .strand_log(its_log_headers, "Opuntia_stricta|D9.1")))
  registry <- rbind(.lib_registry(names(matk), "matK"), .lib_registry(its_log_headers, "ITS"))
  meta <- file.path(tmp, "meta.csv")
  utils::write.csv(data.frame(sid = c("B3.1", "D9.1"), Species_gb = "Opuntia_stricta",
                              Description_gb = c("B3.1 Opuntia stricta desc", "D9.1 Opuntia stricta pseudogene")),
                   meta, row.names = FALSE)

  # rbcL did not enter the screening and has no curated alignment: the library must not look for it.
  dirs <- .lib_step_files(tmp, registry,
                          data.frame(locus = c("matK", "ITS", "rbcL"), enters = c(TRUE, TRUE, FALSE),
                                     stringsAsFactors = FALSE))

  lib <- suppressMessages(finalize_barcoding_library(
    assembly_dir = dirs$assembly_dir, curated_dir = cur, screening_dir = dirs$screening_dir,
    output_dir = file.path(tmp, "4_library"), metadata_file = meta, min_species_with_replica = 2L
  ))

  out_dir <- file.path(tmp, "4_library")

  # Exclusion, declared with locus, species, description and reason
  ex <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_excluded_nonhomologous.csv"), stringsAsFactors = FALSE)
  expect_setequal(ex$sid, c("B3.1", "D9.1"))
  expect_true(all(c("locus", "sid", "species", "description", "reason", "in_curated_alignment") %in% names(ex)))
  expect_true(all(ex$reason == "no_match_either_direction"))
  expect_equal(ex$description[ex$sid == "D9.1"], "D9.1 Opuntia stricta pseudogene")
  expect_equal(ex$in_curated_alignment[ex$sid == "B3.1"], TRUE)
  expect_equal(ex$in_curated_alignment[ex$sid == "D9.1"], FALSE)

  # matK library: excluded sequence out, identical copy collapsed onto A1.1
  f_matk <- file.path(out_dir, "LIB_matK.fasta")
  expect_true(file.exists(f_matk))
  expect_setequal(names(Biostrings::readDNAStringSet(f_matk)),
                  c("Opuntia_robusta|A1.1", "Opuntia_robusta|A2.1", "Opuntia_stricta|B1.1", "Opuntia_stricta|B2.1"))

  col <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_collapsed_identical.csv"), stringsAsFactors = FALSE)
  expect_equal(col$representative[col$locus == "matK" & col$sid == "A3.1"], "A1.1")
  expect_false("B3.1" %in% col$sid)

  # Threshold after the collapse: matK keeps 2 species with replica; ITS drops from 2 to 0 and leaves the library
  s <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_library_summary.csv"), stringsAsFactors = FALSE)
  expect_equal(s$species_with_replicate[s$locus == "matK"], 2L)
  expect_equal(s$species_with_replicate_before_collapse[s$locus == "ITS"], 2L)
  expect_equal(s$species_with_replicate[s$locus == "ITS"], 0L)
  expect_equal(s$accessions_excluded_nonhomologous[s$locus == "matK"], 1L)
  expect_equal(s$accessions_collapsed[s$locus == "matK"], 1L)
  expect_equal(s$accessions_collapsed[s$locus == "ITS"], 2L)
  expect_equal(s$enters[s$locus == "matK"], TRUE)
  expect_equal(s$enters[s$locus == "ITS"], FALSE)
  expect_false(file.exists(file.path(out_dir, "LIB_ITS.fasta")))

  expect_identical(lib$summary$locus, s$locus)
  expect_setequal(s$locus, c("matK", "ITS"))

  # Per-locus funnel (decision of 2026-09-22, sec. 10): every locus assembled or screened, with its count
  # at each stage; a locus that did not enter the screening keeps its row, with NA after that point.
  fn <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_funnel.csv"), stringsAsFactors = FALSE)
  expect_setequal(fn$locus, c("matK", "ITS", "rbcL"))
  expect_true(all(c("provisional_locus", "possible_paralog") %in% names(fn)))
  m <- fn[fn$locus == "matK", ]
  expect_equal(c(m$accessions_assembled, m$accessions_curated, m$accessions_excluded_nonhomologous,
                 m$accessions_collapsed, m$accessions_after_collapse, m$species_with_replicate), c(6L, 6L, 1L, 1L, 4L, 2L))
  expect_equal(c(m$enters_screening, m$enters), c(TRUE, TRUE))
  i <- fn[fn$locus == "ITS", ]
  expect_equal(c(i$accessions_assembled, i$accessions_curated, i$accessions_collapsed, i$accessions_after_collapse),
               c(5L, 4L, 2L, 2L))
  expect_equal(c(i$enters_screening, i$enters), c(TRUE, FALSE))
  r <- fn[fn$locus == "rbcL", ]
  expect_equal(r$accessions_assembled, 0L)
  expect_equal(r$enters_screening, FALSE)
  expect_equal(r$enters, FALSE)
  expect_true(is.na(r$accessions_curated))
  expect_true(is.na(r$accessions_after_collapse))
})

test_that("a locus under paralogy surveillance is flagged in every table of the library, and never excluded for it", {
  skip_if_not_installed("Biostrings")
  tmp <- withr::local_tempdir()

  matk <- c(`Opuntia_robusta|A1.1` = "ACGTACGTAC", `Opuntia_robusta|A2.1` = "ACGTACGTAA",
            `Opuntia_stricta|B1.1` = "TCGTACGTAC", `Opuntia_stricta|B2.1` = "TCGTACGTCC")
  pepc <- c(`Opuntia_robusta|P1.1` = "ACGTACGTAC", `Opuntia_robusta|P2.1` = "ACGTACGTAC",
            `Opuntia_robusta|P3.1` = "ACGTACGTAA",
            `Opuntia_stricta|Q1.1` = "TCGTACGTAC", `Opuntia_stricta|Q2.1` = "TCGTACGTCC",
            `Opuntia_stricta|Q9.1` = "GGGGCCCCTT")
  cur <- .lib_fixture(tmp, list(matK = matk, pepC_like = pepc),
                      list(matK = .strand_log(names(matk)),
                           pepC_like = .strand_log(names(pepc), "Opuntia_stricta|Q9.1")))
  registry <- rbind(.lib_registry(names(matk), "matK"), .lib_registry(names(pepc), "pepC_like"))
  dirs <- .lib_step_files(tmp, registry, data.frame(locus = c("matK", "pepC_like"), enters = TRUE))

  suppressMessages(finalize_barcoding_library(
    assembly_dir = dirs$assembly_dir, curated_dir = cur, screening_dir = dirs$screening_dir,
    output_dir = file.path(tmp, "4_library"), metadata_file = NULL, min_species_with_replica = 2L
  ))

  out_dir <- file.path(tmp, "4_library")
  for (f in c("TABLE_barcoding_library_summary.csv", "TABLE_barcoding_excluded_nonhomologous.csv",
              "TABLE_barcoding_collapsed_identical.csv")) {
    t <- utils::read.csv(file.path(out_dir, f), stringsAsFactors = FALSE)
    expect_true("possible_paralog" %in% names(t), info = f)
    expect_true(all(t$possible_paralog[t$locus == "pepC_like"]), info = f)
    expect_false(any(t$possible_paralog[t$locus == "matK"]), info = f)
  }
  s <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_library_summary.csv"), stringsAsFactors = FALSE)
  expect_equal(s$enters[s$locus == "pepC_like"], TRUE)
  expect_true(file.exists(file.path(out_dir, "LIB_pepC_like.fasta")))
})

test_that("pepC_like is the locus under paralogy surveillance by default", {
  expect_identical(eval(formals(finalize_barcoding_library)$paralog_loci), "pepC_like")
  expect_identical(eval(formals(screen_barcoding_markers)$paralog_loci), "pepC_like")
})

test_that("the paralogy flag marks only the named loci and leaves an empty table valid", {
  flagged <- .bc_flag_paralog(data.frame(locus = c("matK", "pepC_like", "ITS"), stringsAsFactors = FALSE), "pepC_like")
  expect_equal(flagged$possible_paralog, c(FALSE, TRUE, FALSE))
  empty <- .bc_flag_paralog(data.frame(locus = character(0)), "pepC_like")
  expect_equal(nrow(empty), 0L)
  expect_true("possible_paralog" %in% names(empty))
})

test_that("the final library refuses to write inside a directory of the phylogeny", {
  expect_error(
    finalize_barcoding_library(assembly_dir = ".", curated_dir = ".", screening_dir = ".",
                               output_dir = file.path("4_Cleaned", "lib")),
    "phylogeny"
  )
})

test_that("the final library says which earlier step is missing when an input file is absent", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "1_assembly"))
  dir.create(file.path(tmp, "3_screening"))
  utils::write.csv(.lib_registry("Opuntia_robusta|A1.1", "matK"),
                   file.path(tmp, "1_assembly", "TABLE_barcoding_accession_registry.csv"), row.names = FALSE)

  expect_error(
    finalize_barcoding_library(assembly_dir = file.path(tmp, "1_assembly"), curated_dir = file.path(tmp, "2_curated"),
                               screening_dir = file.path(tmp, "3_screening"), output_dir = file.path(tmp, "4_library")),
    "screen_barcoding_markers"
  )
  expect_error(
    finalize_barcoding_library(assembly_dir = file.path(tmp, "none"), curated_dir = file.path(tmp, "2_curated"),
                               screening_dir = file.path(tmp, "3_screening"), output_dir = file.path(tmp, "4_library")),
    "assemble_barcoding_dataset"
  )
})
