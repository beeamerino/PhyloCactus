test_that("run_joint_realignment() cleans up its temporary per-marker files even when the underlying MAFFT run fails (real function call)", {
  skip_if_not_installed("Biostrings")

  tmp_dir <- withr::local_tempdir()
  input_dir <- file.path(tmp_dir, "raw")
  output_fasta_dir <- file.path(tmp_dir, "fasta_out")
  output_aln_dir <- file.path(tmp_dir, "aln_out")
  dir.create(input_dir, recursive = TRUE)

  seqs <- Biostrings::DNAStringSet(c(sp1 = "ACGTACGTACGT", sp2 = "ACGTACGTACGA"))
  Biostrings::writeXStringSet(seqs, file.path(input_dir, "marker1.fasta"))

  # This test targets the tryCatch/finally cleanup block itself, so it does not require (or
  # skip based on) a real `mafft` installation: whether MAFFT is present and succeeds, or is
  # absent and the run fails with a non-zero exit status, the TEMP_* files staged in
  # output_fasta_dir must never be left behind afterward.
  suppressMessages(try(
    run_joint_realignment(input_dir, output_fasta_dir, output_aln_dir),
    silent = TRUE
  ))

  leftover_temp_files <- list.files(output_fasta_dir, pattern = "^TEMP_", full.names = TRUE)
  expect_length(leftover_temp_files, 0)
})

test_that("clean_taxonomic_names() matches the header format the package itself writes", {
  skip_if_not_installed("Biostrings")

  # The function compared underscore-joined binomials taken from the FASTA header against
  # space-separated names taken from the checklist, so keep_idx was empty for every header this
  # package writes, and the exported file held no sequences. No run called it, so the defect never
  # showed. The two sides are normalised before they are compared.
  tmp <- withr::local_tempdir()

  fasta <- file.path(tmp, "matK.fasta")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(stats::setNames(
      c(paste(rep("A", 120), collapse = ""),
        paste(rep("C", 120), collapse = ""),
        paste(rep("G", 120), collapse = "")),
      c("Opuntia_ficus-indica", "Mammillaria_polyedra|AY015284.1|12", "Nothing_here"))),
    fasta
  )

  checklist <- file.path(tmp, "checklist.csv")
  utils::write.csv(data.frame(pureName = c("Opuntia ficus-indica", "Mammillaria polyedra")),
                   checklist, row.names = FALSE)

  out <- clean_taxonomic_names(raw_input_fasta = fasta, checklist_path = checklist,
                               output_clean_dir = file.path(tmp, "clean"))
  kept <- Biostrings::readDNAStringSet(out)

  expect_equal(length(kept), 2L)
  expect_setequal(names(kept), c("Opuntia ficus-indica", "Mammillaria polyedra"))
  expect_false("Nothing here" %in% names(kept))
})

test_that("clean_taxonomic_names() keeps hyphenated epithets whole", {
  skip_if_not_installed("Biostrings")

  # The separator inside a hyphenated epithet is written as a hyphen in the checklist, and as a
  # hyphen or an underscore in the headers, depending on which stage wrote them. All three spellings
  # have to reach the same accepted name, and the epithet has to survive: reducing the binomial to
  # its first two fragments merged every "Opuntia ficus-*" into one species.
  tmp <- withr::local_tempdir()

  fasta <- file.path(tmp, "rbcL.fasta")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(stats::setNames(
      rep(paste(rep("A", 120), collapse = ""), 3L),
      c("Opuntia_ficus-indica", "Opuntia_ficus_indica", "Opuntia_ficus-barbara"))),
    fasta
  )

  checklist <- file.path(tmp, "checklist.csv")
  utils::write.csv(data.frame(pureName = c("Opuntia ficus-indica", "Opuntia ficus-barbara")),
                   checklist, row.names = FALSE)

  out <- clean_taxonomic_names(raw_input_fasta = fasta, checklist_path = checklist,
                               output_clean_dir = file.path(tmp, "clean"))
  kept <- Biostrings::readDNAStringSet(out)

  # The two spellings of ficus-indica are one species, deduplicated to one record; ficus-barbara is
  # a different one and survives as itself.
  expect_setequal(names(kept), c("Opuntia ficus-indica", "Opuntia ficus-barbara"))
  expect_equal(length(kept), 2L)
})
