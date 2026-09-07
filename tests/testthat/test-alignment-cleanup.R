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
