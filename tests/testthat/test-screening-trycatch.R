test_that("a short alignment is skipped without aborting run_marker_screening (real function call)", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")

  # The skip payload is built inside a tryCatch. If it were wrapped in return() it would
  # exit the enclosing function instead of being captured as the value of the tryCatch,
  # aborting the whole run on the first short marker. This test exercises that path
  # through the real function: a 10-column alignment is far below min_cols_to_evaluate,
  # and a second, usable marker must still be processed afterwards.
  tmp_dir <- withr::local_tempdir()
  fasta_dir <- file.path(tmp_dir, "alignments")
  dir.create(fasta_dir)

  # run_marker_screening() only reads files matching ^ALN_masked_final_.*\\.fasta$,
  # the naming that run_alignment_pipeline() emits, and strips that prefix to derive
  # the marker name.
  short <- Biostrings::DNAStringSet(c(sp1 = "ACGTACGTAC", sp2 = "ACGTACGTAC",
                                      sp3 = "ACGTACGTAT", sp4 = "ACGTACGTAG"))
  Biostrings::writeXStringSet(short, file.path(fasta_dir, "ALN_masked_final_short_marker.fasta"))

  # The usable marker is built by mutating a common ancestor rather than by drawing
  # independent random sequences. Independent draws are saturated by construction, which
  # makes phangorn::dist.ml() return NaN for every pair and floods the run with warnings
  # that have nothing to do with what is being tested. A 12 percent divergence keeps the
  # distances computable, which is also what a real marker looks like.
  set.seed(1014)
  bases <- c("A", "C", "G", "T")
  ancestor <- sample(bases, 300, replace = TRUE)
  long_seqs <- vapply(1:6, function(i) {
    child <- ancestor
    mut_pos <- sample(300, size = round(0.12 * 300))
    child[mut_pos] <- sample(bases, length(mut_pos), replace = TRUE)
    paste(child, collapse = "")
  }, character(1))
  names(long_seqs) <- paste0("sp", 1:6)
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(long_seqs),
                              file.path(fasta_dir, "ALN_masked_final_long_marker.fasta"))

  summary_tbl <- run_marker_screening(
    fasta_folder = fasta_dir,
    out_base = file.path(tmp_dir, "screening_out"),
    min_cols_to_evaluate = 50L,
    min_aln_len_to_retain = 200L,
    min_nseq_to_retain = 2L
  )

  expect_s3_class(summary_tbl, "data.frame")

  # Both markers reached the summary: the short one did not abort the run.
  expect_true(all(c("short_marker", "long_marker") %in% summary_tbl$marker))

  short_row <- summary_tbl[summary_tbl$marker == "short_marker", ]
  expect_equal(short_row$status, "skipped_short_alignment")
  expect_equal(short_row$decision, "NO")
  expect_match(short_row$decision_reason, "alignment_shorter_than_50_columns")
})

test_that("run_marker_screening fails on a directory with no FASTA files (real function call)", {
  tmp_dir <- withr::local_tempdir()
  empty_dir <- file.path(tmp_dir, "no_fastas")
  dir.create(empty_dir)

  # An empty input directory is a configuration error, not something to process silently.
  expect_error(
    run_marker_screening(fasta_folder = empty_dir, out_base = file.path(tmp_dir, "out"))
  )
})
