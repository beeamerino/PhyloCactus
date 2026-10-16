# Decision M1 of BMM (13-10, Phase 12): when one marker fails in Module 2, run_alignment_pipeline()
# reports it in the manifest with its message and goes on. The row of a failed marker had fewer
# columns than the row of a processed marker (no ambiguity columns), so one failure stopped the whole
# module at do.call(rbind, ...) and hid its own message. Written before the code.

test_that("the row of a failed marker has the columns of a processed marker, in the same order (M1)", {
  err <- PhyloCactus:::.alignment_error_row("x", TRUE, "Input FASTA contains zero sequences.")
  ok <- PhyloCactus:::build_marker_consolidated_row("x", data.frame(stage = character(0)), NULL, FALSE, TRUE, "OK")
  expect_identical(names(err), names(ok))
  expect_equal(nrow(err), 1L)
  expect_equal(err$status, "ERROR: Input FASTA contains zero sequences.")
  expect_false(err$final_alignment_written)
  expect_true(err$masking_applied)
})

test_that("an empty marker FASTA is reported in the manifest and does not stop Module 2 (M1)", {
  skip_if_not_installed("Biostrings")
  skip_if(!nzchar(Sys.which("mafft")), "MAFFT not installed")
  dir <- withr::local_tempdir()
  input <- file.path(dir, "in"); dir.create(input)
  withr::with_seed(1, {
    anc <- sample(c("A", "C", "G", "T"), 300, replace = TRUE)
    seqs <- vapply(1:6, function(i) { s <- anc; p <- sample(300, 20); s[p] <- sample(c("A", "C", "G", "T"), 20, TRUE)
                                      paste(s, collapse = "") }, "")
  })
  names(seqs) <- paste0("Genus_sp", 1:6)
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(seqs), file.path(input, "good.fasta"))
  file.create(file.path(input, "empty.fasta"))
  man <- suppressWarnings(suppressMessages(
    run_alignment_pipeline(input_folder = input, output_dir = file.path(dir, "out"), mask_alignment_regions = FALSE)))
  expect_setequal(man$marker, c("good", "empty"))
  expect_match(man$status[man$marker == "empty"], "^ERROR:")
  expect_false(grepl("^ERROR", man$status[man$marker == "good"]))
  expect_true(file.exists(file.path(dir, "out", "tables", "TABLE_alignment_summary_all_markers.csv")))
})
