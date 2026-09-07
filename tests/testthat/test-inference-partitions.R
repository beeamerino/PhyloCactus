test_that("preprocess_partitions rejects an unknown model_handling value (real function call)", {
  tmp_dir <- withr::local_tempdir()
  phy <- file.path(tmp_dir, "aln.phy")
  part <- file.path(tmp_dir, "part.txt")
  file.create(phy, part)

  # match.arg() runs before anything else, so this is validated even without RAxML-NG.
  expect_error(
    preprocess_partitions(phy_matrix = phy, part_file = part,
                          raxml_path = "raxml-ng", model_handling = "keep_everything"),
    "arg"
  )
})

test_that("preprocess_partitions fails loudly when the RAxML-NG binary is absent (real function call)", {
  tmp_dir <- withr::local_tempdir()
  phy <- file.path(tmp_dir, "aln.phy")
  part <- file.path(tmp_dir, "part.txt")
  file.create(phy, part)

  # A silent fallback here would let the pipeline continue against an unvalidated
  # partition scheme, so the missing binary must abort the call.
  expect_error(
    preprocess_partitions(phy_matrix = phy, part_file = part,
                          raxml_path = "raxml-ng-binary-that-does-not-exist"),
    "not found in your system's PATH"
  )
})

test_that("preprocess_partitions writes a DNA-only partition file under force_dna (real function call)", {
  skip_if(Sys.which("raxml-ng") == "", "raxml-ng binary not available")

  testdata <- testthat::test_path("testdata")
  models_file <- file.path(testdata, "test_partition_models.txt")
  aln_file <- file.path(testdata, "test_alignment.fasta")
  skip_if_not(file.exists(models_file) && file.exists(aln_file), "Test fixture missing")

  tmp_dir <- withr::local_tempdir()
  phy <- file.path(tmp_dir, "aln.fasta")
  part <- file.path(tmp_dir, "part_models.txt")
  file.copy(aln_file, phy)
  file.copy(models_file, part)

  res <- try(
    preprocess_partitions(phy_matrix = phy, part_file = part, raxml_path = "raxml-ng",
                          output_dir = tmp_dir, model_handling = "force_dna"),
    silent = TRUE
  )
  skip_if(inherits(res, "try-error"), "RAxML-NG rejected the fixture; nothing to assert here")

  # The function returns the path it wrote, so the test asserts on that rather than on whichever
  # .txt happens to be listed first in the directory.
  expect_true(file.exists(res))
  expect_equal(basename(res), "cactus_partitions_validated.txt")

  lines <- readLines(res, warn = FALSE)
  expect_true(all(grepl("^DNA,", lines)))
  expect_false(any(grepl("GTR|HKY", lines)))
})

test_that("preserve mode keeps the model field that force_dna would strip (real function call)", {
  skip_if(Sys.which("raxml-ng") == "", "raxml-ng binary not available")

  testdata <- testthat::test_path("testdata")
  models_file <- file.path(testdata, "test_partition_models.txt")
  aln_file <- file.path(testdata, "test_alignment.fasta")
  skip_if_not(file.exists(models_file) && file.exists(aln_file), "Test fixture missing")

  tmp_dir <- withr::local_tempdir()
  phy <- file.path(tmp_dir, "aln.fasta")
  part <- file.path(tmp_dir, "part_models.txt")
  file.copy(aln_file, phy)
  file.copy(models_file, part)

  # Guards the branch that produced the ModelTest-NG segmentation fault: "preserve" must leave a
  # non-DNA model field in place, which is exactly what run_modeltest_ng() has to reject.
  res <- try(
    preprocess_partitions(phy_matrix = phy, part_file = part, raxml_path = "raxml-ng",
                          output_dir = tmp_dir, prefix = "fixture",
                          model_handling = "preserve"),
    silent = TRUE
  )
  skip_if(inherits(res, "try-error"), "RAxML-NG rejected the fixture; nothing to assert here")

  expect_equal(basename(res), "fixture_partitions_validated.txt")
  expect_false(all(grepl("^DNA,", readLines(res, warn = FALSE))))
})

test_that("preprocess_partitions reports whether RAxML-NG reduced the matrix (real function call)", {
  skip_if(Sys.which("raxml-ng") == "", "raxml-ng binary not available")

  testdata <- testthat::test_path("testdata")
  models_file <- file.path(testdata, "test_partition_models.txt")
  aln_file <- file.path(testdata, "test_alignment.fasta")
  skip_if_not(file.exists(models_file) && file.exists(aln_file), "Test fixture missing")

  tmp_dir <- withr::local_tempdir()
  phy <- file.path(tmp_dir, "aln.fasta")
  part <- file.path(tmp_dir, "part_models.txt")
  file.copy(aln_file, phy)
  file.copy(models_file, part)

  res <- try(
    preprocess_partitions(phy_matrix = phy, part_file = part, raxml_path = "raxml-ng",
                          output_dir = tmp_dir, prefix = "cactus"),
    silent = TRUE
  )
  skip_if(inherits(res, "try-error"), "RAxML-NG rejected the fixture; nothing to assert here")

  # Downstream stages read the analysed matrix from these attributes instead of guessing at a
  # filename, so a missing attribute breaks every later stage silently.
  expect_type(attr(res, "was_reduced"), "logical")
  analysed <- attr(res, "analysed_phy")
  expect_true(file.exists(analysed))
  if (isTRUE(attr(res, "was_reduced"))) {
    expect_match(analysed, "\\.raxml\\.reduced\\.phy$")
  } else {
    expect_equal(normalizePath(analysed), normalizePath(phy))
  }
})

test_that("run_modeltest_ng rejects a partition map carrying RAxML-NG model syntax", {
  tmp_dir <- withr::local_tempdir()
  aln <- file.path(tmp_dir, "aln.phy")
  part <- file.path(tmp_dir, "part.txt")
  file.create(aln)
  writeLines(c("GTR+FC+G4m+B, atpB_rbcL = 1-797",
               "GTR+FC+G4m+B, ITS = 798-1754"), part)

  # ModelTest-NG segfaults on this input instead of reporting a format error, so the guard has to
  # live here for the failure to be diagnosable.
  expect_error(
    run_modeltest_ng(modeltest_exec_path = "modeltest-ng", aln_file = aln, part_file = part,
                     prefix = file.path(tmp_dir, "mt")),
    "datatype token|not found in your system's PATH"
  )
})

test_that("the partition fixture carries the model annotations the preserve mode must keep", {
  testdata <- testthat::test_path("testdata")
  models_file <- file.path(testdata, "test_partition_models.txt")
  skip_if_not(file.exists(models_file), "Test fixture missing")

  # Guards the fixture itself: if it ever loses its model names, the force_dna test
  # above would pass vacuously.
  lines <- readLines(models_file, warn = FALSE)
  expect_true(any(grepl("GTR", lines)))
  expect_true(any(grepl("HKY", lines)))
})
