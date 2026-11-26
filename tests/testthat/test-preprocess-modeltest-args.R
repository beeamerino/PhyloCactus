# AP-05 and AP-23 (closure review of 0.5.4): preprocess_partitions() reported as reduced a matrix
# that only an earlier check had reduced, and fixed --threads at 4; run_modeltest_ng() had no
# output_dir and no seed.

test_that("a new check does not report the reduced matrix of an earlier one", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  phy <- file.path(dir, "m.phy"); writeLines(c("2 4", "A ACGT", "B ACGA"), phy)
  part <- file.path(dir, "p.txt"); writeLines("DNA, l1 = 1-4", part)
  args_file <- file.path(dir, "args.txt")
  fake <- file.path(dir, "raxml-ng")
  writeLines(c("#!/bin/sh", paste0("echo \"$@\" > ", shQuote(args_file)), "exit 0"), fake)
  Sys.chmod(fake, "0755")
  out <- file.path(dir, "out"); dir.create(out)
  writeLines("stale", file.path(out, "cactus.raxml.reduced.phy"))
  writeLines("DNA, l1 = 1-2", file.path(out, "cactus.raxml.reduced.partition"))

  res <- suppressMessages(preprocess_partitions(phy, part, fake, output_dir = out, force_check = TRUE,
                                                threads = 2))
  expect_false(attr(res, "was_reduced"))
  expect_identical(attr(res, "analysed_phy"), phy)
  expect_false(file.exists(file.path(out, "cactus.raxml.reduced.phy")))
  expect_match(readLines(args_file), "--threads 2", fixed = TRUE)
})

test_that("run_modeltest_ng() writes into output_dir and passes the seed", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  part <- file.path(dir, "p.txt"); writeLines("DNA, l1 = 1-4", part)
  args_file <- file.path(dir, "args.txt")
  fake <- file.path(dir, "modeltest-ng")
  # Writes <output>.part.aicc, as ModelTest-NG does, at the path given after --output.
  writeLines(c("#!/bin/sh", paste0("echo \"$@\" > ", shQuote(args_file)),
               "while [ $# -gt 0 ]; do if [ \"$1\" = --output ]; then o=$2; fi; shift; done",
               "o=$(echo $o | tr -d \"'\")", "echo DNA > \"$o.part.aicc\"", "exit 0"), fake)
  Sys.chmod(fake, "0755")
  out <- file.path(dir, "mt")
  res <- suppressMessages(run_modeltest_ng(fake, "aln.phy", part, prefix = "x", output_dir = out, seed = 11))
  expect_identical(res, file.path(out, "x.part.aicc"))
  expect_true(file.exists(res))
  expect_match(readLines(args_file), "--rngseed 11", fixed = TRUE)
})
