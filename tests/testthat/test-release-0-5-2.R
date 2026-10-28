# Version 0.5.2: the supervisor of step 7 in parallel stops when a worker ends without its status
# (tests in test-barcoding-workers.R). The package carries the version and NEWS opens with it.
# Written before the code.

test_that("NEWS keeps the section of 0.5.2 (the version of the release is tested in test-release-0-5-3.R)", {
  f <- testthat::test_path("..", "..", "NEWS.md")
  skip_if(!file.exists(f), "NEWS.md is not in this tree")
  news <- readLines(f, warn = FALSE)
  expect_length(grep("^# PhyloCactus 0\\.5\\.2$", news), 1L)
  expect_length(grep("^# PhyloCactus 0\\.5\\.1$", news), 1L)
})
