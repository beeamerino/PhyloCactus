# Version 0.5.2: the supervisor of step 7 in parallel stops when a worker ends without its status
# (tests in test-barcoding-workers.R). The package carries the version and NEWS opens with it.
# Written before the code.

test_that("the package is version 0.5.2 and NEWS opens with it", {
  expect_identical(as.character(utils::packageVersion("PhyloCactus")), "0.5.2")
  f <- testthat::test_path("..", "..", "NEWS.md")
  skip_if(!file.exists(f), "NEWS.md is not in this tree")
  news <- readLines(f, warn = FALSE)
  expect_identical(news[1], "# PhyloCactus 0.5.2")
  expect_length(grep("^# PhyloCactus 0\\.5\\.1$", news), 1L)
})
