# AP-08 (closure review of 0.5.4): the row of cactus_constraints.csv in PROVENANCE.csv described an
# earlier version of the file (md5 3c0618f0... against a201603f... distributed). Every row of
# PROVENANCE.csv must describe the file inst/extdata distributes, byte for byte.

test_that("every md5 and size in PROVENANCE.csv is that of the distributed file", {
  prov_path <- system.file("extdata", "PROVENANCE.csv", package = "PhyloCactus")
  skip_if(!nzchar(prov_path), "PROVENANCE.csv not installed")
  prov <- utils::read.csv(prov_path, stringsAsFactors = FALSE, colClasses = "character")
  files <- vapply(prov$canonical, function(f) system.file("extdata", f, package = "PhyloCactus"), "")
  expect_true(all(nzchar(files)), info = paste(prov$canonical[!nzchar(files)], collapse = ", "))
  ok <- nzchar(files)
  md5 <- unname(tools::md5sum(files[ok]))
  mismatch <- prov$canonical[ok][md5 != prov$md5[ok]]
  expect_identical(mismatch, character(0))
  expect_identical(as.numeric(prov$bytes[ok]), as.numeric(file.size(files[ok])))
})
