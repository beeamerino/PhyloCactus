# The manifest is the declared inventory of everything the package distributes in inst/extdata.
# These tests keep it honest in both directions: nothing ships undeclared, and nothing the code
# asks for is missing. That closes the loop the 2026-08-28 audit found open, where a distributed
# tree had drifted from what the documentation described and nothing compared the two.

manifest_path <- function() {
  path <- system.file("extdata", "MANIFEST.csv", package = "PhyloCactus")
  skip_if(!nzchar(path), "MANIFEST.csv not installed")
  path
}

extdata_dir <- function() {
  path <- system.file("extdata", package = "PhyloCactus")
  skip_if(!nzchar(path), "extdata not installed")
  path
}

test_that("the manifest declares the columns the publication script depends on", {
  m <- utils::read.csv(manifest_path(), stringsAsFactors = FALSE)

  expect_true(all(c("canonical", "legacy_name", "role", "produced_by", "source_path",
                    "validator", "loaded_by_code", "description") %in% names(m)))

  expect_true(all(m$role %in% c("input", "run_output", "download_only", "manifest")))
  expect_true(all(m$validator %in% c("none", "table", "tree_rooted_clade", "tree_unrooted_clade")))

  # A duplicated canonical name would make the publication step ambiguous.
  expect_equal(anyDuplicated(m$canonical), 0L)

  # Every row needs a description: the manifest doubles as the inventory documentation.
  expect_true(all(nzchar(m$description)))
})

test_that("every run output declares where it comes from and how it is validated", {
  m <- utils::read.csv(manifest_path(), stringsAsFactors = FALSE)
  outputs <- m[m$role %in% c("run_output", "download_only"), ]

  expect_gt(nrow(outputs), 0L)
  # Without a source path the publication script cannot find the file in the run directory.
  expect_true(all(nzchar(outputs$source_path)))
  expect_true(all(nzchar(outputs$produced_by)))

  # Curated inputs are never copied from a run, so they must not declare a source path.
  inputs <- m[m$role == "input", ]
  expect_true(all(!nzchar(inputs$source_path)))
})

test_that("nothing is distributed in extdata without being declared", {
  m <- utils::read.csv(manifest_path(), stringsAsFactors = FALSE)
  present <- setdiff(list.files(extdata_dir()), ".DS_Store")

  # legacy_name is accepted so the test passes both before and after the rename step.
  accepted <- c(m$canonical, m$legacy_name[nzchar(m$legacy_name)])
  expect_setequal(setdiff(present, accepted), character(0))
})

test_that("every file the code loads by name is declared in the manifest", {
  m <- utils::read.csv(manifest_path(), stringsAsFactors = FALSE)

  # Read from the installed package: this is what a user actually gets.
  code_dirs <- c(system.file("scripts", package = "PhyloCactus"),
                 system.file("doc", package = "PhyloCactus"))
  code_files <- unlist(lapply(code_dirs[nzchar(code_dirs)], list.files,
                              pattern = "\\.(R|Rmd)$", full.names = TRUE))
  skip_if(length(code_files) == 0, "no installed scripts or vignettes to scan")

  cited <- unlist(lapply(code_files, function(f) {
    lines <- readLines(f, warn = FALSE)
    hits <- regmatches(lines, gregexpr('system\\.file\\("extdata", *"[^"]+"', lines))
    sub('.*"extdata", *"', "", unlist(hits))
  }))
  cited <- unique(sub('"$', "", cited))
  skip_if(length(cited) == 0, "no extdata references found in the installed sources")

  accepted <- c(m$canonical, m$legacy_name[nzchar(m$legacy_name)])
  expect_setequal(setdiff(cited, accepted), character(0))
})

test_that("every declared file the code loads is actually present", {
  m <- utils::read.csv(manifest_path(), stringsAsFactors = FALSE)
  present <- list.files(extdata_dir())

  # PROVENANCE.csv only exists after the first publication, so it is exempt.
  required <- m[m$loaded_by_code == TRUE & m$canonical != "PROVENANCE.csv", ]

  missing <- vapply(seq_len(nrow(required)), function(i) {
    row <- required[i, ]
    !(row$canonical %in% present || (nzchar(row$legacy_name) && row$legacy_name %in% present))
  }, logical(1))

  expect_setequal(required$canonical[missing], character(0))
})
