# AP-02 (closure review of 0.5.4): automate_treePL() writes one seed into the configuration, and
# run_treePL_cv() ran the same priming configuration n_prime times, so every repeat returned the same
# parameters and prime_rule could not change the selection. From 0.5.5 each repeat has its own seed,
# s to s + n_prime - 1, recorded in prime_<label>.csv.

test_that("each priming repeat runs with its own seed, starting at the seed of the configuration", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  tree <- file.path(dir, "t.tre")
  writeLines("((A:1,B:1):1,(C:1,D:1):1);", tree)
  cfg <- file.path(dir, "base.cfg")
  writeLines(c("numsites = 1000", "mrca = AB A B", "min = AB 1", "max = AB 1", "seed = 100"), cfg)

  # A fake treePL: in a priming configuration it prints the block treePL prints, with opt equal to
  # the seed it was given minus 98, so that each seed gives a different opt; any other run fails.
  fake <- file.path(dir, "treePL")
  writeLines(c(
    "#!/bin/sh",
    "if grep -q '^prime$' \"$1\"; then",
    "  s=$(sed -n 's/^seed = //p' \"$1\")",
    "  echo 'PLACE THE LINES BELOW IN THE CONFIGURATION FILE'",
    "  echo \"opt = $((s - 98))\"; echo 'optad = 1'; echo 'optcvad = 1'",
    "  exit 0",
    "fi",
    "exit 1"), fake)
  Sys.chmod(fake, "0755")

  expect_error(
    run_treePL_cv(cfg, tree, label = "t", n_prime = 4L, treepl_bin = fake, work_dir = dir, quiet = TRUE),
    "cross-validation"
  )
  primes <- utils::read.csv(file.path(dir, "prime_t.csv"))
  expect_equal(primes$seed, 100:103)
  expect_equal(primes$opt, 2:5)
  # The configuration of the last repeat carries its seed and only one seed line.
  last_cfg <- readLines(file.path(dir, "configure_prime_t"))
  expect_equal(grep("^seed", last_cfg, value = TRUE), "seed = 103")
})
