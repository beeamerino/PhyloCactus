test_that("start tree specifications are counted as RAxML-NG counts them", {
  expect_equal(PhyloCactus:::.parse_n_start_trees("rand{25},pars{25}"), 50L)
  expect_equal(PhyloCactus:::.parse_n_start_trees("pars{10}"), 10L)
  # A term without an explicit count contributes exactly one tree, so a mixed specification must
  # not silently drop it.
  expect_equal(PhyloCactus:::.parse_n_start_trees("rand{25},pars"), 26L)
  expect_equal(PhyloCactus:::.parse_n_start_trees("rand"), 1L)
  expect_equal(PhyloCactus:::.parse_n_start_trees(50), 50L)
})

test_that("the worker plan never leaves a worker idle in the final round", {
  # A worker count that does not divide the tree count buys nothing: it reaches the same round
  # count while holding more cores. The plan must therefore pick a divisor.
  for (threads in c(8L, 16L, 64L, 120L, 128L)) {
    plan <- PhyloCactus:::.plan_ml_workers(n_trees = 50L, threads = threads,
                                           min_threads_per_worker = 4L)
    expect_equal(50L %% plan$workers, 0L)
    expect_equal(plan$rounds, 50L / plan$workers)
    expect_gte(plan$threads_per_worker, 4L)
    expect_lte(plan$threads_used, threads)
  }
})

test_that("the worker plan reproduces the Leftraru and workstation allocations", {
  # 120 cores, 50 trees, floor of 4 threads per worker: 25 workers is the largest divisor of 50
  # that keeps the floor, giving two full rounds on 100 of the 120 cores.
  epu <- PhyloCactus:::.plan_ml_workers(50L, 120L, 4L)
  expect_equal(epu$workers, 25L)
  expect_equal(epu$threads_per_worker, 4L)
  expect_equal(epu$rounds, 2L)

  # The 8-thread workstation can only afford two workers at the same floor, which still halves the
  # 50 sequential rounds that n_workers = 1 would produce.
  mac <- PhyloCactus:::.plan_ml_workers(50L, 8L, 4L)
  expect_equal(mac$workers, 2L)
  expect_equal(mac$rounds, 25L)
})

test_that("the worker plan degrades to a single worker rather than violating the thread floor", {
  plan <- PhyloCactus:::.plan_ml_workers(n_trees = 50L, threads = 2L,
                                         min_threads_per_worker = 8L)
  expect_equal(plan$workers, 1L)
  expect_gte(plan$threads_per_worker, 1L)
  expect_equal(plan$rounds, 50L)
})

test_that("generate_ml_search_script writes a submittable single job, not an array", {
  tmp_dir <- withr::local_tempdir()

  script <- generate_ml_search_script(
    alignment_file = file.path(tmp_dir, "aln.phy"),
    partition_file = file.path(tmp_dir, "part.aicc"),
    constraint_file = file.path(tmp_dir, "constraints.tree"),
    outgroup = c("Portulaca_grandiflora", "Portulaca_oleracea"),
    n_init_trees = "rand{25},pars{25}",
    seed = 1111,
    threads = 120,
    output_dir = tmp_dir
  )

  expect_true(file.exists(script))
  lines <- readLines(script, warn = FALSE)

  expect_match(lines[1], "^#!/bin/bash")
  expect_true(any(grepl("^#SBATCH --cpus-per-task=120", lines)))

  # A job array here would produce several best trees with no step to compare their likelihoods.
  expect_false(any(grepl("#SBATCH --array", lines)))

  # The plan must reach the generated command, not just the console message.
  expect_true(any(grepl("NUM_WORKERS=25", lines)))
  expect_true(any(grepl("--workers \\$\\{NUM_WORKERS\\}", lines)))
  expect_true(any(grepl("--tree-constraint", lines)))
  # RAxML-NG takes --outgroup as one comma-separated list. A vector that reached the command as
  # several arguments, or as only its first element, would root the output on the wrong terminals.
  expect_true(any(grepl("--outgroup ['\"]Portulaca_grandiflora,Portulaca_oleracea['\"]", lines)))
  expect_true(any(grepl("--seed \\$\\{SEED\\}", lines)))

  # The module name is not the version, so the job has to report what it actually ran.
  expect_true(any(grepl("--version", lines)))

  plan <- attr(script, "ml_plan")
  expect_equal(plan$workers, 25L)
  expect_equal(plan$rounds, 2L)
})

test_that("generate_ml_search_script honours an explicit worker count", {
  tmp_dir <- withr::local_tempdir()

  script <- generate_ml_search_script(
    alignment_file = file.path(tmp_dir, "aln.phy"),
    partition_file = file.path(tmp_dir, "part.aicc"),
    constraint_file = file.path(tmp_dir, "constraints.tree"),
    n_init_trees = "rand{25},pars{25}",
    seed = 7L,
    threads = 64,
    workers = 10,
    output_dir = tmp_dir
  )

  lines <- readLines(script, warn = FALSE)
  expect_true(any(grepl("NUM_WORKERS=10", lines)))
  expect_equal(attr(script, "ml_plan")$rounds, 5L)

  # No outgroup was declared, so the command must not carry an empty --outgroup flag.
  expect_false(any(grepl("--outgroup", lines)))
})
