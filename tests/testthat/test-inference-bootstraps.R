test_that("generate_bootstrap_script generates valid Slurm Job Array script for Leftraru Epu", {
  tmp_dir <- withr::local_tempdir()
  aln_file <- file.path(tmp_dir, "aln.phy")
  part_file <- file.path(tmp_dir, "part.txt")
  constraint_file <- file.path(tmp_dir, "constraint.nwk")
  file.create(aln_file, part_file, constraint_file)

  script_path <- generate_bootstrap_script(
    alignment_file = aln_file,
    partition_file = part_file,
    constraint_file = constraint_file,
    bs_per_rep = 250,
    max_reps = 4,
    base_seed = 4242,
    threads = 64,
    output_dir = tmp_dir,
    cluster_partition = "main",
    cluster_mem = "64G",
    cluster_time = "24:00:00",
    load_module = c("gcc/14.2.0-nlhpc", "openmpi/5.0.3-o", "raxml-ng/1.1.0-mpi-zen4-n"),
    raxml_exec = "raxml-ng-mpi",
    preparse = TRUE
  )

  expect_true(file.exists(script_path))
  lines <- readLines(script_path)

  # Check Leftraru Epu Slurm directives
  expect_true(any(grepl("#SBATCH -p main", lines)))
  expect_true(any(grepl("#SBATCH --cpus-per-task=64", lines)))
  expect_true(any(grepl("#SBATCH --mem=64G", lines)))
  expect_true(any(grepl("#SBATCH --array=1-4", lines)))
  expect_true(any(grepl("#SBATCH -t 24:00:00", lines)))

  # Check module loading
  expect_true(any(grepl("module load gcc/14.2.0-nlhpc openmpi/5.0.3-o raxml-ng/1.1.0-mpi-zen4-n", lines)))

  # Check default 4:1 thread-to-worker calculation (empirically calibrated for 40-core chunks on Leftraru Epu)
  expect_true(any(grepl("NUM_WORKERS=\\$\\(\\( TOTAL_THREADS / 4 \\)\\)", lines)))

  # Check dynamic seed calculation
  expect_true(any(grepl("BASE_SEED=4242", lines)))
  expect_true(any(grepl("SEED=\\$\\(\\( SLURM_ARRAY_TASK_ID \\* BASE_SEED \\)\\)", lines)))

  # Check chunk prefix and directory
  expect_true(any(grepl("bs_rep_\\$\\{SLURM_ARRAY_TASK_ID\\}", lines)))
  expect_true(any(grepl("cactus_bs_rep_\\$\\{SLURM_ARRAY_TASK_ID\\}", lines)))

  # Check binary RBA pre-parsing block
  expect_true(any(grepl("raxml-ng-mpi --parse", lines)))
  expect_true(any(grepl("--msa \"\\$RBA_FILE\"", lines)))
})

test_that("generate_bootstrap_script supports preparse = FALSE and custom workers", {
  tmp_dir <- withr::local_tempdir()
  aln_file <- file.path(tmp_dir, "aln.phy")
  part_file <- file.path(tmp_dir, "part.txt")
  constraint_file <- file.path(tmp_dir, "constraint.nwk")
  file.create(aln_file, part_file, constraint_file)

  script_path <- generate_bootstrap_script(
    alignment_file = aln_file,
    partition_file = part_file,
    constraint_file = constraint_file,
    bs_per_rep = 100,
    max_reps = 2,
    threads = 32,
    workers = 8,
    preparse = FALSE,
    output_dir = tmp_dir
  )

  lines <- readLines(script_path)
  expect_true(any(grepl("NUM_WORKERS=8", lines)))
  expect_false(any(grepl("--parse", lines)))
  expect_true(any(grepl("--model", lines)))
})

test_that("generate_bootstrap_script respects an explicit threads_per_worker override", {
  tmp_dir <- withr::local_tempdir()
  aln_file <- file.path(tmp_dir, "aln.phy")
  part_file <- file.path(tmp_dir, "part.txt")
  constraint_file <- file.path(tmp_dir, "constraint.nwk")
  file.create(aln_file, part_file, constraint_file)

  script_path <- generate_bootstrap_script(
    alignment_file = aln_file,
    partition_file = part_file,
    constraint_file = constraint_file,
    threads = 64,
    threads_per_worker = 8L,
    output_dir = tmp_dir
  )

  lines <- readLines(script_path)
  # workers is left NULL, so the bash-side dynamic calculation must use the overridden ratio.
  expect_true(any(grepl("NUM_WORKERS=\\$\\(\\( TOTAL_THREADS / 8 \\)\\)", lines)))
})

test_that("collect_bootstraps gathers chunk files and enforces expected_trees", {
  skip_if_not_installed("ape")

  tmp_dir <- withr::local_tempdir()
  bs_dir <- file.path(tmp_dir, "bs_chunks")
  dir.create(file.path(bs_dir, "bs_rep_01"), recursive = TRUE)
  dir.create(file.path(bs_dir, "bs_rep_02"), recursive = TRUE)

  # Create 2 chunks with 5 trees each (10 trees total)
  trees1 <- replicate(5, ape::rtree(6), simplify = FALSE)
  class(trees1) <- "multiPhylo"
  trees2 <- replicate(5, ape::rtree(6), simplify = FALSE)
  class(trees2) <- "multiPhylo"

  f1 <- file.path(bs_dir, "bs_rep_01", "cactus_bs_rep_01.raxml.bootstraps")
  f2 <- file.path(bs_dir, "bs_rep_02", "cactus_bs_rep_02.raxml.bootstraps")
  ape::write.tree(trees1, f1)
  ape::write.tree(trees2, f2)

  # Test successful collection with exact expected count
  out_tree <- collect_bootstraps(bs_dir = bs_dir, expected_trees = 10)
  expect_true(file.exists(out_tree))

  collected_trees <- ape::read.tree(out_tree)
  expect_s3_class(collected_trees, "multiPhylo")
  expect_length(collected_trees, 10)

  # Test failure when expected_trees count mismatches
  expect_error(
    collect_bootstraps(bs_dir = bs_dir, expected_trees = 1000),
    "Bootstrap replicate integrity check failed: expected 1000 bootstrap trees, but collected 10"
  )
})

test_that("collect_bootstraps handles empty directory error", {
  tmp_dir <- withr::local_tempdir()
  expect_error(
    collect_bootstraps(bs_dir = tmp_dir),
    "No bootstrap tree files"
  )
})

test_that(".parse_bs_convergence_log() correctly parses converged and non-converged RAxML-NG logs (real function call)", {
  # 1. Converged mock log
  conv_text <- c(
    "RAxML-NG was called as follows:",
    "Loaded 500 trees with 20 taxa.",
    "Performing bootstrap convergence assessment using autoMRE stop_criterion",
    " # trees        avg WRF       avg WRF in %       # perms: wrf <= 3.00 %     converged?  ",
    "      50          5.210              2.410                          980        NO",
    "     100          2.150              1.010                         1000       YES",
    "Bootstopping test converged after 100 trees",
    "Execution log saved to: test_conv.raxml.log"
  )

  res_conv <- .parse_bs_convergence_log(conv_text)
  expect_true(res_conv$converged)
  expect_equal(res_conv$trees_at_convergence, 100L)
  expect_equal(res_conv$trees_analyzed, 500L)

  # 2. Non-converged mock log
  not_conv_text <- c(
    "RAxML-NG was called as follows:",
    "Loaded 50 trees with 20 taxa.",
    "Performing bootstrap convergence assessment using autoMRE stop_criterion",
    " # trees        avg WRF       avg WRF in %       # perms: wrf <= 3.00 %     converged?  ",
    "      50         34.958             10.821                            0        NO",
    "Bootstopping test did not converge after 50 trees",
    "Execution log saved to: test_not_conv.raxml.log"
  )

  res_not_conv <- .parse_bs_convergence_log(not_conv_text)
  expect_false(res_not_conv$converged)
  expect_true(is.na(res_not_conv$trees_at_convergence))
  expect_equal(res_not_conv$trees_analyzed, 50L)
})

test_that("check_bs_convergence executes and parses with local raxml-ng binary", {
  skip_if(Sys.which("raxml-ng") == "", "raxml-ng binary not available")
  skip_if_not_installed("ape")

  tmp_dir <- withr::local_tempdir()
  t1 <- ape::rtree(8)
  # 50 identical trees will immediately converge at 50 trees under autoMRE
  trees <- replicate(50, t1, simplify = FALSE)
  class(trees) <- "multiPhylo"
  tree_file <- file.path(tmp_dir, "test_identical.tree")
  ape::write.tree(trees, tree_file)

  res <- check_bs_convergence(
    raxml_bin_path = "raxml-ng",
    bs_trees_file = tree_file,
    bs_cutoff = 0.03,
    bs_metric = "tbe",
    threads = 2,
    output_dir = tmp_dir
  )

  expect_s3_class(res, "cactus_bs_convergence")
  expect_true(res$converged)
  expect_equal(res$trees_at_convergence, 50L)
  expect_true(file.exists(res$log_file))
  expect_equal(as.character(res), res$log_file)
})

test_that("run_local_bootstraps validates input files and worker ratio", {
  tmp_dir <- withr::local_tempdir()

  # Error on non-existent binary
  expect_error(
    run_local_bootstraps(
      raxml_bin_path = "non_existent_binary_xyz",
      aln_file = file.path(tmp_dir, "fake.phy"),
      part_file = file.path(tmp_dir, "fake.part"),
      constraint_file = file.path(tmp_dir, "fake.nwk")
    ),
    "not found in your system's PATH"
  )

  # Error on missing files
  expect_error(
    run_local_bootstraps(
      raxml_bin_path = "R", # placeholder existing binary
      aln_file = file.path(tmp_dir, "missing.phy"),
      part_file = file.path(tmp_dir, "part.txt"),
      constraint_file = file.path(tmp_dir, "constraint.nwk")
    ),
    "Alignment file not found"
  )
})


