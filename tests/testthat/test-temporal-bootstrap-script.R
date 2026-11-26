# AP-04 and AP-22 (closure review of 0.5.4): the SLURM script of the temporal bootstrap wrote no
# --blopt, so the cluster route used the RAxML-NG default (NR-FAST) while the local route used
# nr_safe; it fixed TOTAL_THREADS where the other generators read SLURM_CPUS_PER_TASK, and it was not
# made executable.

test_that("the temporal bootstrap script writes --blopt, takes its threads from SLURM and is executable", {
  dir <- withr::local_tempdir()
  script <- suppressMessages(generate_temporal_bootstrap_script(
    alignment_file = "aln.phy", partition_file = "models.part", best_tree_file = "best.tre",
    base_seed = 1, threads = 8, output_dir = dir, cluster_mail_user = ""
  ))
  x <- readLines(script)
  expect_true(any(grepl("--blopt nr_safe", x, fixed = TRUE)))
  expect_true(any(grepl("TOTAL_THREADS=${SLURM_CPUS_PER_TASK:-8}", x, fixed = TRUE)))
  skip_on_os("windows")
  expect_true(file.access(script, 1L) == 0L)

  script2 <- suppressMessages(generate_temporal_bootstrap_script(
    alignment_file = "aln.phy", partition_file = "models.part", best_tree_file = "best.tre",
    base_seed = 1, output_dir = dir, script_name = "fast.sh", cluster_mail_user = "", blopt = "nr_fast"
  ))
  expect_true(any(grepl("--blopt nr_fast", readLines(script2), fixed = TRUE)))
  expect_identical(eval(formals(calculate_temporal_bootstraps)$blopt),
                   eval(formals(generate_temporal_bootstrap_script)$blopt))
})
