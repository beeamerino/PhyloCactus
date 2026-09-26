# Tests of the monophyly per locus of the molecular diagnostic branch (11_barcoding/12_monophyly/).
# Written before the code. Phase 6E of PhyloCactus 0.5.0, decisions N1 to N10 of BMM (26-09).
#
# A group (species or genus) is monophyletic when its tips are one side of a split of the unrooted
# gene tree; it is rejected only against a split of support 70 or more; grades are classified on a
# root taken from the earliest lineage of the reference present in the locus. The trees are written
# by hand, so every status can be checked by hand. No test runs RAxML-NG.

# The gene tree of the fixture, unrooted, with the bootstrap support of each internal branch:
#   A_x      its two tips are a split (95): monophyletic whatever the root
#   A_y      one tip: not evaluable
#   B_z      one tip sits with A_y (90) and the other with B_w (40): rejected
#   B_w      broken only by a branch of 40: not rejected
#   Leu_a    the earliest lineage of the reference, a split (100): the root
#   genus A  broken by A_y|1 + B_z|1 (90): rejected; once rooted, B_z|1 alone is nested in it (a grade)
#   genus B  rejected; nested in it are A_x and A_y, which do not form one lineage (polyphyletic)
.mono_tree <- paste0("((Leu_a|L1,Leu_a|L2)100,((A_x|X1,A_x|X2)95,(A_y|Y1,B_z|Z1)90)80,",
                     "((B_z|Z2,B_w|W1)40,B_w|W2)60);")
# The same unrooted tree written from another node
.mono_tree_other <- paste0("((A_x|X1,A_x|X2)95,(A_y|Y1,B_z|Z1)90,((Leu_a|L1,Leu_a|L2)100,",
                           "((B_z|Z2,B_w|W1)40,B_w|W2)60)80);")
# Leu_a scattered: the tree cannot be rooted on it
.mono_tree_scattered <- paste0("((Leu_a|L1,(Leu_a|L2,A_y|Y1)85)90,((A_x|X1,A_x|X2)95,B_z|Z1)80,",
                               "((B_z|Z2,B_w|W1)40,B_w|W2)60);")
.mono_reference <- "(Out_a,Out_b,(Leu_a,(A_x,(A_y,(B_z,B_w)))));"

.mono_fixture <- function(tmp, tree = .mono_tree, locus = "matK") {
  trees <- file.path(tmp, "trees"); dir.create(trees, recursive = TRUE, showWarnings = FALSE)
  writeLines(tree, file.path(trees, paste0(locus, ".raxml.support")))
  ref <- file.path(tmp, "reference.tree"); writeLines(.mono_reference, ref)
  list(trees_dir = trees, reference_tree = ref, output_dir = file.path(tmp, "12_monophyly"))
}

.mono_run <- function(f, ...) {
  suppressMessages(utils::capture.output(
    res <- assess_barcoding_monophyly(trees_dir = f$trees_dir, output_dir = f$output_dir,
                                      reference_tree = f$reference_tree,
                                      reference_outgroup = c("Out_a", "Out_b"), ...)))
  res
}

.mono_read <- function(f, what) {
  utils::read.csv(file.path(f$output_dir, paste0("TABLE_barcoding_monophyly_", what, ".csv")),
                  stringsAsFactors = FALSE)
}

test_that("a species whose tips are a split is monophyletic, whatever node the tree is written from", {
  for (tr in c(.mono_tree, .mono_tree_other)) {
    tmp <- withr::local_tempdir(); f <- .mono_fixture(tmp, tr)
    .mono_run(f)
    s <- .mono_read(f, "species")
    ax <- s[s$group == "A_x", ]
    expect_identical(ax$status, "monophyletic")
    expect_equal(ax$support, 95)
    expect_identical(s$status[s$group == "Leu_a"], "monophyletic")
  }
})

test_that("non-monophyly is rejected only against a split of support 70 or more", {
  tmp <- withr::local_tempdir(); f <- .mono_fixture(tmp)
  .mono_run(f)
  s <- .mono_read(f, "species")
  expect_identical(s$status[s$group == "B_z"], "rejected")
  expect_identical(s$status[s$group == "B_w"], "not_rejected")
  expect_true(is.na(s$support[s$group == "B_w"]))
  # With the cutoff at 95 the branch of 90 no longer rejects B_z
  tmp2 <- withr::local_tempdir(); f2 <- .mono_fixture(tmp2)
  .mono_run(f2, support_cutoff = 95)
  expect_identical(.mono_read(f2, "species")$status[.mono_read(f2, "species")$group == "B_z"], "not_rejected")
})

test_that("a species with one tip is not evaluable, and is counted", {
  tmp <- withr::local_tempdir(); f <- .mono_fixture(tmp)
  .mono_run(f)
  s <- .mono_read(f, "species")
  expect_identical(s$status[s$group == "A_y"], "not_evaluable")
  sm <- .mono_read(f, "summary")
  expect_identical(c(sm$species_monophyletic, sm$species_rejected, sm$species_not_rejected,
                     sm$species_not_evaluable), c(2L, 1L, 1L, 1L))
  expect_equal(sm$species_prop_monophyletic, 2 / 4)
})

test_that("genera with two or more species are evaluated; a genus of one species is not", {
  tmp <- withr::local_tempdir(); f <- .mono_fixture(tmp)
  .mono_run(f)
  g <- .mono_read(f, "genus")
  expect_identical(g$status[g$group == "A"], "rejected")
  expect_identical(g$status[g$group == "B"], "rejected")
  expect_identical(g$status[g$group == "Leu"], "not_evaluable")
  expect_identical(g$species_in_locus[g$group == "A"], 2L)
})

test_that("grades: one lineage nested is paraphyletic, anything else polyphyletic", {
  tmp <- withr::local_tempdir(); f <- .mono_fixture(tmp)
  .mono_run(f)
  g <- .mono_read(f, "genus")
  expect_identical(g$grade[g$group == "A"], "paraphyletic")
  expect_identical(g$grade[g$group == "B"], "polyphyletic")
  s <- .mono_read(f, "species")
  expect_identical(s$grade[s$group == "B_z"], "polyphyletic")
  expect_true(is.na(s$grade[s$group == "A_x"]))
})

test_that("the root is the earliest lineage of the reference present in the locus, and a scattered one is flagged", {
  tmp <- withr::local_tempdir(); f <- .mono_fixture(tmp)
  .mono_run(f)
  sm <- .mono_read(f, "summary")
  expect_true(sm$rooted)
  expect_identical(sm$root_lineage, "Leu")
  tmp2 <- withr::local_tempdir(); f2 <- .mono_fixture(tmp2, .mono_tree_scattered)
  .mono_run(f2)
  sm2 <- .mono_read(f2, "summary")
  expect_false(sm2$rooted)
  expect_true(is.na(sm2$root_lineage))
  # Without a root, no grade is classified; the statuses do not need one
  g2 <- .mono_read(f2, "genus")
  expect_true(all(is.na(g2$grade)))
  expect_identical(.mono_read(f2, "species")$status[.mono_read(f2, "species")$group == "Leu_a"], "rejected")
})

test_that("each genus is set beside the reference and its compartment", {
  tmp <- withr::local_tempdir(); f <- .mono_fixture(tmp)
  .mono_run(f)
  r <- .mono_read(f, "reference")
  a <- r[r$genus == "A" & r$locus == "matK", ]
  expect_identical(a$compartment, "plastid")
  expect_identical(a$gene_tree_status, "rejected")
  # In the reference, A_x and A_y are not sister: genus A is not monophyletic there either
  expect_false(a$reference_monophyletic)
  expect_identical(.mono_read(f, "genus")$compartment[1], "plastid")
})

test_that("the job script of a locus carries the locus, its seed, the model and nothing unasked", {
  tmp <- withr::local_tempdir()
  lib <- file.path(tmp, "4_library"); dir.create(lib)
  writeLines(c(">Leu_a|L1", "ACGTACGTAC", ">A_x|X1", "ACGTACGTAA", ">A_x|X2", "ACGTACGTAG", ">Zz_q|Q1", "ACGTACGTTT"),
             file.path(lib, "LIB_matK.fasta"))
  trees <- file.path(tmp, "trees")
  suppressMessages(generate_barcoding_gene_tree_scripts(library_dir = lib, trees_dir = trees, seed = 1L,
                                                         cluster_mail_user = ""))
  job <- readLines(file.path(trees, "job", "job_matK.sh"))
  txt <- paste(job, collapse = "\n")
  # The conventions of the ML search of the phylogeny on Leftraru (run_ml_search.sh of 17-09):
  # the MPI build and its modules, the alignment parsed to RBA first, workers over the starting
  # trees, perf_threads and thread-nopin, thread binding off, the tree specification quoted
  expect_match(txt, "module load gcc/14.2.0-nlhpc openmpi/5.0.3-o raxml-ng/1.1.0-mpi-zen4-n", fixed = TRUE)
  expect_match(txt, "raxml-ng-mpi --parse --msa", fixed = TRUE)
  expect_match(txt, "raxml-ng-mpi --all --msa", fixed = TRUE)
  expect_match(txt, ".raxml.rba", fixed = TRUE)
  expect_match(txt, "--model GTR+G4", fixed = TRUE)
  expect_match(txt, "--tree 'pars{10},rand{10}'", fixed = TRUE)
  expect_match(txt, "--bs-trees 'autoMRE{1000}'", fixed = TRUE)
  expect_match(txt, "--workers ", fixed = TRUE)
  expect_match(txt, "--force perf_threads --extra thread-nopin", fixed = TRUE)
  expect_match(txt, "export OMP_PROC_BIND=false", fixed = TRUE)
  expect_match(txt, "--bs-metric fbp", fixed = TRUE)
  expect_match(txt, paste0("--seed ", .bc_query_seed(1L, "matK", "gene_tree", 0L, "<raxml>")), fixed = TRUE)
  expect_match(txt, "LIB_matK.fasta", fixed = TRUE)
  expect_false(grepl("tree-constraint", txt, fixed = TRUE))
  expect_false(grepl("mail-user|#SBATCH -q ", txt))
  expect_true(file.exists(file.path(trees, "job", "submit.sh")))
})

test_that("with the constraint of Module 9, each tip joins the clade of its species and absent species stay free", {
  tmp <- withr::local_tempdir()
  lib <- file.path(tmp, "4_library"); dir.create(lib)
  writeLines(c(">Leu_a|L1", "ACGTACGTAC", ">Leu_a|L2", "ACGTACGTAC", ">Opu_b|O1", "ACGTACGTAA",
               ">Opu_c|O2", "ACGTACGTAG", ">Cac_d|C1", "ACGTACGTTT", ">Cac_d|C2", "ACGTACGTTA",
               ">Zz_q|Q1", "ACGTACGTTG"), file.path(lib, "LIB_matK.fasta"))
  csv <- file.path(tmp, "constraints.csv")
  utils::write.csv(data.frame(Genus = c("Leu", "Opu", "Opu", "Cac"),
                              Specie_name = c("Leu_a", "Opu_b", "Opu_c", "Cac_d"),
                              Clade = c("Leuenbergeria", "Opuntieae", "Opuntieae", "Cacteae")),
                   csv, row.names = FALSE)
  trees <- file.path(tmp, "trees")
  suppressMessages(generate_barcoding_gene_tree_scripts(library_dir = lib, trees_dir = trees,
                                                         constraint = "module9", constraints_csv = csv,
                                                         cluster_mail_user = ""))
  job <- paste(readLines(file.path(trees, "job", "job_matK.sh")), collapse = "\n")
  expect_match(job, "--tree-constraint", fixed = TRUE)
  ct <- ape::read.tree(file.path(trees, "job", "constraint_matK.tree"))
  expect_setequal(ct$tip.label, c("Leu_a|L1", "Leu_a|L2", "Opu_b|O1", "Opu_c|O2", "Cac_d|C1", "Cac_d|C2"))
  u <- ape::unroot(ct)
  for (cl in list(c("Leu_a|L1", "Leu_a|L2"), c("Opu_b|O1", "Opu_c|O2"), c("Cac_d|C1", "Cac_d|C2"))) {
    expect_true(ape::is.monophyletic(u, cl, reroot = TRUE))
  }
})
