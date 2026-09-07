# A node whose age equals one of its own bounds was not estimated: penalized likelihood returned
# the constraint. These tests hold the part of that check that a single chronogram cannot see.
#
# On 2026-09-02 the maximum-likelihood tree placed ACP_root at 52.96 Ma, 0.41 Ma inside its upper
# bound of 53.37, and report_bound_adherence() called it interior. Across the 100 bootstrap
# replicates the interval was 53.29 to 53.37 and 96 returned the bound exactly. The point estimate
# was one realisation of a node whose age the data cannot identify, and it landed just inside by
# chance. The fixtures below reproduce that situation deliberately.

fixture_calibs <- function(min = 1, max = 5) {
  data.frame(mrca = "root", column = "Family", value = "F1;F2",
             min = min, max = max, used_in_analysis = TRUE, stringsAsFactors = FALSE)
}
fixture_constraints <- function() {
  data.frame(Specie_name = c("A", "B", "C", "D"), Family = c("F1", "F1", "F2", "F2"),
             stringsAsFactors = FALSE)
}
# Ultrametric by construction, with the root age written into the branch lengths.
tree_rooted_at <- function(age) {
  a <- age - 3; c <- age - 2
  ape::read.tree(text = sprintf("((A:%s,B:%s):3,(C:%s,D:%s):2);", a, a, c, c))
}

test_that("a point estimate just inside its bound is reported as interior, and said to be a point", {
  skip_if_not_installed("ape")
  ml <- tree_rooted_at(4.9)          # bound at 5, tolerance 0.05: interior by 0.05 to spare

  expect_message(
    out <- report_bound_adherence(ml, fixture_calibs(), fixture_constraints()),
    "single point estimate"
  )
  expect_equal(out$status, "interior")
  expect_true(is.na(out$n_bs))
})

test_that("the replicates overrule the point when most of them sit on the bound", {
  skip_if_not_installed("ape")
  ml <- tree_rooted_at(4.9)
  # Nine replicates on the bound, one inside: the ACP_root pattern in miniature.
  bs <- c(replicate(9, tree_rooted_at(5.0), simplify = FALSE), list(tree_rooted_at(4.9)))
  class(bs) <- "multiPhylo"

  out <- suppressWarnings(
    report_bound_adherence(ml, fixture_calibs(), fixture_constraints(), bootstraps = bs))

  expect_equal(out$n_bs, 10L)
  expect_equal(out$bs_pct_at_bound, 90)
  expect_equal(out$status, "interior")          # el punto sigue diciendo lo mismo
  expect_equal(round(out$bs_median, 2), 5)

  # Two warnings, and the second is the one that matters: it names the disagreement.
  w <- character(0)
  withCallingHandlers(
    report_bound_adherence(ml, fixture_calibs(), fixture_constraints(), bootstraps = bs),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("most bootstrap replicates", w)))
  expect_true(any(grepl("would have said otherwise", w)))
})

test_that("a node interior in the replicates as well raises nothing", {
  skip_if_not_installed("ape")
  ml <- tree_rooted_at(3.0)
  bs <- lapply(c(2.8, 2.9, 3.0, 3.1, 3.2), tree_rooted_at)
  class(bs) <- "multiPhylo"

  w <- character(0)
  out <- withCallingHandlers(
    report_bound_adherence(ml, fixture_calibs(), fixture_constraints(), bootstraps = bs),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })

  expect_equal(out$bs_pct_at_bound, 0)
  expect_length(w, 0L)
  expect_equal(out$bs_lo, 2.8)
  expect_equal(out$bs_hi, 3.2)
})

test_that("bootstraps can be given as the Newick file the pipeline writes", {
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "bsTree_treePL.tree")
  bs <- c(replicate(8, tree_rooted_at(5.0), simplify = FALSE), list(tree_rooted_at(4.5)))
  class(bs) <- "multiPhylo"
  ape::write.tree(bs, f)

  out <- suppressWarnings(report_bound_adherence(
    tree_rooted_at(4.9), fixture_calibs(), fixture_constraints(), bootstraps = f))
  expect_equal(out$n_bs, 9L)
  expect_true(out$bs_pct_at_bound > 50)

  expect_error(report_bound_adherence(
    tree_rooted_at(4.9), fixture_calibs(), fixture_constraints(),
    bootstraps = file.path(tmp, "absent.tree")), "does not exist")
})

test_that("a bound hit in the maximum-likelihood tree is still reported without replicates", {
  skip_if_not_installed("ape")
  # The original behaviour has to survive: this is the August 2026 case, where every calibrated
  # node came back on a bound and the chronogram gave no sign of it.
  expect_warning(
    suppressMessages(report_bound_adherence(
      tree_rooted_at(5.0), fixture_calibs(), fixture_constraints())),
    "returned their own bound"
  )
})
