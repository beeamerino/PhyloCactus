# AP-24 (closure review of 0.5.4): validate_phylogenies() fixed the family at "Cactaceae" and rooted
# with ape::root() directly, which leaves a tree rooted inside the genus as it was.

test_that("validate_phylogenies() takes the family and roots through root_on_clade()", {
  expect_identical(eval(formals(validate_phylogenies)$family), "Cactaceae")
  body_txt <- paste(deparse(body(validate_phylogenies)), collapse = " ")
  expect_match(body_txt, "root_on_clade(tr_family, genus_tips)", fixed = TRUE)
  expect_false(grepl("ape::root(tr_family", body_txt, fixed = TRUE))
  # root_on_clade() roots on the stem of a genus where ape::root() keeps a root inside it.
  tr <- ape::read.tree(text = "((L_a,(X_b,Y_c)),L_d);")
  r <- root_on_clade(tr, c("L_a", "L_d"))
  expect_true(ape::is.monophyletic(r, c("L_a", "L_d")))
})
