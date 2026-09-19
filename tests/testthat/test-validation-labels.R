test_that(".format_reference_labels() italicises et al. and leaves other labels as text", {
  lbl <- .format_reference_labels(c(a = "Zuntini et al. (2024)", b = "This study",
                                    c = "de Vos et al. (2025)", d = 'Tree "B"'))
  expect_true(is.expression(lbl))
  expect_length(lbl, 4L)
  expect_match(deparse(lbl[[1]]), 'italic("et al.")', fixed = TRUE)
  expect_match(deparse(lbl[[1]]), "Zuntini", fixed = TRUE)
  expect_match(deparse(lbl[[1]]), "(2024)", fixed = TRUE)
  expect_identical(eval(lbl[[2]]), "This study")
  expect_match(deparse(lbl[[3]]), "de Vos", fixed = TRUE)
  expect_identical(eval(lbl[[4]]), 'Tree "B"')
})
