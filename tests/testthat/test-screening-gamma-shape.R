# AP-06 (closure review of 0.5.4): the saturation distance used ape::dist.dna(gamma = TRUE), which ape
# takes as a gamma shape of 1, without declaring it. From 0.5.5 the shape is the argument gamma_shape,
# default 1, and the default gives the slopes of earlier versions.

test_that("gamma_shape = 1 gives the slope of gamma = TRUE, and a smaller shape lowers it", {
  data(woodmouse, package = "ape", envir = environment())
  one <- PhyloCactus:::.test_saturation_proxy(woodmouse)
  old_d <- as.numeric(ape::dist.dna(woodmouse, model = "K80", pairwise.deletion = TRUE, gamma = TRUE))
  raw_d <- as.numeric(ape::dist.dna(woodmouse, model = "raw", pairwise.deletion = TRUE))
  expect_equal(one$slope, unname(stats::coef(stats::lm(raw_d ~ old_d))[2]))
  half <- PhyloCactus:::.test_saturation_proxy(woodmouse, gamma_shape = 0.2)
  expect_lt(half$slope, one$slope)
  expect_identical(eval(formals(run_marker_screening)$gamma_shape), 1)
  expect_error(run_marker_screening(tempdir(), tempdir(), gamma_shape = 0), "gamma_shape")
})
