context("genomic inflation factor")

test_that("positive genomic inflation factors are retained", {
  zscores <- matrix(c(-2, -1, 1, 2), ncol = 1)

  mahalanobis <- pcadapt:::get_statistics(
    zscores,
    method = "mahalanobis",
    pass = rep(TRUE, nrow(zscores))
  )
  componentwise <- pcadapt:::get_statistics(
    zscores,
    method = "componentwise",
    pass = rep(TRUE, nrow(zscores))
  )

  expect_gt(mahalanobis$gif, 0)
  expect_gt(componentwise$gif, 0)
  expect_true(all(is.finite(mahalanobis$pvalues)))
  expect_true(all(is.finite(componentwise$pvalues)))
})

test_that("zero genomic inflation factors stop the analysis", {
  zscores <- matrix(1, nrow = 10, ncol = 1)

  expect_error(
    pcadapt:::get_statistics(
      zscores,
      method = "mahalanobis",
      pass = rep(TRUE, nrow(zscores))
    ),
    "Mahalanobis statistics.*strictly positive"
  )
  expect_error(
    pcadapt:::get_statistics(
      matrix(0, nrow = 10, ncol = 2),
      method = "componentwise",
      pass = rep(TRUE, 10)
    ),
    "component 1.*strictly positive"
  )
})

test_that("non-finite statistics stop the analysis", {
  zscores <- matrix(NA_real_, nrow = 10, ncol = 1)

  expect_error(
    pcadapt:::get_statistics(
      zscores,
      method = "componentwise",
      pass = rep(TRUE, nrow(zscores))
    ),
    "component 1.*no finite test statistics"
  )
})
