context("GET BARYCENTRIC COORDINATES")

testthat::test_that("cart2bary_cpp behaves similarly to cart2bary", {
  X <- rbind(c(0, 0), c(0, 1), c(1, 0))
  P <- rbind(c(0.5, 0.5), c(0.1, 0.8))
  testthat::expect_equal(geometry::cart2bary(X, P), cart2bary_cpp(X, P), tolerance = 0.01)
})
