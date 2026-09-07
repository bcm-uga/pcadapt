################################################################################

context("POOLSEQ_PCA")

################################################################################

# The two-pool pathway is a rank-one PCA and should return genuine SVD scores
# and loadings rather than placeholders.
pool <- matrix(
  c(
    0.10, 0.20, 0.70, 0.80, 0.35, 0.65,
    0.30, 0.40, 0.50, 0.60, 0.55, 0.45
  ),
  nrow = 2,
  byrow = TRUE
)
class(pool) <- c("pcadapt_pool", "matrix", "array")

result <- pcadapt(pool, K = 1, min.maf = 0)
centred <- scale(unclass(pool), center = TRUE, scale = FALSE)
expected <- svd(centred)

expect_true(all(is.finite(result$scores)))
expect_equal(abs(result$scores[, 1]), abs(expected$u[, 1]))
expect_equal(abs(result$loadings[, 1]), abs(expected$v[, 1]))
expect_equal(result$singular.values, 1)
expect_equal(sum(result$loadings[, 1]^2), 1)

# PCA signs are arbitrary, but scores and loadings must reconstruct the same
# rank-one centred matrix when their signs are considered together.
observed.reconstruction <-
  result$scores[, 1, drop = FALSE] %*%
  t(result$loadings[, 1, drop = FALSE]) * expected$d[1]
expect_equal(
  observed.reconstruction,
  unname(centred),
  tolerance = 1e-12,
  check.attributes = FALSE
)

################################################################################
