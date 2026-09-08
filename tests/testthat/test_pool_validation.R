################################################################################

context("POOLSEQ_VALIDATION")

################################################################################

set.seed(1)
pool <- matrix(stats::runif(4 * 20, min = 0.1, max = 0.9), nrow = 4)
class(pool) <- c("pcadapt_pool", "matrix", "array")

expect_s3_class(pcadapt(pool, K = 2), "pcadapt")

# Frequency-matrix validation
invalid.low <- pool
invalid.low[1, 1] <- -0.1
expect_error(pcadapt(invalid.low, K = 2), "between 0 and 1", fixed = TRUE)

invalid.high <- pool
invalid.high[1, 1] <- 1.1
expect_error(pcadapt(invalid.high, K = 2), "between 0 and 1", fixed = TRUE)

invalid.infinite <- pool
invalid.infinite[1, 1] <- Inf
expect_error(pcadapt(invalid.infinite, K = 2), "finite values", fixed = TRUE)

invalid.nan <- pool
invalid.nan[1, 1] <- NaN
expect_error(pcadapt(invalid.nan, K = 2), "use NA instead", fixed = TRUE)

missing.marker <- pool
missing.marker[, 1] <- NA_real_
expect_error(pcadapt(missing.marker, K = 2),
             "marker must have at least one observed frequency", fixed = TRUE)

missing.pool <- pool
missing.pool[1, ] <- NA_real_
expect_error(pcadapt(missing.pool, K = 2),
             "pool must have at least one observed marker", fixed = TRUE)

# Pool-seq argument validation
expect_error(pcadapt(pool, K = 1.5), "one positive integer", fixed = TRUE)
expect_error(pcadapt(pool, K = c(1, 2)), "one positive integer", fixed = TRUE)
expect_error(pcadapt(pool, K = NA_real_), "one positive integer", fixed = TRUE)
expect_error(pcadapt(pool, K = 2, method = "unknown"),
             "mahalanobis", fixed = TRUE)
expect_error(pcadapt(pool, K = 2, min.maf = NA_real_),
             "one finite number", fixed = TRUE)
expect_error(pcadapt(pool, K = 2, ploidy = 2),
             "must be NULL", fixed = TRUE)
expect_error(pcadapt(pool, K = 2, LD.clumping = list(size = 10, thr = 0.2)),
             "not implemented", fixed = TRUE)
expect_error(pcadapt(pool, K = 2, pca.only = NA),
             "TRUE or FALSE", fixed = TRUE)
expect_error(pcadapt(pool, K = 2, pca.only = TRUE),
             "not implemented", fixed = TRUE)
expect_error(pcadapt(pool, K = 2, tol = 1e-4),
             "not used for Pool-seq", fixed = TRUE)

# Rank is limited by centring and the retained marker count.
three.pools <- pool[1:3, , drop = FALSE]
class(three.pools) <- c("pcadapt_pool", "matrix", "array")
expect_error(pcadapt(three.pools, K = 3), "K cannot exceed 2", fixed = TRUE)

two.pools <- pool[1:2, , drop = FALSE]
class(two.pools) <- c("pcadapt_pool", "matrix", "array")
expect_error(pcadapt(two.pools, K = 2), "K cannot exceed 1", fixed = TRUE)

no.markers <- pool
no.markers[] <- 0
expect_error(pcadapt(no.markers, K = 1), "No Pool-seq markers remain",
             fixed = TRUE)

# The MAF threshold is inclusive, as it is for genotype input.
boundary <- matrix(
  c(0, 0, 0.2, 0.2, rep(c(0.2, 0.4, 0.6, 0.8), 2)),
  nrow = 4
)
class(boundary) <- c("pcadapt_pool", "matrix", "array")
boundary.result <- pcadapt(boundary, K = 1, min.maf = 0.1)
expect_true(1L %in% boundary.result$pass)

################################################################################
