################################################################################

context("VALIDATION")

################################################################################

set.seed(1)
geno <- matrix(rbinom(20 * 30, size = 2, prob = 0.3), nrow = 20)
class(geno) <- c("pcadapt_matrix", "matrix", "array")

# Scalar analysis arguments
expect_error(pcadapt(geno, K = 1.5), "one positive integer", fixed = TRUE)
expect_error(pcadapt(geno, K = c(1, 2)), "one positive integer", fixed = TRUE)
expect_error(pcadapt(geno, K = NA_real_), "one positive integer", fixed = TRUE)
expect_error(pcadapt(geno, K = Inf), "one positive integer", fixed = TRUE)
expect_error(pcadapt(geno, K = 0), "one positive integer", fixed = TRUE)

expect_error(pcadapt(geno, K = 2, min.maf = NA_real_),
             "one finite number", fixed = TRUE)
expect_error(pcadapt(geno, K = 2, ploidy = 0),
             "either 1 (haploid) or 2 (diploid)", fixed = TRUE)
expect_error(pcadapt(geno, K = 2, ploidy = 1.5),
             "either 1 (haploid) or 2 (diploid)", fixed = TRUE)
expect_error(pcadapt(geno, K = 2, ploidy = 3),
             "either 1 (haploid) or 2 (diploid)", fixed = TRUE)
expect_error(pcadapt(geno, K = 2, pca.only = NA),
             "TRUE or FALSE", fixed = TRUE)
expect_error(pcadapt(geno, K = 2, tol = 0),
             "positive finite number", fixed = TRUE)

# LD-clumping arguments
expect_error(
  pcadapt(geno, K = 2, LD.clumping = list(size = 10)),
  "containing size and thr",
  fixed = TRUE
)
expect_error(
  pcadapt(geno, K = 2, LD.clumping = list(size = 1.5, thr = 0.2)),
  "size must be one positive integer",
  fixed = TRUE
)
expect_error(
  pcadapt(geno, K = 2, LD.clumping = list(size = 10, thr = 1.1)),
  "thr must be one finite number between 0 and 1",
  fixed = TRUE
)

# Rank after marker filtering
few.variable <- matrix(0L, nrow = 20, ncol = 10)
few.variable[, 1] <- rep(c(0L, 1L), each = 10)
few.variable[, 2] <- rep(c(1L, 2L), each = 10)
class(few.variable) <- c("pcadapt_matrix", "matrix", "array")

expect_error(
  pcadapt(few.variable, K = 2),
  "2 retained markers",
  fixed = TRUE
)
expect_error(
  pcadapt(few.variable, K = 1, min.maf = 0.45),
  "No markers remain",
  fixed = TRUE
)

# Individuals without retained genotype calls
missing.individual <- geno
missing.individual[1, ] <- NA_integer_
expect_error(
  pcadapt(missing.individual, K = 2),
  "individual has no called genotypes",
  fixed = TRUE
)

################################################################################
