context("linear regression")

current_regression_reference <- function(genotypes, scores, ploidy = 2) {
  observed <- !is.na(genotypes)
  scores <- scores[observed, , drop = FALSE]
  genotypes <- genotypes[observed]

  allele.frequency <- sum(genotypes) / (ploidy * length(genotypes))
  scaled.genotypes <-
    (genotypes - ploidy * allele.frequency) /
    sqrt(ploidy * allele.frequency * (1 - allele.frequency))

  coefficients <- crossprod(scores, scaled.genotypes)
  fitted <- scores %*% coefficients
  residual.sum.squares <- sum((scaled.genotypes - fitted)^2)
  residual.variance <-
    residual.sum.squares / (length(genotypes) - ncol(scores))

  drop(coefficients) /
    sqrt(colSums(scores^2) * residual.variance)
}

test_that("observed standardized values equal to three are retained", {
  genotypes <- c(2L, 2L, rep(0L, 9))
  scores <- qr.Q(qr(cbind(seq_along(genotypes), seq_along(genotypes)^2)))
  allele.frequency <- sum(genotypes) / (2 * length(genotypes))

  expect_identical(allele.frequency, 2 / 11)
  expect_identical(
    (2 - 2 * allele.frequency) /
      sqrt(2 * allele.frequency * (1 - allele.frequency)),
    3
  )

  observed <- pcadapt:::multLinReg(
    matrix(genotypes, ncol = 1),
    1L,
    allele.frequency,
    2,
    scores
  )
  expected <- current_regression_reference(genotypes, scores)

  expect_equal(drop(observed), expected)
})

test_that("genuinely missing genotypes remain excluded", {
  genotypes <- c(2L, NA_integer_, 1L, rep(0L, 8))
  scores <- qr.Q(qr(cbind(seq_along(genotypes), seq_along(genotypes)^2)))
  allele.frequency <- sum(genotypes, na.rm = TRUE) /
    (2 * sum(!is.na(genotypes)))

  observed <- pcadapt:::multLinReg(
    matrix(genotypes, ncol = 1),
    1L,
    allele.frequency,
    2,
    scores
  )
  expected <- current_regression_reference(genotypes, scores)

  expect_equal(drop(observed), expected)
})
