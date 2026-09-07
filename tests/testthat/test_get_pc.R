context("get.pc")

test_that("get.pc returns one deterministic component per marker", {
  x <- list(
    zscores = rbind(
      c(2, -2, 1),
      c(NA, 3, 1),
      c(NA, NA, NA),
      c(1, 2, 3)
    )
  )

  result <- get.pc(x, 1:4)

  expect_identical(result$SNP, 1:4)
  expect_identical(result$PC, c(1L, 2L, NA_integer_, 3L))
})

test_that("get.pc preserves marker order and duplicate indices", {
  x <- list(zscores = rbind(c(1, 3), c(4, 2)))

  result <- get.pc(x, c(2, 1, 2))

  expect_identical(result$SNP, c(2, 1, 2))
  expect_identical(result$PC, c(1L, 2L, 1L))
})

test_that("get.pc handles an empty marker selection", {
  x <- list(zscores = matrix(numeric(), nrow = 0, ncol = 2))

  result <- get.pc(x, integer())

  expect_identical(result, data.frame(SNP = integer(), PC = integer()))
})
