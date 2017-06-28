context("GET BARYCENTRIC COORDINATES")

testthat::test_that("cart2bary_cpp behaves similarly to cart2bary", {
  X <- rbind(c(0, 0), c(0, 1), c(1, 0))
  P <- matrix(runif(4), nrow = 2)
  testthat::expect_equal(geometry::cart2bary(X, P), cart2bary_cpp(X, P), tolerance = 0.01)
})

testthat::test_that("centroids are exact", {
  n.ind <- 6
  K <- 3
  u <- matrix(runif(n = n.ind * K, min = -1, max = 1), ncol = K)
  pop <- c("POP1", "POP2", "POP3", "POP2", "POP1", "POP1")
  res <- matrix(0, nrow = 3, ncol = K)
  res[1, ] = (u[1, ] + u[5, ] + u[6, ]) / 3
  res[2, ] = (u[2, ] + u[4, ]) / 2
  res[3, ] = u[3, ]
  
  testthat::expect_equal(res, scores_centroids(u, pop), tolerance = 0.001)
  testthat::expect_equal(res[2:3, ], centroids_to_simplex(res, pop, "POP1"), tolerance = 0.001)
  testthat::expect_equal(res[c(1, 3), ], centroids_to_simplex(res, pop, "POP2"), tolerance = 0.001)
  testthat::expect_equal(res[1:2, ], centroids_to_simplex(res, pop, "POP3"), tolerance = 0.001)
})

