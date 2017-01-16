# Binomial sampling test routine

context("SAMPLE_GENO_MATRIX")

test_that("sampled frequencies match with observed frequencies", {
  nPOOL <- 2
  nSNP <- 2
  nINDperPOOL <- 10000
  sample_size <- rep(nINDperPOOL, nPOOL)
  
  nIND <- sum(sample_size)
  p <- runif(nPOOL * nSNP)
  maf <- pmin(p, 1 - p)
  dt <- matrix(maf, nrow = nPOOL, ncol = nSNP)
  ploidy <- 2
  
  x <- t(sample_geno_matrix(dt, ploidy, sample_size))
  
  freq <- array(0, dim = c(nPOOL, nSNP))
  sd <- array(0, dim = c(nPOOL, nSNP))
  beg.pool <- 1
  end.pool <- 0
  for (k in 1:nPOOL){
    end.pool <- end.pool + sample_size[k]
    freq[k,] <- apply(x[beg.pool:end.pool, ], MARGIN = 2, FUN = function(h){mean(h, na.rm = TRUE)}) / ploidy
    beg.pool <- end.pool + 1
  }
  
  for (pool.id in 1:nPOOL){
    testthat::expect_equal(freq[pool.id, ], dt[pool.id, ], tolerance = 0.01)
  }
  
})

test_that("missing frequencies raise missing genotypes encoded with 9", {
  nPOOL <- 2
  nSNP <- 2
  nINDperPOOL <- 10
  sample_size <- rep(nINDperPOOL, nPOOL)
  
  nIND <- sum(sample_size)
  p <- runif(nPOOL * nSNP)
  maf <- pmin(p, 1 - p)
  dt <- matrix(maf, nrow = nPOOL, ncol = nSNP)
  dt[1, 1] <- 9
  dtna <- dt
  dtna[1,1] <- NA
  ploidy <- 2
  
  x <- sample_geno_matrix(dt, ploidy, sample_size)
  y <- sample_geno_matrix(dtna, ploidy, sample_size)
  
  for (ind in 1:nINDperPOOL){
    testthat::expect_equal(x[1, ind], 9)
    testthat::expect_equal(y[1, ind], 9)
  }
  
})