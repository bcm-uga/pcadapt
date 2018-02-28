################################################################################

context("MISSING_VALUES")

res <- list()
tol <- 1e-5

for (ext in c("bed", "lfmm", "pcadapt")) {
  file_path <- system.file("extdata", 
                           paste0("missing.", ext), 
                           package = "pcadapt") 
  x <- read.pcadapt(file_path, type = ext)
  res[[ext]] <- pcadapt(x, K = 5)
}

testthat::expect_equal(res[["bed"]]["scores"], res[["lfmm"]]["scores"], 
                       tolerance = tol)
testthat::expect_equal(res[["bed"]]["scores"], res[["pcadapt"]]["scores"], 
                       tolerance = tol)
testthat::expect_equal(res[["lfmm"]]["scores"], res[["pcadapt"]]["scores"], 
                       tolerance = tol)
