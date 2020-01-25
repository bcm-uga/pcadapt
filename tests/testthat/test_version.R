################################################################################

context("NEW VERSION VS OLD VERSION")

################################################################################

G <- readRDS(system.file("testdata", "full.rds", package = "pcadapt")) 

old.pcadapt <- readRDS(
  system.file("testdata", "pcadapt_v3.rds", package = "pcadapt")) 

mat <- read.pcadapt(G)
expect_s3_class(mat, "pcadapt_matrix")
new.pcadapt <- pcadapt(mat, K = 5)

# Same results as before when no missing values
expect_equal(cor(old.pcadapt$scores, new.pcadapt$scores) ** 2, diag(5))
expect_equal(as.vector(old.pcadapt$stat), new.pcadapt$stat, tol = 1e-6)

# New 'tol' parameters is used
expect_equal(pcadapt(mat, K = 5), new.pcadapt)
expect_failure(expect_equal(pcadapt(mat, K = 5, tol = 1e-1), new.pcadapt))

################################################################################
