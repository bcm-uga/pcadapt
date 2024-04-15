################################################################################

context("VAR_EXPL")

################################################################################

path_to_file <- system.file("extdata", "SSMPG2017.rds", package = "pcadapt")
mat <- read.pcadapt(readRDS(path_to_file), type = "pcadapt")

af <- colMeans(mat) / 2
ind_col <- which(pmin(af, 1 - af) > runif(1, 0.01, 0.1))
ind_col <- sample(ind_col)
var1 <- pcadapt:::total_var_scaled(mat, ind_col, af, 2)
var2 <- norm(scale(mat, 2 * af, sqrt(2 * af * (1 - af)))[, ind_col], "F")^2
expect_equal(var1, var2, tolerance = 1e-5)

ind_na <- sample(length(mat), length(mat) * 0.1)
mat[ind_na] <- NA
var3 <- pcadapt:::total_var_scaled(mat, ind_col, af, 2)
expect_equal(var1, var3, tolerance = 1e-2)

################################################################################

path_to_file <- system.file("extdata", "geno3pops.bed", package = "pcadapt")
input <- read.pcadapt(path_to_file, type = "bed")
mat2 <- pcadapt:::bed2matrix(input)
af2 <- colMeans(mat2) / 2
ind_col2 <- which(pmin(af2, 1 - af2) > runif(1, 0.01, 0.1))
ind_col2 <- sample(ind_col2)
var4 <- norm(scale(mat2, 2 * af2, sqrt(2 * af2 * (1 - af2)))[, ind_col2], "F")^2

n <- attr(input, "n")
p <- attr(input, "p")
xptr <- structure(pcadapt:::bedXPtr(input, n, p), n = n, p = p, class = "xptr_bed")
var5 <- pcadapt:::total_var_scaled(mat2, ind_col2, af2, 2)
expect_equal(var4, var5, tolerance = 1e-5)

################################################################################

x <- pcadapt(input = input, K = n - 1)
expect_equal(sum(x$singular.values^2), 1, tolerance = 1e-5)

x2 <- pcadapt(input = input, K = n - 1, LD.clumping = list(size = 200, thr = 0.2)) 
expect_lt(sum(x2$singular.values^2), 0.99)
expect_gt(sum(x2$singular.values^2), 0.90)

################################################################################
