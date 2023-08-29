library(pcadapt)
library(bigsnpr)
N <- 150; M <- 1500
test <- snp_fake(N, M)
G <- test$genotypes
G[] <- sample(as.raw(0:3), size = length(G), replace = TRUE)
G[1:5, 1:5]
input <- snp_writeBed(test, bedfile = tempfile(fileext = ".bed"))

min.maf <- 0.05
LD.clumping <- TRUE
size <- 200
thr <- 0.2
K <- 5

if (is.character(input)) {
  path_to_bed <- normalizePath(input)
  p <- nrow(data.table::fread(sub("\\.bed$", ".bim", path_to_bed)))
  n <- nrow(data.table::fread(sub("\\.bed$", ".fam", path_to_bed)))
  
  # File mapping
  xptr <- pcadapt:::bedXPtr(path_to_bed, n, p)
  
} else if (is.matrix(input)) {
  # an input matrix has nIND rows and nSNP columns
  xptr <- input
  n <- nrow(xptr)
  p <- ncol(xptr)
}

# Get allele frequencies
# Uses a non-scaled lookup table
af <- pcadapt:::get_af(xptr)
stopifnot(all.equal(af, colMeans(G[], na.rm = TRUE) / 2))


# Create a logical vector to locate SNPs with mAF >= min.maf
maf <- pmin(af, 1 - af)
ind.pass.af <- which(maf >= min.maf)

# Create a logical vector to locate SNPs that have been clumped
if (LD.clumping) {
  # Take also number of NAs into account?
  ord <- order(maf[ind.pass.af], decreasing = TRUE)
  pass <- pcadapt:::clumping(xptr, ind.pass.af, 
                   ord, rep(TRUE, length(ord)), 
                   size, thr)
  ind.pass <- ind.pass.af[pass]
} else {
  ind.pass <- ind.pass.af
}
p2 <- length(ind.pass)

no.pass <- setdiff(ind.pass.af, ind.pass)
corr <- cor(G[], use = "pairwise.complete.obs")
tmp <- which(corr^2 > 0.2, arr.ind = TRUE)
tmp[tmp[, 1] > tmp[, 2] & (tmp[, 1] - size) <= tmp[, 2], ]
no.pass
stopifnot(all(no.pass %in% tmp[tmp[, 1] > tmp[, 2] & (tmp[, 1] - size) <= tmp[, 2], ]))

# Get number of non-missing values per row and per column
nb_nona <- pcadapt:::nb_nona(xptr, ind.pass)
stopifnot(all.equal(nb_nona$n, n - big_counts(G, ind.col = ind.pass)[4, ]))
stopifnot(all.equal(nb_nona$p, p2 - big_counts(G, ind.col = ind.pass, byrow = TRUE)[4, ]))

### SVD using RSpectra
obj.svd <- RSpectra::svds(
  A = function(x, args) {
    # When filtering, the actual number of SNPs that we have is actually
    # sum(pass) and not p anymore
    pcadapt:::pMatVec4(xptr, ind.pass, af, x) / nb_nona$p * p2
  }, 
  Atrans = function(x, args) {
    # NB: nb_nona$n depends on 'pass' as well
    pcadapt:::cpMatVec4(xptr, ind.pass, af, x) / nb_nona$n * n
  },
  k = K, 
  dim = c(n, p2),
  opts = list(tol = 1e-4, maxitr = 100)
)

# Multiple Linear Regression is performed also on SNPs that have been clumped,
# that is why we recompute the lookup table

ind.pass.af <- which(pass.af)
obj.svd$zscores <- pcadapt:::multLinReg(xptr, ind.pass.af, af, obj.svd$u)

# test2 <- matrix(NA_real_, length(ind.pass.af), K) 
# test3 <- matrix(NA_real_, length(ind.pass.af), K) 
# for (j in seq_along(ind.pass.af)) {
#   print(j)
#   coefs <- summary(lm(G[, ind.pass.af[j]] ~ obj.svd$u - 1))$coefficients
#   test2[j, ] <- coefs[, 3]
#   test3[j, ] <- coefs[, 4]
# }
# j <- 3
# y <- G[, ind.pass.af[j]]
# y <- lookup_scale2[y + 1, 1]
# coefs <- summary(lm(y ~ obj.svd$u - 1))$coefficients
# nona <- !is.na(y)
# y2 <- y[nona]
# x <- obj.svd$u[nona, ]   
# 
# crossprod(x)             ## not orthogonal anymore :-(
# sd <- solve(crossprod(x))
# coefs[, 2] / sqrt(diag(sd))         ## OK if no intercept
# # sd2 <- solve(crossprod(cbind(1, x)))
# # coefs[-1, 2] / sqrt(diag(sd2)[-1])    ## OK
# # coefs[-1, 2] / sqrt(diag(sd2)[-1]) / sd(y2)
# 
# (beta <- sd %*% crossprod(x, y2))
# coefs[, 1]
# eps <- y2 - x %*% beta
# tmp <- sqrt(sum(eps * eps) * diag(sd) / (length(eps) - K))
# coefs[, 2] / tmp
# 
# cbind(res <- beta / tmp, coefs[, 3])
# cbind(pval <- 2 * pt(abs(res), df = length(eps) - K, lower.tail = FALSE),
#       coefs[, 4])
# 
# res2 <- beta / sqrt(diag(sd))
# plot(res, res2); abline(0, 1, col = "red")



obj.svd$pass <- pass.af
obj.svd$d <- obj.svd$d^2 / sum(pass)
obj.svd$maf <- pmin(af, 1 - af)
obj.svd$nona1 <- nb_nona$p
obj.svd$nona2 <- nb_nona$n

obj.svd
