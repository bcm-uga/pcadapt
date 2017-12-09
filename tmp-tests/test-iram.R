N <- 150; M <- 1500
library(bigsnpr)
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

if (class(input) == "character") {
  path_to_bed <- normalizePath(input)
  p <- nrow(data.table::fread(sub("\\.bed$", ".bim", path_to_bed)))
  n <- nrow(data.table::fread(sub("\\.bed$", ".fam", path_to_bed)))
  
  # File mapping
  xptr <- pcadapt:::bedXPtr(path_to_bed, n, p)
  
} else if (class(input) == "matrix") {
  # an input matrix has nIND rows and nSNP columns
  xptr <- input
  n <- nrow(xptr)
  p <- ncol(xptr)
}

lookup_byte <- pcadapt:::getCode()
lookup_geno <- rbind(rep(0, p), 1, 2, 3)

# Get allele frequencies
# Uses a non-scaled lookup table
af <- pcadapt:::get_af(xptr, lookup_geno, lookup_byte, 1:p)
stopifnot(all.equal(af, colMeans(G[], na.rm = TRUE) / 2))


# Create a logical vector to locate SNPs with mAF >= min.maf
maf <- pmin(af, 1 - af)
pass.af <- (maf >= min.maf)

# Create a logical vector to locate SNPs that have been clumped
if (LD.clumping) {
  pass <- pcadapt:::clumping_r(xptr = xptr, 
                     lookup_geno = lookup_geno, 
                     lookup_byte = lookup_byte, 
                     ind_col = which(pass.af),
                     maf = maf[pass.af], 
                     size = size, 
                     thr = thr)
} else {
  pass <- pass.af
}

no.pass <- which(!pass)
corr <- cor(G[], use = "pairwise.complete.obs")
tmp <- which(corr^2 > 0.2, arr.ind = TRUE)
tmp[tmp[, 1] > tmp[, 2] & (tmp[, 1] - size) <= tmp[, 2], ]
no.pass

ind.pass <- which(pass)
p2 <- length(ind.pass)

# Get number of non-missing values per row and per column
# Uses a non-scaled lookup table
nb_nona <- pcadapt:::nb_nona(xptr, lookup_geno, lookup_byte, ind.pass)
stopifnot(all.equal(nb_nona$n, n - big_counts(G, ind.col = ind.pass)[4, ]))
stopifnot(all.equal(nb_nona$p, p2 - big_counts(G, ind.col = ind.pass, byrow = TRUE)[4, ]))

# Scaled lookup table 
lookup_scale <- rbind(
  outer(0:2, af[ind.pass], function(g, p) {
    (g - 2 * p) / sqrt(2 * p * (1 - p))
  }), 
  0
)

### SVD using RSpectra
obj.svd <- RSpectra::svds(
  A = function(x, args) {
    # When filtering, the actual number of SNPs that we have is actually
    # sum(pass) and not p anymore
    pcadapt:::pMatVec4(xptr, x, lookup_scale, lookup_byte, ind.pass) / 
      nb_nona$p * p2
  }, 
  Atrans = function(x, args) {
    # NB: nb_nona$n depends on 'pass' as well
    pcadapt:::cpMatVec4(xptr, x, lookup_scale, lookup_byte, ind.pass) / 
      nb_nona$n * n
  },
  k = K, 
  dim = c(n, p2),
  opts = list(tol = 1e-4, maxitr = 100)
)

# Lookup table
lookup_scale2 <- rbind(
  outer(0:2, af[pass.af], function(g, p) {
    (g - 2 * p) / sqrt(2 * p * (1 - p))
  }), 
  0
)

# Multiple Linear Regression is performed also on SNPs that have been clumped,
# that is why we recompute the lookup table

ind.pass.af <- which(pass.af)
obj.svd$zscores <- pcadapt:::multLinReg(xptr,
                              lookup_scale2,
                              lookup_byte,
                              ind.pass.af,
                              obj.svd$u,
                              obj.svd$d,
                              obj.svd$v)


test2 <- matrix(NA_real_, length(ind.pass.af), K) 
for (j in seq_along(ind.pass.af)) {
  print(j)
  for (k in seq_len(K)) {
    test2[j, k] <- 
      summary(lm(G[, ind.pass.af[j]] ~ obj.svd$u[, k]))$coefficients[2, 3]
  }
}

obj.svd$pass <- pass.af
obj.svd$d <- obj.svd$d^2 / sum(pass)
obj.svd$maf <- pmin(af, 1 - af)
obj.svd$nona1 <- nb_nona$p
obj.svd$nona2 <- nb_nona$n

obj.svd