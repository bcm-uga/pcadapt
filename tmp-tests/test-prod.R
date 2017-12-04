popres_bed <- "tmp-data/testfile.bed"
# popres_bed <- "../POPRES_data/POPRES_allchr_QC_norel.bed"
popres_bed <- normalizePath("~/Downloads/plink_mac/POPRES_filtered_0.0001.bed")
p <- nrow(data.table::fread(sub("\\.bed$", ".bim", popres_bed)))
n <- nrow(data.table::fread(sub("\\.bed$", ".fam", popres_bed)))
Rcpp::sourceCpp('src/bed-acc-xptr.cpp')
test <- bedXPtr(popres_bed, n, p)

getCode <- function(NA.VAL = 3L) {
  geno.raw <- as.logical(rawToBits(as.raw(0:255)))
  s <- c(TRUE, FALSE)
  geno1 <- geno.raw[s]
  geno2 <- geno.raw[!s]
  geno <- geno1 + geno2
  geno[geno1 & !geno2] <- NA.VAL
  dim(geno) <- c(4, 256)
  geno
}

Rcpp::sourceCpp('src/af.cpp')
tmp <- af(test, rbind(rep(0, p), 1, 2, 3), getCode())


lookup_scale <- rbind(outer(0:2, tmp, function(g, p) {
  (g - 2 * p) / sqrt(2 * p * (1 - p))
}), 0)

x <- rnorm(p)
Rcpp::sourceCpp('src/prodMat.cpp')
lookup_byte <- getCode()
tmp2 <- pMatVec4(test, x, lookup_scale, lookup_byte)


# library(pcadapt)
# test2 <- bedadapt(popres_bed)
# tmp3 <- pcadapt:::prodMatVec(test2$address, x, 
#                              m = 2 * tmp, 
#                              s = sqrt(2 * tmp * (1 - tmp)), 
#                              pass = rep(TRUE, p)) 

Rcpp::sourceCpp('src/nb-na.cpp')
system.time({
  nb_nona <- nb_nona(test, rbind(rep(0, p), 1, 2, 3), lookup_byte)
  tmp3 <- RSpectra::svds(
    A = function(x, args) {
      cat(".")
      pcadapt:::pMatVec4(test, x, lookup_scale, lookup_byte) / nb_nona[[1]] * n
    }, 
    k = 5, 
    Atrans = function(x, args) {
      pcadapt:::cpMatVec4(test, x, lookup_scale, lookup_byte) / nb_nona[[2]] * p
    },
    dim = c(n, p),
    opts = list(tol = 1e-4)
  )
})

# plot(tmp3$u)

# mat <- as.matrix(read.table("tmp-data/testfile.pcadapt"))
# tmp4 <- pcadapt::pcadapt(mat, K = 10, min.maf = 0)
# plot(tmp4$scores)


system.time(
  tmp5 <- flashpcaR::flashpca(sub("\\.bed$", "", popres_bed), ndim = 5)
)
sqrt(tmp5$values / n) / tmp3$d

# plot(tmp5$vectors)
plot(tmp3$u[, 1], tmp5$vectors[, 1])
abline(0, 1, col = "red")

# library(bigsnpr)
# G <- snp_attach("../paper-packages/backingfiles/POPRESQC.rds")$genotypes
# system.time(
#   tmp6 <- big_randomSVD(G, snp_scaleBinom(), k = 5, ncores = 2)
# )
# plot(tmp6$u[, 1], tmp3$u[, 1])

G <- as.matrix(read.table("~/Documents/thesis/datasets/mv_nunif_0.2.pcadapt"))
G[G == 9] <- NA
X <- runif(nrow(G))

tmp <- apply(G, MARGIN = 1, FUN = function(x) {mean(x, na.rm = TRUE) / 2})

lookup_byte <- pcadapt:::getCode()

lookup_geno <- rbind(outer(0:2, tmp, function(g, p) {
  if (p > 0 || p < 1) {
    return((g - 2 * p) / sqrt(2 * p * (1 - p)))
  } else {
    return(0)
  }
}), 0)

pass <- rep(TRUE, nrow(G))
nb_nona <- pcadapt:::nb_nona(t(G), lookup_geno, lookup_byte, pass, sum(pass))

x <- pcadapt:::pMatVec4(t(G), X, lookup_scale = lookup_geno, lookup_byte = lookup_byte) / nb_nona[[1]]
y <- pcadapt:::prodtGx(G, X, tmp)
all.equal(x, y)

xx <- pcadapt:::prodGx(G, X[1:150], tmp)
yy <- pcadapt:::cpMatVec4(t(G), X[1:150], lookup_scale = lookup_geno,lookup_byte = lookup_byte) / nb_nona[[2]]
all.equal(xx, yy)
