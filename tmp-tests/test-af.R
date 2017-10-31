popres_bed <- "../Dubois2010_data/FinnuncorrNLITUK3hap550.bed"
popres_bed <- "../POPRES_data/POPRES_allchr_QC_norel.bed"
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

lookup_byte <- getCode()
lookup_scale <- rbind(rep(0, p), 1, 2, 3)

Rcpp::sourceCpp('src/af.cpp')
tmp <- af(test, lookup_scale, lookup_byte)

library(BEDMatrix)
x <- BEDMatrix(popres_bed)
colMeans(x[, 1:10], na.rm = TRUE) / 2
pmin(tmp[1:10], 1 - tmp[1:10])


X <- x[, 1:10]
typeof(X)
tmp2 <- af(X, lookup_scale, lookup_byte)

library(pcadapt)
test2 <- bedadapt(popres_bed)
tmp3 <- pcadapt:::cmpt_af(test2$address) 

microbenchmark::microbenchmark(
  tmp <- af(test, lookup_scale, lookup_byte),
  tmp3 <- pcadapt:::cmpt_af(test2$address),
  times = 5
)
all.equal(pmin(tmp, 1 - tmp), pmin(tmp3, 1 - tmp3))


test4 <- gaston::read.bed.matrix(sub("\\.bed$", "", popres_bed))
test4@p
