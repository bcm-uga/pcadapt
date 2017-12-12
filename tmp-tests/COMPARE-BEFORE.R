set.seed(1)

FILE <- "inst/testdata/to-compare.pcadapt"

N <- 300; P <- 1000
mat <- matrix(0, N, P)
for (j in 1:P) mat[, j] <- rbinom(N, size = 2, prob = runif(1))
mat[sample(length(mat), size = length(mat) / 5)] <- 9  ## NAs
write.table(t(mat), file = FILE, quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

# Version 3.0.4
library(pcadapt)
res <- pcadapt(FILE, K = 10, min.maf = 0.05)
str(res)
saveRDS(res, file = sub("\\.pcadapt$", ".rds", FILE))

# New version
bedfile <- writeBed(mmapcharr::mmapchar(FILE, mmapcharr:::CODE_012), 
                    is.pcadapt = TRUE)
obj.svd <- pcadapt:::iram(bedfile, K = 10, min.maf = 0.05)
str(obj.svd)

af <- obj.svd$af
maf <- pmin(af, 1 - af)
all.equal(maf, res$maf, tolerance = 1e-5)
all.equal(ind.NA <- which(is.na(obj.svd$v[, 1])), which(is.na(res$loadings[, 1])))

round(corr <- cor(obj.svd$v[-ind.NA, ], res$loadings[-ind.NA, ]), 2)
mean(diag(abs(corr)))

round(corr <- cor(res$scores, obj.svd$u), 2)
mean(diag(abs(corr)))


# Linear regression
zscores <- res$zscores[ind.pass.af, ]
zscores[1:5, 1:5]
xptr <- pcadapt:::bedXPtr(bedfile, N, P)
ind.pass.af <- which(maf >= 0.05)
af2 <- sapply(seq_along(maf), function(i) {
  maf <- res$maf[i]
  `if`(which.min(abs(maf - c(af[i], 1 - af[i]))) == 1, maf, 1 - maf)
})
plot(af, af2)
zscores1 <- pcadapt:::multLinReg(xptr, ind.pass.af, af2, res$scores)
zscores1[1:5, 1:5]
Rcpp::sourceCpp('tmp-save/linear-regression-noapprox.cpp')
zscores2 <- multLinReg(xptr, ind.pass.af, af, res$scores)
zscores2[1:5, 1:5]

plot(zscores1, zscores)
zscores1 / zscores
round(cor(zscores2, res$zscores[ind.pass.af, ]), 2)
plot(zscores2, zscores)
hist(zscores[, 1])
unname(apply(zscores, 2, function(x) shapiro.test(x)$p.value))
apply(zscores2, 2, function(x) shapiro.test(x)$p.value)


# CovRob
d <- pcadapt:::covRob_rcpp(zscores)$dist
gif <- median(d) / qchisq(0.5, df = 10)
all.equal(gif, res$gif)
all.equal(d, as.numeric(res$stat[ind.pass.af]))
all.equal(d / gif, as.numeric(res$chi2.stat[ind.pass.af]))
all.equal(pchisq(d / gif, df = 10, lower.tail = FALSE), res$pvalues[ind.pass.af])

plot(bigsnpr:::getD(zscores), d)
all.equal(robust::covRob(zscores, estim="pairwiseGK")$dist, d)
all.equal(bigsnpr:::getD(zscores), d)
microbenchmark::microbenchmark(
  pcadapt:::covRob_rcpp(zscores)$dist,
  robust::covRob(zscores, estim="pairwiseGK")$dist,
  bigsnpr:::getD(zscores)
) # 16 vs 41 ms

dim(zscores_rep <- zscores[rep(seq_len(nrow(zscores)), 100), ])
microbenchmark::microbenchmark(
  pcadapt:::covRob_rcpp(zscores_rep)$dist,
  robust::covRob(zscores_rep, estim="pairwiseGK")$dist,
  bigsnpr:::getD(zscores_rep),
  times = 10
)

## CLEAN
unlink(list.files("inst/testdata", "to-compare.pcadapt.", full.names = TRUE))
