input <- readRDS("tmp-data/poolData_Robject_14_pools_AngelaFuentes.rds")
dim(input)
class(input)

library(pcadapt)
system.time(
  obj.pcadapt <- pcadapt(input, K = 11, LD.clumping = list(size = 10000, thr = 0.1))
)

str(obj.pcadapt)
plot(obj.pcadapt, option = "stat.distribution")
hist(obj.pcadapt$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(obj.pcadapt$loadings[, 1], pch = 20)

S2 <- bigutilsr::rollmean(sqrt(obj.pcadapt$stat[obj.pcadapt$pass]), 100)
S2.thr <- bigutilsr::tukey_mc_up(S2, alpha = 0.05)
hist(S2); abline(v = S2.thr, col = "red")
ind.col.excl <- obj.pcadapt$pass[S2 > S2.thr]
input[, ind.col.excl] <- NA

system.time(
  obj.pcadapt2 <- pcadapt(input, K = 11, LD.clumping = list(size = 10000, thr = 0.1))
)
hist(obj.pcadapt2$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(obj.pcadapt2$loadings[, 1], pch = 20)
