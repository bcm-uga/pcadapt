tmp <- "../thesis/POPRES_data/POPRES_allchr.bed"

library(pcadapt)
tmp2 <- read.pcadapt(tmp, "bed")
system.time(
  tmp3 <- pcadapt(tmp2, K = 10, LD.clumping = list(size = 200, thr = 0.1))
)

plot(tmp3)

sum(tmp3$maf > 0.05)
sum(is.na(tmp3$loadings[, 1]))


plot(tmp3$loadings[, 7])
pcadapt:::plot.pcadapt
