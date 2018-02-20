tmp <- "../thesis/POPRES_data/POPRES_allchr.bed"

library(pcadapt)
tmp2 <- read.pcadapt(tmp, "bed")
system.time(
  tmp3 <- pcadapt(tmp2, K = 2)
)

plot(tmp3)
