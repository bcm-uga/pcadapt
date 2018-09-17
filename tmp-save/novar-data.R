library(bigsnpr)

tmp <- snp_fake(100, 1000)
G <- tmp$genotypes
G[] <- sample(as.raw(0:3), length(G), TRUE)
G[, 1] <- 0
G[, 2] <- 1
G[, 3] <- 2
snp_writeBed(tmp, "inst/testdata/novar.bed")

library(pcadapt)

bed <- read.pcadapt("inst/testdata/novar.bed", type = "bed")
# debugonce(pcadapt)
obj.pcadapt <- pcadapt(bed, K = 5)
plot(obj.pcadapt)
head(obj.pcadapt$pvalues)


M <- read.pcadapt(G[])
class(M)
obj.pcadapt <- pcadapt(M, K = 5)
M[1, ] <- NA
obj.pcadapt <- pcadapt(M, K = 5)
