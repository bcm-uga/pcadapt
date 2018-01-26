lfmm <- "inst/extdata/geno3pops.lfmm"
file.copy(lfmm, tmp <- "test.lfmm")
library(pcadapt)

bed <- pcadapt::read.pcadapt(tmp, type = "lfmm")

str(test <- pcadapt(bed, LD.clumping = TRUE))
sum(is.na(test$loadings[, 1]))
sum(is.na(test$zscores))

plot(test)
plot(test, option = "scores")

