lfmm <- system.file("extdata", "geno3pops.lfmm", package = "pcadapt")
file.copy(lfmm, tmp <- tempfile(fileext = ".lfmm"))
library(pcadapt)

bed <- pcadapt::read.pcadapt(tmp, type = "lfmm")
# xptr <- pcadapt:::bedXPtr(bed, 150, 1500)

obj.pcadapt <- pcadapt::pcadapt(bed)

str(test <- pcadapt(bed, LD.clumping = TRUE))
sum(is.na(test$loadings[, 1]))
sum(is.na(test$zscores))

plot(test)
plot(test, option = "scores")

