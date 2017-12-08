lfmm <- "inst/extdata/geno3pops.lfmm"
lfmm.xptr <- mmapcharr::mmapchar(lfmm, code = mmapcharr:::CODE_012)

lfmm.xptr[, 1:5]
bedfile <- writeBed(lfmm.xptr)

library(bigsnpr)
test <- snp_attach(snp_readBed(bedfile, tempfile()))
test$genotypes[, 1:5]
all.equal(test$genotypes[], lfmm.xptr[])


pcadapt <- "inst/extdata/geno3pops.pcadapt"
pcadapt.xptr <- mmapcharr::mmapchar(pcadapt, code = mmapcharr:::CODE_012)

pcadapt.xptr[1:5, ]
bedfile2 <- writeBed(pcadapt.xptr)

test2 <- snp_attach(snp_readBed(bedfile2, tempfile()))
test2$genotypes[, 1:5]
all.equal(test$genotypes[], t(pcadapt.xptr[]))
