#POPRES example
#Generate a file with 1 percent of missing values
file.pcadapt <- "/Users/mblum/home/mblum/Courant/Data/POPRES/POPRES_allchr.pcadapt"
# file_POPRES<-pcadapt::read.pcadapt(,type="pcadapt")
tmp <- mmapcharr::mmapchar(file.pcadapt, code = mmapcharr:::CODE_012)
# devtools::install_github("privefl/bigsnpr")
require(bigsnpr)
tmp2 <- snp_fake(n = ncol(tmp), m = nrow(tmp))
str(tmp2)
options(bigstatsr.check.args = FALSE)

big_apply(tmp2$genotypes, function(X, ind) {
  print(min(ind))
  X[, ind] <- t(tmp[ind, ])
  NULL
})

G <- tmp2$genotypes
ind <- sort(sample(length(G), 3*length(G) / 10))
pryr::object_size(ind)
G[ind] <- as.raw(3)
counts <- big_counts(G)

bed <- snp_writeBed(tmp2, "/Users/mblum/home/mblum/Courant/Data/POPRES/POPRES_allchr_NA30p100.bed")

file_POPRES<-pcadapt::read.pcadapt("/Users/mblum/home/mblum/Courant/Data/POPRES/POPRES_allchr.pcadapt",type="pcadapt")
y <- pcadapt::pcadapt(file_POPRES,K=2)

file_POPRES2<-pcadapt::read.pcadapt("/Users/mblum/home/mblum/Courant/Data/POPRES/POPRES_allchr_NA30p100.bed",type="bed")
y2 <- pcadapt::pcadapt(file_POPRES2,K=2)


plot(-log10(y$pvalues),-log10(y2$pvalues),pch=19,cex=.5)

yc<-pcadapt::pcadapt(file_POPRES,K=5,LD.clumping=list(size = 200, thr = 0.1))
plot(yc)
yc2 <- pcadapt::pcadapt(file_POPRES2,K=10,LD.clumping=list(size = 200, thr = 0.1))

plot(-log10(yc$pvalues),-log10(yc2$pvalues),pch=19,cex=.5);abline(0,1)
plot(counts[4, ],-log10(yc2$pvalues),pch=19,cex=.5)

keep<-(-log10(yc$pvalues)>10)
plot(yc2$maf[keep],-log10(yc2$pvalues)[keep],pch=19,cex=.5)


qval <- qvalue(yc$pvalues)$qvalues
plot(-log10(yc$pvalues),col=(qval < alpha)+1,ylim=c(0,10),pch=19,cex=.5)
length(yc$pvalues)

names<-read.table("/Users/mblum/home/mblum/Courant/Data/POPRES/POPRES_allchr.map")
