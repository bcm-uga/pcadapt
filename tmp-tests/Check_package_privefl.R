#Functions to check privefl update of pcadapt

##First, we check that results are the same when comparing with CRAN version 3.0.4 of pcadapt

#3 pop example of the R vignette
install.packages("pcadapt")
path_to_file <- system.file("extdata","geno3pops",package="pcadapt")
filename <- pcadapt::read.pcadapt(path_to_file,type="lfmm")
x <- pcadapt::pcadapt(filename,K=2)
saveRDS(object = x$pvalues,file="3pops_pvalues.rds")

#POPRES example
file_POPRES<-pcadapt::read.pcadapt("/Users/mblum/home/mblum/Courant/Data/POPRES/POPRES_allchr.pcadapt",type="pcadapt")
y <- pcadapt::pcadapt(file_POPRES,K=2)
saveRDS(object = y$pvalues,file="POPRES_pvalues.rds")

#Case 1, 2 of SSMPG 2015
install.packages("pcadapt")
filename<-pcadapt::read.pcadapt("tmp-tests/Case1.lfmm",type="lfmm")
z1<-pcadapt::pcadapt(filename,K=18)
plot(-log10(z1$pvalues),pch=19)
saveRDS(object = z1$pvalues,file="Case1_pvalues.rds")

filename<-pcadapt::read.pcadapt("tmp-tests/Case2.lfmm",type="lfmm")
groundtruth<-readRDS("tmp-tests/Case2.rds")$ground.truth
z2<-pcadapt::pcadapt(filename,K=3)
p<-length(z2$pvalues)
plot(1:p,-log10(z2$pvalues),pch=19,cex=0.5);points((1:p)[groundtruth],-log10(z2$pvalues)[groundtruth],pch=19,col=2,cex=1)
saveRDS(object = z2$pvalues,file="Case2_pvalues.rds")


##Run Florian package and compare
devtools::install_github("privefl/pcadapt")


filename <- pcadapt::read.pcadapt("inst/extdata/geno3pops.pcadapt",type="pcadapt")
x<-pcadapt::pcadapt(filename,K=2)
cat("3 pops: Correlation between CRAN 3.0.4 version and privefl version (pcadapt) is equal to",cor(-log10(readRDS("3pops_pvalues.rds")),-log10(x$pvalues)))


filename<-pcadapt::read.pcadapt("tmp-tests/Case1.lfmm",type="lfmm")
x<-pcadapt::pcadapt(filename,K=18)
cat("Case1: Correlation between CRAN 3.0.4 version and privefl version (lfmm) is equal to",cor(-log10(readRDS("Case1_pvalues.rds")),-log10(x$pvalues)))

filename<-pcadapt::read.pcadapt("tmp-tests/Case1.bed",type="bed")
x<-pcadapt::pcadapt(filename,K=18)
cat("Case1: Correlation between CRAN 3.0.4 version and privefl version (bed) is equal to",cor(-log10(readRDS("Case1_pvalues.rds")),-log10(x$pvalues)))


filename<-pcadapt::read.pcadapt("tmp-tests/Case2.lfmm",type="lfmm")
x<-pcadapt::pcadapt(filename,K=3)
cat("Case2: Correlation between CRAN 3.0.4 version and privefl version (lfmm) is equal to",cor(-log10(readRDS("Case2_pvalues.rds")),-log10(x$pvalues),use="pairwise.complete.obs"))

filename<-pcadapt::read.pcadapt("tmp-tests/Case2.bed",type="bed")
x<-pcadapt::pcadapt(filename,K=3)
cat("Case2: Correlation between CRAN 3.0.4 version and privefl version (bed) is equal to",cor(-log10(readRDS("Case2_pvalues.rds")),-log10(x$pvalues),use="pairwise.complete.obs"))

file_POPRES<-pcadapt::read.pcadapt("/Users/mblum/home/mblum/Courant/Data/POPRES/POPRES_allchr.pcadapt",type="pcadapt")
x <- pcadapt::pcadapt(file_POPRES,K=2)
cat("POPRES: Correlation between CRAN 3.0.4 version and privefl version (pcadapt) is equal to",cor(-log10(readRDS("POPRES_pvalues.rds")),-log10(x$pvalues),use="pairwise.complete.obs"))
