#Test Keurcien version based on loadings for pool-seq data

require(pcadapt)
pool.data <- system.file("extdata", "pool3pops", package = "pcadapt")
filename <- read.pcadapt(pool.data, type = "pool")

#Keurcien version with K=2
res2 <- pcadapt(filename,K=2)
plot(res2)
plot(res2,option="scores")
plot(res2,option="screeplot")

#Michael version that does not use PCA
aux<-scale(filename,center=TRUE,scale=FALSE)
aux<-t(aux)
aux<-aux[,-1]

cc<-robust::covRob(aux)
gif<-median(cc$dist)/qchisq(0.5,df=2)
stat2<-cc$dist/gif

#Comparison between Keurcien and M Blum version
cat("Correlation between PCA (K=2) and frequency version (no PCA) of pooled-seq stat: ",cor(stat2,res2$chi2.stat),"\n")
plot(stat2,res2$chi2.stat,pch=19,cex=.5);abline(0,1)

#Keurcien version with K=1
res1 <- pcadapt(filename,K=1)
plot(res1)
plot(res1,option="scores")
plot(res1,option="screeplot")

#Comparison between K=1 and K=2
cat("Correlation between PCA (K=2) and frequency version (no PCA) of pooled-seq stat: ",cor(res1$chi2.stat,res2$chi2.stat),"\n")
plot(stat2,res$chi2.stat,pch=19,cex=.5);abline(0,1)

#Test of get.pc
get.pc(res2,1:150)->aux
plot(aux)

#Do not provide K
res0<-pcadapt(filename)
all.equal(res0,res2)