#With 2 populations only, pcadapt Pooled-seq does not work
filename<-filename[1:2,]
file2<- read.pcadapt(filename, type = "pool")

res2 <- pcadapt(file2)

#M Blum solution we can implement in the package
aux<-scale(file2,center=TRUE,scale=FALSE)
aux<-t(aux)
aux<-aux[,-1]

#Dos not work with n=2 pops
##cc<-robust::covRob(aux)
#Compute pval and stat
cc<-NULL
cc$dist<-aux^2
gif<-median(cc$dist)/qchisq(0.5,df=1)
stat2<-cc$dist/gif
pval<-pchisq(stat2,lower.tail=F,df=1)
hist(pval,breaks=20)

#List of outliers
library(qvalue)
qval <- qvalue(pval)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)
length(outliers)
plot(-log10(pval))
