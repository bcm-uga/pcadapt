#2 pools
pool.data <- read.table(system.file("extdata", "pool3pops", package = "pcadapt"))
filename <- read.pcadapt(pool.data, type = "pool")
filename<-filename[1:2,]

file2<- read.pcadapt(filename, type = "pool")

#M Blum solution we can implement in the package
aux<-scale(file2,center=TRUE,scale=FALSE)
aux<-t(aux)
aux<-aux[,-1]

#Dos not work with n=2 pops
##cc<-robust::covRob(aux)
#Compute pval and stat
cc<-NULL
cc$dist<-(aux-median(aux))^2
gif<-median(cc$dist)/qchisq(0.5,df=1)
stat2<-cc$dist/gif
pval<-pchisq(stat2,lower.tail=F,df=1)
hist(pval,breaks=20)

#Comparison betwxwen Blum solution and Keurcien's solution

res2 <- pcadapt(file2)
plot(stat2,res2$chi2.stat);abline(0,1)
plot(-log10(pval),-log10(res2$pvalues),pch=19);abline(0,1)


#List of outliers
library(qvalue)
qval <- qvalue(res2$pval)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)
length(outliers)
plot(-log10(pval))
