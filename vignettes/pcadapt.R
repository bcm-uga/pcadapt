## ----eval=FALSE----------------------------------------------------------
#  install.packages("pcadapt")
#  library(pcadapt)

## ----echo=TRUE,include=FALSE---------------------------------------------
library(pcadapt)#Comment

## ---- results="hide"-----------------------------------------------------
path_to_file <- system.file("extdata","geno3pops",package="pcadapt")
filename <- read.pcadapt(path_to_file,type="lfmm")

## ----  results="hide"----------------------------------------------------
pool.data <- read.table(system.file("extdata","pool3pops",package="pcadapt"))
filename <- read.pcadapt(pool.data,type="pool",local.env = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  x <- pcadapt(filename,K=20)

## ----echo=FALSE----------------------------------------------------------
output.filename <- paste0(path_to_file,20)
K <- 20
method = "mahalanobis"
data.type <- "genotype"
min.maf <- 0.05
x <- create.pcadapt(output.filename,K,method,data.type,min.maf)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="screeplot")

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="screeplot",K=10)

## ------------------------------------------------------------------------
# With integers
poplist.int <- c(rep(1,50),rep(2,50),rep(3,50))
# With names
poplist.names <- c(rep("POP1",50),rep("POP2",50),rep("POP3",50))
print(poplist.int)
print(poplist.names)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="scores",pop=poplist.int)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="scores",pop=poplist.names)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="scores",i=3,j=4,pop=poplist.names)

## ----eval=FALSE----------------------------------------------------------
#  x <- pcadapt(filename,K=2)

## ----echo=FALSE----------------------------------------------------------
output.filename <- path_to_file
K <- 2
method = "mahalanobis"
data.type <- "genotype"
min.maf <- 0.05
x <- create.pcadapt(output.filename,K,method,data.type,min.maf)

## ----eval=FALSE----------------------------------------------------------
#  summary(x)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="manhattan")

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="qqplot",threshold=0.1)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
hist(x$pvalues,xlab="p-values",main=NULL,breaks=50)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="stat.distribution")

## ----echo=FALSE----------------------------------------------------------
output.filename <- path_to_file
K <- 2
method = "communality"
data.type <- "genotype"
min.maf <- 0.05
x_com <- create.pcadapt(output.filename,K,method,data.type,min.maf)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x_com,option="stat.distribution")

## ----echo=FALSE----------------------------------------------------------
output.filename <- path_to_file
K <- 2
method = "componentwise"
data.type <- "genotype"
min.maf <- 0.05
x_cw <- create.pcadapt(output.filename,K,method,data.type,min.maf)
summary(x_cw$pvalues)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x_cw,option="stat.distribution",K=2)

