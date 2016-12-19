## ----eval=FALSE----------------------------------------------------------
#  install.packages("pcadapt")
#  library(pcadapt)

## ----echo=TRUE,include=FALSE---------------------------------------------
library(pcadapt)

## ----eval=FALSE----------------------------------------------------------
#  install.packages("pcadapt")
#  library(pcadapt)

## ----echo=TRUE,include=FALSE---------------------------------------------
library(pcadapt)#Comment

## ----eval=FALSE----------------------------------------------------------
#  pool.data <- read.table("path_to_directory/foo")
#  filename <- read.pcadapt(pool.data,type="pool",local.env = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  pool.data <- read.table(system.file("extdata","pool3pops",package="pcadapt"))
#  filename <- read.pcadapt(pool.data,type="pool",local.env = TRUE)

## ---- include=FALSE------------------------------------------------------
pool.data <- read.table(system.file("extdata","pool3pops",package="pcadapt"))
path_to_file <- system.file("extdata","pool3pops",package="pcadapt")
filename <- read.pcadapt(pool.data,type="pool",local.env = TRUE,pop.sizes = c(200,200,200))

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

