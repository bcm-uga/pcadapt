## ----echo=FALSE,include=FALSE--------------------------------------------
library(pcadapt)

## ----eval=FALSE----------------------------------------------------------
#  install.packages("pcadapt")
#  library(pcadapt)

## ----echo=TRUE,include=FALSE---------------------------------------------
library(pcadapt)#Comment

## ------------------------------------------------------------------------
pooldata <- system.file("extdata","pool3pops",package="pcadapt")

## ------------------------------------------------------------------------
x.pool <- pcadapt(pooldata,data.type="pool")

## ------------------------------------------------------------------------
plot(x.pool,option="manhattan")

