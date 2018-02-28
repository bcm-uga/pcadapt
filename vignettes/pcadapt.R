## ---- eval = FALSE-------------------------------------------------------
#  install.packages("pcadapt")
#  library(pcadapt)

## ---- include = FALSE----------------------------------------------------
library(pcadapt)#Comment

## ------------------------------------------------------------------------
path_to_file <- system.file("extdata", "geno3pops.lfmm", package = "pcadapt")
filename <- read.pcadapt(path_to_file, type = "lfmm")

## ------------------------------------------------------------------------
x <- pcadapt(input = filename, K = 20) 

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
plot(x, option = "screeplot")

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
plot(x, option = "screeplot", K = 10)

## ------------------------------------------------------------------------
# With integers
poplist.int <- c(rep(1, 50), rep(2, 50), rep(3, 50))
# With names
poplist.names <- c(rep("POP1", 50),rep("POP2", 50),rep("POP3", 50))
print(poplist.int)
print(poplist.names)

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
plot(x, option = "scores", pop = poplist.int)

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
plot(x, option = "scores", pop = poplist.names)

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)

## ------------------------------------------------------------------------
x <- pcadapt(filename, K = 2)

## ----eval=FALSE----------------------------------------------------------
#  summary(x)

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
plot(x , option = "manhattan")

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
plot(x, option = "qqplot")

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
plot(x, option = "stat.distribution")

## ------------------------------------------------------------------------
path_to_file <- system.file("extdata", "SSMPG2017.rds", package = "pcadapt")
genotypes<-readRDS(path_to_file)
matrix<- read.pcadapt(genotypes, type = "pcadapt")
res<-pcadapt(matrix,K=20)
plot(res,option="screeplot")

## ------------------------------------------------------------------------
res<-pcadapt(matrix,K=4)
plot(res)

## ------------------------------------------------------------------------
par(mfrow = c(2, 2))
for (i in 1:4)
  plot(res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

## ------------------------------------------------------------------------
res <- pcadapt(matrix, K = 20, LD.clumping = list(size = 200, thr = 0.1))
plot(res, option = "screeplot")

## ------------------------------------------------------------------------
res <- pcadapt(matrix, K = 2, LD.clumping = list(size = 200, thr = 0.1))
par(mfrow = c(1, 2))
for (i in 1:2)
  plot(res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

## ------------------------------------------------------------------------
plot(res)

## ----  eval = TRUE-------------------------------------------------------
pool.data <- system.file("extdata", "pool3pops", package = "pcadapt")
filename <- read.pcadapt(pool.data, type = "pool")

## ----  eval = TRUE-------------------------------------------------------
res<-pcadapt(filename)
summary(res)

## ----  eval = TRUE-------------------------------------------------------
plot(-log10(res$pvalues),pch=19,cex=.5)
padj <- p.adjust(res$pvalues,method="BH")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)

## ---- eval=FALSE---------------------------------------------------------
#  path_to_file <- system.file("extdata", "SSMPG2017.rds", package = "pcadapt")
#  genotypes<-readRDS(path_to_file)
#  print(dim(genotypes))
#  matrix<- read.pcadapt(genotypes, type = "pcadapt")
#  res<-pcadapt(matrix,K=20)
#  plot(res,option="screeplot")

## ------------------------------------------------------------------------
path_to_file <- system.file("extdata", "geno3pops.lfmm", package = "pcadapt")
filename <- read.pcadapt(path_to_file, type = "lfmm")
x_cw <- pcadapt(filename, K = 2, method = "componentwise")
summary(x_cw$pvalues)

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
plot(x_cw, option = "stat.distribution", K = 2)

