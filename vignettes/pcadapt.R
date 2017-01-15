## ----eval=FALSE----------------------------------------------------------
#  install.packages("pcadapt")
#  library(pcadapt)

## ----echo=TRUE,include=FALSE---------------------------------------------
library(pcadapt)#Comment

## ---- eval = FALSE-------------------------------------------------------
#  path_to_file <- system.file("extdata", "geno3pops.lfmm", package = "pcadapt")
#  filename <- read.pcadapt(path_to_file, type = "lfmm")

## ---- echo = FALSE-------------------------------------------------------
filename <- system.file("extdata", "geno3pops.pcadapt", package = "pcadapt")

## ----  eval = FALSE------------------------------------------------------
#  pool.data <- read.table(system.file("extdata", "pool3pops", package = "pcadapt"))
#  filename <- read.pcadapt(pool.data, type = "pool")

## ------------------------------------------------------------------------
x <- pcadapt(input = filename, K = 20) #or x <- pcadapt(input = matrix, K = 20)

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
plot(x, option = "qqplot", threshold = 0.1)

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50)

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
plot(x, option = "stat.distribution")

## ------------------------------------------------------------------------
x_com <- pcadapt(filename, K = 2, method = "communality")

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
plot(x_com, option = "stat.distribution")

## ------------------------------------------------------------------------
x_cw <- pcadapt(filename, K = 2, method = "componentwise")
summary(x_cw$pvalues)

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'----------------
plot(x_cw, option = "stat.distribution", K = 2)

