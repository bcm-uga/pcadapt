context("CONVERT")


N <- 150; P <- 400
mat <- matrix(0, N, P)
for (j in 1:P) mat[, j] <- rbinom(N, size = 2, prob = runif(1))
mat[sample(length(mat), size = length(mat) / 5)] <- 9  ## NAs
pcadapt:::write.table2(mat, tmpfile <- tempfile(), sep = " ")
# readLines(tmpfile)

obj <- mmapcharr::mmapchar(tmpfile, code = mmapcharr:::CODE_012)
bedfile <- writeBed(obj, is.pcadapt = FALSE)
expect_error(writeBed(obj, is.pcadapt = FALSE), 
             "The bed file already exists!", fixed = TRUE)


xptr <- pcadapt:::bedXPtr(bedfile, N, P)
pcadapt:::get_af()