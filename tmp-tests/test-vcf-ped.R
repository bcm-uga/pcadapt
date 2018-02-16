# bed
library(pcadapt)
file <- "inst/extdata/geno3pops"
bed0 <- read.pcadapt(paste0(file, ".bed"), type = "bed")
mat0 <- bed2matrix(unclass(bed0))
p0 <- pcadapt(bed0)



plink <- bigsnpr::download_plink()

tmp <- tempfile()


# ped
system(glue::glue(
  "{plink} --bfile {file} --recode --out {tmp}"
))

bed1 <- read.pcadapt(paste0(tmp, ".ped"), type = "ped")
mat1 <- bed2matrix(unclass(bed1), 150, 1500)
all(rowSums(cbind(colSums(mat1 == mat0), colSums(mat1 == (2 - mat0))) == 150) > 0)

p1 <- pcadapt(bed1)
all.equal(p1$maf, p0$maf)
plot(p1$loadings, p0$loadings)
plot(p1$pvalues, p0$pvalues)


# vcf
system(glue::glue(
  "{plink} --bfile {file} --recode vcf --out {tmp}"
))

bed2 <- read.pcadapt(paste0(tmp, ".vcf"), type = "vcf")
mat2 <- bed2matrix(unclass(bed2), 150, 1500)
all(rowSums(cbind(colSums(mat2 == mat0), colSums(mat2 == (2 - mat0))) == 150) > 0)

p2 <- pcadapt(bed2)
all.equal(p2$maf, p2$maf)
plot(p2$loadings, p0$loadings)
plot(p2$pvalues, p0$pvalues)
ind2 <- order(abs(p2$pvalues - p0$pvalues), decreasing = TRUE)[1:22]
p2$af[ind2]
