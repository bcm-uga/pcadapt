# Conversion and quality control using PLINK
plink <- bigsnpr::download_plink("tmp-data")
prefix <- "tmp-data/cornell_canine"
system(glue::glue(
  "{plink} --bfile {prefix} --dog",
  " --maf 0.05 --mind 0.1 --geno 0.1 --hwe 1e-10 --autosome",
  " --make-bed --out {prefix}_qc"
))

library(pcadapt)
library(data.table)
bedfile <- paste0(prefix, "_qc.bed")
mat0 <- bed2matrix(bedfile)  # 4342 x 144474 -> 2.3 GB
n <- nrow(mat0)
m <- ncol(mat0)

for (i in 1:5) {

  data_file <- paste0("tmp-data/missing-", i, ".pcadapt")
  if (!file.exists(data_file)) {

    set.seed(i)

    nbNA <- VGAM::rbetabinom.ab(m, size = n, shape1 = 0.8, shape2 = 3)
    print(sum(nbNA) / (n * m)) # 21%
    print(range(nbNA))

    # Generate indices of missing values (as a two-column matrix)
    indNA <- cbind(
      unlist(lapply(nbNA, function(nb) {
        `if`(nb > 0, sample(n, size = nb), NULL)
      })),
      rep(1:m, nbNA)
    )
    # Fill a copy of the matrix with NAs (represented as 9s)
    matNA <- mat0; matNA[indNA] <- 9L
    # Write matrix as '.pcadapt' file
    fwrite(as.data.table(t(matNA)), file = data_file,
           col.names = FALSE, sep = " ", na = "9")
  }
}
