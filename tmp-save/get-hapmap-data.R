# https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.qc.poly.tar.bz2

plink <- bigsnpr::download_plink()
dir <- "tmp-data/hapmap3_pop"

# Merge all files
write(rev(list.files(dir, full.names = TRUE)), 
      file = (tmp <- tempfile()),
      ncolumns = 2)

system(glue::glue(
  "{plink} --merge-list {tmp} --make-bed --out {dir}/hapmap_all"
))

system(sprintf(
  "awk '{ if ($1 == 8) print $2 }' %s/hapmap_all.bim > %s/hapmap_all_chr8.txt",
  dir, dir
))

dir.create("hapmap_data")
system(glue::glue(
  "{plink} --bfile {dir}/hapmap_all --make-bed --out hapmap_data/hapmap_chr8",
  " --maf 0.01",
  " --geno 0.2",
  " --mind 0.2",
  " --hwe 1e-50",
  " --extract {dir}/hapmap_all_chr8.txt"
))


library(pcadapt)

bed <- read.pcadapt("hapmap_data/hapmap_chr8.bed", type = "bed")

tmp <- pcadapt(bed, K = 10)
plot(tmp, option = "screeplot")
tmp2 <- pcadapt(bed, K = 2)
plot(tmp2)
plot(tmp2, option = "scores")

tmp3 <- pcadapt(bed, K = 6)
plot(tmp3)
plot(tmp3, option = "scores", i = 3, j = 4)
plot(tmp3$loadings[, 6])

tmp4 <- pcadapt(bed, K = 10, LD.clumping = list(size = 200, thr = 0.1))
plot(tmp4, option = "screeplot")
tmp5 <- pcadapt(bed, K = 2, LD.clumping = list(size = 200, thr = 0.1))
plot(tmp5)
plot(tmp5, option = "qqplot")
hist(tmp5$pvalues)
plot(tmp5, option = "scores", i = 5, j = 4)
plot(tmp4$loadings[, 6])


bed <- read.pcadapt("tmp-data/hapmap3_pop/hapmap_all.bed", type = "bed")
tmp <- pcadapt(bed, K = 10, LD.clumping = list(size = 200, thr = 0.1))
plot(tmp, option = "screeplot")
plot(tmp)

tmp2 <- pcadapt(bed, K = 6, LD.clumping = list(size = 200, thr = 0.1))
plot(tmp2)
