# Install pcadapt v3.0.4
install.packages("https://cran.r-project.org/src/contrib/Archive/pcadapt/pcadapt_3.0.4.tar.gz",
                 repos = NULL, type = "source")

# Download dataset from https://datadryad.org/stash/dataset/doi:10.5061/dryad.266k4
unzip("tmp-data/doi_10.5061_dryad.266k4__v1.zip", exdir = "tmp-data")

# Conversion and quality control using PLINK
plink <- bigsnpr::download_plink("tmp-data")
prefix <- "tmp-data/cornell_canine"
system(glue::glue(
  "{plink} --bfile {prefix} --dog",
  " --maf 0.05 --mind 0.1 --geno 0.1 --hwe 1e-10 --autosome",
  " --recode A --out {prefix}"
))
