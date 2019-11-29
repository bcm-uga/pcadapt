# Install pcadapt v3.0.4 -- *old* version
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
  " --recode --out {prefix}"
))

library(pcadapt)
input <- read.pcadapt(paste0(prefix, ".ped"), type = "ped")
obj.pcadapt <- runonce::save_run({
  pcadapt(input, K = 20)
}, file = "tmp-data/canine-pcadapt-old.rds") # 2111 sec

pcadapt::scree.plotting(obj.pcadapt, K = 20)

str(obj.pcadapt)
plot(obj.pcadapt) +
  geom_hline(yintercept = -log10(0.05 / length(obj.pcadapt$pvalues)), linetype = 2)

pcadapt::score.plotting(obj.pcadapt, 19, 20)

################################################################################

# Install pcadapt v4.1.0 -- *new* version
install.packages("pcadapt")

# Conversion and quality control using PLINK
plink <- bigsnpr::download_plink("tmp-data")
prefix <- "tmp-data/cornell_canine"
system(glue::glue(
  "{plink} --bfile {prefix} --dog",
  " --maf 0.05 --mind 0.1 --geno 0.1 --hwe 1e-10 --autosome",
  " --make-bed --out {prefix}_qc"
))

library(pcadapt)
input <- read.pcadapt(paste0(prefix, "_qc.bed"), type = "bed")
obj.pcadapt <- runonce::save_run({
  pcadapt(input, K = 20)
}, file = "tmp-data/canine-pcadapt-new.rds")  # 102 sec
system.time(pcadapt(input, K = 5))            #  35 sec
system.time(pcadapt(input, K = 10))           #  60 sec

pcadapt::scree_plot(obj.pcadapt, K = 20)
pcadapt::score_plot(obj.pcadapt)

obj.pcadapt.old <- readRDS("tmp-data/canine-pcadapt-old.rds")
round(100 * cor(obj.pcadapt$scores, obj.pcadapt.old$scores), 1)
