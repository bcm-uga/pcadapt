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
  " --recode --out {prefix}_qc"
))

library(pcadapt)
input <- read.pcadapt(paste0(prefix, "_qc.ped"), type = "ped")
(nlines <- bigreadr::nlines(input)) # 144474
file.size(input) / 1024^2  # 1197 MB
file.size(input) / nlines # 8685 (~ 4342 * 2)

obj.pcadapt <- runonce::save_run({
  pcadapt(input, K = 20)
}, file = "tmp-data/canine-pcadapt-old.rds") # 2111 sec

scree.plotting(obj.pcadapt, K = 20)

str(obj.pcadapt)
plot(obj.pcadapt) +
  geom_hline(yintercept = -log10(0.05 / length(obj.pcadapt$pvalues)), linetype = 2)

score.plotting(obj.pcadapt, 19, 20)

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
file.size(input) / 1024^2  # 150 MB (~ 1197 / 8)

obj.pcadapt <- runonce::save_run({
  pcadapt(input, K = 20)
}, file = "tmp-data/canine-pcadapt-new.rds")  # 102 sec
system.time(pcadapt(input, K = 5))            #  35 sec
system.time(pcadapt(input, K = 10))           #  60 sec

scree_plot(obj.pcadapt, K = 20)
score_plot(obj.pcadapt)

obj.pcadapt.old <- readRDS("tmp-data/canine-pcadapt-old.rds")
diag(round(100 * cor(obj.pcadapt$scores, obj.pcadapt.old$scores), 2))
#  [1]  100 -100  100 -100  100 -100  100  100  100  100 -100  100  100
# [14]  100 -100 -100 -100 -100  100  100

################################################################################

library(ggplot2)
ggplot(data.frame(K = c(5, 10, 20), t = c(35, 60, 102) / 60), aes(K, t)) +
  bigstatsr::theme_bigstatsr() +
  geom_point(size = 3) +
  geom_line(size = 1.5) +
  geom_hline(yintercept = 2111 / 60, linetype = 2, color = "blue", size = 1.5) +
  labs(x = "Number of PCs used (K)", y = "Computation time (in minutes)") +
  scale_y_continuous(breaks = seq(0, 60, by = 10),
                     minor_breaks = seq(0, 60, by = 2))

ggsave("timings.pdf", width = 6, height = 5)
