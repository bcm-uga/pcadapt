################################################################################

sapply(1:5, function(i) {
  cor(readRDS(paste0("tmp-data/missing-res-new-", i, ".rds"))$chi2.stat,
      readRDS(paste0("tmp-data/missing-res-old-", i, ".rds"))$chi2.stat,
      use = "pairwise.complete.obs")
})
# 0.9975642 0.9975758 0.9976246 0.9973503 0.9975562

################################################################################

library(ggplot2)

qplot(readRDS(paste0("tmp-data/missing-res-old-1.rds"))$chi2.stat,
      readRDS(paste0("tmp-data/missing-res-new-1.rds"))$chi2.stat) +
  bigstatsr::theme_bigstatsr() +
  geom_abline(color = "red") +
  labs(x = "Statistic (old version)", y = "Statistic (new version)")

ggsave("compare-with-old.png", width = 8, height = 6)

qplot(readRDS(paste0("tmp-data/canine-pcadapt-new.rds"))$chi2.stat,
      readRDS(paste0("tmp-data/missing-res-new-1.rds"))$chi2.stat) +
  bigstatsr::theme_bigstatsr() +
  geom_abline(color = "red") +
  labs(x = "Statistic (no missing values)", y = "Statistic (many missing values)")

ggsave("compare-with-missing.png", width = 8, height = 6)

################################################################################
