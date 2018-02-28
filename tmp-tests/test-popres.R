tmp <- "../thesis/POPRES_data/POPRES_allchr.bed"

library(pcadapt)
tmp2 <- read.pcadapt(tmp, "bed")
system.time(
  tmp3 <- pcadapt(tmp2, K = 10, LD.clumping = list(size = 200, thr = 0.1))
)

plot(tmp3)

sum(tmp3$maf > 0.05)
sum(is.na(tmp3$loadings[, 1]))


plot(tmp3$loadings[, 7])
pcadapt:::plot.pcadapt

library(bigsnpr)
popres <- snp_attach(snp_readBed(tmp, tempfile())) 
G <- popres$genotypes
CHR <- popres$map$chromosome
POS <- popres$map$physical.pos
maf <- snp_MAF(G)
NCORES <- nb_cores()
obj.svd <- snp_autoSVD(G, CHR, POS, ind.col = which(maf > 0.05), ncores = NCORES)

obj.gwas <- snp_pcadapt(G, obj.svd$u, ind.col = which(maf > 0.05))
obj.gwas2 <- snp_pcadapt(G, obj.svd$u[, 1:5], ind.col = which(maf > 0.05))

plot(predict(snp_gc(obj.gwas)), predict(snp_gc(obj.gwas2)), pch = 20, cex = .5)
abline(0, 1, col = "red")

qval <- qvalue::qvalue(predict(snp_gc(obj.gwas), log10 = FALSE))$qvalues
ind <- which(maf > 0.05)[qval < 0.1]
CHR[ind]
tmp <- data.table::fread("../thesis/POPRES_data/POPRES_Snps_QC2.txt",
                         data.table = FALSE, header = FALSE)

library(dplyr)
pval <- predict(snp_gc(obj.gwas3), log10 = FALSE)
ind <- which(maf > 0.05)[pval < 5e-8]
data.frame(ID = tmp$V2[match(popres$map$marker.ID[ind], tmp$V1)],
           pvalue = pval[pval < 5e-8],
           chr = CHR[ind]) %>%
  arrange(pvalue)


obj.svd2 <- snp_autoSVD(G, CHR, POS, ind.col = which(maf > 0.05), 
                        thr.r2 = 0.1, ncores = NCORES)
plot(obj.svd2)
obj.gwas3 <- snp_pcadapt(G, obj.svd2$u[, 1:5], ind.col = which(maf > 0.05))
plot(predict(snp_gc(obj.gwas3)), predict(snp_gc(obj.gwas2)), pch = 20, cex = .5)
abline(0, 1, col = "red")

plot(obj.svd2, type = "loadings", loadings = 1:10)
plot(snp_gc(obj.gwas3)) +
  ggplot2::geom_hline()