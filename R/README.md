# Introgression

# Example with populus: (https://www.dropbox.com/sh/85cgelc1af1s7az/AABvFT5W7NrDd_cXQHTff2Pma?dl=0)

filename <- "populus3pops.pcadapt" 
popfile <- "populus3pops.pop" 
lab <- read.table(popfile)[, 1] 

geno <- impute.pcadapt(filename, ploidy = 2)
stat <- scan.intro(geno, K = 2, lab = lab, min.maf = 0.05, ploidy = 2, window.size = 15000, direction = 0, ancstrl.1 = 1, ancstrl.2 = 3, admxd = 4)
