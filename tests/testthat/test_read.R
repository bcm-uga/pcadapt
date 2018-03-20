################################################################################

context("READ_PCADAPT")

# Ped and Vcf are not maintained anymore

################################################################################

# Bed
bed <- system.file("extdata", "geno3pops.bed", package = "pcadapt")
bed2 <- read.pcadapt(bed, type = "bed")
expect_equal(unclass(bed2), normalizePath(bed), check.attributes = FALSE)
expect_equal(attr(bed2, "n"), 150)
expect_equal(attr(bed2, "p"), 1500)
expect_s3_class(bed2, "pcadapt_bed")

################################################################################

lfmm <- system.file("extdata", "geno3pops.lfmm", package = "pcadapt")
input <- file.copy(lfmm, tmp <- tempfile(fileext = ".lfmm"))

################################################################################


################################################################################


################################################################################


################################################################################