################################################################################

context("READ_PCADAPT")

# Ped and Vcf are not maintained any more

################################################################################

# Bed
bed <- system.file("extdata", "geno3pops.bed", package = "pcadapt")
bed2 <- read.pcadapt(bed, type = "bed")
expect_equal(unclass(bed2), bed, check.attributes = FALSE)
expect_equal(attr(bed2, "n"), 517)
expect_equal(attr(bed2, "p"), 4542)
expect_s3_class(bed2, "pcadapt_bed")

################################################################################

lfmm <- system.file("extdata", "geno3pops.lfmm", package = "pcadapt")
input <- file.copy(lfmm, tmp <- tempfile(fileext = ".lfmm"))

################################################################################

lfmm <- system.file("extdata", "geno3pops.lfmm", package = "pcadapt")

################################################################################

lfmm <- system.file("extdata", "geno3pops.lfmm", package = "pcadapt")

################################################################################


################################################################################


################################################################################

################################################################################