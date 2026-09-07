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

# Matrix hard-call validation
valid <- matrix(c(0, 1, 2, NA), nrow = 2)
valid.input <- read.pcadapt(valid, type = "pcadapt")
expect_s3_class(valid.input, "pcadapt_matrix")
expect_type(valid.input, "integer")
expect_equal(unclass(valid.input), t(valid), check.attributes = FALSE)

expect_error(
  read.pcadapt(matrix(c(0, 1.9, 2, NA), nrow = 2), type = "pcadapt"),
  "hard calls coded as 0, 1, or 2",
  fixed = TRUE
)
expect_error(
  read.pcadapt(matrix(c(0, 1, 3, NA), nrow = 2), type = "pcadapt"),
  "hard calls coded as 0, 1, or 2",
  fixed = TRUE
)
expect_error(
  read.pcadapt(matrix(c(0, 1, Inf, NA), nrow = 2), type = "pcadapt"),
  "hard calls coded as 0, 1, or 2",
  fixed = TRUE
)
expect_error(
  read.pcadapt(matrix(c(0, 1, NaN, NA), nrow = 2), type = "pcadapt"),
  "NaN is not a valid missing genotype",
  fixed = TRUE
)
expect_error(
  read.pcadapt(matrix(c("0", "1", "2", NA), nrow = 2), type = "pcadapt"),
  "numeric hard calls",
  fixed = TRUE
)

################################################################################
