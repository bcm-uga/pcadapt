################################################################################

context("PLOT")

################################################################################

path_to_file <- system.file("extdata", "geno3pops.bed", package = "pcadapt")
filename <- read.pcadapt(path_to_file, type = "bed")
x <- pcadapt(input = filename, K = 20) 

expect_error(plot(x, option = "hist"), "'arg' should be one of")

expect_s3_class(plot(x, option = "screeplot"), "ggplot")

expect_s3_class(plot(x, option = "scores"), "ggplot")
expect_s3_class(plot(x, option = "scores", pop = rep(1:3, each = 50)), "ggplot")
expect_s3_class(plot(x, option = "scores", pop = rep(1:3, each = 50), 
                     col = c("black", "chartreuse3", "orange")),
                "ggplot")
expect_s3_class(plot(x, option = "scores", plt.pkg = "plotly", pop = rep(1:3, each = 50),
                     col = c("black", "chartreuse3", "orange")),
                "plotly")
expect_error(plot(x, option = "scores", plt.pkg = "ggplot2"), "should be either")

expect_s3_class(plot(x, option = "manhattan"), "ggplot")
expect_s3_class(plot(x, option = "manhattan", chr.info = rep(1:15, each = 100)), "ggplot")
expect_s3_class(plot(x, option = "manhattan", chr.info = rep(1:15, each = 100),
                     plt.pkg = "plotly"), "plotly")

expect_s3_class(plot(x, option = "stat.distribution"), "ggplot")

expect_s3_class(plot(x, option = "qqplot"), "ggplot")

################################################################################

# Componentwise plots use the selected component and one degree of freedom.
componentwise <- structure(
  list(
    chi2.stat = cbind(
      c(1, NA, 4, 9),
      c(16, 25, NA, 36)
    ),
    pvalues = cbind(
      stats::pchisq(c(1, NA, 4, 9), df = 1, lower.tail = FALSE),
      stats::pchisq(c(16, 25, NA, 36), df = 1, lower.tail = FALSE)
    ),
    maf = rep(0.25, 4)
  ),
  K = 2,
  method = "componentwise",
  min.maf = 0.05,
  class = "pcadapt"
)

expected <- -stats::pchisq(
  c(16, 25, 36),
  df = 1,
  lower.tail = FALSE,
  log.p = TRUE
) / log(10)

manhattan <- pcadapt:::manhattan_plot(
  componentwise,
  chr.info = NULL,
  snp.info = NULL,
  K = 2
)
expect_equal(manhattan$data$x, c(1L, 2L, 4L))
expect_equal(manhattan$data$y, expected)

qq <- pcadapt:::qq_plot(componentwise, K = 2)
expect_equal(sort(qq$data$y), sort(expected))

################################################################################
