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
