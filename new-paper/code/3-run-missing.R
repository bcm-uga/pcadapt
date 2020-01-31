# Install pcadapt v3.0.4 -- *old* version
# install.packages("https://cran.r-project.org/src/contrib/Archive/pcadapt/pcadapt_3.0.4.tar.gz",
#                  repos = NULL, type = "source")

# Install pcadapt v4.1.0 -- *new* version
# install.packages("pcadapt")

version <- if (packageVersion("pcadapt") == "4.1.0") {
  "new"
} else if (packageVersion("pcadapt") == "3.0.4") {
  "old"
} else {
  stop("Wrong version.")
}
cat(version)

library(pcadapt)

for (i in 1:5) {

  data_file <- paste0("tmp-data/missing-", i, ".pcadapt")
  res_file <- sub(
    pattern = "-([0-9]+)\\.pcadapt$",
    replacement = paste0("-res-", version, "-\\1.rds"),
    x = data_file)

  if (!file.exists(res_file)) {
    input <- read.pcadapt(data_file, type = "pcadapt")
    obj.pcadapt <- pcadapt(input, K = 10)
    saveRDS(obj.pcadapt, file = res_file)
  }
}
