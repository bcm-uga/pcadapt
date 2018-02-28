################################################################################

context("NEW VERSION VS OLD VERSION")

res <- list()
tol <- 1e-5

G <- readRDS(system.file("testdata", 
                         paste0("full", ".rds"), 
                         package = "pcadapt")) 

old.pcadapt <- readRDS(system.file("testdata", 
                                   paste0("pcadapt_v3", ".rds"), 
                                   package = "pcadapt")) 

tmp <- read.pcadapt(G)
new.pcadapt <- pcadapt(tmp, K = 5)

testthat::expect_equal(apply(old.pcadapt[["scores"]], 2, mean),
                       apply(new.pcadapt[["scores"]], 2, mean), 
                       tolerance = tol)

