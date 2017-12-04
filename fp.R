# pMatVec4 and cpMatVec4 are implemented without unrolling optimization at the
# moment. The test below allows to control that any modification to these 
# functions do not affect the results. 

# devtools::install_github("privefl/pcadapt")
library(pcadapt)

G <- as.matrix(read.table("mv_nunif_0.2.pcadapt"))
G[G == 9] <- NA # pcadapt files had missing values encoded by 9

X <- runif(nrow(G)) # random vector 
Y <- runif(ncol(G)) # random vector

tmp <- apply(G, MARGIN = 1, FUN = function(x) {mean(x, na.rm = TRUE) / 2})

lookup_byte <- pcadapt:::getCode()

lookup_scale <- rbind(outer(0:2, tmp, function(g, p) {
  if (p > 0 || p < 1) {
    return((g - 2 * p) / sqrt(2 * p * (1 - p)))
  } else {
    return(0)
  }
}), 0)

lookup_geno <- outer(0:3, seq_len(nrow(G)), function(g, p) g)

pass <- rep(TRUE, nrow(G))
nb_nona <- pcadapt:::nb_nona(t(G), lookup_geno, lookup_byte, pass)

x <- pcadapt:::pMatVec4(t(G), X, lookup_scale = lookup_scale, 
                        lookup_byte = lookup_byte) / nb_nona[[1]]
y <- pcadapt:::prodtGx(G, X, tmp) # pcadapt used the transposed matrix
all.equal(x, y)

x <- pcadapt:::cpMatVec4(t(G), Y, lookup_scale = lookup_scale,
                         lookup_byte = lookup_byte) / nb_nona[[2]]
y <- pcadapt:::prodGx(G, Y, tmp) # pcadapt used the transposed matrix

all.equal(x, y)
