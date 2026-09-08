library(pcadapt)

get_statistics <- getFromNamespace("get_statistics", "pcadapt")

check_scaling <- function(z) {
  result <- get_statistics(z, method = "componentwise",
                           pass = rep(TRUE, nrow(z)))
  expected_gif <- apply(z^2, 2, median, na.rm = TRUE) / qchisq(.5, 1)
  expected <- z^2
  expected_p <- z^2
  for (k in seq_len(ncol(z))) {
    expected[, k] <- z[, k]^2 / expected_gif[k]
    expected_p[, k] <- pchisq(expected[, k], 1, lower.tail = FALSE)
  }
  stopifnot(isTRUE(all.equal(unname(result$gif), unname(expected_gif))),
            identical(dim(result$chi2.stat), dim(z)),
            isTRUE(all.equal(unname(result$chi2.stat), unname(expected))),
            isTRUE(all.equal(unname(result$pvalues), unname(expected_p))))
  result
}

# Columns differing only by scale should have identical calibrated statistics.
z <- cbind(PC1 = c(-2, -1, 1, 2), PC2 = c(-20, -10, 10, 20))
result <- check_scaling(z)
stopifnot(isTRUE(all.equal(result$chi2.stat[, 1], result$chi2.stat[, 2])))

# Single component; row count not divisible by component count; missing values.
invisible(check_scaling(z[, 1, drop = FALSE]))
invisible(check_scaling(cbind(c(-2, -1, 1, 2, 3),
                             c(-20, -10, 10, 20, 30))))
z[1, ] <- NA_real_
invisible(check_scaling(z))
