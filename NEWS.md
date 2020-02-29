|Date       | Version          | Changes|
|--------   | ----   | -------------------------------------------------------|
|02-21-2018 | 4.0              | - `read.pcadapt` generates `bed` files instead of `pcadapt` files.|
|           |                  | - Computation of PCA is based on the package `RSpectra`.|
|           |                  | - Missing values are handled by specifying vector-matrix operation in `RSpectra` that accounts for missing values.|
|           |                  | - Includes LD thinning to compute PCs.|
|           |                  | - No more dependency to RcppArmadillo|
|           |                  | - For Pooled-seq data, use Mahalanobis distances based on PCA loadings, no more simulations of individuals genotypes|
|01-15-2017 | 3.1.0            | - Switch from C/Lapack to Rcpp/RcppArmadillo.|
|           |                  | - `pcadapt` can take genotype matrices as input.|
|           |                  | - Modified code for binomial sampling.|
|           |                  | - `pcadapt` argument `clean.files` is now deprecated.|
|           |                  | - `pcadapt` argument `output.filename` is now deprecated.|
|           |                  | - `read.pcadapt` argument `local.env` is now deprecated.|
|           |                  | - Latest update of `vcfR` taken into account.|
|12-20-2016 | 3.0.4            | - Method based on sampling genotypes added to handle pooled-sequencing.|
|10-06-2016 | 3.0.3            | - Option `type="vcfR"` has been added to `read.pcadapt` to overcome some conversion issues occurring with VCF files.|
|           |                   | - Argument `transpose` is now deprecated. Read section A for more details.|
| | |
|06-13-2016 | 3.0.2             | - The function `get.pc` has been added. For each SNP, it returns the most correlated principal component.|
| | |
|06-01-2016 | 3.0.1             | - The `read4pcadapt` is now deprecated, it is now called `read.pcadapt`.|
|           |                   | - Using the `pop` option when plotting scores now provides the color legend.|
| | |
|04-04-2016 | 3.0               | - All analyses are now included in the `R` package. Users should not use the `C` software PCAdapt fast anymore.|
|           |                   | - Big datasets can be handled directly within the `R` session.|
|           |                   | - `read4pcadapt` now converts files to the `pcadapt` format.|
|           |                   | - The first argument of `pcadapt` can be either a small genotype matrix or the output of `read4pcadapt`.|
| | |
|02-12-2016 | 2.2               | - The Mahalanobis distance is now estimated from the z-scores rather than the loadings.|
|           |                   | - Make sure you have downloaded the latest version of the C software PCAdapt (last updated on February 11, 2016).|
| | |
|01-05-2016 | 2.1.1             | - Bug fix: vignette header added.|
| | |
|12-18-2015 | 2.1               | - The scaling of the SNP before computing PCA has been changed. Instead of using standard deviation, we now use the square root of $p(1-p)$ (haploid data) or of $2p(1-p)$ (diploid data) where $p$ is the minimum allele frequency.|
|           |                   | - Bug fix: the genomic inflation factor has been corrected when `K=1`.|
|           |                   | - Bug fix: a problem due to high proportion of missing data slowing the program has been fixed.|
|           |                   | - Argument `"minmaf"` has been replaced with `"min.maf"`.|
| | |
|           | 2.0.1             | - Vignette corrected: when reading output from the software PCAdapt, we do not mention the deprecated argument `PCAdapt` in the function `pcadapt`.|
|           |                   | - Bug fix: an issue occurring when reading outputs from PCAdapt has been fixed.|
|           |                   | - We mention in the vignette how to use Pool-seq data with the C software PCAdapt fast.|
| | |
|           | 2.0               | - The default test statistic is not the communality statistic anymore but the Mahalanobis distance.|
|           |                   | - Test statistic for Pool-seq data.|
