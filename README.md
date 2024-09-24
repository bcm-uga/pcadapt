<!-- badges: start -->
[![DOI](https://zenodo.org/badge/doi/10.1093/molbev/msaa053.svg)](https://doi.org/10.1093/molbev/msaa053)
[![CRAN status](https://www.r-pkg.org/badges/version/pcadapt)](https://CRAN.R-project.org/package=pcadapt)
[![R-CMD-check](https://github.com/bcm-uga/pcadapt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bcm-uga/pcadapt/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# pcadapt

R package {pcadapt} was developed to detect genetic markers involved in biological adaptation. 
This package uses statistical tools for outlier detection based on Principal Component Analysis (PCA).

A tutorial for {pcadapt} is available on [pcadapt's website](https://bcm-uga.github.io/pcadapt/articles/pcadapt.html).

**News:** Our paper describing version 4 of pcadapt (which has been available for quite some time now) is out. Please cite the newest paper (see Refs below).


## Installation

```r
# From CRAN
install.packages("pcadapt")

# From GitHub
remotes::install_github("bcm-uga/pcadapt")
```


### References

- Privé, F., Luu, K., Vilhjálmsson, B. J., & Blum, M. G.B. (2020). [Performing highly efficient genome scans for local adaptation with R package pcadapt version 4.](https://doi.org/10.1093/molbev/msaa053) *Molecular Biology and Evolution* 37.7 (2020): 2153–2154.

- Luu, K., Bazin, E., & Blum, M. G.B. (2017). [pcadapt: an R package to perform genome scans for selection based on principal component analysis.](https://doi.org/10.1111/1755-0998.12592) *Molecular Ecology Resources*, 17(1), 67–77.

- Duforet-Frebourg, N., Luu, K., Laval, G., Bazin, E., & Blum, M. G.B. (2015). [Detecting genomic signatures of natural selection with principal component analysis: application to the 1000 Genomes data.](https://doi.org/10.1093/molbev/msv334) *Molecular biology and evolution*, 33(4), 1082–1093.
