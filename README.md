# pcadapt

[![DOI](https://zenodo.org/badge/doi/10.1093/molbev/msaa053.svg)](https://doi.org/10.1093/molbev/msaa053)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/pcadapt)](https://cran.r-project.org/package=pcadapt)
[![R build status](https://github.com/bcm-uga/pcadapt/workflows/R-CMD-check/badge.svg)](https://github.com/bcm-uga/pcadapt)

R package {pcadapt} has been developed to detect genetic markers involved in biological adaptation. 
This package uses statistical tools for outlier detection based on Principal Component Analysis (PCA).

A tutorial for {pcadapt} is available on [pcadapt's website](https://bcm-uga.github.io/pcadapt/articles/pcadapt.html).

**News:** A new paper is finally out describing version 4 of pcadapt (which has been available for quite some time now). Please cite the new paper.


## Installation

```r
# From CRAN
install.packages("pcadapt")

# From GitHub
remotes::install_github("bcm-uga/pcadapt")
```


### References

- Privé, F., Luu, K., Vilhjálmsson, B. J., & Blum, M. G.B. (2020). [Performing highly efficient genome scans for local adaptation with R package pcadapt version 4.](https://doi.org/10.1093/molbev/msaa053) Molecular Biology and Evolution.

- Luu, K., Bazin, E., & Blum, M. G.B. (2017). [pcadapt: an R package to perform genome scans for selection based on principal component analysis.](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12592/full) Molecular Ecology Resources, 17(1), 67-77.

- Duforet-Frebourg, N., Luu, K., Laval, G., Bazin, E., & Blum, M. G.B. (2015). [Detecting genomic signatures of natural selection with principal component analysis: application to the 1000 Genomes data.](http://mbe.oxfordjournals.org/content/33/4/1082) Molecular biology and evolution, 33(4), 1082-1093.
