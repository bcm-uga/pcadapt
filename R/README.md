# Introgression

```{r}
require(pcadapt)
require(data.table)
```

# Introgression

### Download the populus dataset here and set the working directory to where the data have been downloaded: https://www.dropbox.com/sh/85cgelc1af1s7az/AABvFT5W7NrDd_cXQHTff2Pma?dl=0


```{r}
setwd("/Users/Keurcien/Documents/thesis/git/Introgression/populus/Original_data_set_ch6_12_15/")
filename <- "populus3pops.pcadapt" 
popfile <- "populus3pops.pop" 
pop <- as.character(read.table(popfile)[, 1]) 
```

## Getting a complete dataset

Our method requires complete datasets, in case your data contain missing values, the package provides a function to impute missing values in two different ways, depending on whether individual labels are known or not.
In our example, `pop` contains the population information, i.e `pop[i]` returns the name of the population which the `i`-th individual belongs to. 

```{r}
head(pop)
```

```{r, echo=FALSE}
geno <- impute.pcadapt(filename, pop = pop)$x
```

If this information is lacking, leave the `pop` argument to blank.

## Running the scan

Provide the function with the following arguments:

- input a genotype matrix or a character string specifying the name of the file to be imputed.
- K a vector of integers specifying the components along which local ancestries may vary.
- pop a vector of integers or strings specifying which subpopulation the individuals belong to.
- ancstrl.1 a string specifying the label of the ancestral population genetically closer to the hybrid population.
- ancstrl.2 a string specifying the label of the ancestral population genetically further from the hybrid population.
- admxd a string specifying the label of thehybrid population.
- ploidy an integer specifying the ploidy of the individuals.
- window.size an integer specifying the window size.

```{r, echo=FALSE}
stat <- scan.intro( input = geno, 
                    K = 1, 
                    pop = pop, 
                    ancsrtrl.1 = "Trichocarpa", 
                    ancstr.2 = "Balsamifera", 
                    admxd = "Hybrid", 
                    ploidy = 2, 
                    window.size = 15000)
```

## Analyzing the output

High peaks correspond to excess of local ancestry from `ancstrl.2`.

```{r, echo=FALSE}
subset <- seq(1, length(stat), by = 10) #to display one out of ten points, thus reducing plotting time
plot(stat[subset], cex = 0.1, col = "red", ylab = "p-values")
```