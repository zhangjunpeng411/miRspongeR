# miRsponge R package

# Introduction
This package provides several functions to study miRNA sponge (also called ceRNA or miRNA decoy), including popular methods for identifying miRNA sponge interactions, and the integrative method to integrate miRNA sponge interactions from different methods, as well as the functions to validate miRNA sponge interactions, and infer miRNA sponge modules, conduct enrichment analysis of modules, and conduct survival analysis of modules.

# Installation
```{r echo=FALSE, results='hide', message=FALSE}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("miRsponge")
```

# Usage
```{r echo=FALSE, results='hide', message=FALSE}
library(miRsponge)
```

# License
GPL-3
