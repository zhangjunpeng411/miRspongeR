# miRspongeR R package

# Introduction
This package provides several functions to explore miRNA sponge (also called ceRNA or miRNA decoy) regulation from putative miRNA-target interactions or/and transcriptomics data (including bulk, single-cell and spatial gene expression data). It provides eight popular methods for identifying miRNA sponge interactions, and an integrative method to integrate miRNA sponge interactions from different methods, as well as the functions to validate miRNA sponge interactions, and infer miRNA sponge modules, conduct enrichment analysis of miRNA sponge modules, and conduct survival analysis of miRNA sponge modules. By using a sample control variable strategy, it provides a function to infer sample-specific miRNA sponge interactions. In terms of sample-specific miRNA sponge interactions, it implements three similarity methods to construct sample-sample correlation network.

# Installation
```{r echo=FALSE, results='hide', message=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("miRspongeR")
```

# A quick example to use miRspongeR package
```{r echo=FALSE, results='hide', message=FALSE}
# Load the package
library(miRspongeR)

# Identifying miRNA sponge interactions using the miRHomology method
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRspongeR")
miRTarget <- read.csv(miR2Target, header=TRUE, sep=",")
miRHomologyceRInt <- spongeMethod(miRTarget, method = "miRHomology")

# Validation of the identified miRNA sponge interactions
Groundtruthcsv <- system.file("extdata", "Groundtruth.csv", package="miRspongeR")
Groundtruth <- read.csv(Groundtruthcsv, header=TRUE, sep=",")
spongenetwork_validated <- spongeValidate(miRHomologyceRInt[, 1:2], directed = FALSE, Groundtruth)

# Module identification from miRNA sponge interaction network
spongenetwork_Cluster <- netModule(miRHomologyceRInt[, 1:2], modulesize = 2)

# Disease and functional enrichment analysis of miRNA sponge modules
sponge_Module_DEA <- moduleDEA(spongenetwork_Cluster)
sponge_Module_FEA <- moduleFEA(spongenetwork_Cluster)

# Survival analysis of miRNA sponge modules
ExpDatacsv <- system.file("extdata", "ExpData.csv", package="miRspongeR")
ExpData <- read.csv(ExpDatacsv, header=TRUE, sep=",")
SurvDatacsv <- system.file("extdata", "SurvData.csv", package="miRspongeR")
SurvData <- read.csv(SurvDatacsv, header=TRUE, sep=",")
sponge_Module_Survival <- moduleSurvival(spongenetwork_Cluster, 
    ExpData, SurvData, devidePercentage=.5)
```

# License
GPL-3
