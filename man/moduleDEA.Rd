\name{moduleDEA}
\alias{moduleDEA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Disease enrichment analysis of modules
}
\description{
Disease enrichment analysis of modules. The disease ontology
databases have three types including DO: Disease Ontology database
(\url{http://disease-ontology.org/}), DGN: DisGeNET database
(\url{http://www.disgenet.org/}), and NCG: Network of Cancer
Genes database (\url{http://ncg.kcl.ac.uk/}).
}
\usage{
moduleDEA(Modulelist, OrgDb = "org.Hs.eg.db", 
padjustvaluecutoff = 0.05, padjustedmethod = "BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Modulelist}{
A list of miRNA sponge modules.
}
    \item{OrgDb}{
OrgDb
}
    \item{padjustvaluecutoff}{
A cutoff value of adjusted p-values.
}
    \item{padjustedmethod}{
Adjusted method of p-values, can select one of "holm", "hochberg",
"hommel", "bonferroni", "BH", "BY", "fdr", "none".
}
}

\value{
A list of disease enrichment analysis results.
}
\references{
Yu G, Wang L, Yan G, et al. DOSE: an R/Bioconductor package for
Disease Ontology Semantic and Enrichment analysis. Bioinformatics, 2015,
31(4):608-609.
}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}
\examples{
# Obtain miRNA-target interaction data file "miR2Target.csv" in csv format
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRspongeR")
miRTarget <- read.csv(miR2Target, header=TRUE, sep=",")
miRHomologyceRInt <- spongeMethod(miRTarget, method = "miRHomology")
spongenetwork_Cluster <- netModule(miRHomologyceRInt[, 1:2])
sponge_Module_DEA <- moduleDEA(spongenetwork_Cluster)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
