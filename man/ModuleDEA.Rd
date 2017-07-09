\name{ModuleDEA}
\alias{ModuleDEA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Disease enrichment analysis of modules
}
\description{
Disease enrichment analysis of modules. The disease ontology databases have three types including DO: Disease Ontology database (\url{http://disease-ontology.org/}), DGN: DisGeNET database (\url{http://www.disgenet.org/}), and NCG: Network of Cancer Genes database (\url{http://ncg.kcl.ac.uk/}).
}
\usage{
ModuleDEA(Modulelist, OrgDb = "org.Hs.eg.db", ont = "DO",
padjustvaluecutoff = 0.05, padjustedmethod = "BH",plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Modulelist}{
A list of miRNA sponge modules.
}
  \item{OrgDb}{
OrgDb
}
  \item{ont}{
One of "DO", and "DOLite" subontologies.
}

  \item{padjustvaluecutoff}{
A cutoff value of adjusted p-values.
}

  \item{padjustedmethod}{
Adjusted method of p-values, can select one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
}

  \item{plot}{
A logical value, plot or not.
}
}

\value{
A list of disease enrichment analysis results.
}
\references{
Yu G, Wang L, Yan G, et al. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics, 2015, 31(4):608-609.
}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}
\examples{
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
miRHomologyceRInt <- miRHomology(miR2Target)
spongenetwork_Cluster <- NetModule(miRHomologyceRInt[, 1:2])
sponge_Module_DEA <- ModuleDEA(spongenetwork_Cluster)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
