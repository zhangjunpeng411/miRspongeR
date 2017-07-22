\name{moduleFEA}
\alias{moduleFEA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functional GO, KEGG and Reactome enrichment analysis of modules
}
\description{
Functional GO, KEGG and Reactome enrichment analysis of modules.
GO: Gene Ontology database (\url{http://www.geneontology.org/}),
KEGG: Kyoto Encyclopedia of Genes and Genomes Pathway Database
(\url{http://www.genome.jp/kegg/}) and Reactome: Reactome
Pathway Database (\url{http://reactome.org/}).
}
\usage{
moduleFEA(Modulelist, ont = "BP", KEGGorganism = "hsa",
Reactomeorganism = "human", OrgDb = "org.Hs.eg.db",
padjustvaluecutoff = 0.05, padjustedmethod = "BH", plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Modulelist}{
A list of miRNA sponge modules.
}
    \item{ont}{
One of "MF", "BP", and "CC" subontologies.
}
    \item{KEGGorganism}{
Organism, supported organism listed in
\url{http://www.genome.jp/kegg/catalog/org_list.html}
}
    \item{Reactomeorganism}{
Organism, one of "human", "rat", "mouse", "celegans", "yeast",
"zebrafish", "fly".
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
    \item{plot}{
A logical value, plot or not.
}
}

\value{
A list of functional GO, KEGG and Reactome enrichment analysis results.
}
\references{
1. Yu G, Wang L, Han Y, et al. clusterProfiler: an R package for
comparing biological themes among gene clusters. OMICS: A Journal
of Integrative Biology, 2012, 16(5):284-287.

2. Yu G and He Q. ReactomePA: an R/Bioconductor package for reactome
pathway analysis and visualization. Molecular BioSystems, 2016,
12(12), pp. 477-479.
}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}
\examples{
\dontrun{
ExpDatacsv <- system.file("extdata","ExpData.csv",package="miRsponge")
ExpData <- read.csv(ExpDatacsv, header=FALSE, sep=",")
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
miRTarget <- read.csv(miR2Target, header=TRUE, sep=",")
pcceRInt <- spongeMethod(miRTarget, ExpData, method = "pc")
spongenetwork_Cluster <- netModule(pcceRInt[, 1:2])
sponge_Module_FEA <- moduleFEA(spongenetwork_Cluster)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
