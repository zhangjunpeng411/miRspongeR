\name{sample_cor_network}
\alias{sample_cor_network}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying sample-sample correlation network
}
\description{
Identifying sample-sample correlation network in terms of sample-specific miRNA sponge networks. 
}
\usage{
sample_cor_network(ceRNet, genes_num, method = "Simpson", 
padjustvaluecutoff = 0.01, padjustmethod = "BH", num.cores = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ceRNet}{
List object, igraph input sample-specific miRNA sponge networks.
}
  \item{genes_num}{
Number of possible miRNA sponges.
}
  \item{method}{
Methods for calculating similatiry between two sample-specific miRNA sponge networks, select one of three methods (Simpson, Jaccard and Lin). Default method is Simpson
}
  \item{padjustvaluecutoff}{
A cutoff value of adjusted p-values. 
}
  \item{padjustmethod}{
Adjusted method of p-values, can select one of "holm", "hochberg",
"hommel", "bonferroni", "BH", "BY", "fdr", "none". 
}
  \item{num.cores}{
The number of CPU cores to be selected. 
}
}
\value{
A list of sample-sample correlations.
}
\references{
1. Tucker CM, Cadotte MW, Carvalho SB, et al. A guide to phylogenetic metrics for conservation, community ecology and macroecology. Biol Rev Camb Philos Soc. 2017;92(2):698-715.

2. Jaccard P. The Distribution of the Flora in the Alpine Zone. The New Phytologist 11, no. 2 (1912): 37–50.

3. Lin D, et al. An information-theoretic definition of similarity. In: Icml, vol.
98. 1998. p. 296–304.
}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}
\examples{
# Obtain expression data file "ExpData.csv" in csv format
ExpDatacsv <- system.file("extdata","ExpData.csv",package="miRspongeR")
ExpData <- read.csv(ExpDatacsv, header=TRUE, sep=",")

# Obtain miRNA-target interaction data file "miR2Target.csv" in csv format
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRspongeR")
miRTarget <- read.csv(miR2Target, header=TRUE, sep=",")

# Identifying sample-specific miRNA sponge interactions, the sppc method 
# is used to identify miRNA sponge interactions
sponge_sample_specific_net <- sponge_sample_specific(miRTarget, ExpData, senscorcutoff = 0.1, method = "sppc")

sample_cor_network_res <- sample_cor_network(sponge_sample_specific_net, 
genes_num = 31, padjustvaluecutoff = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

