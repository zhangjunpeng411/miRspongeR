\name{netModule}
\alias{netModule}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying miRNA sponge modules from network
}
\description{
Identifying miRNA sponge modules from network. Possible methods include
FN, MCL, LINKCOMM, MCODE, betweenness, infomap, prop, eigen, louvain, walktrap.
}
\usage{
netModule(spongenetwork, method = "MCL", directed = FALSE, modulesize = 3, save = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{spongenetwork}{
Input miRNA sponge interaction network.
}
    \item{method}{
Cluster method, can select one of FN, MCL, LINKCOMM, MCODE, betweenness, infomap, prop, eigen, louvain, walktrap.
}
    \item{directed}{
A logical value, the network is directed or not.
}
    \item{modulesize}{
The size cutoff of the identified modules.
}
    \item{save}{
A logical value, save the identified modules or not.
}
}

\value{
A list of miRNA sponge modules.
}
\references{

1. Clauset A, Newman ME, Moore C. Finding community structure in very
large networks. Phys Rev E Stat Nonlin Soft Matter Phys.,
2004, 70(6 Pt 2):066111.

2. Enright AJ, Van Dongen S, Ouzounis CA. An efficient algorithm for
large-scale detection of protein families. Nucleic Acids Res., 2002,
30(7):1575-84.

3. Kalinka AT, Tomancak P. linkcomm: an R package for the generation,
visualization, and analysis of link communities in networks of arbitrary
size and type. Bioinformatics, 2011, 27(14):2011-2.

4. Bader GD, Hogue CW. An automated method for finding molecular complexes
in large protein interaction networks. BMC Bioinformatics, 2003, 4:2.

5. Newman ME, Girvan M. Finding and evaluating community structure in networks. Phys Rev E Stat Nonlin Soft Matter Phys. 2004;69(2 Pt 2):026113. 

6. Rosvall M, Bergstrom CT. Maps of random walks on complex networks reveal community structure. Proc Natl Acad Sci U S A. 2008;105(4):1118-1123.

7. Raghavan UN, Albert R, Kumara S. Near linear time algorithm to detect community structures in large-scale networks. Phys Rev E Stat Nonlin Soft Matter Phys. 2007;76(3 Pt 2):036106.

8. Newman ME. Finding community structure in networks using the eigenvectors of matrices. Phys Rev E Stat Nonlin Soft Matter Phys. 2006;74(3 Pt 2):036104.

9. Blondel VD, Guillaume JL, Lambiotte R, Lefebvre E. Fast unfolding of communities in large networks. Journal of statistical mechanics: theory and experiment, 2008, 2008(10): P10008.

10. Pons P, Latapy M. Computing communities in large networks using random walks. Graph Algorithms Appl. 2006.

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
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
