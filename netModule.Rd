\name{netModule}
\alias{netModule}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying miRNA sponge modules from network
}
\description{
Identifying miRNA sponge modules from network. Possible methods include
FN, MCL, LINKCOMM and MCODE.
}
\usage{
netModule(spongenetwork, method = "MCL", directed = FALSE, save = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{spongenetwork}{
Input miRNA sponge interaction network.
}
    \item{method}{
Cluster method, can select one of FN, MCL, LINKCOMM and MCODE.
}
    \item{directed}{
A logical value, the network is directed or not.
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

}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}

\examples{
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
miRTarget <- read.csv(miR2Target, header=TRUE, sep=",")
miRHomologyceRInt <- spongeMethod(miRTarget, method = "miRHomology")
spongenetwork_Cluster <- netModule(miRHomologyceRInt[, 1:2])
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
