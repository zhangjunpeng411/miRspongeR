\name{spongeValidate}
\alias{spongeValidate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Validation of computationally predicted miRNA sponge interactions
}
\description{
Validation of computationally predicted miRNA sponge interactions.
The built-in groundtruth of miRNA sponge interactions 
are from miRSponge (<http://bio-bigdata.hrbmu.edu.cn/miRSponge/>), 
lncACTdb (<http://bio-bigdata.hrbmu.edu.cn/LncACTdb/>), 
LncCeRBase (<http://www.insect-genome.com/LncCeRBase/front/>).
}
\usage{
spongeValidate(spongenetwork, directed = FALSE, Groundtruth)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{spongenetwork}{
Input miRNA sponge interaction network.
}
    \item{directed}{
A logical value, the network is directed or not.
}
    \item{Groundtruth}{
The groundtruth of miRNA sponge interactions.
}
}

\value{
A list of experimentally validated miRNA sponge interactions.
}

\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}

\examples{
# Obtain miRNA-target interaction data file "miR2Target.csv" in csv format
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRspongeR")
miRTarget <- read.csv(miR2Target, header=TRUE, sep=",")

# Obtain experimentally validated miRNA sponge interaction data file "Groundtruth.csv" in csv format
Groundtruthcsv <- system.file("extdata", "Groundtruth.csv", package="miRspongeR")
Groundtruth <- read.csv(Groundtruthcsv, header=TRUE, sep=",")
miRHomologyceRInt <- spongeMethod(miRTarget, method = "miRHomology")
spongenetwork_validated <- spongeValidate(miRHomologyceRInt[, 1:2], directed = FALSE, Groundtruth)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
