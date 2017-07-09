\name{SpongeValidate}
\alias{SpongeValidate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Validation of computationally predicted miRNA sponge interactions
}
\description{
Validation of computationally predicted miRNA sponge interactions. The groundtruth of miRNA sponge interactions are from miRSponge v2.0 (\url{http://www.bio-bigdata.net/miRSponge/}).
}
\usage{
SpongeValidate(spongenetwork, Groundtruthcsv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spongenetwork}{
Input miRNA sponge interaction network.
}
  \item{Groundtruthcsv}{
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
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
miRHomologyceRInt <- miRHomology(miR2Target)
Groundtruth <- system.file("extdata", "Groundtruth.csv", package="miRsponge")
spongenetwork_validated <- SpongeValidate(miRHomologyceRInt[, 1:2], Groundtruth)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
