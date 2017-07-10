\name{Integrate}
\alias{Integrate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrate method for identifying miRNA sponge interactions
by integrating different methods
}
\description{
Integrate method for identifying miRNA sponge interactions
by integrating different methods.
}
\usage{
Integrate(Interlist, Intersect_num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Interlist}{
List object, a list of miRNA sponge interactions from different methods.
}
    \item{Intersect_num}{
The least number of different methods intersected for integration.
The value of 1 means the union of miRNA sponge interactions
from different methods.
}
}

\value{
A list of integrated miRNA sponge interactions.
}

\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}

\examples{
ExpData <- system.file("extdata","ExpData.csv",package="miRsponge")
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
miRHomologyceRInt <- miRHomology(miR2Target)
PCceRInt <- PC(miR2Target, ExpData)
SPPCceRInt <- SPPC(miR2Target, ExpData)
Interlist=list(miRHomologyceRInt[, 1:2],PCceRInt[, 1:2],SPPCceRInt[, 1:2])
IntegrateceRInt <- Integrate(Interlist, 3)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
