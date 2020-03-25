\name{integrateMethod}
\alias{integrateMethod}
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
integrateMethod(Interlist, Intersect_num)
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
# Obtain expression data file "ExpData.csv" in csv format
ExpDatacsv <- system.file("extdata","ExpData.csv",package="miRspongeR")
ExpData <- read.csv(ExpDatacsv, header=FALSE, sep=",", stringsAsFactors = TRUE)

# Obtain miRNA-target interaction data file "miR2Target.csv" in csv format
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRspongeR")
miRTarget <- read.csv(miR2Target, header=TRUE, sep=",")
miRHomologyceRInt <- spongeMethod(miRTarget, method = "miRHomology")
pcceRInt <- spongeMethod(miRTarget, ExpData, method = "pc")
sppcceRInt <- spongeMethod(miRTarget, ExpData, method = "sppc")
Interlist <- list(miRHomologyceRInt[, 1:2], pcceRInt[, 1:2], sppcceRInt[, 1:2])
IntegrateceRInt <- integrateMethod(Interlist, 2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
