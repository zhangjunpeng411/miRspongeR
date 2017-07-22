\name{querymiRTargetbinding}
\alias{querymiRTargetbinding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Query miRNA-target interactions by combining expression data
and putative miRNA-target interactions
%%  ~~function to do ... ~~
}
\description{
Query miRNA-target interactions by combining expression data
and putative miRNA-target interactions.
}
\usage{
querymiRTargetbinding(ExpData, miRTarget)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{ExpData}{
An input gene expression data frame, the rows are samples
and the columns are genes.
}
    \item{miRTarget}{
An input miRNA-target interaction data frame, the first
column is miRNA with name "mir" and the second column is target
with name "gene".
}
}

\value{
A list of queried miRNA-target interactions
}

\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}

\examples{
ExpDatacsv <- system.file("extdata", "ExpData.csv", package="miRsponge")
ExpData <- read.csv(ExpDatacsv, header=FALSE, sep=",")
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
miRTarget <- read.csv(miR2Target, header=TRUE, sep=",")
miR2Target_queried <- querymiRTargetbinding(ExpData, miRTarget)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
