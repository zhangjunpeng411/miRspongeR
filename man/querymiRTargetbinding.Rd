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
querymiRTargetbinding(ExpDatacsv, miRTargetbinding)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{ExpDatacsv}{
An input gene expression dataset in .csv format, the rows are samples
and the columns are genes.
}
    \item{miRTargetbinding}{
An input miRNA-target interaction dataset in .csv format, the first
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
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
miR2Target_queried <- querymiRTargetbinding(ExpDatacsv, miR2Target)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
