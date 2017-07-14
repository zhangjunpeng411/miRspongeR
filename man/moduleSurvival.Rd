\name{moduleSurvival}
\alias{moduleSurvival}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Survival analysis of modules
}
\description{
Survival analysis of modules.
}
\usage{
moduleSurvival(Modulelist, ExpDatacsv, SurvDatacsv,
devidePercentage=.5, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Modulelist}{
A list of miRNA sponge modules.
}
    \item{ExpDatacsv}{
An input expression dataset, the columns are genes and the rows are samples.
}
    \item{SurvDatacsv}{
An input survival dataset, three columns contain the information of sample,
time (Months) and status, respectively.
}
    \item{devidePercentage}{
A percentage value, the percentage of high risk group.
}
    \item{plot}{
A logical value, plot or not.
}
}

\value{
Survival analysis result of modules.
}
\references{
Terry M. Therneau and Patricia M. Grambsch. Modeling Survival Data:
Extending the Cox Model. Springer, New York. ISBN 0-387-98784-3, 2000.
}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}
\examples{
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
ExpData <- system.file("extdata", "ExpData.csv", package="miRsponge")
SurvData <- system.file("extdata", "SurvData.csv", package="miRsponge")
pcceRInt <- pc(miR2Target, ExpData)
spongenetwork_Cluster <- netModule(pcceRInt[, 1:2])
sponge_Module_Survival <- moduleSurvival(spongenetwork_Cluster,
ExpData, SurvData)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
