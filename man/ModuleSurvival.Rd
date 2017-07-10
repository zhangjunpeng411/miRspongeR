\name{ModuleSurvival}
\alias{ModuleSurvival}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Survival analysis of modules
}
\description{
Survival analysis of modules.
}
\usage{
ModuleSurvival(Modulelist, ExpDatacsv, SurvDatacsv, plot = FALSE)
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
time and status, respectively.
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
PCceRInt <- PC(miR2Target, ExpData)
spongenetwork_Cluster <- NetModule(PCceRInt[, 1:2])
sponge_Module_Survival <- ModuleSurvival(spongenetwork_Cluster,
ExpData, SurvData)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
