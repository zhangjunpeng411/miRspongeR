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
moduleSurvival(Modulelist, ExpData, SurvData,
devidePercentage=.5, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Modulelist}{
A list of miRNA sponge modules.
}
    \item{ExpData}{
An input expression data, the columns are genes and the rows are samples.
}
    \item{SurvData}{
An input survival data, three columns contain the information of sample,
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
# Obtain expression data file "ExpData.csv" in csv format
ExpDatacsv <- system.file("extdata","ExpData.csv",package="miRspongeR")
ExpData <- read.csv(ExpDatacsv, header=FALSE, sep=",")

# Obtain miRNA-target interaction data file "miR2Target.csv" in csv format
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRspongeR")
miRTarget <- read.csv(miR2Target, header=TRUE, sep=",")

# Obtain survival data file "SurvData.csv" in csv format
SurvDatacsv <- system.file("extdata", "SurvData.csv", package="miRspongeR")
SurvData <- read.csv(SurvDatacsv, header=TRUE, sep=",")
pcceRInt <- spongeMethod(miRTarget, ExpData, method = "pc")
spongenetwork_Cluster <- netModule(pcceRInt[, 1:2])
sponge_Module_Survival <- moduleSurvival(spongenetwork_Cluster,
ExpData, SurvData)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
