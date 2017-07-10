\name{SPPC}
\alias{SPPC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying miRNA sponge interactions based on Sensitivity
Partial Pearson Correlation (SPPC) method
}
\description{
Identifying miRNA sponge interactions based on Sensitivity
Partial Pearson Correlation (SPPC) method.
}
\usage{
SPPC(miR2Target, ExpDatacsv, minSharedmiR = 3, poscorcutoff = 0,
padjustvaluecutoff = 0.01, padjustmethod = "BH", senscorcutoff = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{miR2Target}{
Putative miRNA-target interactions.
}
    \item{ExpDatacsv}{
An input expression dataset, the columns are genes and the rows are samples.
}
    \item{minSharedmiR}{
The minimum number of shared miRNAs between targets.
}
    \item{poscorcutoff}{
A cutoff value of positive correlation.
}
    \item{padjustvaluecutoff}{
A cutoff value of adjusted p-values.
}
    \item{padjustmethod}{
Adjusted method of p-values, can select one of "holm", "hochberg",
"hommel", "bonferroni", "BH", "BY", "fdr", "none".
}
    \item{senscorcutoff}{
A cutoff value of sensitivity partial pearson correlation.
}
}

\value{
A list of identified miRNA sponge interactions.
}

\references{

1. Paci P, Colombo T, Farina L, Computational analysis identifies a
sponge interaction network between long non-coding RNAs and messenger
RNAs in human breast cancer. BMC Syst Biol., 2014, 8:83.

2. Le TD, Zhang J, Liu L, et al. Computational methods for identifying
miRNA sponge interactions. Brief Bioinform., 2016,
doi: 10.1093/bib/bbw042.

}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}

\examples{
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
ExpData <- system.file("extdata", "ExpData.csv", package="miRsponge")
SPPCceRInt <- SPPC(miR2Target, ExpData)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
