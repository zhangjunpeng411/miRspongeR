\name{CERNIA}
\alias{CERNIA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying miRNA sponge interactions using CERNIA method

}
\description{
Identifying miRNA sponge interactions based on using CERNIA method.

}
\usage{
CERNIA(miR2Target, MREs, ExpDatacsv, minSharedmiR = 3, poscorcutoff = 0,
padjustvaluecutoff = 0.01, padjustmethod = "BH", scorecutoff = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{miR2Target}{
Putative miRNA-target interactions.
}
    \item{MREs}{
Putative MiRNA Response Elements (MREs) dataset, each row
contains five elements: Mirna, Target, energy, gap_l, gap_r.
}
    \item{ExpDatacsv}{
An input expression dataset, the columns are genes and
the rows are samples.
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
Adjusted method of p-values, can select one of "holm",
"hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
}
    \item{scorecutoff}{
A cutoff value of normalized score (range from 0 to 1).
}
}

\value{
A list of identified miRNA sponge interactions.
}
\references{
Sardina DS, Alaimo S, Ferro A, et al. A novel computational method
for inferring competing endogenous interactions. Brief Bioinform.,
2016, DOI: 10.1093/bib/bbw084.
}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}
\examples{
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
MREs <- system.file("extdata", "MREs.csv", package="miRsponge")
ExpData <- system.file("extdata", "ExpData.csv", package="miRsponge")
CERNIAceRInt <- CERNIA(miR2Target, MREs,ExpData)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
