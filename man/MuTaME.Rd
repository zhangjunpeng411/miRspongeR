\name{MuTaME}
\alias{MuTaME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying miRNA sponge interactions using MuTaME method

}
\description{
Identifying miRNA sponge interactions based on using MuTaME method.

}
\usage{
MuTaME(miR2Target, MREs, minSharedmiR = 3, padjustvaluecutoff
= 0.01,padjustmethod = "BH", scorecutoff = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{miR2Target}{
Putative miRNA-target interactions.
}
    \item{MREs}{
Putative MiRNA Response Elements (MREs) dataset, each row contains
five elements: Mirna, Target, energy, gap_l, gap_r.
}
    \item{minSharedmiR}{
The minimum number of shared miRNAs between targets.
}
    \item{padjustvaluecutoff}{
A cutoff value of adjusted p-values.
}
    \item{padjustmethod}{
Adjusted method of p-values, can select one of "holm", "hochberg",
"hommel", "bonferroni", "BH", "BY", "fdr", "none".
}
    \item{scorecutoff}{
A cutoff value of normalized score (range from 0 to 1).
}
}

\value{
A list of identified miRNA sponge interactions.
}

\references{
Tay Y, Kats L, Salmena L, et al. Coding-independent regulation of
the tumor suppressor PTEN by competing endogenous mRNAs. Cell,
2011, 147(2):344-57.
}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}
\examples{
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
MREs <- system.file("extdata", "MREs.csv", package="miRsponge")
MuTaMEceRInt <- MuTaME(miR2Target,MREs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
