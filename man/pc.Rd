\name{pc}
\alias{pc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying miRNA sponge interactions based on the homology
of miRNA binding sites present in the 3'UTR sequence and
positive correlation
}
\description{
Identifying miRNA sponge interactions based on the homology
of miRNA binding sites present in the 3'UTR sequence and
positive correlation.
}
\usage{
pc(miR2Target, ExpDatacsv, minSharedmiR = 3, poscorcutoff = 0,
padjustvaluecutoff = 0.01, padjustmethod = "BH")
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
}

\value{
A list of identified miRNA sponge interactions.
}

\references{
1. Le TD, Zhang J, Liu L, et al. Computational methods for
identifying miRNA sponge interactions. Brief Bioinform. 2017,
18(4):577-590.

2. Zhou X, Liu J, Wang W, Construction and investigation of
breast-cancer-specific ceRNA network based on the mRNA and
miRNA expression data. IET Syst Biol., 2014, 8(3):96-103.

3. Xu J, Li Y, Lu J, et al. The mRNA related ceRNA-ceRNA landscape
and significance across 20 major cancer types. Nucleic Acids Res.,
2015, 43(17):8169-82.
}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}

\examples{
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
ExpData <- system.file("extdata", "ExpData.csv", package="miRsponge")
PCceRInt <- pc(miR2Target, ExpData)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
