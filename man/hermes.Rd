\name{hermes}
\alias{hermes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying miRNA sponge interactions using Hermes method

}
\description{
Identifying miRNA sponge interactions using Hermes method
(\url{http://califano.c2b2.columbia.edu/hermes}).
}
\usage{
hermes(miR2Target, ExpDatacsv, minSharedmiR = 3,
num_perm = 100, padjustvaluecutoff = 0.01, padjustmethod = "BH")
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
    \item{num_perm}{
The number of permutations.
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
Sumazin P, Yang X, Chiu HS, et al. An extensive microRNA-mediated
network of RNA-RNA interactions regulates established oncogenic
pathways in glioblastoma. Cell, 2011, 147(2):370-81.
}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}
\examples{

miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
ExpData <- system.file("extdata", "ExpData.csv", package="miRsponge")
HermesceRInt <- hermes(miR2Target, ExpData, num_perm =5)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
