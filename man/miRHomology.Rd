\name{miRHomology}
\alias{miRHomology}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying miRNA sponge interactions based on the homology of miRNA binding sites present in the 3'UTR sequence
}
\description{
Identifying miRNA sponge interactions based on the homology of miRNA binding sites present in the 3'UTR sequence.
}
\usage{
miRHomology(miR2Target, minSharedmiR = 3, padjustvaluecutoff = 0.01, padjustmethod = "BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{miR2Target}{
Putative miRNA-target interactions.
}
  \item{minSharedmiR}{
The minimum number of shared miRNAs between targets.
}
  \item{padjustvaluecutoff}{
A cutoff value of adjusted p-values.
}
  \item{padjustmethod}{
Adjusted method of p-values, can select one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
}
}

\value{
A list of identified miRNA sponge interactions.
}

\references{
1. Le TD, Zhang J, Liu L, et al. Computational methods for identifying miRNA sponge interactions. Brief Bioinform. 2016, DOI: 10.1093/bib/bbw042.

2. Li JH, Liu S, Zhou H, et al. starBase v2.0: decoding miRNA-ceRNA, miRNA-ncRNA and protein-RNA interaction networks from large-scale CLIP-Seq data. Nucleic Acids Res., 2014, 42(Database issue):D92-7.

3. Sarver AL, Subramanian S. Competing endogenous RNA database. Bioinformation,2012, 8(15):731-3.

}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}

\examples{
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRsponge")
miRHomologyceRInt <- miRHomology(miR2Target)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
