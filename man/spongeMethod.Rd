\name{spongeMethod}
\alias{spongeMethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying miRNA sponge interactions using spongeMethod
}
\description{
Identifying miRNA sponge interactions using spongeMethod. We
implement seven popular methods (miRHomology, pc, sppc, ppc,
hermes, muTaME, and cernia) to identify miRNA sponge
interactions.
}
\usage{
spongeMethod(miRTarget, ExpData = NULL, mres = NULL, minSharedmiR = 3,
poscorcutoff = 0, num_perm = 100, padjustvaluecutoff = 0.01,
padjustmethod = "BH", senscorcutoff = 0.3, scorecutoff = 0.5,
method = c("miRHomology", "pc", "sppc", "ppc", "hermes", "muTaME", "cernia"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{miRTarget}{
Putative miRNA-target interactions. Required option for method "miRHomology", "pc", "sppc", "ppc", "hermes", "muTaME" and "cernia".
}
  \item{ExpData}{
An input expression data frame, the columns are genes and
the rows are samples. Required option for method "pc", "sppc", "ppc", "hermes" and "cernia".
}
  \item{mres}{
Putative MiRNA Response Elements (mres) data frame, each row
contains five elements: Mirna, Target, energy, gap_l, gap_r. Required option for method "muTaME" and "cernia".
}
  \item{minSharedmiR}{
The minimum number of shared miRNAs between targets. Required option for method "miRHomology", "pc", "sppc", "ppc", "hermes", "muTaME" and "cernia".
}
  \item{poscorcutoff}{
A cutoff value of positive correlation. Required option for method "pc", "sppc" and "cernia".
}
  \item{num_perm}{
The number of permutations. Required option for method "ppc" and "hermes".
}
  \item{padjustvaluecutoff}{
A cutoff value of adjusted p-values. Required option for method "miRHomology", "pc", "sppc", "ppc", "hermes", "muTaME" and "cernia".
}
  \item{padjustmethod}{
Adjusted method of p-values, can select one of "holm", "hochberg",
"hommel", "bonferroni", "BH", "BY", "fdr", "none". Required option for method "miRHomology", "pc", "sppc", "ppc", "hermes", "muTaME" and "cernia".
}
  \item{senscorcutoff}{
A cutoff value of sensitivity partial pearson correlation. Required option for method "sppc".
}
  \item{scorecutoff}{
A cutoff value of normalized score (range from 0 to 1). Required option for method "muTaME" and "cernia".
}
  \item{method}{
Select a method for identifying miRNA sponge interactions, can select
one of "miRHomology", "pc", "sppc", "ppc", "hermes", "muTaME", "cernia".
}
}
\value{
A list of identified miRNA sponge interactions.
}
\references{
1. Le TD, Zhang J, Liu L, et al. Computational methods for identifying
miRNA sponge interactions. Brief Bioinform., 2017, 18(4):577-590.

2. Li JH, Liu S, Zhou H, et al. starBase v2.0: decoding miRNA-ceRNA,
miRNA-ncRNA and protein-RNA interaction networks from large-scale
CLIP-Seq data. Nucleic Acids Res., 2014, 42(Database issue):D92-7.

3. Sarver AL, Subramanian S. Competing endogenous RNA database.
Bioinformation, 2012, 8(15):731-3.

4. Zhou X, Liu J, Wang W, Construction and investigation of
breast-cancer-specific ceRNA network based on the mRNA and
miRNA expression data. IET Syst Biol., 2014, 8(3):96-103.

5. Xu J, Li Y, Lu J, et al. The mRNA related ceRNA-ceRNA landscape
and significance across 20 major cancer types. Nucleic Acids Res.,
2015, 43(17):8169-82.

6. Paci P, Colombo T, Farina L, Computational analysis identifies a
sponge interaction network between long non-coding RNAs and messenger
RNAs in human breast cancer. BMC Syst Biol., 2014, 8:83.

7. Sumazin P, Yang X, Chiu HS, et al. An extensive microRNA-mediated
network of RNA-RNA interactions regulates established oncogenic pathways
in glioblastoma. Cell, 2011, 147(2):370-81.

8. Tay Y, Kats L, Salmena L, et al. Coding-independent regulation of
the tumor suppressor PTEN by competing endogenous mRNAs. Cell,
2011, 147(2):344-57.

9. Sardina DS, Alaimo S, Ferro A, et al. A novel computational method
for inferring competing endogenous interactions. Brief Bioinform.,
2016, DOI: 10.1093/bib/bbw084.
}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}
\examples{
# Obtain expression data file "ExpData.csv" in csv format
ExpDatacsv <- system.file("extdata","ExpData.csv",package="miRspongeR")
ExpData <- read.csv(ExpDatacsv, header=FALSE, sep=",", stringsAsFactors = TRUE)

# Obtain miRNA-target interaction data file "miR2Target.csv" in csv format
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRspongeR")
miRTarget <- read.csv(miR2Target, header=TRUE, sep=",")

# miRHomology method
miRHomologyceRInt <- spongeMethod(miRTarget, method = "miRHomology")

# pc method
pcceRInt <- spongeMethod(miRTarget, ExpData, method = "pc")

# sppc method
sppcceRInt <- spongeMethod(miRTarget, ExpData, senscorcutoff = 0.1, method = "sppc")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

