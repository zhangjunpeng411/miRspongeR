\name{sponge_sample_specific}
\alias{sponge_sample_specific}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying sample-specific miRNA sponge interactions
}
\description{
A sample control variable strategy is used to identify sample-specific miRNA sponge interactions. In the strategy, seven popular methods (pc, sppc, ppc,
hermes, muTaME, cernia, and SPONGE) to identify miRNA sponge interactions.
}
\usage{
sponge_sample_specific(miRTarget, ExpData = NULL, mres = NULL, consider.miR.expr = "TRUE", minSharedmiR = 3, poscorcutoff = 0, num_perm = 100, padjustvaluecutoff = 0.01, padjustmethod = "BH", senscorcutoff = 0.3, scorecutoff = 0.5, null_model, method = c("pc", "pc_parallel", "sppc", "sppc_parallel", "ppc", "ppc_parallel", "hermes", "hermes_parallel", "cernia", "cernia_parallel", "sponge_parallel"), num.cores = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{miRTarget}{
Putative miRNA-target interactions. Required option for method "pc", "pc_parallel", "sppc", "sppc_parallel", "ppc", "ppc_parallel", "hermes", "hermes_parallel", "muTaME", "muTaME_parallel", "cernia", "cernia_parallel", and "sponge_parallel".
}
  \item{ExpData}{
An input expression data frame, the columns are genes and
the rows are samples. Required option for method "pc", "pc_parallel", "sppc", "sppc_parallel", "ppc", "ppc_parallel", "hermes" "hermes_parallel", "cernia", "cernia_parallel", and "sponge_parallel".
}
  \item{mres}{
Putative MiRNA Response Elements (mres) data frame, each row
contains five elements: Mirna, Target, energy, gap_l, gap_r. Required option for method "muTaME", "muTaME_parallel", "cernia", and "cernia_parallel".
}
  \item{consider.miR.expr}{
Logical value, TRUE for considering miRNA expression data and FALSE for ignoring miRNA expression data
}
  \item{minSharedmiR}{
The minimum number of shared miRNAs between targets. Required option for method "pc", "pc_parallel", "sppc", "sppc_parallel", "ppc", "ppc_parallel", "hermes", "hermes_parallel", "muTaME", "muTaME_parallel", "cernia", "cernia_parallel", and "sponge_parallel".
}
  \item{poscorcutoff}{
A cutoff value of positive correlation. Required option for method "pc", "pc_parallel", "sppc", "sppc_parallel", "cernia", "cernia_parallel", and "sponge_parallel".
}
  \item{num_perm}{
The number of permutations. Required option for method "ppc", "ppc_parallel", "hermes", "hermes_parallel".
}
  \item{padjustvaluecutoff}{
A cutoff value of adjusted p-values. Required option for method "pc", "pc_parallel", "sppc", "sppc_parallel", "ppc", "ppc_parallel", "hermes", "hermes_parallel", "muTaME", "muTaME_parallel", "cernia", "cernia_parallel", and "sponge_parallel".
}
  \item{padjustmethod}{
Adjusted method of p-values, can select one of "holm", "hochberg",
"hommel", "bonferroni", "BH", "BY", "fdr", "none". Required option for method "miRHomology", "miRHomology_parallel", "pc", "pc_parallel", "sppc", "sppc_parallel", "ppc", "ppc_parallel", "hermes", "hermes_parallel", "muTaME", "muTaME_parallel", "cernia", "cernia_parallel", and "sponge_parallel".
}
  \item{senscorcutoff}{
A cutoff value of sensitivity partial pearson correlation. Required option for method "sppc", "sppc_parallel", and "sponge_parallel".
}
  \item{scorecutoff}{
A cutoff value of normalized score (range from 0 to 1). Required option for method "muTaME", "muTaME_parallel", "cernia", and "cernia_parallel".
}
  \item{null_model}{
Optional, pre-computed null model. Users can also build null model using "sponge_build_null_model" function in SPONGE R package (https://bioconductor.org/packages/SPONGE/). Required option for method "sponge_parallel".  
}
  \item{method}{
Select a method for identifying miRNA sponge interactions, can select
one of "pc", "pc_parallel", "sppc", "sppc_parallel", "ppc", "ppc_parallel", "hermes", "hermes_parallel", "muTaME", "muTaME_parallel", "cernia", "cernia_parallel", "sponge_parallel". The seven methods ("miRHomology_parallel", "pc_parallel", "sppc_parallel", "ppc_parallel", "hermes_parallel", "muTaME_parallel", "cernia_parallel") are the parallel versions of the seven original methods ("miRHomology", "pc", "sppc", "ppc", "hermes", "muTaME", "cernia").
}
  \item{num.cores}{
The number of CPU cores to be selected. Required option for method "pc_parallel", "sppc_parallel", "ppc_parallel",  "hermes_parallel", "muTaME_parallel", "cernia_parallel", and "sponge_parallel".
}
}
\value{
A list of sample-specific miRNA sponge interactions.
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

9. Sardina DS, Alaimo S, Ferro A, Pulvirenti A, Giugno R. A novel computational method for inferring competing endogenous interactions. Brief Bioinform. 2017;18(6):1071-1081.

10. List M, Dehghani Amirabad A, Kostka D, Schulz MH. Large-scale inference of competing endogenous RNA networks with sparse partial correlation. Bioinformatics. 2019;35(14):i596-i604.
}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}
\examples{
# Obtain expression data file "ExpData.csv" in csv format
ExpDatacsv <- system.file("extdata","ExpData.csv",package="miRspongeR")
ExpData <- read.csv(ExpDatacsv, header=TRUE, sep=",")

# Obtain miRNA-target interaction data file "miR2Target.csv" in csv format
miR2Target <- system.file("extdata", "miR2Target.csv", package="miRspongeR")
miRTarget <- read.csv(miR2Target, header=TRUE, sep=",")

# Identifying sample-specific miRNA sponge interactions, the sppc method is used to # identify miRNA sponge interactions
sponge_sample_specific_net <- sponge_sample_specific(miRTarget, ExpData, senscorcutoff = 0.1, method = "sppc")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

