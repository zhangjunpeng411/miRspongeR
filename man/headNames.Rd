\name{headNames}
\alias{headNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get data header from .csv dataset
}
\description{
Get data header from .csv dataset
}
\usage{
headNames(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{dataset}{
    the input dataset in .csv format.
}
}

\value{
A vector of header names in .csv dataset
}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
}
\examples{
ExpDatacsv <- system.file("extdata", "ExpData.csv", package="miRsponge")
Header_names <- headNames(ExpDatacsv)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


