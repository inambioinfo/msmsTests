\name{pval.by.fc}
\alias{pval.by.fc}

\title{ Table of cumulative frequencies of p-values by log fold change bins
}

\description{
Given the set of p-values and log fold changes that result from a test,
computes a table of cumulative frequencies of features by p-values in bins
of log fold changes.
}

\usage{
pval.by.fc(pvals,lfc)
}

\arguments{
  \item{lfc}{
The log fold changes estimated from the tests.
}
  \item{pvals}{
The p-values, adjusted or not, obtained from the tests.
}
}

\value{
A matrix of cumulated frequencies with descriptive row and column names.
}

\references{
 Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva (2013). An Effect Size Filter Improves the Reproducibility in Spectral Counting-based Comparative Proteomics. Journal of Proteomics,
 DOI \url{http://dx.doi.org/10.1016/j.jprot.2013.05.030}
}

\author{
Josep Gregori i Font
}

\seealso{
\code{\link{test.results}}
}

\examples{
library(msmsTests)
data(msms.spk)
# Subset
treat <- pData(msms.spk)
jdx <- which(treat=="U200" | treat=="U600")
e <- msms.spk[,jdx]
pData(e)$treat <- treat[jdx,1,drop=TRUE]
# Pre-process expression matrix
e <- pp.msms.data(e)
# Models and normalizing condition
null.f <- "y~1"
alt.f <- "y~treat"
div <- apply(exprs(e),2,sum)
#Test
res <- msms.glm.pois(e,alt.f,null.f,div=div)
# Post-test filter
lst <- test.results(res,e,pData(e)$treat,"U600","U200",div,
                    alpha=0.05,minSpC=2,minLFC=1,
                    method="BH")

##  On all features, with multitest adjusted p-values
pval.by.fc(lst$tres$adjp, lst$tres$LogFC)

### On all features deemed significant and biologically relevant
flags <- lst$tres$DEP=="T"
pval.by.fc(lst$tres$adjp[flags], lst$tres$LogFC[flags])
}
\keyword{ htest }
\keyword{ univar }

