\name{msms.edgeR}
\alias{msms.edgeR}
\title{ Spectral counts differential expression by edgeR }
\description{
Given a null and an alternative model, with a two level treatment factor 
as the two conditions to compare, executes the negative binomial test by 
edgeR functions to discover differentially expressed proteins between the two conditions. The null and alternative models may include blocking
factors.The reference level of the main factor is considered to be the control condition
}

\usage{
msms.edgeR(msnset,form1,form0,facs=NULL,div=NULL,fnm=NULL)
}
\arguments{
  \item{msnset}{A MSnSet object with spectral counts in the expression matrix.}
  \item{form1}{ The alternative hypothesis model as an standard R formula,
  with the treatment factor of interest, and eventual blocking factors. }
  \item{form0}{ The null hypothesis model as an standard R formula.It may be 
  the standard null model (y~.) or contain one or multiple blocking factors. }
  \item{facs}{ NULL or a data frame with the factors in its columns. }
  \item{div}{ NULL or a vector with the divisors used to compute the offsets. }
  \item{fnm}{ NULL or a character string with the treatment factor name, as used in the column names of the factors data frame, and in the formula. }
}

\details{ The right hand site of the formulas is expected to be "y~", with
the combination of factors after the tilde. If \code{facs} is NULL the factors 
are taken as default from \code{pData(msnset)}. If \code{div} is NULL all 
divisors are taken equal to one. If \code{fnm} is NULL it is taken to be the 
first factor in \code{facs}.
}

\value{ A data frame with column names 'LogFC', 'LR', 'p.value', with the
estimated log fold changes, likelihood ratio statistic and corresponding
p-value as obtaimed from a call to glmLRT() from the edgeR package. }

\references{
Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor
  package for differential expression analysis of digital gene
  expression data. Bioinformatics 26, 139-140
  
Robinson MD and Smyth GK (2007). Moderated statistical tests for
  assessing differences in tag abundance. Bioinformatics 23, 2881-2887

Robinson MD and Smyth GK (2008). Small-sample estimation of negative
  binomial dispersion, with applications to SAGE data. Biostatistics,
  9, 321-332
  
Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva (2013).
 An Effect Size Filter Improves the Reproducibility in Spectral Counting-based
 Comparative Proteomics. Journal of Proteomics, 
 DOI \url{http://dx.doi.org/10.1016/j.jprot.2013.05.030}
}

\author{
Josep Gregori i Font
}

\seealso{
\code{\link{MSnSet}}, \code{\link{edgeR}}, \code{\link{glmLRT}}, \code{\link{msmsEDA}} 
}

\examples{
## Example
library(msmsTests)
data(msms.dataset)
e <- pp.msms.data(msms.dataset)
e

null.f <- "y~batch"
alt.f <- "y~treat+batch"
div <- apply(exprs(e),2,sum)
res <- msms.edgeR(e,alt.f,null.f,div=div,fnm="treat")

str(res)
head(res)
}

\keyword{ univar }
\keyword{ models }
\keyword{ design }
