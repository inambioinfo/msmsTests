\name{msms.glm.qlll}
\alias{msms.glm.qlll}
\title{  Spectral counts differential expression by quasi-likelihood GLM}

\description{
Given a null and an alternative model, with a two level treatment factor
as the two conditions to compare, executes a quasi-likelihood based 
GLM regression to discover differentially expressed proteins between the two
conditions.
The null and alternative models may include blocking factors.The reference 
level of the main factor is considered to be the control condition.
}
\usage{
msms.glm.qlll(msnset,form1,form0,facs=NULL,div=NULL)
}

\arguments{
  \item{msnset}{A MSnSet object with spectral counts in the expression matrix.}
  \item{form1}{ The alternative hypothesis model as an standard R formula,
  with the treatment factor of interest, and eventual blocking factors. }
  \item{form0}{ The null hypothesis model as an standard R formula.It may be 
  the standard null model (y~.) or contain one or multiple blocking factors. }
  \item{facs}{ NULL or a data frame with the factors in its columns. }
  \item{div}{ NULL or a vector with the divisors used to compute the offsets. }
}

\details{ The right hand site of the formulas is expected to be "y~", with
the combination of factors after the tilde. If \code{facs} is NULL the factors 
are taken as default from \code{pData(msnset)}. If \code{div} is NULL all 
divisors are taken equal to one. 
}

\value{
A data frame with the following columns:
\item{LogFC}{Log fold change estimated from the model parameters.}
\item{D}{Residual deviance as statistic of the test.}
\item{p.value}{The p-values obtained from the test.}
}

\references{
Agresti, A. (2002) \emph{Categorical Data Analysis}, 2nd Edition, John Wiley & Sons, Inc., Hoboken, New Jersey

Thompson L.A. (2009) \emph{R (and S-PLUS) Manual to Accompany Agresti s 
Categorical Data Analysis}, 2nd edition
\url{https://home.comcast.net/~lthompson221/Splusdiscrete2.pdf}

Dobson, A.J. (2002) \emph{An Introduction to Generalized Linear Models}, 2nd
Edition, Chapman & Hall/CRC, New York

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S},
New York: Springer

Li, M.; Gray, W.; Zhang, H.; Chung, C. H.; Billheimer, D.; Yarbrough, W. G.; Liebler, D. C.; Shyr, Y.; Slebos, R. J. C. \emph{Comparative shotgun proteomics
using spectral count data and quasi-likelihood modeling}, J Proteome Res 2010,
9, 4295-4305.

Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva 
(2013). \emph{An Effect Size Filter Improves the Reproducibility in Spectral
Counting-based Comparative Proteomics}, Journal of Proteomics, 
http://dx.doi.org/10.1016/j.jprot.2013.05.030
}

\author{
Josep Gregori i Font
}

\seealso{
\code{\link{MSnSet}}, \code{\link{glm}}
}
\examples{
library(msmsTests)
data(msms.dataset)
# Pre-process expression matrix
e <- pp.msms.data(msms.dataset)
# Factors
pData(e)
# Control condition
levels(pData(e)$treat)[1]
# Treatment condition
levels(pData(e)$treat)[2]

# Models and normalizing condition
null.f <- "y~batch"
alt.f <- "y~treat+batch"
div <- apply(exprs(e),2,sum)

#Test
res <- msms.glm.qlll(e,alt.f,null.f,div=div)

str(res)
head(res)
}

\keyword{ univar }
\keyword{ models }
\keyword{ design }
