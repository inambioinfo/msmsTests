\name{test.results}
\alias{test.results}

\title{ Multitest p-value adjustment and post-test filter  }

\description{
Operates on the statistic tests results obtained from \code{msms.glm.pois()},
\code{msms.glm.qlll()} or \code{msms.edgeR()}. The following variables are
computed: Raw expression mean values for each condition (control and treatment),
log fold change based on these expression levels and taking into account the
normalizing divisors (\code{div}), multitest adjusted p-values with FDR control,
and a post test filter based on minimum spectral counts and minimum absolute log
fold change as estimated by the statistic test. According to the results of this
post-test filter, features are flagged as T or F depending on whether they
result relevant or not, beyond their statistic signicance.  
}
\usage{
test.results(test, msnset, gpf, gp1, gp2, div, alpha = 0.05, 
             minSpC = 2, minLFC = 1, method = "BH")
}

\arguments{
  \item{test}{
  The dataframe obtained from either \code{msms.glm.pois()},
  \code{msms.glm.qlll()} or \code{msms.edgeR()}
  }
  \item{msnset}{A MSnSet object with spectral counts in the expression matrix.}
  \item{gpf}{The factor used in the tests.}
  \item{gp1}{The treatment level name.}
  \item{gp2}{
  The control level name. Should be the factor's reference level.
  See R function \code{relevel}.
}
  \item{div}{
 The weights used as divisors (offsets) in the GLM model. Usually the sum of
 spectral counts of each sample.
}
  \item{alpha}{
 The multi test adjusted p-value significance threshold.
}
  \item{minSpC}{
 The minimum spectral counts considered as relevant in the most abundant
 condition. This filter aims at reaching good reproducibility.
}
  \item{minLFC}{
 The minimum absolute log fold change considered both, relevant and biologically
 significant. This filter aims at assuring enough biological effect size and at
 reaching good reproducibility.
}
  \item{method}{
 One among \code{BH} or \code{qval}. The p-values are FDR ajusted by the 
 Benjamini-Hochberg method (\code{BH}) or by qvalue (\code{qval}).
}
}
\details{
No feature is removed in the filter, but instead they are flagged as T or F
depending on whether they are considered as differentially expressed or not, in the \code{DEP} column, taking into account statistic significance and 
reproducibility metrics.
}

\value{
A data frame with the following columns:
\item{first column}{ Column named as the treatment level with the mean raw
 spectral counts observed for this condition}
\item{second column}{ Column named as the control level with the mean raw
 spectral counts observed for this condition}
\item{lFC.Av}{ Log fold change computed from the mean expression levesl taking
 into account the given normalization factors.}
\item{logFC}{ Log fold change estimated by fitting the given GLM model. 
The reference level of the main factor is taken as control. }
\item{D or LR}{ The statistic obtained from the tests. The residual deviance
 \code{D} for Poisson and quasi-likelihood, or the likelihood ratio \code{LR}
 for edgeR.}
\item{p.val}{ The unadjusted p-values obtained from the tests. }
\item{adjp}{ The multitest adjusted p-values with FDR control. }
\item{DEP}{ A character type column with T or F flagging the features
 considered both as statistically significant and relevant for reproducibility.}

}

\references{

 Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva 
 (2013). An Effect Size Filter Improves the Reproducibility in Spectral
 Counting-based Comparative Proteomics. Journal of Proteomics, 
 DOI \url{http://dx.doi.org/10.1016/j.jprot.2013.05.030}
 
 Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery 
 rate: a  practical and powerful approach to multiple testing. Journal of 
 the Royal Statistical Society Series B, 57, 289-300. 

 Alan Dabney, John D. Storey and with assistance from Gregory R.
 Warnes. qvalue: Q-value estimation for false discovery rate
 control. R package version 1.30.0.
}

\author{
Josep Gregori i Font
}

\seealso{
 \code{\link{pval.by.fc}},
 \code{\link{p.adjust}},
 \code{\link{qvalue}},
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
# Post-test filter
lst <- test.results(res,e,pData(e)$treat,"U600","U200",div,
                    alpha=0.05,minSpC=2,minLFC=1,
                    method="BH")
str(lst)
lst$cond
head(lst$tres)
rownames(lst$tres)[which(lst$tres$DEP=='T')]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
