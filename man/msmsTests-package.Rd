\name{msmsTests-package}
\alias{msmsTests-package}
\alias{msmsTests}
\docType{package}
\title{
LC-MS/MS Differential Expression Tests
}
\description{
Statistical tests for label-free LC-MS/MS data by spectral counts, to discover
differentially expressed proteins between two biological conditions. Three tests are available: Poisson GLM regression, quasi-likelihood GLM regression, and the negative binomial of the edgeR package.
The three models admit blocking factors to control for nuissance variables.
To assure a good level of reproducibility a post-test filter is available, where we may set the minimum effect size considered biologicaly relevant, and the minimum expression of the most abundant condition.
}
\details{
\tabular{ll}{
Package: \tab msmsTests\cr
Type: \tab Package\cr
Version: \tab 0.99.1\cr
Date: \tab 2013-07-26\cr
License: \tab GPL-2\cr
}

\tabular{ll}{
\code{msms.glm.pois}: \tab { Poisson based GLM regression} \cr
\code{msms.glm.qlll}: \tab { Quasi-likelihood GLMregression} \cr
\code{msms.edgeR}: \tab {The binomial negative of edgeR} \cr
\code{pval.by.fc}: \tab {Table of cumulative frequencies of features by p-values
    in bins of log fold change} \cr
\code{test.results}: \tab {Multitest p-value adjustement and post-test filter} \cr
\code{res.volcanoplot}: \tab {Volcanplot of the results} 
}
}

\author{
Josep Gregori, Alex Sanchez, and Josep Villanueva \cr
Maintainer: Josep Gregori <josep.gregori@gmail.com>
}
\references{
 Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva 
 (2013). An Effect Size Filter Improves the Reproducibility in Spectral
 Counting-based Comparative Proteomics. Journal of Proteomics,
 DOI \url{http://dx.doi.org/10.1016/j.jprot.2013.05.030}
 }
\keyword{ htest }
\keyword{ hplot }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}
