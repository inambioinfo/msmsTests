\name{res.volcanoplot}
\alias{res.volcanoplot}

\title{ Volcanoplot }

\description{
Given the data frame obtained from \code{test.results()} a volcano
plot is drawn.The features are colored according to significance and
relevance.
}
\usage{
res.volcanoplot(tres, max.pval=0.05, min.LFC=1, maxx=3, maxy=10,
                ylbls=20)
}

\arguments{
  \item{tres}{
 The dataframe with test results as obtained from \code{test.results()}.
 Or a data frame with, at least, the following columns: \code{LogFC} with
 log fold changes, \code{adjp} with multitest adjusted p-values, and
 \code{DEP} with TRUE or FALSE as post test filter results, being the TRUE
 features both statistically significant and relevant for reproducibility.
}
  \item{max.pval}{
 The maximum adjusted p-value considered as statistically significant.
}
  \item{min.LFC}{
 The minimum absolute log fold change considered as biologically relevant.
}
  \item{maxx}{
 The maximum value in abcissas (i.e. log2(fold change)).
}
  \item{maxy}{
 The maximum value in ordinates (i.e. -log10(p.val))
}
  \item{ylbls}{
 All features with -log10(p.val) above this value will be ploted with feature labels.
}
}

\details{
Abscissas and ordinates may be limited giving a value other than NULL to the
parameters \code{maxx} and \code{maxy}. All features deemed significant and
relevant are ploted by a blue dot, all features deemed significant but not
passing the post test filter are plotted by a red dot. The non-significant
features are plotted as smaller black dots. All features deemed significant
and relevant and with a -log10 p-value above \code{ylbls} are plotted with a
label showing their row index in the test results dataframe.
The borders limiting the values given by \code{max.pval} and \code{min.LFC}
are ploted as dash-and-dot red lines.
}

\value{
No return value.
}
\references{
 Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva (2013). An Effect Size Filter Improves the Reproducibility in Spectral Counting-based Comparative Proteomics. Journal of Proteomics,
 DOI \url{http://dx.doi.org/10.1016/j.jprot.2013.05.030}
}
\author{
Josep Gregori i Font
}

\seealso{
\code{\link{test.results}}, \code{\link[limma]{volcanoplot}}
}

\examples{
library(msmsTests)
data(msms.dataset)
# Pre-process expression matrix
e <- pp.msms.data(msms.dataset)
# Models and normalizing condition
null.f <- "y~batch"
alt.f <- "y~treat+batch"
div <- apply(exprs(e),2,sum)
#Test
res <- msms.glm.qlll(e,alt.f,null.f,div=div)
lst <- test.results(res,e,pData(e)$treat,"U600","U200",div,
                    alpha=0.05,minSpC=2,minLFC=log2(1.8),
                    method="BH")
# Plot
res.volcanoplot(lst$tres, max.pval=0.05, min.LFC=1, maxx=3, maxy=NULL,
                        ylbls=4)
}

\keyword{ hplot }
\keyword{ univar }
\keyword{ htest }

