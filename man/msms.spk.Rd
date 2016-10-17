\name{msms.spk}
\alias{msms.spk}
\docType{data}
\title{
Yeast lisate samples spiked with human proteins
}
\description{
A MSnSet with a spectral counts in the expression matrix and a treatment factor
in the phenoData slot. \cr
The spectral counts matrix has samples in the columns, and proteins in the rows.
Each sample consists in 500ng of standard yeast lisate spiked with 100, 200, 400 and 600fm of a mix of 48 equimolar human proteins (UPS1, Sigma-Aldrich). The dataset contains a different number of technical replicates of each sample.
}
\usage{
data(msms.spk)
}
\format{
  A MSnSet
}
\references{
Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva (2013). An Effect Size Filter Improves the Reproducibility in Spectral Counting-based Comparative Proteomics. Journal of Proteomics, 
DOI \url{http://dx.doi.org/10.1016/j.jprot.2013.05.030}

Laurent Gatto and Kathryn S. Lilley, MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation, Bioinformatics 28(2), 288-289 (2012).
}

\seealso{
See \code{\link{MSnSet}} for detail on the class, and the \code{exprs} and \code{pData} accessors.
}

\examples{
data(msms.spk)
msms.spk
dim(msms.spk)
table(pData(msms.spk))
head(exprs(msms.spk))
}
\keyword{datasets}
