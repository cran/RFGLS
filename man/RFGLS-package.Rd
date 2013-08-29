\name{RFGLS-package}
\alias{RFGLS-package}
\alias{RFGLS}
\docType{package}
\title{
Rapid Feasible Generalized Least Squares
}

\description{
\emph{RFGLS} uses a generalized least-squares method to perform single-marker association analysis, in datasets of nuclear families containing parents, twins, and/or adoptees.  It is designed for families of no greater than four members.  When conducting association analysis with a large number of markers, as in GWAS, \emph{RFGLS} uses \emph{rapid} feasible generalized least-squares, an approximation to feasible generalized least-squares (FGLS) that considerably reduces computation time with minimal bias in \emph{p}-values, and with negligible loss in power.
  
The package includes four functions.  Function \code{\link{gls.batch}()} actually conducts GWAS using the rapid feasible generalized least-squares approximation, under which the residual variance-covariance matrix is estimated once from a regression of the phenotype onto covariates only, and is subsequently "plugged in" for use in all subsequent single-SNP analyses.  Function \code{\link{fgls}()} is called by \code{\link{gls.batch}()}, and conducts a single FGLS regression.  It can be used to simultaneously estimate fixed-effects regression coefficients and the residual covariance matrix.  Function \code{\link{gls.batch.get}()} is useful to restructure data, for use with \code{\link{fgls}()}. Function \code{\link{FSV.frompedi}()} creates family-structure variables based upon available information in a pedigree file.  Functions \code{\link{gls.batch}()} and \code{\link{gls.batch.get}()} use these family-structure variables, which represent the type of family to which each participant belongs and how s/he fits into that family.
%Reference: Li X, Basu S, Miller MB, Iacono WG, McGue M:
%A Rapid Generalized Least Squares Model for a Genome-Wide Quantitative Trait Association Analysis in Families.
%Hum Hered 2011;71:67-82 (DOI: 10.1159/000324839)
}

\references{
Li X, Basu S, Miller MB, Iacono WG, McGue M:
A Rapid Generalized Least Squares Model for a Genome-Wide Quantitative Trait Association Analysis in Families.
\emph{Human Heredity} 2011;71:67-82 (DOI: 10.1159/000324839)
}

\details{
\tabular{ll}{
Package: \tab RFGLS\cr
Version: \tab 1.1\cr
Date: \tab 2013/8/29\cr
Depends: \tab R (>= 2.15.0), stats, bdsmatrix, Matrix\cr
License: \tab GPL (>= 2)\cr
}}

\author{
Robert M. Kirkpatrick \email{rkirkpatrick2@vcu.edu}, Xiang Li \email{lixxx554@umn.edu}, and Saonli Basu \email{saonli@umn.edu} .

Maintainer: Saonli Basu <saonli@umn.edu>
}

\keyword{ package }
