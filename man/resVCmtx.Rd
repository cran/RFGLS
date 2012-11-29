\name{resVCmtx}
\alias{resVCmtx}
\docType{data}
\title{
 Residual variance-covariance matrix.
}
\description{
The residual covariance matrix used in package examples.  
}

\usage{data(resVCmtx)}

\format{
  An object of class \code{\link[bdsmatrix:bdsmatrix-class]{bdsmatrix}} (from package \emph{bdsmatrix}).  It is a 4050-by-4050 block-diagonal sparse matrix, with off-diagonal elements of 0.  It has six slots, the descriptions of which may be found in the documentation for \code{\link[bdsmatrix:bdsmatrix-class]{bdsmatrix-class}}.
}
  
\details{
It is perhaps most instructive to see the syntax by which \code{resVCmatrix} can be reproduced from datasets \code{pheno}, \code{geno}, and \code{pedigree}:
\preformatted{
data(pheno)
data(geno)
data(pedigree)
foo <- gls.batch.get(
  phenfile=pheno,
  genfile=data.frame(t(geno)),
  pedifile=pedigree, 
  outfile="example_output.txt", 
  covmtxfile.in=NULL,covmtxfile.out=paste(phen,"_cov_matrix.txt",sep=""),
  phen="Zscore", covars = "IsFemale",
  med = "rfgls", sizeLab = "OOPP", Mz = TRUE, Bo = TRUE, Ad = TRUE, Mix = TRUE,
  indobs = TRUE, col.names = TRUE, pediheader = FALSE,
  pedicolname=c("FAMID","ID","PID","MID","SEX"),
  sep.phe = " ", sep.gen = " ", sep.ped = " ")
bar <- fgls(
  fixed=Zscore ~ IsFemale, data=foo$test.dat, sizelist=foo$sizelist, 
  tlist=foo$tlist, 
  sizeLab="OOPP",Mz = TRUE, Bo = TRUE, Ad = TRUE, Mix = TRUE, indobs = TRUE,
  get.hessian=FALSE,vmat = NULL, subset=NULL, weights=NULL, na.action=NULL)
}
Then, \code{bar$sigma} is identical to \code{resVCmtx}.
}

\examples{
data(resVCmtx)
str(resVCmtx)
##Also see examples for functions fgls() and gls.batch().
}

\keyword{datasets}
