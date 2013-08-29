\name{rescovmtx}
\alias{rescovmtx}
\alias{resVCmtx}
\docType{data}
\title{
 Residual variance-covariance matrix.
}
\description{
The residual covariance matrix used in package examples.  In the previous package version, its name was "resVCmtx".
}

\usage{data(rescovmtx)}

\format{
  An object of class \code{\link[bdsmatrix:bdsmatrix-class]{bdsmatrix}} (from package \emph{bdsmatrix}).  It is a 4050-by-4050 block-diagonal sparse matrix, with off-diagonal elements of 0.  It has six slots, the descriptions of which may be found in the documentation for \code{\link[bdsmatrix:bdsmatrix-class]{bdsmatrix-class}}.
}
  
\details{
It is perhaps most instructive to see the syntax by which \code{rescovmtx} can be reproduced from datasets \code{pheno}, \code{geno}, and \code{pedigree}:
\preformatted{
data(pheno)
data(geno)
data(pedigree)
foo <- gls.batch.get(
  phenfile=pheno,genfile=data.frame(t(geno)),pedifile=pedigree,
  covmtxfile.in=NULL,theta=NULL,snp.names=NULL,input.mode=c(1,2,3),
  pediheader=FALSE,pedicolname=c("FAMID","ID","PID","MID","SEX"),
  sep.phe=" ",sep.gen=" ",sep.ped=" ",
  phen="Zscore",covars="IsFemale",med=c("UN","VC"),
  outfile,col.names=TRUE,return.value=FALSE,
  covmtxfile.out=NULL,
  covmtxparams.out=NULL,
  sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL)
bar <- fgls(
  Zscore ~ IsFemale, data=foo$test.dat, tlist=foo$tlist,
  sizelist=foo$sizelist,med=c("UN","VC"), 
  vmat=NULL, 
  start=NULL, theta=NULL, drop=NULL, get.hessian=FALSE, 
  optim.method="BFGS", control=list(), weights=NULL,
  sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL)
}
Then, \code{bar$sigma} is identical to \code{rescovmtx}.
}

\examples{
data(rescovmtx)
str(rescovmtx)
##Also see examples for functions fgls() and gls.batch().
}

\keyword{datasets}
