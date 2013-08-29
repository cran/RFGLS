\name{pedigree}
\alias{pedigree}
\docType{data}
\title{
 Pedigree table
}
\description{
A pedigree table for the same simulees as in dataset \code{\link{pheno}}.
}
\usage{data(pedigree)}
\format{
  A data frame with 4050 observations on the following 6 integer-valued variables:
  \describe{
    \item{\code{FAMID}}{"Family ID."  Each family in the dataset is uniquely identified by a value of \code{FAMID}, which are all multiples of 10.}
    \item{\code{ID}}{\emph{Individual} ID.  Each subject in the dataset is uniquely identified by a value of \code{ID}.}
    \item{\code{PID}}{"Paternal ID."  Coded 0 for founders (parents, adoptees, and "independent observations.")}
    \item{\code{MID}}{"Maternal ID."  Coded 0 for founders (parents, adoptees, and "independent observations.")}
    \item{\code{SEX}}{Coded 1 for male and 2 for female.}
  }}
\details{
 Merely a pedigree table in a commonly used format.  Note that its column names are the default names that \code{\link{gls.batch}()} or \code{\link{gls.batch.get}()} assign to the pedigree file.  However, the only column that those two functions \emph{strictly} require is \code{ID}; see examples below.
}
\examples{
data(pedigree)
data(pheno)
data(geno)
data(map)

foo <- gls.batch.get(
  phenfile=pheno,genfile=data.frame(t(geno)),pedifile=pedigree,
  covmtxfile.in=NULL,theta=NULL,snp.names=map[,2],
  input.mode=1,
  pediheader=FALSE,pedicolname=c("FAMID","ID","PID","MID","SEX"),
  sep.phe=" ",sep.gen=" ",sep.ped=" ",
  phen="Zscore",covars="IsFemale",med=c("UN","VC"),
  outfile,col.names=TRUE,return.value=FALSE,
  covmtxfile.out=NULL,
  covmtxparams.out=NULL,
  sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL)
str(foo)

##Also works, since phenfile provides 'FTYPE' and 'INDIV',
##and input.mode=1:
pedigree2 <- pedigree
pedigree2[,-2] <- NA    ##<--Change all but column 'ID' to NA.
foo2 <- gls.batch.get(
  phenfile=pheno,genfile=data.frame(t(geno)),
  pedifile=pedigree2,                       ##<--Note change.
  covmtxfile.in=NULL,theta=NULL,snp.names=map[,2],
  input.mode=1, ##<-- =2 or =3 would need more pedifile columns
  pediheader=FALSE,pedicolname=c("FAMID","ID","PID","MID","SEX"),
  sep.phe=" ",sep.gen=" ",sep.ped=" ",
  phen="Zscore",covars="IsFemale",med=c("UN","VC"),
  outfile,col.names=TRUE,return.value=FALSE,
  covmtxfile.out=NULL,
  covmtxparams.out=NULL,
  sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL)
str(foo2)

}
\keyword{datasets}
