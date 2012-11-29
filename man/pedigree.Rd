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
    \item{\code{SEX}}{Binary indicator; a value of 1 indicates female sex.}
  }}
\details{
 Merely a pedigree table in a commonly used format.  Note that its column names are the default names that \code{\link{gls.batch}()} or \code{\link{gls.batch.get}()} assign to the pedigree file.  However, the only column that those two functions require is \code{ID}; see examples below.
}
\examples{
data(pedigree)
data(pheno)
data(geno)

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
str(foo)

##Also works:
pedigree2 <- pedigree
pedigree2[,-2] <- NA    ##<--Change all but column 'ID' to NA.
foo2 <- gls.batch.get(
  phenfile=pheno,
  genfile=data.frame(t(geno)),
  pedifile=pedigree2,             ##<--Note change.
  outfile="example_output.txt", 
  covmtxfile.in=NULL,covmtxfile.out=paste(phen,"_cov_matrix.txt",sep=""),
  phen="Zscore", covars = "IsFemale",
  med = "rfgls", sizeLab = "OOPP", Mz = TRUE, Bo = TRUE, Ad = TRUE, Mix = TRUE,
  indobs = TRUE, col.names = TRUE, pediheader = FALSE,
  pedicolname=c("FAMID","ID","PID","MID","SEX"),
  sep.phe = " ", sep.gen = " ", sep.ped = " ")

}
\keyword{datasets}
