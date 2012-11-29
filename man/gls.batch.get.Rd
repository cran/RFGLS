\name{gls.batch.get}
\alias{gls.batch.get}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data restructuring for \code{\link{fgls}()}.
}
\description{
     Carries out the data restructuring performed by \code{\link{gls.batch}()}, before it estimates the residual covariance matrix.  Useful if calling \code{\link{fgls}()} directly.
}
\usage{
gls.batch.get(phenfile,genfile,pedifile,outfile,covmtxfile.in=NULL,
  covmtxfile.out=paste(phen,"_cov_matrix.txt",sep=""),phen,covars=NULL,
  med="rfgls",sizeLab="OOPP",Mz=TRUE,Bo=TRUE,Ad=TRUE,Mix=TRUE,
  indobs=TRUE,col.names=TRUE,pediheader=FALSE,
  pedicolname=c("FAMID","ID","PID","MID","SEX"),
  sep.phe=" ",sep.gen=" ",sep.ped=" ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phenfile}{
 This can be either (1) a character string specifying a phenotype file on disk which includes the phenotypes and other covariates, or (2) a data frame object containing the same data.  In either case, the data must be appropriately structured.    See below under "Details."
}
  \item{genfile}{
 This can be either (1) a character string specifying a genotype file of genotype scores (such as 0,1,2, for the additive genetic model) to be read from disk, or (2) a data frame object containing them.  In such a file, each row must represent a SNP, each column must represent a subject, and there should NOT be column headers or row numbers.  In such a data frame, the reverse holds: each row must represent a subject, and each column, a SNP.  If the data frame--say, \code{geno}--need be transposed, then use \code{genfile=data.frame(t(geno))}.  Using a matrix instead of a data frame is not recommended, because it makes the process of merging data very memory-intensive, unless the sample size is quite small.
}
  \item{pedifile}{
This can be either (1) a character string specifying the pedigree file corresponding to \option{genfile}, to be read from disk, or (2) a data frame object containing this pedigree information.  At minimum, \option{pedifile} must have a column of subject IDs, named \code{'ID'}, ordered in the same order as subjects' genotypic data in \option{genfile}.  Every row in \option{pedifile} is matched to a SNP in \option{genfile}.  That is, if reading files from disk (which is recommended), each row \emph{i} of the pedigree file, which has \emph{n} rows, matches column \emph{i} of the genotype file, which has \emph{n} columns.  This is how the program matches subjects in the phenotype file to their genotypic data.

 The pedigree file or data frame  can also include other columns of pedigree information, like father's ID, mother's ID, etc.  Argument \option{pediheader} (see below) is an indicator of whether the pedigree file on disk has a header or not, with default as \code{FALSE}.  Argument \option{pedicolnames} (see below) gives the names that \command{gls.batch.get()} will assign to the columns of \option{pedifile}, and the default, \code{c("FAMID","ID","PID","MID","SEX")}, is the familiar "pedigree table" format.  In any event, the user's input must somehow provide the program with a column of IDs, labeled as \code{'ID'}.
}
  \item{phen}{
A character string specifying the phenotype (column name) in the phenotype file to be analyzed.
}
  \item{covars}{
A character string or character vector that holds the (column) names of the covariates, in the phenotype file, to be used in the regression model.
}
  \item{pediheader}{
  A logical indicator specifying whether the pedigree file to be read from disk has a header row, to ensure it is read in correctly.  Even if \code{TRUE}, \command{gls.batch()} assigns the values in \option{pedicolname} to the columns after it has been read in.  Defaults to \code{FALSE}.  Also see \option{pedifile} above and under "Details" below.
}
  \item{pedicolname}{
 A vector of character strings giving the column names that \command{gls.batch()} will assign to the columns of the pedigree file. The default,\cr \code{c("FAMID","ID","PID","MID","SEX")}, is the familiar "pedigree table" format.  The two criteria this vector must have are that it must (1) assign the name "ID" to the column of subject IDs in the pedigree file, and (2) its length must equal the number of columns of the pedigree file.  Also see \option{pedifile} above, and under "Details" below.
}
  \item{sep.phe}{
 Separator character of the phenotype file to be read from disk.  Defaults to a single space.
}
  \item{sep.gen}{
 Separator character of the genotype file to be read from disk.  Defaults to a single space.
}
  \item{sep.ped}{
 Separator character of the pedigree file to be read from disk.  Defaults to a single space.
}
  \item{covmtxfile.in, covmtxfile.out, med, outfile}{
These arguments are accepted but not used, in order for \code{\link{gls.batch.get}()} to parallel \code{\link{gls.batch}()} as closely as possible.
}
  \item{sizeLab, Mz, Bo, Ad, Mix, indobs, col.names}{
These arguments are likewise accepted but not used.  
}
}
\details{
  Though originally used for debugging purposes, \command{gls.batch.get()} was included because it facilitates directly invoking \code{\link{fgls}()} when the need arises.  This function first reads in the files and merges the files into a data frame with columns of pedigree information, phenotypes, covariates, and genotypes.  It then creates a \option{tlist} vector and a \option{sizelist} vector, which comprise the family labels and family sizes in the data.
 It returns a list containing the merged data frame, and the \option{tlist} and \option{sizelist} vectors.
 
The phenotype file must conform to the following guidelines:
\itemize{
  \item It must have the following four named columns: \code{'FAMID'} (family ID), \code{'ID'} (\emph{unique} individual ID), \code{'FTYPE'}  (family type), and \code{'INDIV'} (individual code).  The value of \code{FTYPE} and \code{FAMID} will be the same for all members of a given family.  There are six recognized family types: \code{FTYPE=1} for MZ-twin, \code{FTYPE=2} for DZ-twin, \code{FTYPE=3} for adoptive-offspring, \code{FTYPE=4} for non-twin bio-offspring, \code{FTYPE=5} for "mixed" families with one bio and one adopted offspring, and \code{FTYPE=6} for "independent observations" who do not fit into a four-person nuclear family.  The individual code \code{INDIV} represents how the subject fits into his/her family: \code{INDIV=1} is for "Offspring #1," \code{INDIV=2} is for "Offspring #2," \code{INDIV=3} is for the mother, and \code{INDIV=4} is for the father.  The distinction between "Offspring #1" and "#2" is mostly arbitrary, except that in "mixed" families, the biological offspring MUST have \code{INDIV=1}, and the adopted offspring, \code{INDIV=2}.
  \item Within each family, members must be ordered by \code{INDIV}, as: offspring, mother, father.  For mixed family type, members must be ordered as: bio-offspring, adopted-offspring, mother, father.
  \item The phenotype file has rows as subjects and columns as variables, whereas the genotype file provided to \option{genfile} must have rows as SNPs and columns as subjects.} 
}
\value{
A list with these three components:
  \item{test.dat}{
The merged data frame of pedigree information, phenotypes, covariates, and genotypes.
}
  \item{tlist}{
A vector of family labels, with length equal to the number of families in the data (each "independent observation" is treated as a separate family).
The names of its components are the family IDs.
}
  \item{sizelist}{
A vector of family sizes, with length equal to the number of families in the data (each "independent observation" is treated as a separate family).
The names of its components are the family IDs.
}
}
\author{
Xiang Li <lixxx554@umn.edu>, Robert M. Kirkpatrick <kirk0191@umn.edu>, and Saonli Basu <saonli@umn.edu>.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{fgls}}, \code{\link{gls.batch.get}}
%%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
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
olsmod <- lm(   ##<--OLS regression could be applied to the merged dataset...
    Zscore ~ snp.1 + IsFemale, data=foo$test.dat)
summary(olsmod)  #<--...but the standard errors and t-statistics will not be valid.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
