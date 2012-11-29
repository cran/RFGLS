\name{gls.batch}
\alias{gls.batch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized least-squares batch analysis.
}
\description{
     Fits a generalized least-squares regression model to test association between a quantitative phenotype and all SNPs in a genotype file, one at a time, via Rapid Feasible Generalized Least Squares. For each SNP, genotype is treated as a fixed     effect, and the residual variance-covariance matrix is also estimated.  In each trait-SNP association test, the \code{\link{fgls}()} function is used for parameter estimation.
}
\usage{
gls.batch(phenfile,genfile,pedifile,outfile,covmtxfile.in=NULL,
  covmtxfile.out=paste(phen,"_cov_matrix.txt",sep=""),phen,covars=NULL,
  med="rfgls",  sizeLab="OOPP",Mz=TRUE,Bo=TRUE,Ad=TRUE,Mix=TRUE,
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
 This can be either (1) a character string specifying a genotype file of genotype scores (such as 0,1,2, for the additive genetic model) to be read from disk, or (2) a data frame object containing them.  In such a file, each row must represent a SNP, each column must represent a subject, and there should NOT be column headers or row numbers.  In such a data frame, the reverse holds: each row must represent a subject, and each column, a SNP (e.g. \code{\link{geno}}).  If the data frame--say, \code{geno}--need be transposed, then use \code{genfile=data.frame(t(geno))}.  Using a matrix instead of a data frame is not recommended, because it makes the process of merging data very memory-intensive, and will likely overflow \R's workspace unless the sample size is quite small.
 
 \emph{Warning:} If \option{genfile} is a data frame, \command{gls.batch} will attempt to remove it from \R's workspace after loading it.  Therefore, it should be a standalone object, and not be part of a list.

 Note that genotype scores need not be integers; they merely need to be numeric.  So, \command{gls.batch()} can be used to analyze imputed dosages, etc.
}
  \item{pedifile}{
 This can be either (1) a character string specifying the pedigree file corresponding to \option{genfile}, to be read from disk, or (2) a data frame object containing this pedigree information.  At minimum, \option{pedifile} must have a column of subject IDs, named \code{'ID'}, ordered in the same order as subjects' genotypic data in \option{genfile}.  Every row in \option{pedifile} is matched to a SNP in \option{genfile}.  That is, if reading files from disk (which is recommended), each row \emph{i} of the pedigree file, which has \emph{n} rows, matches column \emph{i} of the genotype file, which has \emph{n} columns.  This is how the program matches subjects in the phenotype file to their genotypic data.

 The pedigree file or data frame  can also include other columns of pedigree information, like father's ID, mother's ID, etc.  Argument \option{pediheader} (see below) is an indicator of whether the pedigree file on disk has a header or not, with default as \code{FALSE}.  Argument \option{pedicolnames} (see below) gives the names that \command{gls.batch()} will assign to the columns of \option{pedifile}, and the default, \code{c("FAMID","ID","PID","MID","SEX")}, is the familiar "pedigree table" format.  In any event, the user's input must somehow provide the program with a column of IDs, labeled as \code{'ID'}.
}
  \item{outfile}{
 A character string specifying the path and filename for the output file to be written.  If a file with the same path and filename already exists, \command{gls.batch()} appends the output to that file, rather than overwriting it. Users are \emph{warned} to be sure that the specified directory exists and can be written to!
}
  \item{covmtxfile.in}{
  Optional; can be either (1) a character string specifying a file on disk from which the residual variance-covariance matrix is to be read, or (2) the matrix itself.  If \code{NULL}, then \command{gls.batch()} will estimate this matrix.  The file to be read in must be a single column, with a header, containing the contents of the 'blocks' of an object of class \code{\link[bdsmatrix:bdsmatrix-class]{bdsmatrix}}; no other file structures are presently compatible.  If \option{covmtxfile.in} is an actual matrix object, then using one of class \code{\link[bdsmatrix:bdsmatrix-class]{bdsmatrix}} is a virtual requirement.  See below under "Details" for more information.
}
  \item{covmtxfile.out}{
 An optional character string specifying the filename and path to which the residual variance-covariance matrix, if it is to be calculated (i.e., \code{covmtxfile.in=NULL}), will be written.  The default is a generic filename that refers to the phenotype (argument \option{phen}, below), to be written to \R's working directory.  If \code{NULL}, then no file is written.  See below under "Details" for more information.
}
  \item{phen}{
 A character string specifying the phenotype (column name) in the phenotype file to be analyzed.
}
  \item{covars}{
 A character string or character vector that holds the (column) names of the covariates, in the phenotype file, to be used in the regression model.
}
  \item{med}{
 "Method."  Presently, only \code{"rfgls"}, the default, is implemented.  
}
  \item{sizeLab}{
 A character string indicating the maximum size of the families in the data.  Must be one of the following strings:
\itemize{
  \item \code{"OOPP"}, if the largest family has two offspring and both parents;
  \item \code{"OPP"}, if the largest family has 1 offspring and both parents (a parent-child trio);
  \item \code{"OO"}, if there are no parents in the data; 
  \item \code{"PP"}, if there are no offspring in the data;
  \item \code{"OP"}, if the largest family has 1 offspring and 1 parent.
  }
  The default is the largest, \code{"OOPP"}.
}
  \item{Mz}{
 Logical (\code{TRUE} or \code{FALSE}).  An indicator of whether Mz-twin families are in the data; must be set to \code{FALSE} if \code{sizeLab="PP"}.  Defaults to \code{TRUE}.
}
  \item{Bo}{
 A logical indicator of whether bio-offspring (including DZ-twin) families are in the data; must be set to \code{FALSE} if \command{sizeLab = }\code{"PP"}.  Defaults to \code{TRUE}.
}
  \item{Ad}{
 A logical indicator of whether adopted-offspring families are in the data; must be set to \code{FALSE} if \command{sizeLab = }\code{"PP"}.  Defaults to \code{TRUE}.
}
  \item{Mix}{
 A logical indicator of whether "mixed" families, with 1 biological and 1 adopted offspring, are in the data; must be set to \code{FALSE} if \command{sizeLab = }\code{"PP"}.  Defaults to \code{TRUE}.
}
  \item{indobs}{
 A logical indicator of whether there are "independent observations" who do not fit into a four-person nuclear family present in the data.  If \code{TRUE}, a separate residual variance parameter will be estimated for those individuals.
}
  \item{col.names}{
 A logical indicator specifying whether to write column names in the output file.  Defaults to \code{TRUE}.
}
  \item{pediheader}{
 A logical indicator specifying whether the pedigree file to be read from disk has a header row, to ensure it is read in correctly.  Even if \code{TRUE}, \command{gls.batch()} assigns the values in \option{pedicolname} to the columns after it has been read in.  Defaults to \code{FALSE}.  Also see \option{pedifile} above and under "Details" below.
}
  \item{pedicolname}{
 A vector of character strings giving the column names that \command{gls.batch()} will assign to the columns of the pedigree file. The default, \cr \code{c("FAMID","ID","PID","MID","SEX")}, is the familiar "pedigree table" format.  The two criteria this vector must have are that it must (1) assign the name "ID" to the column of subject IDs in the pedigree file, and (2) its length must equal the number of columns of the pedigree file.  Also see \option{pedifile} above, and under "Details" below.
}
  \item{sep.phe}{
 Separator character of the phenotype file to be read from disk.  Defaults to a single space.
}
  \item{sep.gen}{
 Separator character of the genotype file to be read from disk.  Defaults to a single space.
}
  \item{sep.ped}{
 Separator character of the pedigree file.  Defaults to a single space.
}
}
\details{
 Reference is frequently made throughout this documentation to the "phenotype file," the "genotype file," and so forth, because \command{gls.batch()} was intended to be used with potentially large datafiles to be read from disk.  This should be evident from the presence of the word "file" in the names of many of this function's arguments, and the fact that all of those arguments may be character strings providing a filename and path.  However, it can also accept the data if the file has already been loaded into \R's workspace as a data frame object, in which case "the [whatever] file" should be taken to refer to such a data frame.    For details specific to each argument, see above.
 
 The function \command{gls.batch()} first reads in the files and merges them into a data frame with columns of pedigree information, phenotypes, covariates, and genotypes.  Then, it creates a \option{tlist} vector and a \option{sizelist} vector, which comprise the family labels and family sizes in the data.  Finally, it carries out single-SNP association analyses for all the SNPs in the genotype file.
 
 The phenotype file must conform to the following guidelines:
\itemize{
  \item It must have the following four named columns: \code{'FAMID'} (family ID), \code{'ID'} (\emph{unique} individual ID), \code{'FTYPE'}  (family type), and \code{'INDIV'} (individual code).  The value of \code{"FTYPE"} and \code{"FAMID"} will be the same for all members of a given family.  There are six recognized family types: \code{FTYPE=1} for MZ-twin, \code{FTYPE=2} for DZ-twin, \code{FTYPE=3} for adoptive-offspring, \code{FTYPE=4} for non-twin bio-offspring, \code{FTYPE=5} for "mixed" families with one bio and one adopted offspring, and \code{FTYPE=6} for "independent observations" who do not fit into a four-person nuclear family.  The individual code \code{"INDIV"} represents how the subject fits into his/her family: \code{INDIV=1} is for "Offspring #1," \code{INDIV=2} is for "Offspring 2," \code{INDIV=3} is for the mother, and \code{INDIV=4} is for the father.  Note that subjects with \code{FTYPE=6} MUST have \code{INDIV=1}.  The distinction between "Offspring #1" and "#2" is mostly arbitrary, except that in "mixed" families, the biological offspring MUST have \code{INDIV=1}, and the adopted offspring, \code{INDIV=2}.
  \item Within each family, members must be ordered by \code{INDIV}, as: offspring, mother, father.  For mixed family type, members must be ordered as: bio-offspring, adopted-offspring, mother, father.  For purposes of ordering the phenotype file, subjects with the same family ID but different values for \code{FTYPE} are treated as being in different family units.
  \item The phenotype file has rows as subjects and columns as variables, whereas the genotype file provided to \option{genfile} must have rows as SNPs and columns as subjects.} 

This function handles the following family structures (see \option{sizeLab}): \code{"OOPP"}, 2 offspring and 2 parents; \code{"OO"}, 2 offspring; \code{"PP"}, 2 parents; \code{"OP"}, 1 offspring and 1 parent; and \code{"OPP"}, 1 offspring and two parents.  For each family structure, it handles any combination of the following family types: Mz-twin family type ("Mz"), non-Mz-twin-bio-offspring family type ("Bo"), adopted-offspring family type ("Ad"), and bio/adopted-offspring ("Mix") family type.

When one is conducting parallel analyses on a computing array, judicious use of arguments \option{covmtxfile.in} and \option{covmtxfile.out} can save time.  For example, suppose one is analyzing different SNP sets in parallel but using a common phenotype file for all.  In this case, one should calculate the residual variance-covariance matrix ahead of time and write it to a file.  Then, use the same filename and path for argument \option{covmtxfile.in}, for all jobs running in parallel.  The matrix can be calculated by using \code{\link{gls.batch.get}()} and then \code{\link{fgls}()}.
}
\value{
 \command{gls.batch()} writes an output file with the following columns: "phen","snp","beta","se","t-stat","df","model","pval","method".  However, the actual value returned by the function is simply \code{NULL}.
}
\references{
Li X, Basu S, Miller MB, Iacono WG, McGue M:
A Rapid Generalized Least Squares Model for a Genome-Wide Quantitative Trait Association Analysis in Families.
Hum Hered 2011;71:67-82 (DOI: 10.1159/000324839) 
}
\author{
Xiang Li <lixxx554@umn.edu>, Robert M. Kirkpatrick <kirk0191@umn.edu>, and Saonli Basu <saonli@umn.edu>.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{fgls}}, \code{\link{pheno}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
setwd(tempdir()); getwd() #<--Temp directory to write to.
data(pheno)
data(geno)
data(pedigree)
data(resVCmtx)
gls.batch(
  phenfile=pheno,
  genfile=data.frame(t(geno)),
  pedifile=pedigree,
  outfile="example_output.txt",
  covmtxfile.in=resVCmtx, #<--Precomputed, to save time.
  covmtxfile.out=NULL,
  phen="Zscore",covars="IsFemale",med="rfgls",sizeLab="OOPP",
  Mz=TRUE,Bo=TRUE,Ad=TRUE,Mix=TRUE,indobs=TRUE,
  col.names=TRUE,pediheader=FALSE,pedicolname=c("FAMID","ID","PID","MID","SEX"),
  sep.phe=" ",sep.gen=" ",sep.ped=" ")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
