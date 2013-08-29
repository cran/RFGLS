\name{gls.batch.get}
\alias{gls.batch.get}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data restructuring for \code{\link{fgls}()}.
}
\description{
Carries out the data restructuring performed by \code{\link{gls.batch}()}.  Useful if calling \code{\link{fgls}()} directly.
     
Several arguments to \command{gls.batch.get()} are accepted only for the sake of parallelism with \code{\link{gls.batch}()}, and are ignored: \option{covmtxfile.in}, \option{theta}, \option{outfile}, \option{col.names}, \option{return.value}, \option{covmtxfile.out}, and \option{covmtxparams.out}.
}
\usage{
gls.batch.get(phenfile,genfile,pedifile,covmtxfile.in=NULL,theta=NULL,
  snp.names=NULL,input.mode=c(1,2,3),pediheader=FALSE, 
  pedicolname=c("FAMID","ID","PID","MID","SEX"),
  sep.phe=" ",sep.gen=" ",sep.ped=" ",
  phen,covars=NULL,med=c("UN","VC"),
  outfile,col.names=TRUE,return.value=FALSE,
  covmtxfile.out=NULL,
  covmtxparams.out=NULL,
  sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{phenfile}{
 This can be either (1) a character string specifying a phenotype file on disk which includes the phenotypes and other covariates, or (2) a data frame object containing the same data.  In either case, the data must be appropriately structured.  See below under "Details."
}

\item{genfile}{
This can be \code{NULL}, in which case no SNP data are loaded.  Otherwise, this argument can be either (1) a character string specifying a genotype file of genotype scores (such as 0,1,2, for the additive genetic model) to be read from disk, or (2) a data frame object containing them.  In such a file, each row must represent a SNP, each column must represent a subject, and there should NOT be column headers or row numbers.  In such a data frame, the reverse holds: each row must represent a subject, and each column, a SNP (e.g. \code{\link{geno}}).  If the data frame--say, \code{geno}--need be transposed, then use \code{genfile=data.frame(t(geno))}.  Using a matrix instead of a data frame is not recommended, because it makes the process of merging data very memory-intensive, and will likely overflow \R's workspace unless the sample size or number of SNPs is quite small.
 
% \emph{Warning:} If \option{genfile} is a data frame, \command{gls.batch} will attempt to remove it from \R's workspace after loading it.  Therefore, it should be a standalone object, and not be part of a list.

 Note that genotype scores need not be integers; they can also be numeric.  So, \command{gls.batch()} can be used to analyze imputed dosages, etc.
}
  \item{pedifile}{
 This can be either (1) a character string specifying the pedigree file corresponding to \option{genfile}, to be read from disk, or (2) a data frame object containing this pedigree information.  At minimum, \option{pedifile} must have a column of subject IDs, named \code{'ID'}, ordered in the same order as subjects' genotypic data in \option{genfile}.  Every row in \option{pedifile} is matched to a participant in \option{genfile}.  That is, if reading files from disk (which is recommended), each row \emph{i} of the pedigree file, which has \emph{n} rows, matches column \emph{i} of the genotype file, which has \emph{n} columns.  This is how the program matches subjects in the phenotype file to their genotypic data.

 The pedigree file or data frame  can also include other columns of pedigree information, like father's ID, mother's ID, etc.  Argument \option{pediheader} (see below) is an indicator of whether the pedigree file on disk has a header or not, with default as \code{FALSE}.  Argument \option{pedicolnames} (see below) gives the names that \command{gls.batch.get()} will assign to the columns of \option{pedifile}, and the default, \code{c("FAMID","ID","PID","MID","SEX")}, is the familiar "pedigree table" format.  In any event, the user's input \emph{must} somehow provide the program with a column of IDs, labeled as \code{"ID"}.
}
\item{covmtxfile.in}{
  Accepted but not used.
}
  \item{theta}{
Accepted but not used.
}
  \item{snp.names}{
An optional character vector with length equal to the number of markers in \option{genfile}, providing names for those markers.  Defaults to \code{NULL}, in which case generic SNP names are used.  Ignored if \option{genfile} is \code{NULL}.
}
  \item{input.mode}{
Either 1 (default), 2, or 3, which tells \command{gls.batch.get()} where to look for the family-structure variables \code{"FTYPE"} and \code{"INDIV"} (see below, under "Details").  By default, \command{gls.batch.get()} first looks in the phenotype file, and if the variables are not found there, then looks in the pedigree file, and if the variables are not there, attempts to create them from information available in the pedigree file, via \code{\link{FSV.frompedi}()}.  If \code{input.mode=2}, then \command{gls.batch.get()} skips looking in the phenotype file, and begins by looking in the pedigree file.  If \code{input.mode=3}, then \command{gls.batch.get()} skips looking in the phenotype file and pedigree file, and goes straight to \code{\link{FSV.frompedi}()}.
}
  \item{pediheader}{
A logical indicator specifying whether the pedigree file to be read from disk has a header row, to ensure it is read in correctly.  Even if \code{TRUE}, \command{gls.batch()} assigns the values in \option{pedicolname} to the columns after it has been read in.  Defaults to \code{FALSE}.  Also see \option{pedifile} above and under "Details" below.
}
  \item{pedicolname}{
A vector of character strings giving the column names that \command{gls.batch.get()} will assign to the columns of the pedigree file (starting with the first column and moving left to right). The default, \code{c("FAMID","ID","PID","MID","SEX")}, is the familiar "pedigree table" format.  The two criteria this vector must have are that it must (1) assign the name "ID" to the column of subject IDs in the pedigree file, and (2) its length must not exceed the number of columns of the pedigree file.  If its length is less than the number of columns, columns to which it does not assign a name are discarded.  Also see \option{pedifile} above, and under "Details" below.
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
  \item{phen}{
 A character string specifying the phenotype (column name) in the phenotype file to be analyzed.
}
  \item{covars}{
 A character string or character vector that holds the (column) names of the covariates, in the phenotype file, to be used in the regression model.  Defaults to \code{NULL}, in which case no covariates are included.
}
  \item{med}{
A character string, either \code{"UN"} or \code{"VC"}, which are the two RFGLS methods described by Li et al. (2011).  If \code{"UN"} (default), which stands for "unstructured," the residual covariance matrix will be constructed from, at most, 12 parameters (8 correlations and 4 variances).  If \code{"VC"}, which stands for "variance components," the residual covariance matrix will be constructed from, at most, 3 variance components (additive-genetic, shared-environmental, and unshared-environmental).
}
  \item{outfile}{
Accepted but not used.
}
  \item{col.names}{
Accepted but not used.
}
  \item{return.value}{
Accepted but not used.
} 
  \item{covmtxfile.out}{
Accepted but not used.
}
  \item{covmtxparams.out}{
Accepted but not used.
}
  \item{sizeLab}{
This is an optional argument, and may be eliminated in future versions of this package.  Defaults to \code{NULL};  otherwise, must be a character string.  If the number of characters in the string is not equal to the size of the largest family in the data, \command{gls.batch.get()} will produce a warning.
}
  \item{Mz, Bo, Ad, Mix}{
These are optional logical indicators that specify whether families containing MZ twins (\option{MZ}; family-type 1), DZ twins or full siblings (\option{Bo}; family-types 2 and 4), two adoptees (\option{Ad}; family-type 3), or 1 biological offspring and 1 adoptee (\option{Mix}; family-type 5) are present in the data.  The values of each are checked against the actual family types present, after loading and merging the data and trimming out incomplete cases, and a warning is generated for each mismatch.  If any of these four arguments is \code{NULL} (which is the default), the check corresponding to that family type is skipped.
}
  \item{indobs}{
An optional logical indicator of whether there are "independent observations" who do not fit into a four-person nuclear family present in the data.  After loading and merging the data and trimming out incomplete cases, the value of \option{indobs} is checked against whether such individuals are actually present, and a warning is generated in case of a mismatch.  If \code{indobs=NULL}, which is the default, this check is skipped.
}
}
%old arguments

\details{
  Though originally used for debugging purposes, \command{gls.batch.get()} was included because it facilitates directly invoking \code{\link{fgls}()} when the need arises.  This function first reads in the files and merges the files into a data frame with columns of family-structure information, phenotype, covariates, and genotypes.  It then creates a \option{tlist} vector and a \option{sizelist} vector, which comprise the family labels and family sizes in the data.
 It returns a list containing the merged data frame, and the \option{tlist} and \option{sizelist} vectors.
 
 At the \emph{bare minimum}, the phenotype file must contain columns named \code{"ID"}, \code{"FAMID"}, and whatever character string is supplied to \option{phen}.  These columns respectively contain individual IDs, family IDs, and phenotype scores; individual IDs must be unique.
 
 At the \emph{bare minimum}, the pedigree file need only contain a column consisting of unique individual IDs, corresponding to the label \code{"ID"} in \option{pedicolname}.  The number of participants in the pedigree file must equal the number of participants in the genotype file, with participants ordered the same way in both files.  However, the default value for argument \option{pedicolname} (see above) assumes five columns, in the familiar "pedigree table" format.
 
 The phenotype file or pedigree file may also contain the two key family-structure variables, \code{"FTYPE"} (family-type) and \code{"INDIV"} (individual code).  If both contain these variables, then by default, they are read from the phenotype file (but see argument \option{input.mode} above).  There are six recognized family types, which are distinguished primarily by how the offspring in the family are related to one another:
 \itemize{
  \item \code{FTYPE=1}, containing MZ twins;
  \item \code{FTYPE=2}, containing DZ twins;
  \item \code{FTYPE=3}, containing adoptees;
  \item \code{FTYPE=4}, containing non-twin full siblings;
  \item \code{FTYPE=5}, "mixed" families containing one biological offspring and one adoptee;
  \item \code{FTYPE=6}, containing "independent observations" who do not fit into a four-person nuclear family.
  }
 It is assumed that all offspring except adoptees are biological children of the parents in the family.  The four individual codes are:
 \itemize{
  \item \code{INDIV=1} is for "Offspring #1;"
  \item \code{INDIV=2} is for "Offspring #2;"
  \item \code{INDIV=3} is for mothers;
  \item \code{INDIV=4} is for fathers.
  }
The distinction between "Offspring #1" and "#2" is mostly arbitrary, except that in "mixed" families(\code{FTYPE=5}), the biological offspring MUST have \code{INDIV=1}, and the adopted offspring, \code{INDIV=2}.  If the phenotype file contains variables \code{"FTYPE"} and \code{"INDIV"}, it should be ordered by family ID (\code{"FAMID"}), and by individual code \code{"INDIV"} within family ID.  Note that \command{gls.batch.get()} treats participants with \code{FTYPE=6} as the sole members of their own family units, and not as part of the family corresponding to their family ID.

If neither the phenotype nor pedigree file contain \code{"FTYPE"} and \code{"INDIV"}, \command{gls.batch()} will construct them via \code{\link{FSV.frompedi}()}.

}
\value{
A list with these three components:
  \item{test.dat}{
The merged data frame of family-structure variables, phenotype, covariates, and genotypes.  Participants of family-type 6 will be moved to the end of the data frame.  There will also be three additional columns:
  \itemize{
    \item \code{famsize} (integer): The size of the family to which each participant belongs.
    \item \code{unisid} (character): Single-character representation of each participants' \code{"FTYPE"} and \code{"INDIV"}.  Adoptees have "a", MZ twins have "c", non-MZ-twin biological offspring have "b", mothers have "m", fathers have "f", and members of family-type 6 have \code{NA}.
    \item \code{famlab} (character): "Family labels;" the \code{unisid}'s of the members of each participant's family, pasted together in order of \code{"INDIV"}.  
  }
}
  \item{tlist}{
A character vector of family labels, with length equal to the number of families in the data (each participant of family-type 6 is treated as a separate family).  The names of its components are the family IDs.
}
  \item{sizelist}{
A vector of family sizes, with length equal to the number of families in the data (each participant of family-type 6 is treated as a separate family).  The names of its components are the family IDs.
}
}
\author{
Xiang Li \email{lixxx554@umn.edu}, Robert M. Kirkpatrick \email{kirk0191@umn.edu}, and Saonli Basu \email{saonli@umn.edu} .
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{fgls}}, \code{\link{gls.batch}}
%%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(pheno)
data(geno)
data(map)
data(pedigree)
foo <- gls.batch.get(
  phenfile=pheno,genfile=data.frame(t(geno)),pedifile=pedigree,
  covmtxfile.in=NULL,theta=NULL,snp.names=map[,2],input.mode=c(1,2,3),
  pediheader=FALSE,pedicolname=c("FAMID","ID","PID","MID","SEX"),
  sep.phe=" ",sep.gen=" ",sep.ped=" ",
  phen="Zscore",covars="IsFemale",med=c("UN","VC"),
  outfile,col.names=TRUE,return.value=FALSE,
  covmtxfile.out=NULL,
  covmtxparams.out=NULL,
  sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL)
olsmod <- lm(   ##<--OLS regression could be applied to the merged dataset...
    Zscore ~ rs3934834 + IsFemale, data=foo$test.dat)
summary(olsmod)  ##<--...but the standard errors and t-statistics will not be valid.

##The 'tlist' vector can be useful for figuring out if any residual-covariance
##parameters are poorly identified in the data:
pheno2 <- subset(pheno, (pheno$INDIV<3 & pheno$FAMID>20) | 
                       (pheno$ID \%in\% c(11,12,13,21,22,23)))
foo2 <- gls.batch.get(
  phenfile=pheno2,
  genfile=data.frame(t(geno)),pedifile=pedigree,
  covmtxfile.in=NULL,theta=NULL,snp.names=map[,2],input.mode=c(1,2,3),
  pediheader=FALSE,pedicolname=c("FAMID","ID","PID","MID","SEX"),
  sep.phe=" ",sep.gen=" ",sep.ped=" ",
  phen="Zscore",covars="IsFemale",med=c("UN","VC"),
  outfile,col.names=TRUE,return.value=FALSE,
  covmtxfile.out=NULL,
  covmtxparams.out=NULL,
  sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL)
table(foo2$tlist)
##Only two families have the label 'ccm', that is, only two have 
##a mother.  So, if calling fgls()
##with med="UN", it would probably be a good idea to drop the
##mother variance [drop=10], or the biological mother-offspring 
##correlation [drop=2], or both [drop=c(2,10)].
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
