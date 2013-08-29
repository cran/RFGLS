\name{gls.batch}
\alias{gls.batch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized least-squares batch analysis.
}
\description{
     Fits a generalized least-squares regression model to test association between a quantitative phenotype and all SNPs in a genotype file, one at a time, via Rapid Feasible Generalized Least Squares. For each SNP, genotype is treated as a fixed effect, and the residual variance-covariance matrix is also estimated.  In each trait-SNP association test, the \code{\link{fgls}()} function is used for parameter estimation.
     
The arguments to \command{gls.batch()} may be regarded as belonging to four groups:
\enumerate{
  \item those concerning how to load the \bold{input} (\option{phenfile}, \option{genfile}, \option{pedifile}, \option{covmtxfile.in}, \option{theta}, \option{snp.names}, \option{input.mode}, \option{pediheader}, \option{pedicolname}, \option{sep.phe}, \option{sep.gen}, \option{sep.ped});
  \item those concerning \bold{what to do} with the input, that is, the actual analysis (\option{phen}, \option{covars}, \option{med});
  \item those concerning how to provide the \bold{output} (\option{outfile}, \option{col.names}, \option{return.value}, \option{covmtxfile.out}, \option{covmtxparams.out});
  \item and those that merely trigger \bold{optional checks} on the input (\option{sizeLab}, \option{Mz}, \option{Bo}, \option{Ad}, \option{Mix}, \option{indobs}).
  }     
}
\usage{
gls.batch(phenfile,genfile,pedifile,covmtxfile.in=NULL,theta=NULL,
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
This can be \code{NULL}, in which case no SNPs are analyzed, and \command{gls.batch()} conducts a single \code{\link{fgls}()} regression of the phenotype onto an intercept and covariates (if any).  Otherwise, this argument can be either (1) a character string specifying a genotype file of genotype scores (such as 0,1,2, for the additive genetic model) to be read from disk, or (2) a data frame object containing them.  In such a file, each row must represent a SNP, each column must represent a subject, and there should NOT be column headers or row numbers.  In such a data frame, the reverse holds: each row must represent a subject, and each column, a SNP (e.g. \code{\link{geno}}).  If the data frame--say, \code{geno}--need be transposed, then use \code{genfile=data.frame(t(geno))}.  Using a matrix instead of a data frame is not recommended, because it makes the process of merging data very memory-intensive, and will likely overflow \R's workspace unless the sample size or number of SNPs is quite small.
 
% \emph{Warning:} If \option{genfile} is a data frame, \command{gls.batch} will attempt to remove it from \R's workspace after loading it.  Therefore, it should be a standalone object, and not be part of a list.

 Note that genotype scores need not be integers; they can also be numeric.  So, \command{gls.batch()} can be used to analyze imputed dosages, etc.
}
  \item{pedifile}{
 This can be either (1) a character string specifying the pedigree file corresponding to \option{genfile}, to be read from disk, or (2) a data frame object containing this pedigree information.  At minimum, \option{pedifile} must have a column of subject IDs, named \code{"ID"}, ordered in the same order as subjects' genotypic data in \option{genfile}.  Every row in \option{pedifile} is matched to a participant in \option{genfile}.  That is, if reading files from disk (which is recommended), each row \emph{i} of the pedigree file, which has \emph{n} rows, matches column \emph{i} of the genotype file, which has \emph{n} columns.  This is how the program matches subjects in the phenotype file to their genotypic data.

 The pedigree file or data frame  can also include other columns of pedigree information, like father's ID, mother's ID, etc.  Argument \option{pediheader} (see below) is an indicator of whether the pedigree file on disk has a header or not, with default as \code{FALSE}.  Argument \option{pedicolnames} (see below) gives the names that \command{gls.batch()} will assign to the columns of \option{pedifile}, and the default, \code{c("FAMID","ID","PID","MID","SEX")}, is the familiar "pedigree table" format.  In any event, the user's input \emph{must} somehow provide the program with a column of IDs, labeled as \code{"ID"}.
}
\item{covmtxfile.in}{
  Optional; can be either (1) a character string specifying a file on disk from which the residual variance-covariance matrix is to be read, or (2) the matrix itself.  If \code{NULL}, then \command{gls.batch()} will estimate this matrix.  The file to be read in must be a single column, with a header, containing the contents of the 'blocks' of an object of class \code{\link[bdsmatrix:bdsmatrix-class]{bdsmatrix}}; no other file structures are presently compatible.  If \option{covmtxfile.in} is an actual matrix object, then using one of class \code{\link[bdsmatrix:bdsmatrix-class]{bdsmatrix}} is a virtual requirement.  See below under "Details" for more information.
}
  \item{theta}{
An optional vector of previously estimated (or known) residual-covariance parameters.  Defaults to \code{NULL}, in which case it is ignored.  Otherwise, it must be a numerical vector of of either length 12 if \code{med="UN"}, or of length 3 if \code{med="VC"}.  Each vector element provides the value for the parameter corresponding to its index (serial position).  Values of \code{NA} are accepted for extraneous parameters.  See \code{\link{fgls}()}, under "Details," for which parameters correspond to which indices.  Note that at least one of \option{covmtxfile.in} and \option{theta} must be \code{NULL}.
}
  \item{snp.names}{
An optional character vector with length equal to the number of markers in \option{genfile}, providing names for those markers.  Defaults to \code{NULL}, in which case generic SNP names are used.  Ignored if \option{genfile} is \code{NULL}.
}
  \item{input.mode}{
Either 1 (default), 2, or 3, which tells \command{gls.batch()} where to look for the family-structure variables \code{"FTYPE"} and \code{"INDIV"} (see below, under "Details").  By default, \command{gls.batch()} first looks in the phenotype file, and if the variables are not found there, then looks in the pedigree file, and if the variables are not there, attempts to create them from information available in the pedigree file, via \code{\link{FSV.frompedi}()}.  If \code{input.mode=2}, then \command{gls.batch()} skips looking in the phenotype file, and begins by looking in the pedigree file.  If \code{input.mode=3}, then \command{gls.batch()} skips looking in the phenotype file and pedigree file, and goes straight to \code{\link{FSV.frompedi}()}.
}
  \item{pediheader}{
A logical indicator specifying whether the pedigree file to be read from disk has a header row, to ensure it is read in correctly.  Even if \code{TRUE}, \command{gls.batch()} assigns the values in \option{pedicolname} to the column names after the pedigree file has been read in.  Defaults to \code{FALSE}.  Also see \option{pedifile} above, and under "Details" below.
}
  \item{pedicolname}{
A vector of character strings giving the column names that \command{gls.batch()} will assign to the columns of the pedigree file (starting with the first column and moving left to right). The default, \cr \code{c("FAMID","ID","PID","MID","SEX")}, is the familiar "pedigree table" format.  The two criteria this vector must have are that it must (1) assign the name "ID" to the column of subject IDs in the pedigree file, and (2) its length must not exceed the number of columns of the pedigree file.  If its length is less than the number of columns, columns to which it does not assign a name are discarded.  Also see \option{pedifile} above, and under "Details" below.
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
A character string, either \code{"UN"} or \code{"VC"}, which are the two RFGLS methods described by Li et al. (2011).  If \code{"UN"} (default), which stands for "unstructured," the residual covariance matrix will be constructed from, at most, 12 parameters (8 correlations and 4 variances).  If \code{"VC"}, which stands for "variance components," the residual covariance matrix will be constructed from, at most, 3 variance components (additive-genetic, shared-environmental, and unshared-environmental).  For more information, see \code{\link{fgls}()}.
}
  \item{outfile}{
Either a character string specifying the path and filename for the output file to be written, or \code{NULL}, in which case no output file is written.  The output file contains the SNP analysis results, so argument \option{outfile} is  ignored if \option{genfile} is \code{NULL}.  Note that \command{gls.batch()} will not simultaneously accept \code{outfile=NULL} and \code{return.value=FALSE}.

\emph{Users are warned} that if a file with the same path and filename already exists, \command{gls.batch()} will overwrite it!
}
  \item{col.names}{
 A logical indicator specifying whether to write column names in the output file to be written to disk.  Defaults to \code{TRUE}.
}
  \item{return.value}{
A logical indicator specifying whether function \command{gls.batch()} should actually return a value.  Defaults to \code{FALSE}, in which case the function merely returns \code{NULL}.  If \code{TRUE} and non-\code{NULL} value was supplied to \option{genfile}, the function returns a data frame containing the results of the SNP analyses(i.e., the output file as a data frame).  If \code{TRUE} and \code{genfile=NULL}, the function returns the \code{\link{fgls}()} output from a regression of the phenotype onto an intercept and covariates (if any).  Note that \command{gls.batch()} will not simultaneously accept \code{outfile=NULL} and \code{return.value=FALSE}.
} 
  \item{covmtxfile.out}{
An optional character string specifying the filename and path to which the residual covariance matrix, if it is to be constructed (i.e., if \code{covmtxfile.in=NULL}), will be written.  The default is \code{NULL}, in which case no such file is written to disk.  See below under "Details" for more information.

\emph{Users are warned} that if a file with the same path and filename already exists, \command{gls.batch()} will overwrite it!
}
  \item{covmtxparams.out}{
An optional character string specifying the filename and path to which the vector of residual-covariance parameters, if they are to be estimated (i.e., if \option{covmtxfile.in} and \option{theta} are both \code{NULL}), will be written.  The default is \code{NULL}, in which case no such file is written to disk.  See below under "Details" for more information.

\emph{Users are warned} that if a file with the same path and filename already exists, \command{gls.batch()} will overwrite it!
}
  \item{sizeLab}{
This is an optional argument, and may be eliminated in future versions of this package.  Defaults to \code{NULL};  otherwise, must be a character string, and if the number of characters in the string is not equal to the size of the largest family in the data, \command{gls.batch()} will produce a warning.
}
  \item{Mz, Bo, Ad, Mix}{
These are optional logical indicators that specify whether families containing MZ twins (\option{MZ}; family-type 1), DZ twins or full siblings (\option{Bo}; family-types 2 and 4), two adoptees (\option{Ad}; family-type 3), or 1 biological offspring and 1 adoptee (\option{Mix}; family-type 5) are present in the data.  The values of each are checked against the actual family types present, after loading and merging the data and trimming out incomplete cases, and a warning is generated for each mismatch.  If any of these four arguments is \code{NULL} (which is the default), the check corresponding to that family type is skipped.
}
  \item{indobs}{
An optional logical indicator of whether there are "independent observations" who do not fit into a four-person nuclear family present in the data.  After loading and merging the data and trimming out incomplete cases, the value of \option{indobs} is checked against whether such individuals are actually present, and a warning is generated in case of a mismatch.  If \code{indobs=NULL}, which is the default, this check is skipped.
}
}
\details{

 Reference is frequently made throughout this documentation to the "phenotype file," the "genotype file," and so forth, because \command{gls.batch()} was intended to be used with potentially large datafiles to be read from disk.  This should be evident from the presence of the word "file" in the names of many of this function's arguments, and the fact that all of those arguments may be character strings providing a filename and path.  However, it can also accept the data if the file has already been loaded into \R's workspace as a data frame object, in which case "the [whatever] file" should be taken to refer to such a data frame.  For details specific to each argument, see above.
 
 The function \command{gls.batch()} first reads in the files and merges them into a data frame with columns of family-structure information, phenotype, covariates, and genotypes.  Then, it creates a \option{tlist} vector and a \option{sizelist} vector, which comprise the family labels and family sizes in the data.  Finally, it carries out single-SNP association analyses for all the SNPs in the genotype file.
 
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
The distinction between "Offspring #1" and "#2" is mostly arbitrary, except that in "mixed" families(\code{FTYPE=5}), the biological offspring MUST have \code{INDIV=1}, and the adopted offspring, \code{INDIV=2}.  If the phenotype file contains variables \code{"FTYPE"} and \code{"INDIV"}, it should be ordered by family ID (\code{"FAMID"}), and by individual code \code{"INDIV"} within family ID.  Note that \command{gls.batch()} treats participants with \code{FTYPE=6} as the sole members of their own family units, and not as part of the family corresponding to their family ID.

If neither the phenotype nor pedigree file contain \code{"FTYPE"} and \code{"INDIV"}, \command{gls.batch()} will construct them via \code{\link{FSV.frompedi}()}.

When one is conducting parallel analyses on a computing array, judicious use of arguments \option{covmtxfile.in},  \option{theta}, \option{covmtxparams.out}, and \option{covmtxfile.out} can save time.  For example, suppose one is analyzing different SNP sets in parallel but using a common phenotype file for all.  In this case, one could calculate the residual covariance matrix ahead of time and write it to a file.  Then, use the same filename and path for argument \option{covmtxfile.in}, for all jobs running in parallel.  The matrix can be calculated by using \code{\link{gls.batch.get}()} and then \code{\link{fgls}()}.  One could similarly obtain the residual-covariance parameters ahead of time, and supply them as a vector to \option{theta} in all jobs running in parallel.

}

\value{
If \code{return.value=FALSE}, then \command{gls.batch()} simply returns \code{NULL}.  If \code{return.value=TRUE} and \code{genfile=NULL}, then \command{gls.batch()} returns the \code{\link{fgls}()} output from a regression of the phenotype onto an intercept and covariates (if any).  If \code{return.value=TRUE} and \option{genfile} is non-\code{NULL}, then \command{gls.batch()} returns a data frame containing the results of the single-SNP analyses, 1 row per SNP.  Specifically, this data frame includes the following named columns:
  \itemize{
    \item \code{snp} (character): the names of the SNPs; equal to \option{snp.names} if any were supplied.
    \item \code{coef} (numeric): the regression coefficients of the SNPs.
    \item \code{se} (numeric): estimated standard errors of SNPs' regression coefficients.
    \item \code{t.stat} (numeric): \emph{t}-statistics, i.e. regression coefficients divided by their estimated standard errors.
    \item \code{df} (integer): degrees-of-freedom (see \code{df.residual}, from \code{\link{fgls}()} output).
    \item \code{pval} (numeric): two-tailed \emph{p}-values, from corresponding \emph{t}-statistics and degrees-of-freedom.
  }

Function \command{gls.batch()} also has optional side effects of writing as many as three files to disk, depending on arguments \option{outfile}, \option{covmtxfile.out}, and \option{covmtxparams.out}.  Note that if a file is written for  \option{outfile}, that file will contain the single-SNP analysis results described above.
}
\references{
Li X, Basu S, Miller MB, Iacono WG, McGue M:
A Rapid Generalized Least Squares Model for a Genome-Wide Quantitative Trait Association Analysis in Families.
\emph{Human Heredity} 2011;71:67-82 (DOI: 10.1159/000324839) 
}
\author{
Xiang Li \email{lixxx554@umn.edu}, Robert M. Kirkpatrick \email{kirk0191@umn.edu}, and Saonli Basu \email{saonli@umn.edu} .
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
data(pheno)
data(geno)
data(map)
data(pedigree)
data(rescovmtx)
minigwas <- gls.batch(
  phenfile=pheno,genfile=data.frame(t(geno)),pedifile=pedigree,
  covmtxfile.in=rescovmtx, #<--Precomputed, to save time.
  theta=NULL,snp.names=map[,2],input.mode=c(1,2,3),pediheader=FALSE, 
  pedicolname=c("FAMID","ID","PID","MID","SEX"),
  sep.phe=" ",sep.gen=" ",sep.ped=" ",
  phen="Zscore",covars="IsFemale",med=c("UN","VC"),
  outfile=NULL,col.names=TRUE,return.value=TRUE,
  covmtxfile.out=NULL,covmtxparams.out=NULL,
  sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL)
minigwas
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
