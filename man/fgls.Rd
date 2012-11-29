\name{fgls}
\alias{fgls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Feasible Generalized Least Squares regression with family GWAS data.
}
\description{
Jointly estimates the fixed-effects coefficients and residual variance-covariance matrix in a generalized least squares model by maximizing the (multivariate-normal) likelihood function, via \code{\link[stats:optim]{optim}()} in the \R base distribution.  The residual variance-covariance matrix is block-diagonal sparse, constructed with \code{\link[bdsmatrix:bdsmatrix]{bdsmatrix}()} from the \command{bdsmatrix} package.
}
\usage{
fgls(fixed, data = parent.frame(), tlist = tlist, sizelist = sizelist,
    sizeLab = "OOPP",Mz = TRUE,Bo = TRUE,Ad = TRUE,Mix = TRUE,indobs = TRUE,
    get.hessian=FALSE,vmat = NULL,subset, weights, na.action)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{
 An object of class 'formula' (or one that can be coerced to that class): a symbolic description of the model
          to be fitted. The RHS of the formula contains the fixed effects of the model.
}
  \item{data}{
    An optional data frame, list or environment (or object
          coercible by \command{as.data.frame()} to a data frame) containing
          the variables in the model.  If not found in \option{data} the
          variables are taken from \command{environment(formula)},
          typically the environment from which \command{fgls()} is called.
}
  \item{tlist}{
The vector of the list of the family labels ("famlab") in the data.  The length of the vector equals the number of family units.  It should be ordered in the same order as the families appear in the data.
Object \option{tlist} is created from the data in the \code{\link{gls.batch}()} function.
}
  \item{sizelist}{
The vector of the list of the family sizes in the data.  The length of the vector equals the number of family units.  It should be ordered as the familin the same order as the families appear in the data.
Object \option{sizelist} is created from the data in the \code{\link{gls.batch}()} function.
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
 Logical (\code{TRUE} or \code{FALSE}).  An indicator of whether Mz-twin families are in the data; must be set to \code{FALSE} if \code{sizeLab="PP"}.
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
  \item{get.hessian}{
Logical; default is \code{FALSE}.  If \code{TRUE}, \command{fgls()} will include the Hessian matrix from \command{optim()} in its output list.  Otherwise, the entry 'hessian' in the list will be \code{NULL}.
}
  \item{vmat}{
The previously estimated (or known) residual variance-covariance matrix (for conducting \emph{Rapid} FGLS).  If it is \code{NULL} (default), the matrix will be jointly estimated with the fixed-effects coefficients (FGLS).  If not \code{NULL}, must be either (1) an object of class \command{bdsmatrix} (from the \code{\link[bdsmatrix]{bdsmatrix}} package), or (2) a character string specifying the filename and path for a single-column text file, with header, containing the "blocks" of a \command{bdsmatrix}.
}
  \item{subset, weights, na.action}{
May be passed through \command{fgls()} as additional arguments to \code{\link[stats:lm]{lm}()}.
}
}
\details{
Function \command{fgls()} was primarily intended to be called automatically, from within \command{gls.batch()}, though it can be called directly if need be.  The difficulty when directly invoking \command{fgls()} is supplying the function with arguments \option{tlist} and \option{sizelist}.  However, these can be obtained easily via \code{\link{gls.batch.get}()}.
%-There are following requirements for the data used with this function:
%-(1) the column of individual id must be named as "iid".
%-(2) within each family, it must be ordered as offspring-1, offspring-2, mother, father.  If it is a mixed family type, it must be: bio-offspring, addopted-offspring, mother, father.
}
\value{
an object of class \command{fgls}.  It includes the following components:
  \item{coefficients}{
The vector of regression coefficients for the fixed effects.  They are ordered as their corresponding terms appear in the RHS of argument \option{fixed}.
}
  \item{estimates}{
The vector of MLEs of the parameters used to construct the residual variance-covariance matrix: correlations and variances.
}
  \item{sigma}{
The estimated residual covariance matrix.  It is of class \command{bdsmatrix}.  One of its slots, \code{sigma@blocks}, can be written to a single-column text file and subsequently read in by \command{gls.batch()}. (Due to its potential size, it is not advised to return \command{sigma} to \R's standard output or print it to the console.)
}
  \item{variance}{
The estimated covariance matrix for the fixed-effects regression coefficients.
}
  \item{ctable}{
Table of coefficients reminiscent of output from \code{\link[stats:summary.lm]{summary.lm}()}.  Each fixed-effect term (including the intercept) has one row of the table,
which are ordered as the terms appear in argument \option{fixed}.  Each row contains a point estimate, an estimated standard error, a t-statistic, and a two-tailed p-value.
}
  \item{residuals}{
Residuals from the feasible generalized-least-squares regression.
}
  \item{fitted.values}{
Predicted phenotype scores from the feasible generalized-least-squares regression.
}
  \item{df.residual}{
Residual degrees of freedom in the feasible generalized-least-squares regression.
}
  \item{loglik}{
The value of the loglikelihood function for the variance-covariance parameters, at its solution.
}
  \item{iter}{
Number of iterations, from the numerical optimization of variance-covariance parameters.  Will be 0 if the residual variance-covariance matrix was provided as argument \code{vmat}.
}
  \item{n}{
Sample size (i.e., number of individual participants).
}
  \item{call}{
Echo of \command{fgls()} function call.
}
  \item{hessian}{
If \code{get.hessian=TRUE}, the Hessian matrix from \code{\link[stats:optim]{optim}()} for the variance-covariance parameters; \code{NULL} otherwise.
}

Function \command{fgls()} also prints to console the vector of variance and correlation parameter estimates it uses to construct the residual variance-covariance matrix.
}
\references{
Li X, Basu S, Miller MB, Iacono WG, McGue M:
A Rapid Generalized Least Squares Model for a Genome-Wide Quantitative Trait Association Analysis in Families.
Hum Hered 2011;71:67-82 (DOI: 10.1159/000324839) 
}
\author{
Xiang Li lixxx554@umn.edu, Robert M. Kirkpatrick <kirk0191@umn.edu>, and Saonli Basu <saonli@umn.edu>.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gls.batch}}
}
\examples{
data(pheno)
data(geno)
data(pedigree)
data(resVCmtx)
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
  fixed=Zscore ~ snp.1 + IsFemale, data=foo$test.dat, sizelist=foo$sizelist, 
  tlist=foo$tlist, 
  sizeLab="OOPP",Mz = TRUE, Bo = TRUE, Ad = TRUE, Mix = TRUE, indobs = TRUE,
  get.hessian=FALSE,
  vmat = resVCmtx, #<--Resid. cov. matrix from fgls onto IsFemale only.
  subset=NULL, weights=NULL, na.action=NULL)
bar$ctable
##  To simultaneously estimate residual covariance matrix
##  and regression coefficients for snp.1 & IsFemale,
##  use the same syntax, except with vmat = NULL .
}
%sizeLab="OOPP";Mz=T;Bo=T;Ad=T;Mix=T;indobs=T
%fgls(test.dat[,phen]~1,data=test.dat,tlist=tlist, famsize=famsize,
%			 sizeLab=sizeLab,Mz=Mz,Bo=Bo,Ad=Ad,Mix=Mix,indobs=indobs,na.action=na.omit)
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
