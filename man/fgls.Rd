\name{fgls}
\alias{fgls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Feasible Generalized Least Squares regression with family GWAS data.
}
\description{
Jointly estimates the fixed-effects coefficients and residual variance-covariance matrix in a generalized least squares model by minimizing the (multivariate-normal) negative loglikelihood function, via \code{\link[stats:optim]{optim}()} in the \R base distribution.  The residual variance-covariance matrix is block-diagonal sparse, constructed with \code{\link[bdsmatrix:bdsmatrix]{bdsmatrix}()} from the \command{bdsmatrix} package.
}
\usage{
fgls(fixed, data=parent.frame(), tlist=tlist, sizelist=sizelist,
  med=c("UN","VC"), vmat=NULL, start=NULL, theta=NULL, drop=NULL,
  get.hessian=FALSE, optim.method="BFGS", control=list(), weights=NULL,
  sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{
An object of class 'formula' (or one that can be coerced to that class): a symbolic description of the regression model
to be fitted. The RHS of the formula contains the fixed effects of the model.
}
  \item{data}{
An optional data frame, list or environment (or object coercible by \code{\link[base:as.data.frame]{as.data.frame()}} to a data frame) containing the variables in the model, as specified in argument \option{fixed}.  If not found in \option{data} the variables are taken from \command{environment(formula)}, typically the environment from which \command{fgls()} is called.  If it contains a column named \code{"ID"}, then that column's values will be the row and column names of output element \code{sigma} (see below, under "Value").
}
  \item{tlist}{
The character vector of the family labels ("famlab") in the data.  The length of the vector equals the number of family units.  It should be ordered in the same order as the families appear in the data.  Object \option{tlist} is created by the \code{\link{gls.batch}()} and \code{\link{gls.batch.get}()} functions.
}
  \item{sizelist}{
The integer vector of the  family sizes in the data.  The length of the vector equals the number of family units.  It should be ordered in the same order as the families appear in the data.  Object \option{sizelist} is created by the \code{\link{gls.batch}()} and \code{\link{gls.batch.get}()} functions.
}
  \item{med}{
A character string, either \code{"UN"} or \code{"VC"}, which are the two RFGLS methods described by Li et al. (2011).  If \code{"UN"} (default), which stands for "unstructured," the residual covariance matrix will be constructed from, at most, 12 parameters (8 correlations and 4 variances).  If \code{"VC"}, which stands for "variance components," the residual covariance matrix will be constructed from, at most, 3 variance components (additive-genetic, shared-environmental, and unshared-environmental).
}
  \item{vmat}{
The previously estimated (or known) residual covariance matrix (for conducting \emph{Rapid} FGLS).  If it is \code{NULL} (default), the matrix will either be jointly estimated with the fixed-effects coefficients (FGLS), or will be constructed values supplied to \option{theta}.  If not \code{NULL}, must be either (1) an object of class 'bdsmatrix' (from the \code{\link[bdsmatrix]{bdsmatrix}} package), or (2) a character string specifying the filename and path for a single-column text file, with header, containing the "blocks" of a \code{bdsmatrix}.  Note that at least one of \option{vmat} and \option{theta} must be \code{NULL}.
}
  \item{start}{
A numeric vector of initial values for the residual-covariance parameters.  If \code{NULL} (default), generic start values are used.  Otherwise, it must be a numerical vector of either length 12 if \code{med="UN"}, or of length 3 if \code{med="VC"}.  Each vector element provides the initial value for the parameter corresponding to its index (serial position).  Values of \code{NA} are accepted, and will be replaced with the generic start value for that parameter.  See below under "Details" for which parameters correspond to which indices.  Users should bear in mind that especially poor start values can cause optimization to fail.  Ignored if no residual-covariance parameters are estimated.
}
  \item{theta}{
A numeric vector of previously estimated (or known) residual-covariance parameters.  Defaults to \code{NULL}, in which case it is ignored.  Otherwise, it must be a numerical vector of of either length 12 if \code{med="UN"}, or of length 3 if \code{med="VC"}.  Each vector element provides the value for the parameter corresponding to its index (serial position).  Values of \code{NA} are accepted for extraneous parameters.  See below under "Details" for which parameters correspond to which indices.  The residual covariance matrix is constructed from the elements of \option{theta}, exactly as-is.  Note that at least one of \option{vmat} and \option{theta} must be \code{NULL}.
}
  \item{drop}{
An integer vector of indices (serial positions) specifying which residual-covariance parameters to drop.  Dropped parameters are not estimated.  In addition to those specified by \option{drop}, \command{fgls()} automatically identifies which parameters are \emph{completely} unidentified from the data (i.e., zero observations in the data are informative about them), and drops them as well.  If \option{drop} is \code{NULL} (default), no user-specified parameters are dropped.  Otherwise, it must be a vector of integers, either between 1 and 12 (inclusive) if \code{med="UN"}, or between 1 and 2 (inclusive) if \code{med="VC"}.  Note that if a user-specified-dropped parameter ends up being needed to construct the residual covariance matrix, its value is taken to be that of its OLS equivalent: zero for correlations (\code{med="UN"}) and for the familial variance components (\code{med="VC"}), and the OLS residual variance for variances (\code{med="UN"}).  It may be prudent to drop parameters when very few observations in the data are informative about them, which can at least save computation time.  See below under "Details" for which parameters correspond to which indices.  Ignored if no residual-covariance parameters are estimated.
}
  \item{get.hessian}{
Logical; default is \code{FALSE}.  If \code{TRUE}, \command{fgls()} will include the Hessian matrix from \command{optim()} in its output list.  Otherwise, the entry 'hessian' in the list will be \code{NULL}.  Ignored if no residual-covariance parameters are estimated.
}
  \item{optim.method}{
Character string, passed as \option{method} to \code{\link[stats:optim]{optim}()}.  Ignored if no residual-covariance parameters are estimated.  The default, \code{"BFGS"}, is usually fast and is recommended for general use.  If method \code{"L-BFGS-B"} is used, \command{fgls()} will supply \code{\link[stats:optim]{optim}()} with reasonable box constraints on the parameters, intended for use with \code{\link[stats:optim]{optim}()}'s default control parameters (see argument \option{control} below).

Method \code{"BFGS"} (the default) may fail when any of the residual-covariance parameters are poorly identified from the data.  In these cases, it may be wise simply to drop the offending parameters.  Other optimization methods, including \code{"L-BFGS-B"}, can succeed where \code{"BFGS"} fails.  Method \code{"SANN"} should \emph{not} generally be relied upon to find the global optimum, but it can sometimes produce reasonable, approximate solutions in instances where no other method works.  As a last-resort diagnostic, one can combine \code{optim.method="SANN"} with \code{hessian=TRUE}, since the resulting Hessian matrix may provide clues as to which parameters are causing problems.
}
  \item{control}{
A list of control parameters passed to \code{\link[stats:optim]{optim}()}, intended for advanced users.  The default is also \code{\link[stats:optim]{optim}()}'s default, which should be adequate for general use.
  }
  \item{weights}{
A numeric vector of weights, with length equal to the number of observations in the data.  Defaults to \code{NULL}.
}
  \item{sizeLab}{
This is an optional argument, and may be eliminated in future versions of this package.  Defaults to \code{NULL};  otherwise, must be a character string.  If the number of characters in the string is not equal to the size of the largest family in the data, \command{fgls()} will produce a warning.
}
  \item{Mz, Bo, Ad, Mix, indobs}{
These arguments are deprecated, and their values are ignored.  They are retained in this package version for legacy reasons, but will be eliminated in future versions.
}
}

\details{
Function \command{fgls()} was originally intended to be called automatically, from within \code{\link{gls.batch}()}.  However, calling it directly is likely to be useful to advanced users.  The difficulty when directly invoking \command{fgls()} is supplying the function with arguments \option{tlist} and \option{sizelist}.  But, these can be obtained easily via \code{\link{gls.batch.get}()}.

When residual-covariance parameters are to be estimated, \command{fgls()} will attempt optimization, at most, two times.  If the initial attempt fails, \command{fgls()} prints a message saying so to the console, and tries a second time.  On the second attempt, before each evaluation of the objective function, the blocks composing the block-diagonal residual covariance matrix are forced to be positive definite.  This uses \code{\link[Matrix:nearPD]{nearPD}()} from the \emph{Matrix} package, which turns each block matrix into its nearest positive-definite approximation (where "nearest" is meant in a least-squares sense).  Forcing positive-definiteness in this way is only used for the second attempt, and not for the initial attempt (which has its own way of ensuring a positive-definite solution), since it slows down optimization and is unnecessary when the parameters are well-identified.  Furthermore, it can have consequences the user might not expect.  For instance, in \command{fgls()}'s output (see below, under "Value"), the elements of the residual covariance matrix \code{sigma} might not correspond to the parameter estimates in \code{estimates}, or covariances that are supposed to be the same across families might not be so in the actual matrix \code{sigma}.  Nevertheless, the second attempt may succeed when the initial attempt fails.

When \code{med="UN"}, the residual covariance matrix is constructed from, at most, 12 parameters--8 correlations and 4 variances.  Below is an enumerated list of those 12 parameters, in which the number of each list entry is the index (serial position) of that parameter, and the quoted text is the element name of each estimated parameter as it appears in \command{fgls()} output:
  \enumerate{
  \item "cor(m,f)", correlation between mothers and fathers.
  \item "cor(c/b,m)", correlation between biological offspring and mothers.
  \item "cor(c/b,f)", correlation between biological offspring and fathers.
  \item "cor(c,c)", MZ-twin correlation.
  \item "cor(b,b)", full-sibling (DZ-twin) correlation.
  \item "cor(a,m)", correlation between adoptees and mothers.
  \item "cor(a,f)", correlation between adoptees and fathers.
  \item "cor(a,a)", adoptive-sibling correlation.
  \item "var(O)", offspring variance.
  \item "var(m)", mother variance.
  \item "var(f)", father variance.
  \item "var(ind)", variance for "independent observations."
  }
  
When \code{med="VC"}, the residual covariance matrix is constructed from, at most, 3 variance components.  Below is an enumerated list of those 3 parameters, in which the number of each list entry is the index (serial position) of that parameter, and the quoted text is the label of each estimated parameter as it appears in \command{fgls()} output:
\enumerate{
  \item "A", additive-genetic variance.
  \item "C", shared-environmental variance (compound-symmetric within families).
  \item "E", unshared-environmental variance (which cannot be dropped).
  }
Additive-genetic variance contributes to covariance between family members commensurately to the expected proportion of segregating alleles they share: 1.0 for MZ twins, 0.5 for first-degree relatives, 0 for spouses and adoptive relatives.  Shared-environmental variance, as defined here, represents covariance between biologically unrelated family members (including spouses).

In package version 1.0, arguments \option{subset} and \option{na.action} were accepted, and passed to \code{\link[stats:lm]{lm}()}.  Neither are accepted any longer.  Subsetting should be done before directly calling \command{fgls()}; the function handles \code{NA}'s in the data by what is (in effect) \code{na.action=na.omit}.
%-There are following requirements for the data used with this function:
%-(1) the column of individual id must be named as "iid".
%-(2) within each family, it must be ordered as offspring-1, offspring-2, mother, father.  If it is a mixed family type, it must be: bio-offspring, addopted-offspring, mother, father.
}
\value{
An object of class 'fgls'.  It includes the following components:
  \item{ctable}{
Table of coefficients reminiscent of output from \code{\link[stats:summary.lm]{summary.lm}()}.  Each fixed-effect term (including the intercept) has one row of the table,
which are ordered as the terms appear in argument \option{fixed}.  Each row contains a point estimate, an estimated standard error, a \emph{t}-statistic, and a two-tailed \emph{p}-value.
}
  \item{Rsqd}{
The generalized-least-squares coefficient of determination, \emph{a la} Buse (1973).
}
  \item{estimates}{
The vector of MLEs of the parameters used to construct the residual covariance matrix, ordered as in the lists above, under "Details."  Dropped parameters are given value \code{NA}.  If no residual-covariance parameters are estimated, will instead be a single \code{NA}.
}
  \item{drop}{
A vector of parameter indices, representing which residual-covariance parameters were dropped (not estimated).  See above, under "Details," for which parameters correspond to which indices.  \code{NULL} if no parameters were dropped or estimated.
}
  \item{iter}{
\code{NULL} if no residual-covariance parameters were estimated.  Otherwise, a single-row data frame, containing miscellaneous output pertaining to the optimization, specifically, the following named columns:
  \enumerate{
  \item \code{iterations} (integer): the number of function iterations, as returned from \code{\link[stats:optim]{optim}()}.
  \item \code{convergence} (integer): convergence code, as returned from \code{\link[stats:optim]{optim}()}; value 0 means that convergence was successful.
  \item \code{message} (character): additional information from the optimizer; a single whitespace means that \code{\link[stats:optim]{optim}()} returned a message of \code{NULL}.
  \item \code{first_try} (logical): Did \command{fgls()}'s first attempt at optimization succeed?  If \code{FALSE}, then during the second attempt, \command{fgls()} had to force the block matrices of the residual covariance matrix to be positive-definite, as described above, under "Details."
  }
}
  \item{loglik}{
The \emph{negative} loglikelihood, at the solution.  If the residual-covariance parameters were estimated, then it equals -1 times the maximized joint loglikelihood of those parameters and the regression coefficients.  If the residual-covariance parameter values were provided with argument \option{vmat} or \option{theta}, then it equals -1 times the maximized joint loglikelihood of the regression coefficients, \emph{conditional} on the values supplied for the residual-covariance parameters.
}
\item{sigma}{
The residual covariance matrix.  It is of class 'bdsmatrix'.  Its row and column names are taken from the column named "ID", if any, in argument \code{data}, otherwise its row and column names are sequential numbers.   One of its slots, \code{sigma@blocks}, can be written to a single-column text file and subsequently read in by \code{\link{gls.batch}()}. Due to its potential size, it is not advised to return \code{sigma} to \R's standard output or print it to the console.
}
  \item{hessian}{
If \code{get.hessian=TRUE} and residual-covariance parameters were estimated, the Hessian matrix from \code{\link[stats:optim]{optim}()} for those parameters; \code{NULL} otherwise.
}
  \item{n}{
Sample size (i.e., number of individual participants), after excluding those with missing data (\code{NA}'s).
}
\item{df.residual}{
Residual degrees of freedom in the feasible generalized-least-squares regression, as returned by \code{\link[stats:lm]{lm}()}.  Note that it only reflects the number of regression coefficients, and not the number of residual-covariance parameters that were estimated.
}
\item{residuals}{
Residuals from the feasible generalized-least-squares regression.  It is a vector of length \code{n}, i.e. it is not padded with \code{NA}'s for participants with missing data.
}
  \item{fitted.values}{
Predicted phenotype scores from the feasible generalized-least-squares regression.  It is a vector of length \code{n}, i.e. it is not padded with \code{NA}'s for participants with missing data.
}
  \item{variance}{
The estimated covariance matrix for (the sampling distribution of) the fixed-effects regression coefficients.
}
\item{call}{
Echo of \command{fgls()} function call.
}
Function \command{fgls()} also prints to console the estimates of non-dropped residual-covariance parameters (if any).
}

\references{
Li X, Basu S, Miller MB, Iacono WG, McGue M:
A Rapid Generalized Least Squares Model for a Genome-Wide Quantitative Trait Association Analysis in Families.
\emph{Human Heredity} 2011;71:67-82 (DOI: 10.1159/000324839)

Buse, A:
Goodness of Fit in Generalized Least Squares Estimation
\emph{The American Statistician} 1973;27:106-108
}
\author{
Xiang Li \email{lixxx554@umn.edu}, Robert M. Kirkpatrick \email{kirk0191@umn.edu}, and Saonli Basu \email{saonli@umn.edu} .
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
data(map)
data(pedigree)
data(rescovmtx)
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

bar <- fgls(
  Zscore ~ rs3934834 + IsFemale, data=foo$test.dat, tlist=foo$tlist,
  sizelist=foo$sizelist,med=c("UN","VC"), 
  vmat=rescovmtx, #<--Resid. cov. matrix from fgls onto IsFemale only.
  start=NULL, theta=NULL, drop=NULL, get.hessian=FALSE, 
  optim.method="BFGS", control=list(), weights=NULL,
  sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL)

bar$ctable
##  To simultaneously estimate residual covariance matrix
##  and regression coefficients for rs3934834 & IsFemale,
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
