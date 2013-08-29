\name{FSV.frompedi}
\alias{FSV.frompedi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Family-Structure Variables from Pedigree File
}
\description{
This function creates the family-structure variables \code{"FTYPE"} (family-type) and \code{"INDIV"} (individual code) from available information in a pedigree file.  Note that \command{FSV.frompedi()} is called internally by \code{\link{gls.batch}()} and \code{\link{gls.batch.get}()} when their argument \option{input.mode} is set to 3.
}
\usage{
FSV.frompedi(pedi.dat,phen.dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pedi.dat}{
A pedigree file, as a data frame, with named columns.  Typically, it will contain at least the following five named columns (which correspond to the default for argument \option{pedicolname} to \code{\link{gls.batch}()}): \code{"FAMID"}, (family IDs), \code{"ID"} (unique individual IDs), \code{"PID"} (paternal ID), \code{"MID"} (maternal ID), and \code{"SEX"} (coded 1 for male, 2 for female).  The paternal and maternal IDs of founders must either be 0 or \code{NA}.

Argument \option{pedi.dat} may also contain any/all of the following three named columns, the effects of which are described below under "Details": \code{"ZYGOSITY"}, \code{"ADOPTED"}, and \code{"INDEP"}.  The \code{"ZYGOSITY"} column must contain a value of 1 for each MZ twin, and a value of 2 for each DZ twin.  The \code{"ADOPTED"} column must be a dummy variable for adoptive status, i.e. with value 1 for adoptees and value 0 otherwise (\code{NA}'s are treated as 0).  The \code{"INDEP"} column must be a dummy variable for whether the individual should be treated as an "independent observation" (family-type 6), with 1 for "yes" and 0 for "no" (\code{NA}'s are treated as 0).
}
  \item{phen.dat}{
A phenotype file, as a data frame with named columns.  At the bare minimum, it must contain a column of unique individual IDs, named \code{"ID"}.  The value returned by \command{FSV.frompedi} is this same data frame, with columns named \code{"FTYPE"} and \code{"INDIV"} appended thereto, unless columns with those names were already present, in which case their contents will be overwritten with new values.  Any other named columns in \option{phen.dat} are ignored.
}
}
\details{
\emph{RFGLS} recognizes six recognized family types, which are distinguished primarily by how the offspring in the family are related to one another:
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
The distinction between "Offspring #1" and "#2" is mostly arbitrary, except that in "mixed" families(\code{FTYPE=5}), the biological offspring MUST have \code{INDIV=1}, and the adopted offspring, \code{INDIV=2}.

The way that \command{FSV.frompedi()} assigns family-types and individual codes to participants depends upon the presence/absence of eight named columns in \option{pedi.dat}: \code{"ID"}, \code{"FAMID"}, \code{"PID"}, \code{"MID"}, \code{"SEX"}, \code{"ZYGOSITY"}, \code{"ADOPTED"}, \code{"INDEP"}.  If any of the first five of these are absent, all participants are assigned \code{FTYPE=6} and \code{INDIV=1}, with a warning.  Assuming that those first five columns are present, what \command{FSV.frompedi()} does depends upon the presence/absence of the other three columns, as follows.

If \code{"INDEP"} is present, then \command{FSV.frompedi()} assigns \code{FTYPE=6}, \code{INDIV=1} to participants with \code{INDEP=1}.  These participants are then disregarded for the rest of the job.  Like the other functions in this package, \command{FSV.frompedi()} treats participants with \code{FTYPE=6} as the sole members of their own family units, and not as part of the family corresponding to their family ID.

If \code{"ZYGOSITY"} and \code{"ADOPTED"} are both absent, then (after first checking for \code{"INDEP"}, as above), all participants are assigned \code{FTYPE=4}.  Non-founders are identified as offspring, and participants whose IDs appear in \code{"MID"} or \code{"PID"} are assigned \code{INDIV=3} or \code{INDIV=4}, respectively.  Offspring individual codes are adjusted so that each family has only one instance each of \code{INDIV=1} and \code{INDIV=2}.  If more than two offspring are identified in a family, or if more than one mother or more than one father are identified in family, these participants are forced to \code{FTYPE=6}, \code{INDIV=1}.  Also, any participant not otherwise assigned an individual code is given \code{FTYPE=6}, \code{INDIV=1}.

If \code{"ZYGOSITY"} is present but \code{"ADOPTED"} is absent, then \command{FSV.frompedi()} behaves similarly, except that (after first checking for \code{"INDEP"}, as above) known twins are identified as offspring, and participants belonging to a family containing at least one twin are assigned \code{FTYPE=1} (for MZ) or \code{FTYPE=2} (for DZ), as the case may be.  Member of families with no twins are assigned \code{FTYPE=4}.  The program then proceeds as described in the immediately preceding paragraph.

If \code{"ADOPTED"} is present, \command{FSV.frompedi()} first makes some simple family-type assignments: if \code{"ZYGOSITY"} is present, to \code{FTYPE=1} and \code{FTYPE=2} as appropriate (see above), and then if \code{"INDEP"} is present, to \code{FTYPE=6}, \code{INDIV=1} as appropriate (see above).  Then, within each family, the program resolves each member in order of ID, from least to greatest.  The first non-founder is assigned \code{INDIV=1}, the second, \code{INDIV=2}, and any thereafter, \code{FTYPE=6}, \code{INDIV=1}.  The first adoptee is assigned \code{INDIV=2}, the second, \code{INDIV=1}, and any thereafter, \code{FTYPE=6}, \code{INDIV=1}.  The first female non-adoptee non-founder is assigned \code{INDIV=3}, and any others are assigned \code{FTYPE=6}, \code{INDIV=1}.  The first male non-adoptee non-founder is assigned \code{INDIV=4}, and any others are assigned \code{FTYPE=6}, \code{INDIV=1}.  If family-type has not yet been assigned, then it is resolved as \code{FTYPE=3} if there are two adoptees, \code{FTYPE=5} if there is one adoptee and one biological offspring, and as \code{FTYPE=4} otherwise.

Function \command{FSV.frompedi()} produces a warning whenever it forces a non-founder to \code{FTYPE=6}, \code{INDIV=1}.

Note that there is definitely a degree of arbitrariness in how ambiguous cases are resolved, in that \command{FSV.frompedi()} scans through the pedigree file from top to bottom \emph{after it has sorted the file} by family ID, and by ID within the same family.  So for example, if two participants in the same family are both provisionally assigned \code{INDIV=3},  then the apparent mother with the smaller ID retains \code{INDIV=3}, and the other is forced to \code{FTYPE=6}, \code{INDIV=1}.
}

\value{
A data frame, containing the same columns as \option{phen.dat}, with the addition of \code{"FTYPE"} and \code{"INDIV"}.  Usually, this data frame will simply be \option{phen.dat} with \code{"FTYPE"} and \code{"INDIV"} appended thereto.  However, if \option{phen.dat} contained columns named \code{"FTYPE"} or \code{"INDIV"}, the values in these columns will be overwritten with the new values produced by \command{FSV.frompedi()}.
}
\author{
Robert M. Kirkpatrick \email{kirk0191@umn.edu}.
}

\seealso{
\code{\link{gls.batch}}, \code{\link{gls.batch.get}}
}
\examples{
data(pheno)
data(pedigree)
table(pheno$FTYPE) ##<--Frequencies of correct family types.


fsvtest1 <- FSV.frompedi(pedi.dat=pedigree,
  phen.dat=data.frame(ID=pheno[,2])) ##<--Bare minimum phenotype file.
table(fsvtest1$FTYPE) ##<--Not correct, because pedigree file
                      ##doesn't have enough additional info
                      ##to recover the actual family-types
                      ##and individual codes.

#Create "ZYGOSITY" column:
pedigree$ZYGOSITY <- NA
pedigree$ZYGOSITY[pheno$FTYPE==1 & pheno$INDIV<3] <- 1
pedigree$ZYGOSITY[pheno$FTYPE==2 & pheno$INDIV<3] <- 2

fsvtest2 <- FSV.frompedi(pedi.dat=pedigree,phen.dat=data.frame(ID=pheno[,2]))
table(fsvtest2$FTYPE) ##<--Still not right, because pedigree file
                      ##lacks info about adoptees.
                      
#Create "ADOPTED" column:
pedigree$ADOPTED <- 0
pedigree$ADOPTED[pheno$FTYPE==3 & pheno$INDIV<3] <- 1
pedigree$ADOPTED[pheno$FTYPE==5 & pheno$INDIV==2] <- 1
fsvtest3 <- FSV.frompedi(pedi.dat=pedigree,phen.dat=data.frame(ID=pheno[,2]))
table(fsvtest3$FTYPE) ##<--Almost there.

#Create "INDEP" column:
pedigree$INDEP <- 0
pedigree$INDEP[pheno$FTYPE==6] <- 1
fsvtest4 <- FSV.frompedi(pedi.dat=pedigree,phen.dat=data.frame(ID=pheno[,2]))
table(fsvtest4$FTYPE) ##<--Correct family types have been recovered.
table(pheno$FTYPE) ##<--Compare.
all(pheno$FTYPE==fsvtest4$FTYPE) ##<--TRUE.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
