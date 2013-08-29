\name{pheno}
\alias{pheno}
\docType{data}
\title{
 Simulated quantitative-trait dataset
}
\description{
A dataset of observations on a normally distributed phenotype, generated with means conditional on genotype at the effect locus.  All six family-types recognized in package \emph{RFGLS} are represented.  
}

\usage{data(pheno)}

\format{
  A data frame with 4050 observations on the following 6 variables:
  \describe{
    \item{\code{FAMID}}{"Family ID."  Each family in the dataset is uniquely identified by a value of \code{FAMID}, which are all multiples of 10.}
    \item{\code{ID}}{\emph{Individual} ID.  Each subject in the dataset is uniquely identified by a value of \code{ID}, which is equal to his/her \code{INDIV} plus his/her \code{FAMID}.}
    \item{\code{FTYPE}}{"Family-type."  \emph{RFGLS} recognizes six different family-types, five of which are, at largest, four-person nuclear families (two parents, two offspring), distinguished by how the two offspring are related to one another:
      \enumerate{
        \item MZ-twin families,
        \item DZ-twin families,
        \item Adoptive-offspring families,
        \item Non-twin bio-offspring families,
        \item "Mixed" families with one bio and one adopted offspring, and
        \item "Independent observations," who do not fit into a four-person nuclear family.
    }
}
  \item{\code{INDIV}}{"Individual code," which represents how a subject fits into his/her family: \code{INDIV=1} is for "Offspring #1," \code{INDIV=2} is for "Offspring 2," \code{INDIV=3} is for the mother, and \code{INDIV=4} is for the father.  In families of \code{FTYPE=5}, the biological offspring has \code{INDIV=1}, and the adopted offspring, \code{INDIV=2}.  All "independent observations" (i.e., \code{FTYPE=6}) have \code{INDIV=1}.  Note that individuals in a given family are ordered by their \code{INDIV}.}
  \item{\code{Zscore}}{The phenotype score.}
  \item{\code{IsFemale}}{Binary indicator; a value of 1 indicates female sex.  All offspring in families of type #2 (DZ twins) happen to be same-sex, but this is not a requirement.}
}}
  
\details{
Each family's phenotype scores were generated from a multivariate normal distribution (\command{mvrnorm} in package \emph{mvtnorm}) with a centroid defined conditional upon the family members' genotypes on the effect locus (rs7681769 in dataset \code{\link{geno}}), and a variance matrix with 1s on its diagonal and covariances (really, correlations) consistent with an additive heritability of 0.5 and a shared-environmentality of 0.2, but zero assortative mating.
}

\examples{
data(pheno)
str(pheno)
qqnorm(pheno$Zscore[pheno$INDIV==1]) ##<--Normally distributed phenotype.
qqline(pheno$Zscore[pheno$INDIV==1])
##Also see examples for functions fgls(), gls.batch(), and gls.batch.get().
}

\keyword{datasets}
