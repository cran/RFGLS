\name{geno}
\alias{geno}
\docType{data}
\title{
 Simulated genotypic dataset
}
\description{
A dataset of simulated genotypes on 10 arbitrary SNPs, for the same simulees in datasets \code{\link{pheno}} and \code{\link{pedigree}}.
}
\usage{data(geno)}
\format{
  A data frame containing only integers, with 1 row per SNP, and 1 column per subject.  The row and column names are rs numbers and individual IDs, respectively.
}

\details{
The genotypes are coded as counts of each SNPs reference allele on the HapMap (\url{http://hapmap.ncbi.nlm.nih.gov/}) positive strand.  First, 10 SNPs, one each from the first 10 human chromosomes, were selected arbitrarily.  Then, genotypes were generated for founders (parents, adoptees, and "independent observations"), under Hardy-Weinberg equilibrium, using the allele frequencies from HapMap's CEU reference data (representing Caucasians of European Ancestry).  After that, genes were "dropped" from parents to offspring.  Subjects' genotypes on the arbitrarily chosen effect locus, rs7681769, were conditioned upon to simulate quantitative phenotype scores (\code{Zscore} in dataset \code{\link{pheno}}).  The true effect size in the data-generating distribution is approximately 0.5\% of phenotypic variance.

Dataset \code{geno} has both row and column names, which is acceptable for a data frame to be provided as argument \option{genfile} to \code{\link{gls.batch}()}.  However, a genotype file to be read from disk should have NEITHER an extra column of row labels nor an extra row of column headers.

Details about each SNP may be found in dataset \code{\link{map}}.
}
\examples{
data(geno)
str(data.frame(t(geno)))
round(cor(t(geno)),3)  ##<--SNPs are on different chromosomes, so no LD.
##Also see examples for functions fgls(), gls.batch(), and gls.batch.get().
}
\keyword{datasets}
