\name{map}
\alias{map}
\docType{data}
\title{
 SNP "map" file
}
\description{
A table of information about the 10 SNPs in dataset \code{\link{geno}}, in a format reminiscent of a PLINK binary map (.bim) file. 
}
\usage{data(map)}
\format{
  A data frame with 10 rows, 1 per SNP.  There are 6 labeled columns, for the following variables:
  \describe{
    \item{\code{chr}}{Chromosome numbers}
    \item{\code{SNP}}{A vector of SNP rs numbers, as character strings.}
    \item{\code{position}}{SNP position, in base-pairs, from HapMap Build 36.}
    \item{\code{Ref_Allele}}{Vector of reference alleles, as character strings, from the HapMap positive strand.  Genotypes on each SNP in dataset \code{\link{geno}} are counts of its reference allele.}
    \item{\code{Other_Allele}}{The "other allele."  Also a vector of character strings.}
    \item{\code{Freq_Ref_Allele}}{The relative frequency of the reference allele, in HapMap CEU reference data.}
  }}

\details{
 See documentation for dataset \code{\link{geno}} for details on how the genotypic data on these SNPs was simulated.  See documentation for dataset \code{\link{pheno}} for details on how the effect locus, rs7681769, was used to simulate phenotypic data.
}
\examples{
data(map)
str(map)
map
}
\keyword{datasets}
