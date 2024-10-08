% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/essential-pro.R
\name{calc_r2}
\alias{calc_r2}
\title{calc_r2: Calculate Linkage Disequilibrium (R²) Between Two Variants}
\usage{
calc_r2(row1, row2)
}
\arguments{
\item{row1}{A vector representing the first variant, where the first 9 columns contain metadata and the genotype data starts from the 10th column onward.}

\item{row2}{A vector representing the second variant, structured similarly to `row1`, where the genotype data starts from the 10th column.}
}
\value{
A numeric value representing the \(R^2\) value between the two variants, which measures the degree of linkage disequilibrium between them.
}
\description{
This function calculates the linkage disequilibrium (LD) measure, \(R^2\), between two genetic variants based on their genotypes. The function uses the allele frequencies of the two variants and their haplotype frequencies to compute \(R^2\), which quantifies the non-random association of alleles at two loci.
}
\details{
The function extracts the genotypes (`GT`) for both variants using the helper function `get.field`, which parses the genotype data from the variant vectors. It then calculates the allele frequencies `pA` and `pB` for both variants using the `allele.freq` function. The haplotype frequencies are determined using the `get.haplotypes` function, and from this, the \(R^2\) value is calculated using the formula:

\deqn{D = p_{AB} - (p_A \times p_B)}
\deqn{R^2 = \frac{D^2}{p_A \times (1 - p_A) \times p_B \times (1 - p_B)}}
}
\note{
The `get.field`, `count.genotypes`, `allele.freq`, and `get.haplotypes` functions must be defined or available in the user's environment for this function to work properly.
}
\examples{
# Example usage:
variant1 <- c("chr1", "12345", "A", "T", "...", "...", "...", "...", "...", "0/0", "0/1", "1/1")
variant2 <- c("chr1", "67890", "G", "C", "...", "...", "...", "...", "...", "0/0", "0/1", "1/1")
r2_value <- calc_r2(variant1, variant2)

}
