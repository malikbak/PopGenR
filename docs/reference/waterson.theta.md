# waterson.theta: Calculate Watterson's Theta from Genetic Data

This function calculates Watterson's Theta \\\theta\\, a measure of
nucleotide diversity, based on the number of segregating sites in a
given dataset. It can be normalized by the number of base pairs.

## Usage

``` r
waterson.theta(data, perBP = TRUE)
```

## Arguments

- data:

  A data frame containing genetic data in VCF format, where the first 9
  columns contain metadata and genotype data starts from the 10th
  column.

- perBP:

  A logical value indicating whether to normalize Watterson's Theta by
  the number of base pairs. If \`TRUE\`, the result is divided by the
  total number of base pairs; if \`FALSE\`, the raw value is returned.

## Value

A numeric value representing Watterson's Theta. If \`perBP\` is
\`TRUE\`, the result is normalized by the number of base pairs.

## Details

The function calculates the number of segregating sites \\S_n\\ in the
data that have a minor allele frequency greater than zero. It then uses
the formula:

\$\$\theta = \frac{S_n}{\sum\_{i=1}^{2N-1} \frac{1}{i}}\$\$

where \\N\\ is the number of individuals in the sample. If \`perBP\` is
set to \`TRUE\`, the result is normalized by the number of base pairs in
the dataset.

## Note

The \`maf\` function must be available in the user's environment for
this function to work correctly.

## Examples

``` r
# Example usage:
if (FALSE) { # \dontrun{
vcf_data <- read.vcf("path/to/sample.vcf")
watterson_theta <- waterson.theta(vcf_data, perBP = TRUE)
} # }

# Output: A numeric value representing Watterson's Theta normalized by base pairs
```
