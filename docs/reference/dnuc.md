# dnuc: Calculate Nucleotide Diversity (Dxy) Between Two Populations

This function calculates the nucleotide diversity \\D\_{xy}\\ between
two populations using allele frequencies from two variant call format
(VCF) datasets. The measure quantifies the average number of nucleotide
substitutions per site between two populations.

## Usage

``` r
dnuc(vcf1, vcf2, perBP = TRUE)
```

## Arguments

- vcf1:

  A data frame representing the VCF data for the first population, where
  the first 9 columns contain metadata, and the genotype data starts
  from the 10th column.

- vcf2:

  A data frame representing the VCF data for the second population,
  structured similarly to \`vcf1\`.

- perBP:

  A logical value indicating whether to normalize \\D\_{xy}\\ per base
  pair (bp). If \`TRUE\`, the function divides the \\D\_{xy}\\ value by
  the total base pairs between the start and end positions in \`vcf1\`.
  Default is \`TRUE\`.

## Value

A numeric value representing the \\D\_{xy}\\ value. If \`perBP\` is
\`TRUE\`, the value is normalized per base pair, otherwise it returns
the raw \\D\_{xy}\\ value.

## Details

The function works by first extracting the genotype data (\`GT\`) from
the VCF files using the \`get.field\` function. It then calculates the
allele frequencies (\`p\`) for each variant in both populations using
\`allele.freq\`. The \\D\_{xy}\\ measure is computed using the formula:

\$\$D\_{xy} = \sum{p_1(1 - p_2) + p_2(1 - p_1)}\$\$

Where \\p_1\\ and \\p_2\\ are the allele frequencies in population 1 and
population 2, respectively. If \`perBP\` is \`TRUE\`, the function
calculates the total number of base pairs covered by the variants and
returns \\D\_{xy}\\ normalized by this number.

## Note

The \`get.field\` and \`allele.freq\` functions must be defined or
available in the user's environment for this function to work.

## Examples

``` r
# Example usage:
if (FALSE) { # \dontrun{
vcf1 <- read.vcf("population1.vcf")
vcf2 <- read.vcf("population2.vcf")
dxy_value <- dnuc(vcf1, vcf2, perBP = TRUE)
} # }
```
