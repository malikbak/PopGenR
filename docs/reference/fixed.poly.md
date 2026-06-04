# fixed.poly: Identify Fixed and Polymorphic Sites Between Two Populations

This function identifies whether a genetic site is fixed or polymorphic
between two populations based on allele frequencies calculated from
their genotype data.

## Usage

``` r
fixed.poly(vcf1, vcf2)
```

## Arguments

- vcf1:

  A data frame representing the VCF data for the first population, where
  the first 9 columns contain metadata and the genotype data starts from
  the 10th column.

- vcf2:

  A data frame representing the VCF data for the second population,
  structured similarly to \`vcf1\`.

## Value

A character vector indicating whether each site is "Fixed" or
"Polymorphic" between the two populations.

## Details

The function first extracts the genotype data (\`GT\`) from both VCF
files using the \`get.field\` function. It then calculates allele
frequencies (\`p\`) for each site in both populations using the
\`allele.freq\` and \`count.genotypes\` functions. A site is considered
"Fixed" if the absolute difference in allele frequencies between the two
populations is exactly 1, and "Polymorphic" otherwise.

## Note

The \`get.field\`, \`allele.freq\`, and \`count.genotypes\` functions
must be available in the user's environment for this function to work.

## Examples

``` r
# Example usage:
if (FALSE) { # \dontrun{
vcf1 <- read.vcf("population1.vcf")
vcf2 <- read.vcf("population2.vcf")
result <- fixed.poly(vcf1, vcf2)
} # }

# Output: A vector indicating "Fixed" or "Polymorphic" for each site
```
