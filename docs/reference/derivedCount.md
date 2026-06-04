# derivedCount: Calculate Derived Allele Count from Genotype Data

This function calculates the derived allele count for a given genetic
variant based on its genotype data. The derived allele count is computed
as twice the number of homozygous alternate genotypes (\`aa\`) plus the
number of heterozygous genotypes (\`Aa\`).

## Usage

``` r
derivedCount(row)
```

## Arguments

- row:

  A vector representing a genetic variant, where the first 9 columns
  contain metadata, and the genotype data starts from the 10th column
  onward. The 9th column should contain the format string for genotype
  fields (e.g., \`"GT"\` for genotype).

## Value

A numeric value representing the derived allele count for the variant.

## Details

The function extracts the genotype data (\`GT\`) from the \`row\` vector
using the \`get.field\` function. It then uses the \`count.genotypes\`
function to count the occurrences of \`AA\`, \`Aa\`, \`aa\`, and missing
genotypes (\`NN\`). The derived allele count is computed as:

\$\$\text{Derived Count} = (2 \times \text{aa count}) + \text{Aa
count}\$\$

## Note

The \`get.field\` and \`count.genotypes\` functions must be defined or
available in the user's environment for this function to work.

## Examples

``` r
# Example usage:
variant <- c("chr1", "12345", ".", "A", "T", ".", "PASS", ".", "GT", "0/0", "0/1", "1/1")
derived_allele_count <- derivedCount(variant)
```
