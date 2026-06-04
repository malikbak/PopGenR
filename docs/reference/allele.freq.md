# allele.freq: Calculate Allele Frequencies from Genotype Counts

This function computes the allele frequencies for a bi-allelic locus
(with alleles \`A\` and \`a\`) based on genotype counts. It calculates
the frequency of the dominant allele (\`p\`) and the recessive allele
(\`q\`).

## Usage

``` r
allele.freq(genotypeCounts)
```

## Arguments

- genotypeCounts:

  A named numeric vector with counts for genotypes \`AA\`, \`Aa\`,
  \`aa\`, and potentially \`NN\` for missing data. The names of the
  vector should be \`AA\`, \`Aa\`, \`aa\`, and optionally \`NN\`.

## Value

A named vector containing the frequencies of the dominant allele \`p\`
and the recessive allele \`q\`.

## Details

The function first calculates the total number of valid genotypes by
subtracting any missing data (\`NN\`) from the total genotype count. It
then computes the frequency of the dominant allele \`A\` (\`p\`) and the
recessive allele \`a\` (\`q\`) using the formula:

\$\$p = \frac{(2 \times \text{AA count}) + \text{Aa count}}{2 \times
\text{total count}}\$\$ \$\$q = 1 - p\$\$

## Note

Ensure that \`genotypeCounts\` contains valid counts for all three
genotypes (\`AA\`, \`Aa\`, \`aa\`) and optionally \`NN\` (if missing
data is present).

## Examples

``` r
# Example usage:
genotypes <- c(AA = 50, Aa = 30, aa = 20, NN = 10) # 10 missing values
alleleFrequencies <- allele.freq(genotypes)

# If no missing values:
genotypes <- c(AA = 50, Aa = 30, aa = 20)
alleleFrequencies <- allele.freq(genotypes)
```
