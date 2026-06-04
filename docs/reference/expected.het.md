# expected.het: Calculate Expected Heterozygosity from Genotype Data

This function calculates the expected heterozygosity (Hexp) for a set of
genotypes at a bi-allelic locus. Expected heterozygosity is a measure of
genetic diversity and is based on allele frequencies.

## Usage

``` r
expected.het(genotypes)
```

## Arguments

- genotypes:

  A character vector containing genotype data, such as \`"0\|0"\`,
  \`"0/1"\`, \`"1\|1"\`, etc. The genotypes can be phased (\`\|\`) or
  unphased (\`/\`).

## Value

A numeric vector containing three elements: the allele frequency \`p\`
of the dominant allele, the number of non-missing genotypes \`n\`, and
the expected heterozygosity \`Hexp\`.

## Details

The function first counts the occurrences of different genotypes using
the \`count.genotypes\` function. It then calculates the total number of
valid (non-missing) genotypes \`n\` and the allele frequencies \`p\` and
\`q\` using the \`allele.freq\` function. The expected heterozygosity is
computed as:

\$\$H\_{exp} = 2 \times p \times q\$\$

Where \`p\` is the frequency of the dominant allele, and \`q = 1 - p\`
is the frequency of the recessive allele.

## Note

The \`count.genotypes\` and \`allele.freq\` functions must be available
in the user's environment.

## Examples

``` r
# Example usage:
genotypes <- c("0|0", "0/1", "1|1", "./.")
result <- expected.het(genotypes)

# Output: [1] p (allele frequency), n (number of genotypes), Hexp (expected heterozygosity)
```
