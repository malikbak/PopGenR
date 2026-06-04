# count.genotypes: Count Genotypes from Genotype Data

This function counts the occurrences of specific genotype patterns
(homozygous reference, heterozygous, homozygous alternate, and missing
data) from a vector of genotype strings.

## Usage

``` r
count.genotypes(genotypes)
```

## Arguments

- genotypes:

  A character vector containing genotype data in the form of \`"0\|0"\`,
  \`"0/1"\`, \`"1\|1"\`, etc. The genotypes can be phased (\`\|\`) or
  unphased (\`/\`).

## Value

A named numeric vector containing the counts of the genotypes \`AA\`
(homozygous reference), \`Aa\` (heterozygous), \`aa\` (homozygous
alternate), and \`NN\` (missing or unknown).

## Details

The function first removes any phase (\`\|\`) or unphase (\`/\`) symbols
in the genotypes, normalizing the genotypes to simple two-character
strings (\`"00"\`, \`"01"\`, \`"10"\`, \`"11"\`, or \`".."\`). It then
counts the occurrences of each pattern. The heterozygous counts for
\`"01"\` and \`"10"\` are combined into a single count for \`Aa\`.

## Note

The function returns a vector with counts for \`AA\`, \`Aa\`, \`aa\`,
and \`NN\`, even if some genotypes are not present in the input.

## Examples

``` r
# Example usage:
genotypes <- c("0|0", "0/1", "1|0", "1|1", "./.")
genotype_counts <- count.genotypes(genotypes)

# Output: named vector with counts for AA, Aa, aa, and NN
```
