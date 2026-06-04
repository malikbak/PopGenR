# get.haplotypes: Generate Haplotype Combinations from Genotype Data

This function generates haplotypes by combining the alleles from two
sets of genotype data. It concatenates the alleles from the two
genotypes to produce haplotype pairs for each sample.

## Usage

``` r
get.haplotypes(genotypes1, genotypes2)
```

## Arguments

- genotypes1:

  A character vector containing the genotype data for the first set of
  samples. The genotypes can be phased (\`\|\`) or unphased (\`/\`).

- genotypes2:

  A character vector containing the genotype data for the second set of
  samples, structured similarly to \`genotypes1\`.

## Value

A character vector of concatenated haplotypes, where each haplotype is
formed by combining alleles from \`genotypes1\` and \`genotypes2\`.

## Details

The function first removes the phase markers (\`\|\`) from the genotype
strings, then splits the combined alleles from both genotype sets into
individual alleles. It finally concatenates the alleles from
\`genotypes1\` and \`genotypes2\` to create haplotypes for each sample.

## Note

The input genotype vectors must have the same length, and corresponding
entries should represent the same samples.

## Examples

``` r
# Example usage:
genotypes1 <- c("0|1", "1|0", "0|0")
genotypes2 <- c("1|0", "0|1", "1|1")
haplotypes <- get.haplotypes(genotypes1, genotypes2)

# Output: A vector of haplotypes: c("01", "10", "00", "10", "01", "11")
```
