# maf: Calculate Minor Allele Frequency (MAF) from VCF Row

This function calculates the minor allele frequency (MAF) for a given
row of variant call format (VCF) data. MAF is the frequency of the least
common allele at a locus and is a measure of genetic variation.

## Usage

``` r
maf(vcf.row)
```

## Arguments

- vcf.row:

  A vector representing a single row of VCF data, where the first 9
  columns contain metadata and genotype information starts from the 10th
  column.

## Value

A numeric value representing the minor allele frequency (MAF) for the
specified VCF row.

## Details

The function extracts genotype data from the VCF row using the
\`get.field\` function. It counts the genotypes with
\`count.genotypes\`, calculates allele frequencies using
\`allele.freq\`, and returns the minimum frequency as the MAF.

## Note

The \`get.field\`, \`count.genotypes\`, and \`allele.freq\` functions
must be available in the user's environment for this function to work
correctly.

## Examples

``` r
# Example usage:
vcf_row <- c("chr1", "123456", ".", "A", "G", ".", "PASS", ".", "GT:DP", "0/1", "1/1", "0/0")
maf_value <- maf(vcf_row)

# Output: Numeric value representing the minor allele frequency
```
