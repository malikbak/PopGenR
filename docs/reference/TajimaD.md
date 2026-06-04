# Calculate Tajima's D by Genomic Window

Calculate Tajima's D and supporting diversity statistics in fixed-width
windows from a VCF file.

## Usage

``` r
TajimaD(vfile, window_size = 150000, nsample = NULL)
```

## Arguments

- vfile:

  Path to a VCF file.

- window_size:

  Window size in base pairs. Defaults to 150,000.

- nsample:

  Optional number of samples. When \`NULL\`, the sample count is
  inferred from the VCF genotype columns.

## Value

A data frame with one row per genomic window and columns for SNP count,
chromosome count, Watterson's theta, nucleotide diversity, variance, and
Tajima's D.

## Examples

``` r
if (FALSE) { # \dontrun{
TajimaD("path/to/sample.vcf", window_size = 150000)
} # }
```
