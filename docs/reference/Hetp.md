# Calculate Pooled Heterozygosity by Genomic Window

Calculate pooled heterozygosity (Hp) from SNP reference-allele counts in
fixed-width windows across a VCF file.

## Usage

``` r
Hetp(vcf_file, window_size = 150000)
```

## Arguments

- vcf_file:

  Path to a VCF file.

- window_size:

  Window size in base pairs. Defaults to 150,000.

## Value

A data frame with chromosome, genomic range, maximum and minimum allele
counts, Hp, and standardized ZHp values.

## Examples

``` r
if (FALSE) { # \dontrun{
Hetp("path/to/sample.vcf", window_size = 150000)
} # }
```
