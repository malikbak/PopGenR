# get.field: Extract Specific Field from Genotype Data

This function extracts a specific field (such as \`GT\`, \`DP\`, etc.)
from the genotype information provided in the VCF format for each
sample.

## Usage

``` r
get.field(samples, format, fieldName)
```

## Arguments

- samples:

  A character vector containing the genotype data for each sample. Each
  element in the vector typically contains multiple fields separated by
  colons (e.g., \`"GT:DP:AD"\`).

- format:

  A character string specifying the format field in the VCF file, which
  describes the fields present in the genotype data (e.g.,
  \`"GT:DP:AD"\`).

- fieldName:

  A character string specifying the name of the field to extract (e.g.,
  \`"GT"\` for genotype, \`"DP"\` for read depth).

## Value

A character vector containing the values of the specified field for each
sample.

## Details

The function works by first splitting the \`samples\` vector and the
\`format\` string based on colons (\`:\`). It identifies the position of
\`fieldName\` within the format and extracts the corresponding values
from the \`samples\` vector. If \`fieldName\` is not found in the
format, the function throws an error.

## Note

This function assumes that all samples follow the same format structure.

## Examples

``` r
# Example usage:
samples <- c("0/1:35:99", "1/1:20:85")
format <- "GT:DP:GQ"
field_values <- get.field(samples, format, "GT")

# Output: A vector with extracted genotype values: "0/1", "1/1"
```
