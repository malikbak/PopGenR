# Annotate Single Genomic Positions with Gene Names

Match SNP or marker positions to gene annotations.

## Usage

``` r
AnnoSin(data, annotation)
```

## Arguments

- data:

  A data frame with \`CHR\` and \`POSITION\` columns.

- annotation:

  A data frame with \`chromosome_name\`, \`start_position\`,
  \`end_position\`, and \`external_gene_name\` columns.

## Value

The input data with an added \`external_gene_name\` column.

## Examples

``` r
positions <- data.frame(CHR = c("1", "2"), POSITION = c(1500, 6000))
genes <- data.frame(
  chromosome_name = c("1", "2"),
  start_position = c(1000, 5000),
  end_position = c(2000, 6500),
  external_gene_name = c("GeneA", "GeneB")
)
AnnoSin(positions, genes)
#>   CHR POSITION external_gene_name
#> 1   1     1500              GeneA
#> 2   2     6000              GeneB
```
