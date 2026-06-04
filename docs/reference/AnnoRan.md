# Annotate Genomic Regions with Gene Names

Match genomic intervals, such as runs of homozygosity, to gene
annotations.

## Usage

``` r
AnnoRan(data, CHR, starP, endP, annotation)
```

## Arguments

- data:

  A data frame containing the regions to annotate.

- CHR:

  Column name or index in \`data\` containing chromosome identifiers.

- starP:

  Column name or index in \`data\` containing region start positions.

- endP:

  Column name or index in \`data\` containing region end positions.

- annotation:

  Either a data frame of gene annotations or a path to a GFF/GTF file.
  Data frames should contain \`chromosome_name\`, \`start_position\`,
  \`end_position\`, and \`external_gene_name\`.

## Value

A data frame with \`ID\`, \`Range\`, and \`Matching_Names\`.
\`Matching_Names\` contains comma-separated gene names fully contained
within each region, or \`"NA"\` when no annotation overlaps.

## Examples

``` r
regions <- data.frame(chr = c("1", "1"), start = c(1000, 5000), end = c(2000, 6500))
genes <- data.frame(
  chromosome_name = c("1", "1"),
  start_position = c(1200, 5500),
  end_position = c(1800, 6000),
  external_gene_name = c("GeneA", "GeneB")
)
AnnoRan(regions, CHR = "chr", starP = "start", endP = "end", annotation = genes)
#>   ID     Range Matching_Names
#> 1  1 1000-2000          GeneA
#> 2  1 5000-6500          GeneB
```
