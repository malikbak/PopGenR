# Create Fixed-Width BED Windows with bedtools

Create genomic windows from a FASTA index or genome file by calling
\`bedtools makewindows\`.

## Usage

``` r
Mbed(Ffile, Window, output = "test.bed", bedtools = "bedtools")
```

## Arguments

- Ffile:

  Path to a genome file accepted by \`bedtools makewindows -g\`.

- Window:

  Window size in base pairs.

- output:

  Output BED file path. Defaults to \`test.bed\`.

- bedtools:

  Path or command name for bedtools.

## Value

Invisibly returns \`output\` after bedtools completes.

## Examples

``` r
if (FALSE) { # \dontrun{
Mbed("genome.fai", Window = 40000, output = "windows.bed")
} # }
```
