# PopGenR

PopGenR is an R package for common population genetics workflows on VCF and
PLINK-style genotype data. It includes helpers for allele frequencies,
heterozygosity, nucleotide diversity, Watterson's theta, Tajima's D, linkage
disequilibrium, PLINK quality control, runs of homozygosity, and annotation of
genomic regions or single positions with gene names.

## Features

- Read VCF files into analysis-ready data frames with `read.vcf()`.
- Summarize genotype counts, allele frequencies, minor allele frequency, and
  expected heterozygosity.
- Calculate diversity statistics such as nucleotide diversity, Dxy,
  Watterson's theta, and Tajima's D.
- Calculate pooled heterozygosity (`Hetp()`) across genomic windows.
- Calculate pairwise linkage disequilibrium (`calc_r2()`).
- Run PLINK QC/PCA (`plQC()`) and runs of homozygosity (`RunHom()`).
- Annotate ranges (`AnnoRan()`) and single positions (`AnnoSin()`) with gene
  annotation data.

## Installation

Install the development version from GitHub:

```r
install.packages("devtools")
devtools::install_git("https://github.com/malikbak/PopGenR.git")
```

PopGenR imports `data.table`, `dplyr`, `ggplot2`, and `vcfR`. Reading GFF/GTF
annotation files with `AnnoRan()` also requires `rtracklayer`, which is a
Bioconductor package:

```r
install.packages("BiocManager")
BiocManager::install("rtracklayer")
```

## Quick Start

```r
library(PopGenR)

genotypes <- c("0/0", "0/1", "1/0", "1/1", "./.")

counts <- count.genotypes(genotypes)
counts

allele.freq(counts)
expected.het(genotypes)
```

Annotate genomic intervals with an annotation data frame:

```r
regions <- data.frame(
  chr = c("1", "1"),
  start = c(1000, 5000),
  end = c(2000, 6500)
)

genes <- data.frame(
  chromosome_name = c("1", "1"),
  start_position = c(1200, 5500),
  end_position = c(1800, 6000),
  external_gene_name = c("GeneA", "GeneB")
)

AnnoRan(regions, CHR = "chr", starP = "start", endP = "end", annotation = genes)
```

Run window-based statistics on a VCF:

```r
tajima <- TajimaD("path/to/sample.vcf", window_size = 150000)
hp <- Hetp("path/to/sample.vcf", window_size = 150000)
```

Run PLINK QC and PCA:

```r
p <- plQC(file = "path/to/plink_prefix", maf = 0.05, chrset = 30)
p
```

The PLINK helpers can use a local development copy of `extdata/plink/plink.exe`
when present, or you can pass `plink_path = "path/to/plink"` to use your own
PLINK installation.

## Developer Workflow

Regenerate function documentation:

```r
roxygen2::roxygenise()
```

Build and check the package:

```r
devtools::document()
devtools::check()
```

Build the package manual:

```powershell
& "C:\Program Files\R\R-4.6.0\bin\R.exe" CMD Rd2pdf . --output=PopGenR-manual.pdf
```

Build the pkgdown site:

```r
pkgdown::build_site()
```

The generated site is written to `docs/`, which can be published with GitHub
Pages.

## License

PopGenR is distributed under the MIT License. See `LICENSE.txt` for details.

## Authors

- Abu Bakar
- Waqar ul Haq
