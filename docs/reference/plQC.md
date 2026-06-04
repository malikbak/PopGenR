# Run PLINK Quality Control and Plot PCA

Filter genotype data with PLINK and plot the first two principal
components.

## Usage

``` r
plQC(
  file,
  maf = 0.05,
  gen = 0.1,
  mind = 0.1,
  hwe = 1e-06,
  chrset = 22,
  plink_path = NULL
)
```

## Arguments

- file:

  Input PLINK file prefix, without \`.ped\` or \`.bed\`.

- maf:

  Minor allele frequency threshold.

- gen:

  Variant missingness threshold passed to PLINK \`–geno\`.

- mind:

  Individual missingness threshold passed to PLINK \`–mind\`.

- hwe:

  Hardy-Weinberg equilibrium p-value threshold.

- chrset:

  Number of chromosomes passed to PLINK \`–chr-set\`.

- plink_path:

  Optional path to a PLINK executable. When omitted, the function looks
  for a package/local development copy and then for \`plink\` on the
  system path.

## Value

A \`ggplot\` object. PLINK output files are written to the current
working directory.

## Examples

``` r
if (FALSE) { # \dontrun{
plQC(file = "data/my_genotypes", maf = 0.05, chrset = 30)
} # }
```
