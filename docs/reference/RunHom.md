# Run PLINK Runs of Homozygosity

Detect runs of homozygosity (ROH) with PLINK.

## Usage

``` r
RunHom(
  file,
  homozygWS = 50,
  homozygS = 50,
  homozygWM = 3,
  homozygKB = 100,
  homozygDN = 1000,
  chrset = 22,
  plink_path = NULL
)
```

## Arguments

- file:

  Input PLINK file prefix, without \`.ped\` or \`.bed\`.

- homozygWS:

  PLINK \`–homozyg-window-snp\` value.

- homozygS:

  PLINK \`–homozyg-snp\` value.

- homozygWM:

  PLINK \`–homozyg-window-missing\` value.

- homozygKB:

  PLINK \`–homozyg-kb\` value.

- homozygDN:

  PLINK \`–homozyg-density\` value.

- chrset:

  Number of chromosomes passed to PLINK \`–chr-set\`.

- plink_path:

  Optional path to a PLINK executable.

## Value

A data frame read from \`roh.hom\`, when PLINK produces it.

## Examples

``` r
if (FALSE) { # \dontrun{
roh <- RunHom(file = "data/my_genotypes", chrset = 30)
} # }
```
