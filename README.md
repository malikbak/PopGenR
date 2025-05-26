# PopGenR

PopGenR is an R package designed for population genetics analysis. It provides a suite of tools to detect and analyze signatures of selection in genomic data, utilizing methods such as Fst, Tajima's D, integrated Haplotype Score (iHS), and pooled heterozygosity (Hp). Additionally, PopGenR offers functionalities for annotating various types of genomic data, including Regions of Homozygosity (ROH) and individual SNP data, to facilitate comprehensive genetic studies.

## Core Functionalities

PopGenR provides a range of tools for population genetic analysis:

*   **Detecting Signatures of Selection:**
    *   Calculation of Fst values to assess population differentiation.
    *   Computation of Tajima's D to infer evolutionary neutrality.
    *   Determination of integrated Haplotype Score (iHS) for detecting recent positive selection.
    *   `Hetp()`: Calculation of pooled heterozygosity (Hp) across genomic windows.
*   **Genomic Annotation:**
    *   `AnnoRan()`: Annotate genomic regions (e.g., Runs of Homozygosity) with corresponding gene information from GFF/GTF files.
    *   `AnnoSin()`: Annotate individual genomic positions (SNPs) with gene information.
*   **Population Genetic Statistics:**
    *   `pi.diversity()`: Calculate nucleotide diversity.
    *   `waterson.theta()`: Estimate Waterson's theta.
    *   `allele.freq()`: Calculate allele frequencies.
    *   `maf()`: Calculate minor allele frequencies.
*   **Data Input/Handling:**
    *   Efficiently processes VCF files for various analyses.

## System Requirements

*   R version 4.0.0 or higher.

## Dependencies

The package imports the following R packages:
*   `vcfR`
*   `ggplot2`
*   `data.table`
*   `dplyr`

These packages will be automatically installed if they are not already present when installing `PopGenR` using `devtools`.

## Installation

You can install the development version of PopGenR from GitHub using the `devtools` package:

```R
devtools::install_git("https://github.com/malikbak/PopGenR.git")
```

If you don't have `devtools` installed, you'll need to install it first:
```R
install.packages("devtools")
```
The installation process should also install any missing dependencies listed in the 'Dependencies' section.

## Basic Examples

Below are some basic examples of how to use PopGenR. For more detailed examples and parameter options, please refer to the documentation for each function (e.g., `?Hetp`, `?AnnoRan`).

### Calculating Pooled Heterozygosity (Hp)

The `Hetp` function calculates pooled heterozygosity across genomic windows from a VCF file.

```R
# Ensure PopGenR is loaded
library(PopGenR)

# Path to your VCF file
vcf_file_path <- "path/to/your/data.vcf" # Replace with your actual VCF file path

# Calculate Hp with a window size of 100kb
# Note: Ensure your VCF file is properly formatted and accessible.
# This is a conceptual example; you may need to load example data or use your own.
# heterozygosity_results <- Hetp(vcf_file = vcf_file_path, window_size = 100000)

# View the results (conceptual)
# print(head(heterozygosity_results))
```
*Note: The example above is conceptual. You'll need to provide a valid path to a VCF file. The `extdata` directory might contain sample data or you might need to guide users to prepare their own.*

### Annotating Genomic Regions (e.g., ROH)

The `AnnoRan` function annotates genomic regions (like Runs of Homozygosity - ROH) with gene information from a GFF/GTF file or a suitable data frame.

```R
# Ensure PopGenR is loaded
library(PopGenR)

# Example ROH data (replace with your actual data or load from your analysis)
roh_data <- data.frame(
  CHR = c("1", "1", "2"), # Chromosome
  STARP = c(10000, 50000, 20000), # Start position
  ENDP = c(30000, 70000, 40000), # End position
  IND = c("Ind1", "Ind1", "Ind2") # Individual ID
)

# Example annotation data (replace with your GFF/GTF file path or a prepared data frame)
# If using a file:
# annotation_source <- "path/to/your/annotation.gff"
# Or, an example annotation data frame (column names should match GFF/GTF structure if converted):
annotation_df <- data.frame(
  chromosome_name = c("1", "1", "2", "X"),
  start_position = c(15000, 55000, 25000, 10000),
  end_position = c(25000, 65000, 35000, 20000),
  feature_type = c("gene", "gene", "gene", "gene"), # GFF/GTF typically has a feature type column
  external_gene_name = c("GeneA", "GeneB", "GeneC", "GeneX"), # Attribute column for gene names
  stringsAsFactors = FALSE
)

# Annotate the ROH data using the example data frame
# The function will look for columns like 'chromosome_name', 'start_position', 'end_position', 
# and an attribute column (often 'external_gene_name' or similar for gene names) in the annotation.
# For `AnnoRan`, the input `data` needs CHR, STARP, ENDP columns.
# conceptual_annotated_roh <- AnnoRan(
#   data = roh_data,
#   annotation = annotation_df,
#   CHR = "CHR", # Name of the chromosome column in roh_data
#   starP = "STARP",  # Name of the start position column in roh_data
#   endP = "ENDP"    # Name of the end position column in roh_data
# )

# View the annotated results (conceptual)
# print(head(conceptual_annotated_roh))
```
*Note: The `AnnoRan` example uses a conceptual data frame for annotation. For actual use, you would provide a path to a GFF/GTF file or a data frame prepared from such a file. Ensure column names in your input data and annotation source are correctly specified to the function if they differ from defaults.*

## License

PopGenR is licensed under the MIT License. See the [LICENSE.txt](LICENSE.txt) file for more details.

## Authors

*   Abu Bakar
*   Waqar ul Haq

Maintainer: The package maintainer <malikabubakar279@gmail.com>

## Contributing

Contributions to PopGenR are welcome! If you have suggestions for improvements, bug fixes, or new features, please feel free to:
*   Open an issue on the GitHub repository.
*   Submit a pull request.
