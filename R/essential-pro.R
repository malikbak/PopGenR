#' allele.freq: Calculate Allele Frequencies from Genotype Counts
#'
#' This function computes the allele frequencies for a bi-allelic locus (with alleles `A` and `a`) based on 
#' genotype counts. It calculates the frequency of the dominant allele (`p`) and the recessive allele (`q`).
#'
#' @param genotypeCounts A named numeric vector with counts for genotypes `AA`, `Aa`, `aa`, and potentially `NN` 
#' for missing data. The names of the vector should be `AA`, `Aa`, `aa`, and optionally `NN`.
#'
#' @return A named vector containing the frequencies of the dominant allele `p` and the recessive allele `q`.
#'
#' @details The function first calculates the total number of valid genotypes by subtracting any missing data 
#' (`NN`) from the total genotype count. It then computes the frequency of the dominant allele `A` (`p`) and the 
#' recessive allele `a` (`q`) using the formula:
#' 
#' \deqn{p = \frac{(2 \times \text{AA count}) + \text{Aa count}}{2 \times \text{total count}}}
#' \deqn{q = 1 - p}
#' 
#' @note Ensure that `genotypeCounts` contains valid counts for all three genotypes (`AA`, `Aa`, `aa`) and 
#' optionally `NN` (if missing data is present).
#'
#' @examples
#' # Example usage:
#' genotypes <- c(AA = 50, Aa = 30, aa = 20, NN = 10) # 10 missing values
#' alleleFrequencies <- allele.freq(genotypes)
#' 
#' # If no missing values:
#' genotypes <- c(AA = 50, Aa = 30, aa = 20)
#' alleleFrequencies <- allele.freq(genotypes)
#'
#' @export
allele.freq <- function(genotypeCounts) {
  n = sum(genotypeCounts) - genotypeCounts["NN"]
  p = ((2*genotypeCounts["AA"]) + genotypeCounts["Aa"])/(2*n)
  q = 1-p
  freqs = c(p,q)
  names(freqs) = c("p", "q")
  return(freqs)
}
#' calc_r2: Calculate Linkage Disequilibrium (R²) Between Two Variants
#'
#' This function calculates the linkage disequilibrium (LD) measure, \(R^2\), between two genetic variants based on their genotypes. The function uses the allele frequencies of the two variants and their haplotype frequencies to compute \(R^2\), which quantifies the non-random association of alleles at two loci.
#'
#' @param row1 A vector representing the first variant, where the first 9 columns contain metadata and the genotype data starts from the 10th column onward.
#' @param row2 A vector representing the second variant, structured similarly to `row1`, where the genotype data starts from the 10th column.
#'
#' @return A numeric value representing the \(R^2\) value between the two variants, which measures the degree of linkage disequilibrium between them.
#'
#' @details The function extracts the genotypes (`GT`) for both variants using the helper function `get.field`, which parses the genotype data from the variant vectors. It then calculates the allele frequencies `pA` and `pB` for both variants using the `allele.freq` function. The haplotype frequencies are determined using the `get.haplotypes` function, and from this, the \(R^2\) value is calculated using the formula:
#' 
#' \deqn{D = p_{AB} - (p_A \times p_B)}
#' \deqn{R^2 = \frac{D^2}{p_A \times (1 - p_A) \times p_B \times (1 - p_B)}}
#'
#' @note The `get.field`, `count.genotypes`, `allele.freq`, and `get.haplotypes` functions must be defined or available in the user's environment for this function to work properly.
#'
#' @examples
#' # Example usage:
#' variant1 <- c("chr1", "12345", "A", "T", "...", "...", "...", "...", "...", "0/0", "0/1", "1/1")
#' variant2 <- c("chr1", "67890", "G", "C", "...", "...", "...", "...", "...", "0/0", "0/1", "1/1")
#' r2_value <- calc_r2(variant1, variant2)
#'
#' @export
calc_r2 <- function(row1, row2) {
  g1 = get.field(row1[10:length(row1)], row1[9], "GT")
  g2 = get.field(row2[10:length(row2)], row2[9], "GT")
  pA = unname(allele.freq(count.genotypes(g1))["p"])
  pB = unname(allele.freq(count.genotypes(g2))["p"])
  h = get.haplotypes(g1, g2)
  pAB = (length(h[h=="00"]))/(length(h))
  D = pAB - (pA * pB)
  rsq = (D**2)/(pA*(1-pA)*pB*(1-pB))
  return(rsq)
}
#' count.genotypes: Count Genotypes from Genotype Data
#'
#' This function counts the occurrences of specific genotype patterns (homozygous reference, heterozygous, homozygous alternate, and missing data) from a vector of genotype strings.
#'
#' @param genotypes A character vector containing genotype data in the form of `"0|0"`, `"0/1"`, `"1|1"`, etc. The genotypes can be phased (`|`) or unphased (`/`).
#'
#' @return A named numeric vector containing the counts of the genotypes `AA` (homozygous reference), `Aa` (heterozygous), `aa` (homozygous alternate), and `NN` (missing or unknown).
#'
#' @details The function first removes any phase (`|`) or unphase (`/`) symbols in the genotypes, normalizing the genotypes to simple two-character strings (`"00"`, `"01"`, `"10"`, `"11"`, or `".."`). It then counts the occurrences of each pattern. The heterozygous counts for `"01"` and `"10"` are combined into a single count for `Aa`.
#'
#' @note The function returns a vector with counts for `AA`, `Aa`, `aa`, and `NN`, even if some genotypes are not present in the input.
#'
#' @examples
#' # Example usage:
#' genotypes <- c("0|0", "0/1", "1|0", "1|1", "./.")
#' genotype_counts <- count.genotypes(genotypes)
#' 
#' # Output: named vector with counts for AA, Aa, aa, and NN
#'
#' @export
count.genotypes <- function(genotypes) {
  genotypes = gsub("(\\||/)", "", genotypes) 
  gen.patterns = c("00", "01", "10", "11", "..") 
  my.counts=table(factor(genotypes, levels=gen.patterns)) 
  final.counts = c(my.counts[1], (my.counts[2]+my.counts[3]), my.counts[4:5]) 
  names(final.counts) = c("AA", "Aa", "aa", "NN") 
  return(final.counts)
}
#' derivedCount: Calculate Derived Allele Count from Genotype Data
#'
#' This function calculates the derived allele count for a given genetic variant based on its genotype data. The derived allele count is computed as twice the number of homozygous alternate genotypes (`aa`) plus the number of heterozygous genotypes (`Aa`).
#'
#' @param row A vector representing a genetic variant, where the first 9 columns contain metadata, and the genotype data starts from the 10th column onward. The 9th column should contain the format string for genotype fields (e.g., `"GT"` for genotype).
#'
#' @return A numeric value representing the derived allele count for the variant.
#'
#' @details The function extracts the genotype data (`GT`) from the `row` vector using the `get.field` function. It then uses the `count.genotypes` function to count the occurrences of `AA`, `Aa`, `aa`, and missing genotypes (`NN`). The derived allele count is computed as:
#' 
#' \deqn{\text{Derived Count} = (2 \times \text{aa count}) + \text{Aa count}}
#'
#' @note The `get.field` and `count.genotypes` functions must be defined or available in the user's environment for this function to work.
#'
#' @examples
#' # Example usage:
#' variant <- c("chr1", "12345", "A", "T", "...", "...", "...", "...", "...", "0/0", "0/1", "1/1")
#' derived_allele_count <- derivedCount(variant)
#'
#' @export
derivedCount <- function(row) {
  row=as.vector(row, mode="character")
  x=count.genotypes(get.field(row[10:length(row)], row[9], "GT"))
  dc=(2*x["aa"])+x["Aa"]
  return(unname(dc))
}
#' dnuc: Calculate Nucleotide Diversity (Dxy) Between Two Populations
#'
#' This function calculates the nucleotide diversity (\(D_{xy}\)) between two populations using allele frequencies from two variant call format (VCF) datasets. The measure quantifies the average number of nucleotide substitutions per site between two populations.
#'
#' @param vcf1 A data frame representing the VCF data for the first population, where the first 9 columns contain metadata, and the genotype data starts from the 10th column.
#' @param vcf2 A data frame representing the VCF data for the second population, structured similarly to `vcf1`.
#' @param perBP A logical value indicating whether to normalize \(D_{xy}\) per base pair (bp). If `TRUE`, the function divides the \(D_{xy}\) value by the total base pairs between the start and end positions in `vcf1`. Default is `TRUE`.
#'
#' @return A numeric value representing the \(D_{xy}\) value. If `perBP` is `TRUE`, the value is normalized per base pair, otherwise it returns the raw \(D_{xy}\) value.
#'
#' @details The function works by first extracting the genotype data (`GT`) from the VCF files using the `get.field` function. It then calculates the allele frequencies (`p`) for each variant in both populations using `allele.freq`. The \(D_{xy}\) measure is computed using the formula:
#' 
#' \deqn{D_{xy} = \sum{p_1(1 - p_2) + p_2(1 - p_1)}}
#' 
#' Where \(p_1\) and \(p_2\) are the allele frequencies in population 1 and population 2, respectively. If `perBP` is `TRUE`, the function calculates the total number of base pairs covered by the variants and returns \(D_{xy}\) normalized by this number.
#'
#' @note The `get.field` and `allele.freq` functions must be defined or available in the user's environment for this function to work.
#'
#' @examples
#' # Example usage:
#' vcf1 <- read.vcf("population1.vcf")
#' vcf2 <- read.vcf("population2.vcf")
#' dxy_value <- dnuc(vcf1, vcf2, perBP = TRUE)
#'
#' @export
dnuc <- function(vcf1, vcf2, perBP=TRUE) {
  g1=t(apply(vcf1, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  g2=t(apply(vcf2, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  af1=apply(g1, 1, function(x) allele.freq(count.genotypes(x))["p"])
  af2=apply(g2, 1, function(x) allele.freq(count.genotypes(x))["p"])
  ## Let x be the allele frequency (p) in pop1 * (1-p) in Pop2
  x = af1 * (1-af2)
  ## Let y be the allele frequency (p) in pop2 * (1-p) in Pop1
  y = af2 * (1-af1)
  dxy=sum((x+y))
  if (perBP) {
    c = unique(vcf1$CHROM)
    s = sapply(c, function(x,y) min(y[which(y$CHROM==x),2]), y=vcf1)
    e = sapply(c, function(x,y) max(y[which(y$CHROM==x),2]), y=vcf1)
    bp=sum(e-s)
    return(dxy/bp)
  } else { 
    return(dxy) 
  }
}
#' expected.het: Calculate Expected Heterozygosity from Genotype Data
#'
#' This function calculates the expected heterozygosity (Hexp) for a set of genotypes at a bi-allelic locus. Expected heterozygosity is a measure of genetic diversity and is based on allele frequencies.
#'
#' @param genotypes A character vector containing genotype data, such as `"0|0"`, `"0/1"`, `"1|1"`, etc. The genotypes can be phased (`|`) or unphased (`/`).
#'
#' @return A numeric vector containing three elements: the allele frequency `p` of the dominant allele, the number of non-missing genotypes `n`, and the expected heterozygosity `Hexp`.
#'
#' @details The function first counts the occurrences of different genotypes using the `count.genotypes` function. It then calculates the total number of valid (non-missing) genotypes `n` and the allele frequencies `p` and `q` using the `allele.freq` function. The expected heterozygosity is computed as:
#' 
#' \deqn{H_{exp} = 2 \times p \times q}
#' 
#' Where `p` is the frequency of the dominant allele, and `q = 1 - p` is the frequency of the recessive allele.
#'
#' @note The `count.genotypes` and `allele.freq` functions must be available in the user's environment.
#'
#' @examples
#' # Example usage:
#' genotypes <- c("0|0", "0/1", "1|1", "./.")
#' result <- expected.het(genotypes)
#' 
#' # Output: [1] p (allele frequency), n (number of genotypes), Hexp (expected heterozygosity)
#'
#' @export
expected.het <- function(genotypes) {
  obs.counts = count.genotypes(genotypes)
  n = sum(obs.counts) - obs.counts["NN"]
  freqs = allele.freq(obs.counts)
  Hexp = 2 * freqs[1] * freqs[2]
  res = c(freqs["p"], n, Hexp)
  res=as.numeric(unname(res))
  return(res)
}
#' fixed.poly: Identify Fixed and Polymorphic Sites Between Two Populations
#'
#' This function identifies whether a genetic site is fixed or polymorphic between two populations based on allele frequencies calculated from their genotype data.
#'
#' @param vcf1 A data frame representing the VCF data for the first population, where the first 9 columns contain metadata and the genotype data starts from the 10th column.
#' @param vcf2 A data frame representing the VCF data for the second population, structured similarly to `vcf1`.
#'
#' @return A character vector indicating whether each site is "Fixed" or "Polymorphic" between the two populations.
#'
#' @details The function first extracts the genotype data (`GT`) from both VCF files using the `get.field` function. It then calculates allele frequencies (`p`) for each site in both populations using the `allele.freq` and `count.genotypes` functions. A site is considered "Fixed" if the absolute difference in allele frequencies between the two populations is exactly 1, and "Polymorphic" otherwise.
#'
#' @note The `get.field`, `allele.freq`, and `count.genotypes` functions must be available in the user's environment for this function to work.
#'
#' @examples
#' # Example usage:
#' vcf1 <- read.vcf("population1.vcf")
#' vcf2 <- read.vcf("population2.vcf")
#' result <- fixed.poly(vcf1, vcf2)
#'
#' # Output: A vector indicating "Fixed" or "Polymorphic" for each site
#'
#' @export
fixed.poly <- function(vcf1, vcf2) {
  g1=t(apply(vcf1, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  g2=t(apply(vcf2, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  af1=apply(g1, 1, function(x) allele.freq(count.genotypes(x))["p"])
  af2=apply(g2, 1, function(x) allele.freq(count.genotypes(x))["p"])
  res = rep("Polymorphic", nrow(g1))
  res[which(abs(af1-af2)==1)] = "Fixed"
  return(res)
}
#' get.field: Extract Specific Field from Genotype Data
#'
#' This function extracts a specific field (such as `GT`, `DP`, etc.) from the genotype information provided in the VCF format for each sample.
#'
#' @param samples A character vector containing the genotype data for each sample. Each element in the vector typically contains multiple fields separated by colons (e.g., `"GT:DP:AD"`).
#' @param format A character string specifying the format field in the VCF file, which describes the fields present in the genotype data (e.g., `"GT:DP:AD"`).
#' @param fieldName A character string specifying the name of the field to extract (e.g., `"GT"` for genotype, `"DP"` for read depth).
#'
#' @return A character vector containing the values of the specified field for each sample.
#'
#' @details The function works by first splitting the `samples` vector and the `format` string based on colons (`:`). It identifies the position of `fieldName` within the format and extracts the corresponding values from the `samples` vector. If `fieldName` is not found in the format, the function throws an error.
#'
#' @note This function assumes that all samples follow the same format structure.
#'
#' @examples
#' # Example usage:
#' samples <- c("0/1:35:99", "1/1:20:85")
#' format <- "GT:DP:GQ"
#' field_values <- get.field(samples, format, "GT")
#'
#' # Output: A vector with extracted genotype values: "0/1", "1/1"
#'
#' @export
get.field <- function(samples, format, fieldName) {
  x=strsplit(samples, split=":")
  fields=unlist(strsplit(format, split=":")) 
  i=which(fields==fieldName)
  if (!(fieldName %in% fields)) stop('fieldName not found in format fields') 
  return(sapply(x, `[[`, i)) 
}
#' get.haplotypes: Generate Haplotype Combinations from Genotype Data
#'
#' This function generates haplotypes by combining the alleles from two sets of genotype data. It concatenates the alleles from the two genotypes to produce haplotype pairs for each sample.
#'
#' @param genotypes1 A character vector containing the genotype data for the first set of samples. The genotypes can be phased (`|`) or unphased (`/`).
#' @param genotypes2 A character vector containing the genotype data for the second set of samples, structured similarly to `genotypes1`.
#'
#' @return A character vector of concatenated haplotypes, where each haplotype is formed by combining alleles from `genotypes1` and `genotypes2`.
#'
#' @details The function first removes the phase markers (`|`) from the genotype strings, then splits the combined alleles from both genotype sets into individual alleles. It finally concatenates the alleles from `genotypes1` and `genotypes2` to create haplotypes for each sample.
#'
#' @note The input genotype vectors must have the same length, and corresponding entries should represent the same samples.
#'
#' @examples
#' # Example usage:
#' genotypes1 <- c("0|1", "1|0", "0|0")
#' genotypes2 <- c("1|0", "0|1", "1|1")
#' haplotypes <- get.haplotypes(genotypes1, genotypes2)
#'
#' # Output: A vector of haplotypes: c("01", "10", "00", "10", "01", "11")
#'
#' @export
get.haplotypes <- function(genotypes1, genotypes2) {
  a1 = gsub("\\|", "", genotypes1) 
  a2 = gsub("\\|", "", genotypes2)
  a1=unlist(strsplit(paste0(a1, collapse=""), split="")) 
  a2=unlist(strsplit(paste0(a2, collapse=""), split=""))
  haps = paste0(a1,a2)
  return(haps)
}
#' maf: Calculate Minor Allele Frequency (MAF) from VCF Row
#'
#' This function calculates the minor allele frequency (MAF) for a given row of variant call format (VCF) data. MAF is the frequency of the least common allele at a locus and is a measure of genetic variation.
#'
#' @param vcf.row A vector representing a single row of VCF data, where the first 9 columns contain metadata and genotype information starts from the 10th column.
#'
#' @return A numeric value representing the minor allele frequency (MAF) for the specified VCF row.
#'
#' @details The function extracts genotype data from the VCF row using the `get.field` function. It counts the genotypes with `count.genotypes`, calculates allele frequencies using `allele.freq`, and returns the minimum frequency as the MAF.
#'
#' @note The `get.field`, `count.genotypes`, and `allele.freq` functions must be available in the user's environment for this function to work correctly.
#'
#' @examples
#' # Example usage:
#' vcf_row <- c("chr1", "123456", ".", "A", "G", ".", "PASS", ".", "GT:DP", "0/1", "1/1", "0/0")
#' maf_value <- maf(vcf_row)
#'
#' # Output: Numeric value representing the minor allele frequency
#'
#' @export
maf <- function(vcf.row) {
  temp=as.vector(vcf.row, mode="character")
  af=allele.freq(count.genotypes(get.field(temp[10:length(temp)], temp[9], "GT")))
  maf=min(unname(af))
  return(maf)
}
#' pi.diversity: Calculate Nucleotide Diversity (Pi) from VCF Data
#'
#' This function calculates the nucleotide diversity (Pi) for a given VCF data frame, which quantifies genetic variation at a particular site or across multiple sites in a population.
#'
#' @param vcf A data frame representing the VCF data, where the first 9 columns contain metadata and the genotype data starts from the 10th column.
#' @param perBP A logical value indicating whether to normalize the nucleotide diversity by the number of base pairs. If `TRUE`, the result is divided by the total number of base pairs; if `FALSE`, the raw nucleotide diversity is returned.
#'
#' @return A numeric value representing the nucleotide diversity (Pi). If `perBP` is `TRUE`, the result is normalized by the number of base pairs.
#'
#' @details The function computes the derived count for each site using the `derivedCount` function. It calculates the total number of alleles (`C`) and uses the formula for nucleotide diversity:
#' 
#' \deqn{\pi = \frac{1}{N(N-1)} \sum_{i=1}^{S} J_i (C - J_i)}
#' 
#' where \( J \) is the derived count, and \( N \) is the total number of alleles. If `perBP` is set to `TRUE`, the function normalizes the result by the total number of base pairs across all chromosomes.
#'
#' @note The `derivedCount`, `count.genotypes`, and other necessary functions must be available in the user's environment for this function to work properly.
#'
#' @examples
#' # Example usage:
#' vcf_data <- read.vcf("sample.vcf")
#' nucleotide_diversity <- pi.diversity(vcf_data, perBP = TRUE)
#'
#' # Output: A numeric value representing nucleotide diversity normalized by base pairs
#'
#' @export
pi.diversity <- function(vcf, perBP=TRUE) {
  J=apply(vcf, 1, derivedCount)
  N=apply(vcf, 1, function(x) sum(count.genotypes(get.field(x[10:length(x)], x[9], "GT"))))
  C=2*N
  pi = sum((2*J*(C-J))/(C*(C-1)))
  if (perBP) {
    c = unique(vcf$CHROM)
    s = sapply(c, function(x,y) min(y[which(y$CHROM==x),2]), y=vcf)
    e = sapply(c, function(x,y) max(y[which(y$CHROM==x),2]), y=vcf)
    bp=sum(e-s)
    return(pi/bp)
  } else { return(pi) }
}	      
#' read.vcf: Read and Clean VCF Files
#'
#' This function reads a Variant Call Format (VCF) file, removes specific comment lines, and returns the data as a data frame.
#'
#' @param file A character string representing the path to the VCF file to be read.
#' @param special.char A character string that indicates the type of comment lines to be removed from the VCF file. The default is `"##"` which is typical for VCF files.
#' @param ... Additional arguments to be passed to the `read.table` function for reading the data.
#'
#' @return A data frame containing the cleaned VCF data, with appropriate column names for the genomic data.
#'
#' @details The function reads all lines from the specified VCF file, removes lines starting with the specified comment character (defaulting to `##`), and replaces the header line containing `#CHROM` with `CHROM`. The cleaned lines are then converted into a data frame using `read.table`.
#'
#' @note The input file must be a valid VCF file. The function assumes that the format of the file is correct and consistent with VCF standards.
#'
#' @examples
#' # Example usage:
#' vcf_data <- read.vcf("path/to/sample.vcf")
#'
#' # Output: A data frame representing the cleaned VCF data
#'
#' @export
read.vcf <- function(file, special.char="##", ...) {
  my.search.term=paste0(special.char, ".*")
  all.lines=readLines(file)
  clean.lines=gsub(my.search.term, "",  all.lines)
  clean.lines=gsub("#CHROM", "CHROM", clean.lines)
  read.table(..., text=paste(clean.lines, collapse="\n"))
}
#' variance.d: Calculate Variance of the Genetic Diversity Index (D)
#'
#' This function calculates the variance of the genetic diversity index (D) using the number of samples and the observed diversity.
#'
#' @param n An integer representing the number of samples.
#' @param S A numeric value representing the observed genetic diversity.
#'
#' @return A numeric value representing the variance of the genetic diversity index (D).
#'
#' @details The function uses the following formulas to compute the variance:
#' 
#' \deqn{a1 = \sum_{i=1}^{n-1} \frac{1}{i}}
#' \deqn{a2 = \sum_{i=1}^{n-1} \frac{1}{i^2}}
#' 
#' Then it calculates coefficients \(b1\), \(b2\), \(c1\), and \(c2\) based on these sums, and finally computes the variance using:
#' 
#' \deqn{var = e1 \cdot S + e2 \cdot S \cdot (S - 1)}
#'
#' where \(e1\) and \(e2\) are derived from \(c1\) and \(c2\).
#'
#' @note This function assumes that the input values are valid and that the calculations adhere to the statistical properties of genetic diversity.
#'
#' @examples
#' # Example usage:
#' n_samples <- 10
#' observed_diversity <- 0.5
#' diversity_variance <- variance.d(n_samples, observed_diversity)
#'
#' # Output: A numeric value representing the variance of the genetic diversity index
#'
#' @export
variance.d <- function(n,S) {
  a1=sum(1/(seq(from=1, to=(n-1), by=1)))
  a2=sum(1/((seq(from=1, to=(n-1), by=1))**2))
  b1=(n+1)/(3*(n-1))
  b2=(2*((n**2)+n+3))/((9*n)*(n-1))
  c1=b1 - (1/a1)
  c2=b2-((n+2)/(a1*n)) + (a2/(a1**2))
  e1=c1/a1
  e2=c2/((a1**2)+a2)
  var=(e1*S) + (e2*S*(S-1))
  return(var)
}
#' waterson.theta: Calculate Watterson's Theta from Genetic Data
#'
#' This function calculates Watterson's Theta (\(\theta\)), a measure of nucleotide diversity, based on the number of segregating sites in a given dataset. It can be normalized by the number of base pairs.
#'
#' @param data A data frame containing genetic data in VCF format, where the first 9 columns contain metadata and genotype data starts from the 10th column.
#' @param perBP A logical value indicating whether to normalize Watterson's Theta by the number of base pairs. If `TRUE`, the result is divided by the total number of base pairs; if `FALSE`, the raw value is returned.
#'
#' @return A numeric value representing Watterson's Theta. If `perBP` is `TRUE`, the result is normalized by the number of base pairs.
#'
#' @details The function calculates the number of segregating sites (\(S_n\)) in the data that have a minor allele frequency greater than zero. It then uses the formula:
#' 
#' \deqn{\theta = \frac{S_n}{\sum_{i=1}^{2N-1} \frac{1}{i}}}
#' 
#' where \(N\) is the number of individuals in the sample. If `perBP` is set to `TRUE`, the result is normalized by the number of base pairs in the dataset.
#'
#' @note The `maf` function must be available in the user's environment for this function to work correctly.
#'
#' @examples
#' # Example usage:
#' vcf_data <- read.vcf("path/to/sample.vcf")
#' watterson_theta <- waterson.theta(vcf_data, perBP = TRUE)
#'
#' # Output: A numeric value representing Watterson's Theta normalized by base pairs
#'
#' @export
waterson.theta <- function(data, perBP=TRUE) {
  num.bp=nrow(data)
  check=apply(data, 1, FUN=maf)
  filter=data[which(check>0),]
  Sn=nrow(filter)
  n.ind=ncol(filter)-9
  i=seq(1, ((2*n.ind)-1))
  theta=Sn/sum(1/i)
  if (perBP) { return(theta/num.bp) }
  else { return(theta) }
}