#' Hetp: Calculate pooled Heterozygosity Per Window from VCF Data
#'
#' This function calculates the heterozygosity per window from a VCF file, analyzing SNPs across specified genomic windows.
#'
#' @param vcf_file A character string representing the path to the VCF file to be analyzed.
#' @param window_size An integer specifying the size of the genomic window (in base pairs) for analysis. The default is set to 150,000.
#'
#' @return A data frame containing the calculated heterozygosity per window, along with maximum and minimum SNP counts for each range.
#'         The returned data frame includes columns for chromosome, range, maximum value, minimum value, heterozygosity (hp), and Z-scores of heterozygosity.
#'
#' @details The function reads a VCF file, processes SNP data to create genomic windows, and counts the number of SNPs within each window. 
#'          It calculates the heterozygosity per window using the formula:
#'          \deqn{hp = \frac{2 \cdot Max\_Value \cdot Min\_Value}{(Max\_Value + Min\_Value)^2}}
#'          It also computes the Z-scores for heterozygosity values.
#'
#' @note The input VCF file must be in a valid format, and the function requires the `data.table` and `dplyr` packages for data manipulation.
#'
#' @examples
#' # Example usage:
#' heterozygosity_results <- Hetp("path/to/sample.vcf", window_size = 150000)
#'
#' # Output: A data frame with heterozygosity calculations for each window
#'
#' @export
Hetp <- function(vcf_file=vcf_file, window_size=150000){
  vcf <- read.vcfR(vcf_file)
  snp <- as.data.frame(vcf@fix)
  snp$CHROM <- as.numeric(snp$CHROM)
  snp <- snp[,c(1,2,4)]
  snp <- na.omit(snp)
  colnames(snp) <- c("chr","position", "Ref")
  snp$position <- as.numeric(snp$position)
  create_windows <- function(data, window_size) {
    unique_chromosomes <- unique(data$CHROM)
    results_list <- list()

    for(chrom in unique_chromosomes) {
      chrom_data <- subset(data, CHROM == chrom)
      win_start <- seq(min(chrom_data$POS), max(chrom_data$POS), by = window_size)
      win_end <- win_start + window_size - 1 # Adjust if you want the window to be exactly window_size wide

      # Create a temporary data frame for this chromosome
      temp_df <- data.frame(Chromosome = chrom, Start = win_start, End = win_end)

      # Count the number of SNPs within each window
      # temp_df$SNPs <- sapply(1:nrow(temp_df), function(i) {
      #   sum(chrom_data$POS >= temp_df$Start[i] & chrom_data$POS <= temp_df$End[i])
      # })

      results_list[[chrom]] <- temp_df
    }

    # Combine all chromosomes' results into one dataframe
    final_results <- do.call(rbind, results_list)
    return(final_results)
  }
  bed <- create_windows(my.data, window_size)
  colnames(bed) <- c("chr","start","end")
  bed$chr <- as.numeric(bed$chr)
  bed <- na.omit(bed)
  setDT(bed)
  setDT(snp)
  setkey(bed, chr)
  setkey(snp, chr)
  # Convbed# Convert to data.table
  # Merge and filter using data.table
  merged_dt <- snp[bed, on = .(chr, position >= start, position <= end)]
  merged_dt$range <- paste0(merged_dt$position,"-",merged_dt$position.1)
  count_data <- merged_dt %>%
    dplyr::group_by(chr,range,Ref)%>%
    dplyr::summarise(count=n())
  count_data <- na.omit(count_data)
  ################
  # Group by 'Factor' column and calculate max and min of 'Value' column
  result <- count_data %>%
    dplyr::group_by(chr,range) %>%
    dplyr::summarize(Max_Value = max(count), Min_Value = min(count))

  hp = (2 * result$Max_Value * result$Min_Value) / ((result$Max_Value + result$Min_Value) ** 2)
  result$hp <- hp
  mean_data <- mean(result$hp)
  sd_data <- sd(result$hp)

  # Calculate Z-scores using the formula
  result$ZHp <- (result$hp - mean_data) / sd_data
  return(result)

}
#results <- Hp(bed_file = "../bosTau9.bed", vcf_file = "../males_out_qc_vcf.vcf")
#
# max(z_scores_formula)
