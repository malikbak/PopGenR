#' Calculate Pooled Heterozygosity by Genomic Window
#'
#' Calculate pooled heterozygosity (Hp) from SNP reference-allele counts in
#' fixed-width windows across a VCF file.
#'
#' @param vcf_file Path to a VCF file.
#' @param window_size Window size in base pairs. Defaults to 150,000.
#'
#' @return A data frame with chromosome, genomic range, maximum and minimum
#'   allele counts, Hp, and standardized ZHp values.
#'
#' @examples
#' \dontrun{
#' Hetp("path/to/sample.vcf", window_size = 150000)
#' }
#'
#' @export
Hetp <- function(vcf_file, window_size = 150000) {
  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)
  snp <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
  snp <- snp[, c("CHROM", "POS", "REF")]
  names(snp) <- c("chr", "position", "Ref")
  snp$position <- as.numeric(snp$position)
  snp <- stats::na.omit(snp)

  create_windows <- function(data, window_size) {
    windows <- lapply(unique(data$chr), function(chrom) {
      chrom_data <- data[data$chr == chrom, ]
      starts <- seq(min(chrom_data$position), max(chrom_data$position), by = window_size)
      data.frame(
        chr = chrom,
        start = starts,
        end = starts + window_size - 1,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, windows)
  }

  bed <- create_windows(snp, window_size)
  data.table::setDT(bed)
  data.table::setDT(snp)

  merged_dt <- snp[bed, on = .(chr, position >= start, position <= end), nomatch = 0]
  if (nrow(merged_dt) == 0) {
    return(data.frame())
  }

  merged_dt$range <- paste0(merged_dt$position, "-", merged_dt$position.1)
  count_data <- dplyr::summarise(
    dplyr::group_by(merged_dt, chr, range, Ref),
    count = dplyr::n(),
    .groups = "drop"
  )

  result <- dplyr::summarise(
    dplyr::group_by(count_data, chr, range),
    Max_Value = max(count),
    Min_Value = min(count),
    .groups = "drop"
  )

  result$hp <- (2 * result$Max_Value * result$Min_Value) /
    ((result$Max_Value + result$Min_Value)^2)

  hp_sd <- stats::sd(result$hp)
  result$ZHp <- if (is.na(hp_sd) || hp_sd == 0) {
    NA_real_
  } else {
    (result$hp - mean(result$hp)) / hp_sd
  }

  as.data.frame(result)
}
