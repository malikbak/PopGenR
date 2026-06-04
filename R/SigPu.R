#' Run PLINK Runs of Homozygosity
#'
#' Detect runs of homozygosity (ROH) with PLINK.
#'
#' @param file Input PLINK file prefix, without `.ped` or `.bed`.
#' @param homozygWS PLINK `--homozyg-window-snp` value.
#' @param homozygS PLINK `--homozyg-snp` value.
#' @param homozygWM PLINK `--homozyg-window-missing` value.
#' @param homozygKB PLINK `--homozyg-kb` value.
#' @param homozygDN PLINK `--homozyg-density` value.
#' @param chrset Number of chromosomes passed to PLINK `--chr-set`.
#' @param plink_path Optional path to a PLINK executable.
#'
#' @return A data frame read from `roh.hom`, when PLINK produces it.
#'
#' @examples
#' \dontrun{
#' roh <- RunHom(file = "data/my_genotypes", chrset = 30)
#' }
#'
#' @export
RunHom <- function(file, homozygWS = 50, homozygS = 50, homozygWM = 3,
                   homozygKB = 100, homozygDN = 1000, chrset = 22,
                   plink_path = NULL) {
  plink <- resolve_plink_path(plink_path)
  input_type <- detect_plink_input(file)

  system2(plink, c(
    input_type$flag, file,
    "--homozyg-group",
    "--homozyg-window-snp", homozygWS,
    "--homozyg-snp", homozygS,
    "--homozyg-window-missing", homozygWM,
    "--homozyg-kb", homozygKB,
    "--homozyg-density", homozygDN,
    "--chr-set", chrset,
    "--allow-extra-chr",
    "--out", "roh"
  ))

  if (!file.exists("roh.hom")) {
    stop("PLINK did not produce roh.hom.", call. = FALSE)
  }

  utils::read.table("roh.hom", header = TRUE)
}

#' Calculate Tajima's D by Genomic Window
#'
#' Calculate Tajima's D and supporting diversity statistics in fixed-width
#' windows from a VCF file.
#'
#' @param vfile Path to a VCF file.
#' @param window_size Window size in base pairs. Defaults to 150,000.
#' @param nsample Optional number of samples. When `NULL`, the sample count is
#'   inferred from the VCF genotype columns.
#'
#' @return A data frame with one row per genomic window and columns for SNP
#'   count, chromosome count, Watterson's theta, nucleotide diversity, variance,
#'   and Tajima's D.
#'
#' @examples
#' \dontrun{
#' TajimaD("path/to/sample.vcf", window_size = 150000)
#' }
#'
#' @export
TajimaD <- function(vfile, window_size = 150000, nsample = NULL) {
  vcf <- vcfR::read.vcfR(vfile, verbose = FALSE)
  snp <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
  my.data <- cbind(snp, as.data.frame(vcf@gt, stringsAsFactors = FALSE))
  my.data$POS <- as.integer(my.data$POS)

  windows <- lapply(unique(my.data$CHROM), function(chrom) {
    chrom_data <- my.data[my.data$CHROM == chrom, ]
    starts <- seq(min(chrom_data$POS), max(chrom_data$POS), by = window_size)
    data.frame(
      chr = chrom,
      start = starts,
      end = starts + window_size - 1,
      stringsAsFactors = FALSE
    )
  })
  bed <- do.call(rbind, windows)

  num.samples <- if (!is.null(nsample)) nsample else ncol(my.data) - 9
  bed$SNPs <- 0L
  bed$Nchr <- 2 * num.samples
  bed$ThetaW <- 0
  bed$Pi <- 0

  harmonic <- sum(1 / seq_len((2 * num.samples) - 1))

  for (i in seq_len(nrow(bed))) {
    window_data <- my.data[
      my.data$CHROM == bed$chr[i] &
        my.data$POS >= bed$start[i] &
        my.data$POS <= bed$end[i],
    ]

    bed$SNPs[i] <- nrow(window_data)
    bed$ThetaW[i] <- bed$SNPs[i] / harmonic

    if (nrow(window_data) > 0) {
      j <- apply(window_data, 1, derivedCount)
      c <- rep(bed$Nchr[i], length(j))
      bed$Pi[i] <- sum((2 * j * (c - j)) / (c * (c - 1)), na.rm = TRUE)
    }
  }

  bed$varD <- vapply(seq_len(nrow(bed)), function(i) {
    variance.d(n = bed$Nchr[i], S = bed$SNPs[i])
  }, numeric(1))

  bed$TajimaD <- ifelse(
    bed$varD > 0,
    (bed$Pi - bed$ThetaW) / sqrt(bed$varD),
    NA_real_
  )

  bed
}
