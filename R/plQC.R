#' Run PLINK Quality Control and Plot PCA
#'
#' Filter genotype data with PLINK and plot the first two principal components.
#'
#' @param file Input PLINK file prefix, without `.ped` or `.bed`.
#' @param maf Minor allele frequency threshold.
#' @param gen Variant missingness threshold passed to PLINK `--geno`.
#' @param mind Individual missingness threshold passed to PLINK `--mind`.
#' @param hwe Hardy-Weinberg equilibrium p-value threshold.
#' @param chrset Number of chromosomes passed to PLINK `--chr-set`.
#' @param plink_path Optional path to a PLINK executable. When omitted, the
#'   function looks for a package/local development copy and then for `plink`
#'   on the system path.
#'
#' @return A `ggplot` object. PLINK output files are written to the current
#'   working directory.
#'
#' @examples
#' \dontrun{
#' plQC(file = "data/my_genotypes", maf = 0.05, chrset = 30)
#' }
#'
#' @export
plQC <- function(file, maf = 0.05, gen = 0.1, mind = 0.1, hwe = 1e-6,
                 chrset = 22, plink_path = NULL) {
  plink <- resolve_plink_path(plink_path)
  input_type <- detect_plink_input(file)

  qc_args <- c(
    input_type$flag, file,
    "--hwe", hwe,
    "--chr-set", chrset,
    "--maf", maf,
    "--mind", mind,
    "--geno", gen,
    "--allow-extra-chr",
    "--make-bed",
    "--out", "qc_filter"
  )
  system2(plink, qc_args)

  system2(plink, c(
    "--bfile", "qc_filter",
    "--pca",
    "--chr-set", chrset,
    "--allow-extra-chr",
    "--out", "pca"
  ))

  if (!file.exists("pca.eigenvec")) {
    stop("PLINK did not produce pca.eigenvec.", call. = FALSE)
  }

  pca <- utils::read.table("pca.eigenvec", sep = " ", header = FALSE)
  pve <- if (file.exists("pca.eigenval")) {
    eigenvalues <- scan("pca.eigenval", quiet = TRUE)
    eigenvalues / sum(eigenvalues) * 100
  } else {
    c(NA_real_, NA_real_)
  }

  pc1_label <- if (is.na(pve[1])) "PC1" else paste0("PC1 (", signif(pve[1], 3), "%)")
  pc2_label <- if (is.na(pve[2])) "PC2" else paste0("PC2 (", signif(pve[2], 3), "%)")

  ggplot2::ggplot(pca, ggplot2::aes(V3, V4)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::coord_equal() +
    ggplot2::theme_light() +
    ggplot2::xlab(pc1_label) +
    ggplot2::ylab(pc2_label)
}

resolve_plink_path <- function(plink_path = NULL) {
  if (!is.null(plink_path)) {
    if (!file.exists(plink_path)) {
      stop("`plink_path` does not exist: ", plink_path, call. = FALSE)
    }
    return(plink_path)
  }

  packaged <- system.file("extdata", "plink", "plink.exe", package = "PopGenR", mustWork = FALSE)
  if (nzchar(packaged) && file.exists(packaged)) {
    return(packaged)
  }

  local <- file.path("extdata", "plink", "plink.exe")
  if (file.exists(local)) {
    return(local)
  }

  path_plink <- Sys.which("plink")
  if (nzchar(path_plink)) {
    return(path_plink)
  }

  stop("Could not find PLINK. Provide `plink_path`.", call. = FALSE)
}

detect_plink_input <- function(file) {
  if (file.exists(paste0(file, ".ped"))) {
    return(list(flag = "--file", type = "ped"))
  }

  if (file.exists(paste0(file, ".bed"))) {
    return(list(flag = "--bfile", type = "bed"))
  }

  stop("Could not find PLINK input files for prefix: ", file, call. = FALSE)
}
