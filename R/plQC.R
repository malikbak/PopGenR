#' Perform Quality Control and PCA Analysis using PLINK
#'
#' This function performs quality control (QC) on genetic data using PLINK and calculates principal components 
#' for the provided genotype data. It generates a PCA plot based on the resulting principal components.
#'
#' @param file A character string specifying the path to the input file (without extension). The function accepts both 
#'             `.ped` and `.bed` file formats.
#' @param maf A numeric value representing the minor allele frequency threshold for filtering. Default is 0.05.
#' @param gen A numeric value representing the genotype missingness threshold. Default is 0.1.
#' @param mind A numeric value representing the individual missingness threshold. Default is 0.1.
#' @param hwe A numeric value representing the p-value threshold for Hardy-Weinberg equilibrium test. Default is 1e-6.
#' @param chrset An integer representing the chromosome set to be analyzed. Default is 22.
#'
#' @return A PCA plot showing the first two principal components (PC1 and PC2) of the filtered genotype data.
#'
#' @details The function uses the PLINK tool to perform QC filtering on the input genetic data and calculates 
#'          the principal components. The PCA plot visualizes the distribution of the samples based on the first 
#'          two principal components, highlighting potential clusters or outliers in the data.
#'
#' @note Ensure that the PLINK executable is accessible at the specified path in the `extdata/plink/` directory.
#'       The input files must be in the correct format, and the function assumes the presence of `ggplot2` for plotting.
#'
#' @examples
#' # Example usage:
#' plQC(file = "path/to/data", maf = 0.05, gen = 0.1, mind = 0.1, hwe = 1e-6, chrset = 22)
#'
#' @export
plQC <- function(file=file, maf = 0.05, gen = 0.1, mind = 0.1 ,hwe=1e-6,chrset = 22){
  if(Sys.info()['sysname'] == "Windows"){
    if(paste0(file,".ped") %in% list.files(pattern = "*.ped")){
      system(glue::glue("extdata/plink/plink.exe --file {file} --hwe {hwe} --chr-set {chrset} --maf {maf} --mind {mind} --geno {gen} --allow-extra-chr --make-bed --out qc_filter"))
      system(glue::glue("extdata/plink/plink.exe --bfile qc_filter --pca --chr-set {chrset} --allow-extra-chr --out pca" ))
      pca <- read.table("pca.eigenvec",sep=" ",header=F)

      # plot pca
      b <- ggplot(pca, aes(V3, V4)) + geom_point(size = 3)
      b <- b + scale_colour_manual(values = c("red", "blue"))
      b <- b + coord_equal() + theme_light()
      b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
    }else if(paste0(file,".bed") %in% list.files(pattern = "*.bed")){
      system(glue::glue("extdata/plink/plink.exe --bfile {file} --hwe {hwe} --chr-set {chrset} --maf {maf} --geno {gen} --mind {mind} --allow-extra-chr --make-bed --out pca"))
      system(glue::glue("extdata/plink/plink.exe --bfile qc_filter --pca --chr-set {chrset} --allow-extra-chr --out pca" ))
      pca <- read.table("pca.eigenvec",sep=" ",header=F)

      # plot pca
      b <- ggplot(pca, aes(V3, V4)) + geom_point(size = 3)
      b <- b + scale_colour_manual(values = c("red", "blue"))
      b <- b + coord_equal() + theme_light()
      b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
    }
  }
}
# plQC(file = "All.GT", maf = 0.01, chrset = 30)
# library(tidyverse)
