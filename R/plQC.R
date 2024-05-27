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
