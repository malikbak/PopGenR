Mbed <- function(Ffile=NULL, Window=NULL){
  file <- c(glue::glue("bedtools makewindows -g {Ffile} -w {Window} >test.bed"))
  writeLines(file, "command.sh", sep = "\n")
  # Your R code goes here
  suppressWarnings(shell("./command.sh", shell = "bash", ignore.stdout = F))
}
#Mbed(Ffile = "GCF_019923935.1_NDDB_SH_1_genomic.fna.fai", Window = 40000)
