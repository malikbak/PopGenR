RunHom <- function(file=file, homozygWS = 50, homozygS = 50, homozygWM = 3,
                 homozygKB= 100 ,homozygDN = 1000 ,chrset = 22){
  if(Sys.info()['sysname'] == "Windows"){
    if(paste0(file,".ped") %in% list.files(pattern = "*.ped")){
      system(glue::glue("E:/thesis/package-pro/plink.exe --file {file} --homozyg-group --homozyg-window-snp {homozygWS} --homozyg-snp {homozygWS} --homozyg-window-missing {homozygWM} --homozyg-kb {homozygKB} --homozyg-density {homozygDN} --chr-set {chrset} --allow-extra-chr --out roh"))
    }else if(paste0(file,".bed") %in% list.files(pattern = "*.bed")){
      system(glue::glue("E:/thesis/package-pro/plink.exe --bfile {file} --homozyg-group --homozyg-window-snp {homozygWS} --homozyg-snp {homozygWS} --homozyg-window-missing {homozygWM} --homozyg-kb {homozygKB} --homozyg-density {homozygDN} --chr-set {chrset} --allow-extra-chr --out roh"))
    }
  }
  if("roh.hom" %in% list.files()){
    ROH <- read.table("roh.hom", header = T)
    return(ROH)
  }else(cat("File not found"))
}
# ROH <- RunHom(file = "All.GT", chrset = 30)
###################### Tajima's D ###################################

TajimaD <- function(vfile=file, window_size=150000, nsample=1){
  vcf <- read.vcfR(vfile)
  snp <- as.data.frame(vcf@fix)
  my.data <- cbind(snp, vcf@gt)
  my.data$POS <- as.integer(my.data$POS)
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
  bed$SNPs <- 0
  for (i in 1:nrow(bed)) { # Notice here that the number of iterations in my loop is equal to the number of rows in my RESULTS, not my input data
    d=subset(my.data, (my.data$CHROM == bed$chr[i] & my.data$POS>=bed$start[i] & my.data$POS<= bed$end[i])) # For each of my windows, I get the subset of data that falls into that window
    bed$SNPs[i]=nrow(d) # Then I count the number of SNPs in that subset and save it in my results table
    message(glue::glue("SNPs for {i} is {nrow(d)} "))
  }
  ### Fill in the Number of Chromosomes for each Window
  if(length(nsample) > 1){
    num.samples = nsample
  }else{
    num.samples=ncol(my.data)-4
  }
  bed$Nchr=2*(num.samples) # Because there is no missing data in this particular file, I can assume the value for 2N is the same for every window

  ### Calculate Waterson's Theta for Each window
  a=seq(from=1, to=((2*num.samples)-1), by=1) # Like our Theta calculations last week, I need to get the sequence of numbers from 1 to 2N-1,
  bed$ThetaW = bed$SNPs/(sum(1/a)) # then theta-w = #SNPs/(1/1 + 1/2 + 1/3 + ... + 1/(2N-1))
  ### Create a function for derived counts
  derivedCount <- function(row) {
    row=as.vector(row, mode="character")
    x=count.genotypes(get.field(row[5:length(row)], row[4], "GT"))
    dc=(2*x["aa"])+x["Aa"]
    return(unname(dc))
  }

  ### Now, calculate Pi for each window
  bed$Pi=rep(0, nrow(bed)) # Here, I add a column of zeros to my results table so I can fill in the Pi calculations
  for (i in 1:100) { # Again, I want a calculation for every WINDOW, not every SNP in the VCF file, so my loop goes through 1:50 (there are 50 windows)
    d=subset(my.data, (my.data$CHROM == bed$chr[i] & my.data$POS>=bed$start[i] & my.data$POS<bed$end[i])) # get the subset in the window
    if(nrow(d) > 1){
      j=apply(d, 1, FUN=derivedCount) # The apply function gets the derived count at EVERY row in my subset of data; so the "j" variable is a VECTOR of counts
      # this is the same as doing a second loop inside of the first loop; it is just a little more efficient
      c=rep(bed$Nchr[i], length(j)) # I don't actually have to do this, but here I am make a vector of my 2N values (same as c in Pi equation) that is the same length as my vector of j values
      bed$Pi[i]=sum((2*j*(c-j))/(c*(c-1))) # Finally, I use the equation for Pi to get the calculation for the whole window
      message("Genotype present...............................")
    }else(message("No Genotype"))
  }
  ### Calculate Var(d) for each window
  bed$varD=rep(0, nrow(bed)) # Set up an empty column to hold the results of the variance function

  ### Create a function to calculate var(d)
  ### n = the number of chromosomes (2 x num. samples)
  ### S = # SNPs
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

  for (i in 1:nrow(bed)) { # Loop through the windows one more time, this time calculate the variance for each window
    bed$varD[i] = variance.d(n=bed$Nchr[i], S=bed$SNPs[i])
  }

  ### Now, calculate Tajima's D
  bed$TajimaD = (bed$Pi - bed$ThetaW)/(sqrt(bed$varD)) # With all of the calculations set up in my table, I can do the last calculation for Tajima's D using R's vector math capabilities (and trust that it understands I want one calculation per row in the table)
  return(my.results)
}

