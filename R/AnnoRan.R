#' AnnoRan: Annotate Regions of Homozygosity (ROH) Based on Genomic Annotations
#'
#' This function annotates genomic ranges within a dataset by matching them with corresponding gene annotations. 
#' It is designed to work with a data frame containing regions of homozygosity (ROH) and a genomic annotation 
#' file or data frame.
#'
#' @param data A data frame containing the regions to be annotated. This should include columns for chromosome 
#' identifier, start position, and end position.
#' @param CHR A character string or numeric value specifying the column in `data` representing the chromosome.
#' @param starP A character string or numeric value specifying the column in `data` representing the start 
#' position of the regions.
#' @param endP A character string or numeric value specifying the column in `data` representing the end position 
#' of the regions.
#' @param annotation A file path to a GFF/GTF file or a data frame containing gene annotations. If a character 
#' path is provided, the function will read the GFF file using `rtracklayer::readGFF`. The data frame or file 
#' should include columns for chromosome name, start position, end position, and gene name.
#'
#' @return A data frame that includes the original regions along with an additional column, `Matching_Names`, 
#' which contains the gene names corresponding to the genomic regions.
#'
#' @details The function works by first creating a new column, `region`, in the input `data` that combines the 
#' start and end positions into a string format. It then matches these ranges with the gene annotations provided 
#' via the `annotation` argument. If a match is found, the gene names are retrieved and associated with the 
#' respective regions. If no match is found, "NA" is assigned.
#'
#' @note The function supports annotations provided either as a file in GFF/GTF format or as a data frame. It 
#' expects the annotation data frame to have specific column names, such as `chromosome_name`, `start_position`, 
#' `end_position`, and `external_gene_name` (or equivalent).
#'
#' @importFrom rtracklayer readGFF
#' @examples
#' # Example usage with data frame and GFF file:
#' ROH_data <- data.frame(chr = c("1", "1"), start = c(1000, 5000), end = c(2000, 6000))
#' annotation_file <- "path/to/annotation.gff"
#' annotated_data <- AnnoRan(data = ROH_data, CHR = "chr", starP = "start", endP = "end", annotation = annotation_file)
#'
#' # Example usage with data frames:
#' annotation_df <- data.frame(chromosome_name = c("1", "1"), start_position = c(1500, 5500), 
#'                             end_position = c(1800, 5800), external_gene_name = c("Gene1", "Gene2"))
#' annotated_data <- AnnoRan(data = ROH_data, CHR = "chr", starP = "start", endP = "end", annotation = annotation_df)
#'
#' @export
AnnoRan <- function(data = data, CHR = chrcol, starP = starP, endP = endP, annotation = Anofile){
  ROH_test <- data
  ROH_test$region <- paste(ROH_test[,starP], ROH_test[,endP], sep = "-")
  # Sample data frame with ranges
  df_ranges <- data.frame(ID = ROH_test[,CHR],
                          Range = ROH_test$region)
  
  # Sample data frame with start, end, and name columns
  if(class(annotation) == "data.frame"){
    if("chromosome_name" %in% colnames(annotation)){
      df_data <- data.frame(ID = annotation[,"chromosome_name"],
                            Start = annotation[,"start_position"],
                            End = annotation[,"end_position"],
                            Name = annotation[,"external_gene_name"])
    }
  }else if(class(annotation) == "character"){
    annotation <- rtracklayer::readGFF(annotation)
    annotation <- as.data.frame(annotation)
    df_data <- data.frame(ID = annotation[,"seqid"],
                          Start = annotation[,"start"],
                          End = annotation[,"end"],
                          Name = annotation[,"Name"])
  }
  
  # Function to match ranges and retrieve corresponding names
  match_ranges <- function(id, range_str, start_values, end_values, names, data_ids) {
    matched_names <- character(length(id))
    
    for (i in seq_along(id)) {
      # Find the index of the matching ID
      match_index <- which(id[i] == data_ids)
      
      # Filter the data frame based on the matching ID
      filtered_data <- df_data[df_data$ID == id[i], ]
      
      # Convert range string to start and end values
      range_parts <- strsplit(range_str[i], "-")
      range_start <- as.numeric(range_parts[[1]][1])
      range_end <- as.numeric(range_parts[[1]][2])
      
      # Find indices of matching ranges within the filtered data frame
      range_indices <- which(filtered_data$Start >= range_start & filtered_data$End <= range_end)
      
      if (length(range_indices) > 0) {
        # Retrieve corresponding names
        matched_names[i] <- paste(filtered_data$Name[range_indices], collapse = ",")
      } else {
        # Assign NA if no match found
        matched_names[i] <- NA_character_
      }
    }
    
    matched_names
  }
  
  # Get the IDs from df_data
  data_ids <- df_data$ID
  
  # Apply the function to each row in df_ranges
  df_ranges$Matching_Names <- match_ranges(df_ranges$ID, df_ranges$Range, df_data$Start, df_data$End, df_data$Name, data_ids)
  
  # Convert NA values to character "NA"
  df_ranges$Matching_Names <- ifelse(is.na(df_ranges$Matching_Names), "NA", df_ranges$Matching_Names)
  return(df_ranges)
}

#' AnnoSin: Single-Position Annotation Based on Genomic Annotations
#'
#' This function annotates individual genomic positions within a dataset by matching them with corresponding gene 
#' annotations. It takes a data frame with genomic positions and annotates them using an external annotation file 
#' or data frame.
#'
#' @param data A data frame containing the genomic positions to be annotated. This should include columns for 
#' chromosome (`CHR`) and position (`POSITION`).
#' @param annotation A data frame containing gene annotations. The annotation data frame must include columns 
#' for `chromosome_name`, `start_position`, `end_position`, and `external_gene_name`.
#'
#' @return A data frame that includes the original positions along with an additional column, `external_gene_name`, 
#' containing the gene name corresponding to each position.
#'
#' @details The function uses the `POSITION` and `CHR` columns from the input `data` to find the corresponding 
#' gene name in the `annotation` data frame. For each position, it checks whether the position lies within the 
#' start and end coordinates of a gene in the specified chromosome. If a match is found, the corresponding gene 
#' name is assigned; otherwise, `NA` is returned.
#'
#' @note The function expects the `annotation` data frame to contain specific columns: `chromosome_name`, 
#' `start_position`, `end_position`, and `external_gene_name`.
#'
#' @importFrom dplyr mutate
#' @importFrom purrr map2_chr
#'
#' @examples
#' # Example usage:
#' positions_data <- data.frame(CHR = c("1", "2"), POSITION = c(1500, 6000))
#' annotation_df <- data.frame(chromosome_name = c("1", "2"), start_position = c(1000, 5000), 
#'                             end_position = c(2000, 6500), external_gene_name = c("GeneA", "GeneB"))
#' annotated_data <- AnnoSin(data = positions_data, annotation = annotation_df)
#'
#' @export
AnnoSin <- function(data=data,annotation=annotation){
  data_gene <- data %>%
    mutate(external_gene_name = map2_chr(POSITION, CHR, function(x, y) {
      inds = x >= annotation$start_position & x <= annotation$end_position & y == annotation$chromosome_name
      if (any(inds)) annotation$external_gene_name[which.max(inds)] else NA
    }))
  return(data_gene)
}
