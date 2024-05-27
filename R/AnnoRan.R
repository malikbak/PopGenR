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

AnnoSin <- function(data=data,annotation=annotation){
  data_gene <- data %>%
    mutate(external_gene_name = map2_chr(POSITION, CHR, function(x, y) {
      inds = x >= annotation$start_position & x <= annotation$end_position & y == annotation$chromosome_name
      if (any(inds)) annotation$external_gene_name[which.max(inds)] else NA
    }))
  return(data_gene)
}
