#' Annotate Genomic Regions with Gene Names
#'
#' Match genomic intervals, such as runs of homozygosity, to gene annotations.
#'
#' @param data A data frame containing the regions to annotate.
#' @param CHR Column name or index in `data` containing chromosome identifiers.
#' @param starP Column name or index in `data` containing region start positions.
#' @param endP Column name or index in `data` containing region end positions.
#' @param annotation Either a data frame of gene annotations or a path to a
#'   GFF/GTF file. Data frames should contain `chromosome_name`,
#'   `start_position`, `end_position`, and `external_gene_name`.
#'
#' @return A data frame with `ID`, `Range`, and `Matching_Names`. `Matching_Names`
#'   contains comma-separated gene names fully contained within each region, or
#'   `"NA"` when no annotation overlaps.
#'
#' @examples
#' regions <- data.frame(chr = c("1", "1"), start = c(1000, 5000), end = c(2000, 6500))
#' genes <- data.frame(
#'   chromosome_name = c("1", "1"),
#'   start_position = c(1200, 5500),
#'   end_position = c(1800, 6000),
#'   external_gene_name = c("GeneA", "GeneB")
#' )
#' AnnoRan(regions, CHR = "chr", starP = "start", endP = "end", annotation = genes)
#'
#' @export
AnnoRan <- function(data, CHR, starP, endP, annotation) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  region_data <- data.frame(
    ID = data[[CHR]],
    Start = as.numeric(data[[starP]]),
    End = as.numeric(data[[endP]])
  )
  region_data$Range <- paste(region_data$Start, region_data$End, sep = "-")

  if (is.data.frame(annotation)) {
    required <- c("chromosome_name", "start_position", "end_position", "external_gene_name")
    missing_cols <- setdiff(required, names(annotation))
    if (length(missing_cols) > 0) {
      stop("`annotation` is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }

    gene_data <- data.frame(
      ID = annotation$chromosome_name,
      Start = as.numeric(annotation$start_position),
      End = as.numeric(annotation$end_position),
      Name = annotation$external_gene_name
    )
  } else if (is.character(annotation) && length(annotation) == 1) {
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
      stop("Install the `rtracklayer` package to read GFF/GTF annotation files.", call. = FALSE)
    }

    gff <- as.data.frame(rtracklayer::readGFF(annotation))
    name_col <- if ("Name" %in% names(gff)) "Name" else if ("gene_name" %in% names(gff)) "gene_name" else NA_character_
    if (is.na(name_col)) {
      stop("The annotation file must contain a `Name` or `gene_name` column.", call. = FALSE)
    }

    gene_data <- data.frame(
      ID = gff$seqid,
      Start = as.numeric(gff$start),
      End = as.numeric(gff$end),
      Name = gff[[name_col]]
    )
  } else {
    stop("`annotation` must be a data frame or a single file path.", call. = FALSE)
  }

  matches <- vapply(seq_len(nrow(region_data)), function(i) {
    hit <- gene_data$ID == region_data$ID[i] &
      gene_data$Start >= region_data$Start[i] &
      gene_data$End <= region_data$End[i]

    if (any(hit)) {
      paste(unique(gene_data$Name[hit]), collapse = ",")
    } else {
      "NA"
    }
  }, character(1))

  data.frame(
    ID = region_data$ID,
    Range = region_data$Range,
    Matching_Names = matches,
    stringsAsFactors = FALSE
  )
}

#' Annotate Single Genomic Positions with Gene Names
#'
#' Match SNP or marker positions to gene annotations.
#'
#' @param data A data frame with `CHR` and `POSITION` columns.
#' @param annotation A data frame with `chromosome_name`, `start_position`,
#'   `end_position`, and `external_gene_name` columns.
#'
#' @return The input data with an added `external_gene_name` column.
#'
#' @examples
#' positions <- data.frame(CHR = c("1", "2"), POSITION = c(1500, 6000))
#' genes <- data.frame(
#'   chromosome_name = c("1", "2"),
#'   start_position = c(1000, 5000),
#'   end_position = c(2000, 6500),
#'   external_gene_name = c("GeneA", "GeneB")
#' )
#' AnnoSin(positions, genes)
#'
#' @export
AnnoSin <- function(data, annotation) {
  required_data <- c("CHR", "POSITION")
  required_annotation <- c("chromosome_name", "start_position", "end_position", "external_gene_name")

  if (!all(required_data %in% names(data))) {
    stop("`data` must contain CHR and POSITION columns.", call. = FALSE)
  }

  missing_cols <- setdiff(required_annotation, names(annotation))
  if (length(missing_cols) > 0) {
    stop("`annotation` is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  data$external_gene_name <- vapply(seq_len(nrow(data)), function(i) {
    hit <- data$POSITION[i] >= annotation$start_position &
      data$POSITION[i] <= annotation$end_position &
      data$CHR[i] == annotation$chromosome_name

    if (any(hit)) {
      annotation$external_gene_name[which(hit)[1]]
    } else {
      NA_character_
    }
  }, character(1))

  data
}
