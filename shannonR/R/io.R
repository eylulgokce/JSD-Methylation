#' @title Input/Output functions for genomic data
#' @description Functions to load and transform genomic data files
#' @name io
NULL

#' Read bedcount files
#'
#' @description bedcount_reader returns a data frame of the bedcount data
#' @param bedcount Path to bedcount file
#' @param compression Compression type (default: NULL)
#' @param chunksize Chunk size for reading (default: 10000)
#' @return Data frame with bedcount data
#' @export
#' @examples
#' \dontrun{
#' data <- bedcount_reader("example.bedcount")
#' }
bedcount_reader <- function(bedcount, compression = NULL, chunksize = 10000) {
  # Read the file with appropriate settings
  data <- read.table(bedcount,
                     header = TRUE,
                     sep = "\t",
                     stringsAsFactors = FALSE,
                     colClasses = c("#chrom" = "character", "start" = "integer"))

  return(data)
}

#' Filter population data based on metadata
#'
#' @description Read metadata and return the population and the quotient set
#' @param metadata Path to metadata file
#' @param subset Subset condition (default: NULL)
#' @param relation Relation for grouping (default: NULL)
#' @return List with reference and qset components
#' @export
#' @examples
#' \dontrun{
#' pop <- population_filter("metadata.txt", subset = "condition == 'treated'")
#' }
population_filter <- function(metadata, subset = NULL, relation = NULL) {
  pop <- list(reference = NULL, qset = NULL)
  meta <- read.table(metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  if (!is.null(subset)) {
    # Simple subset evaluation - in practice you might want more sophisticated parsing
    pop$reference <- meta$sample[eval(parse(text = subset), envir = meta)]
  } else {
    pop$reference <- meta$sample
  }

  if (!is.null(relation)) {
    reference_meta <- meta[meta$sample %in% pop$reference, ]
    grouped <- split(reference_meta$sample, reference_meta[[relation]])
    pop$qset <- as.list(grouped)
  }

  return(pop)
}
