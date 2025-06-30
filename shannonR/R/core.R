#' Compute Within-Group Divergence
#'
#' Computes within-group divergence for population using Jensen-Shannon divergence.
#'
#' @param sample Data frame containing sample information with 'url' and 'label' columns
#' @param chrom Character string specifying the chromosome to analyze
#' @param data_columns List of column specifications for data extraction
#' @param outfile Character string specifying the output file path
#' @param chunksize Integer specifying the expected number of sites per chunk (default: 1000)
#' @return NULL (results are written to file)
#' @export
divergence <- function(sample, chrom = NULL, data_columns = NULL,
                       outfile = NULL, chunksize = 1000) {

  if (is.null(chrom)) {
    stop("Chromosome must be specified")
  }

  if (is.null(outfile)) {
    stop("Output file must be specified")
  }

  if (is.null(data_columns)) {
    stop("Data columns must be specified")
  }

  # Check if sample has required columns
  if (!all(c("url", "label") %in% colnames(sample))) {
    stop("Sample data frame must contain 'url' and 'label' columns")
  }

  # Get regions for processing
  regions_info <- get_regions(sample$url, chrom = chrom, exp_numsites = chunksize)

  if (is.null(regions_info$regions)) {
    message("No regions to process")
    return(NULL)
  }

  regions_pct <- regions_info$progress
  regions <- regions_info$regions

  # Get data for each region
  regions_data <- get_data(sample$url, labels = sample$label,
                           data_columns = data_columns, regions = regions)

  # Process each region
  for (i in seq_along(regions_data)) {
    progress <- regions_pct[i]
    data <- regions_data[[i]]

    if (nrow(data) == 0) {
      cat(sprintf("...%5.1f %% (skipped empty region)\n", progress))
      next
    }

    # Compute JS divergence
    div <- js_divergence(data)

    if (nrow(div) == 0) {
      cat(sprintf("...%5.1f %% (skipped low-quality region)\n", progress))
      next
    }

    # Determine if header should be written
    header <- FALSE
    if (!file.exists(outfile)) {
      header <- TRUE
    } else if (file.info(outfile)$size == 0) {
      header <- TRUE
    }

    # Write results to file
    write.table(div, file = outfile, sep = "\t",
                col.names = header, row.names = TRUE,
                append = !header, quote = FALSE)

    cat(sprintf("...%5.1f %%\n", progress))
  }

  return(NULL)
}
