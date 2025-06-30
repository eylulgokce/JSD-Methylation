#' shannon: Information-Theoretic Analysis of Genomic Data
#'
#' The shannon package provides tools for computing Jensen-Shannon divergence
#' and other information-theoretic quantities for genomic position files (GPFs).
#' It is particularly useful for analyzing DNA methylation diversity and other
#' genomic features using Shannon entropy and related measures.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{divergence}}: Computes within-group divergence for population
#'   \item \code{\link{js_divergence}}: Calculates Jensen-Shannon divergence
#'   \item \code{\link{shannon_entropy}}: Computes Shannon entropy
#'   \item \code{\link{get_regions}}: Gets regions for analysis
#'   \item \code{\link{get_data}}: Retrieves and combines genomic data
#' }
#'
#' @section Data input:
#' The package works with:
#' \itemize{
#'   \item Tabix-indexed genome position files (GPFs)
#'   \item Metadata tables specifying sample information
#'   \item BED-like format files
#' }
#'
#' @keywords package
"_PACKAGE"

#' @importFrom data.table data.table setkey
#' @importFrom Rsamtools TabixFile scanTabix
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics start end
#' @importFrom rtracklayer import
#' @importFrom utils read.table write.table
#' @importFrom stats weighted.mean
NULL
