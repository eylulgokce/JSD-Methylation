#' @title Preprocessing functions
#' @description Functions for preprocessing genomic data
#' @name preprocessing
NULL

#' Impute missing values
#'
#' @description Imputes missing values of a data frame using specified method
#' @param data Data frame with missing values
#' @param method Imputation method (default: 'pseudocount')
#' @return Imputed value or modified data frame
#' @export
#' @examples
#' \dontrun{
#' imputed_value <- impute(data, method = "pseudocount")
#' }
impute <- function(data, method = "pseudocount") {
  if (method == "pseudocount") {
    # Given by the parameters of a uninformative Dirichlet prior on the probabilities
    return(1)
  }
  # Add other imputation methods as needed
  return(1)
}

#' Generate group filename
#'
#' @description Return filename of the subgroup based on grouping criteria
#' @param by Grouping variables
#' @param name Group name
#' @param fname Original filename
#' @return Modified filename
#' @export
#' @examples
#' \dontrun{
#' new_filename <- groupname(by = c("tissue", "condition"),
#'                          name = c("leaf", "control"),
#'                          fname = "results.txt")
#' }
groupname <- function(by = NULL, name = NULL, fname = NULL) {
  if (is.null(by) || is.null(name) || is.null(fname)) {
    return(fname)
  }

  # Ensure that name is a character vector
  name <- as.character(name)

  # Create group identifier
  group_parts <- paste(by, name, sep = "_")
  group <- paste(group_parts, collapse = "_and_")

  # Get file extension
  parts <- strsplit(fname, "\\.")[[1]]
  if (length(parts) > 1) {
    old_suffix <- parts[length(parts)]
    base_name <- paste(parts[-length(parts)], collapse = ".")
    new_suffix <- paste(group, old_suffix, sep = ".")
    new_fname <- paste(base_name, new_suffix, sep = ".")
  } else {
    new_fname <- paste(fname, group, sep = "_")
  }

  return(new_fname)
}
