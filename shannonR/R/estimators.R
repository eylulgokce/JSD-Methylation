#' Shannon Entropy Estimation
#'
#' Computes Shannon entropy (in nats) of the feature frequency profile using the plug-in estimator.
#'
#' @param countmatrix A matrix of counts
#' @param axis Integer indicating which axis to compute entropy along (1 for rows, 2 for columns)
#' @param method Character string indicating the estimation method (currently only "plug-in" is supported)
#' @return A vector of Shannon entropy values
#' @export


shannon_entropy <- function(countmatrix, axis = 1, method = "plug-in") {
  if (method != "plug-in") {
    stop("Only 'plug-in' method is currently supported")
  }

  countmatrix <- as.matrix(countmatrix)

  if (is.null(countmatrix) || length(dim(countmatrix)) != 2 || any(dim(countmatrix) == 0)) {
    return(rep(NA_real_, max(1, nrow(countmatrix))))
  }

  countmatrix <- apply(countmatrix, 2, as.numeric)

  if (axis == 1) {
    count_distribution <- rowSums(countmatrix, na.rm = TRUE)
    prob <- countmatrix / count_distribution
  } else {
    count_distribution <- colSums(countmatrix, na.rm = TRUE)
    prob <- t(t(countmatrix) / count_distribution)
  }

  prob[is.nan(prob) | is.infinite(prob)] <- 0
  entropy_terms <- ifelse(prob > 0, -prob * log(prob), 0)

  if (axis == 1) {
    entropy <- rowSums(entropy_terms, na.rm = TRUE)
  } else {
    entropy <- colSums(entropy_terms, na.rm = TRUE)
  }

  return(entropy)
}



#' Jensen-Shannon Divergence (Fixed Version)
#'
#' Computes Jensen-Shannon divergence for genomic data with quality control filters.
#'
#' @param indata A data frame with hierarchical column structure (sampling_unit, feature)
#' @param weights Optional weights for computing weighted averages
#' @return A data frame with JSD values and associated statistics
#' @export
js_divergence <- function(indata, weights = NULL) {

  if (nrow(indata) == 0) {
    return(indata)
  }

  # Extract sampling unit and feature information from column names
  col_names <- colnames(indata)

  # Parse hierarchical column structure
  if (any(grepl("\\.", col_names))) {
    # Split column names to extract sampling units and features
    split_names <- strsplit(col_names, "\\.")
    sampling_units <- sapply(split_names, function(x) x[1])
    features <- sapply(split_names, function(x) if(length(x) > 1) x[2] else x[1])
  } else {
    # If no hierarchical structure, treat each column as a separate unit
    sampling_units <- col_names
    features <- rep("feature", length(col_names))
  }

  # Create a mapping for sampling units and features
  unique_units <- unique(sampling_units)
  unique_features <- unique(features)

  # Aggregate data by sampling unit
  count_per_unit <- matrix(0, nrow = nrow(indata), ncol = length(unique_units))
  colnames(count_per_unit) <- unique_units
  rownames(count_per_unit) <- rownames(indata)

  for (i in seq_along(unique_units)) {
    unit_cols <- which(sampling_units == unique_units[i])
    if (length(unit_cols) > 1) {
      count_per_unit[, i] <- rowSums(indata[, unit_cols], na.rm = TRUE)
    } else if (length(unit_cols) == 1) {
      count_per_unit[, i] <- indata[, unit_cols]
    }
  }

  # Calculate sample size (number of non-null units per position)
  samplesize <- rowSums(!is.na(count_per_unit) & count_per_unit > 0)

  # Quality control filters
  min_samplesize <- 2
  min_count <- 3

  count_filter <- apply(count_per_unit >= min_count, 1, any, na.rm = TRUE)
  samplesize_filter <- (samplesize >= min_samplesize)
  combined_filter <- (count_filter & samplesize_filter)

  # Apply filters
  data <- indata[combined_filter, , drop = FALSE]

  if (nrow(data) == 0) {
    return(data)
  }

  data_unit <- count_per_unit[combined_filter, , drop = FALSE]

  # Aggregate by feature
  data_feature <- matrix(0, nrow = nrow(data), ncol = length(unique_features))
  colnames(data_feature) <- unique_features
  rownames(data_feature) <- rownames(data)

  for (i in seq_along(unique_features)) {
    feature_cols <- which(features == unique_features[i])
    if (length(feature_cols) > 1) {
      data_feature[, i] <- rowSums(data[, feature_cols], na.rm = TRUE)
    } else if (length(feature_cols) == 1) {
      data_feature[, i] <- data[, feature_cols]
    }
  }

  # Convert to integer
  data_feature <- matrix(as.integer(data_feature),
                         nrow = nrow(data_feature),
                         ncol = ncol(data_feature))
  colnames(data_feature) <- unique_features
  rownames(data_feature) <- rownames(data)

  # Calculate mixture entropy 
  # For each position, calculate entropy across features
  mix_entropy <- shannon_entropy(data_feature, axis = 1)  # Changed from axis = 2 to axis = 1

  # For average entropy, we need to compute entropy for each sampling unit
  # Initialize matrix to store unit entropies
  unit_entropies <- matrix(NA, nrow = nrow(data_feature), ncol = ncol(data_unit))

  for (i in 1:ncol(data_unit)) {
    # For each sampling unit, calculate the feature distribution
    unit_data <- matrix(0, nrow = nrow(data_feature), ncol = length(unique_features))
    colnames(unit_data) <- unique_features
    rownames(unit_data) <- rownames(data_feature)

    for (j in seq_along(unique_features)) {
      feature_cols <- which(features == unique_features[j] & sampling_units == unique_units[i])
      if (length(feature_cols) > 0) {
        if (length(feature_cols) > 1) {
          unit_data[, j] <- rowSums(data[, feature_cols], na.rm = TRUE)
        } else {
          unit_data[, j] <- data[, feature_cols]
        }
      }
    }

    # Calculate entropy for this unit across all positions
    # This should return one entropy value per position (row)
    unit_entropies[, i] <- shannon_entropy(unit_data, axis = 1)  
  }

  # Calculate weighted average entropy
  # Replace NA values in data_unit with 0 for weighting
  weights_matrix <- data_unit
  weights_matrix[is.na(weights_matrix)] <- 0

  # Calculate weighted average across units for each position
  # Ensure both matrices are numeric
  unit_entropies <- as.matrix(unit_entropies)
  unit_entropies <- apply(unit_entropies, 2, as.numeric)

  weights_matrix <- as.matrix(weights_matrix)
  weights_matrix <- apply(weights_matrix, 2, as.numeric)

  # Compute weighted average across units
  weighted_sum <- rowSums(unit_entropies * weights_matrix, na.rm = TRUE)
  weight_totals <- rowSums(weights_matrix, na.rm = TRUE)
  avg_entropy <- weighted_sum / weight_totals

  # Handle division by zero
  avg_entropy[weight_totals == 0 | is.na(weight_totals)] <- 0

  # Create result data frame
  div <- data.frame(
    JSD_bit_ = LOG2E * (mix_entropy - avg_entropy),
    sample_size = samplesize[combined_filter],
    HMIX_bit_ = LOG2E * mix_entropy,
    stringsAsFactors = FALSE
  )

  # Add feature columns
  div <- cbind(div, data_feature)

  # Round specific columns
  div$JSD_bit_ <- round(div$JSD_bit_, 3)
  div$HMIX_bit_ <- round(div$HMIX_bit_, 3)

  # Set row names
  rownames(div) <- rownames(data)

  return(div)
}





