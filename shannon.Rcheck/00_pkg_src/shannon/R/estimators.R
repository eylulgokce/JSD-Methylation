#' Shannon Entropy Estimation
#'
#' Computes Shannon entropy (in nat) of the feature frequency profile using the plug-in estimator.
#'
#' @param countmatrix A matrix of counts
#' @param axis Integer indicating which axis to compute entropy along (1 for rows, 2 for columns)
#' @param method Character string indicating the estimation method
#' @return A vector of Shannon entropy values
#' @export

library(dplyr)
library(tidyr)
library(tibble)

shannon_entropy <- function(countmatrix, axis = 1, method = "plug-in") {
  if (method != "plug-in") stop("Only 'plug-in' method is currently supported")
  if (is.null(countmatrix)) stop("Count matrix cannot be NULL")

  # Convert to matrix and handle NA
  countmatrix <- as.matrix(countmatrix)
  countmatrix[is.na(countmatrix)] <- 0

  # Handle 1D or invalid input
  if (is.vector(countmatrix) || length(dim(countmatrix)) < 2 || min(dim(countmatrix)) < 2) {
    n <- if (is.vector(countmatrix)) length(countmatrix) else if (axis == 1) nrow(countmatrix) else ncol(countmatrix)
    if (n == 0) return(numeric(0))
    return(rep(0, n)) # Return zeros for invalid cases
  }

  # Ensure numeric
  countmatrix <- apply(countmatrix, 2, as.numeric)

  # Compute count distribution
  count_distribution <- if (axis == 1) rowSums(countmatrix) else colSums(countmatrix)
  n <- length(count_distribution)
  entropy <- rep(0, n)

  # Avoid division by zero
  count_distribution[count_distribution == 0] <- 1

  # Normalize to probabilities
  prob <- if (axis == 1) countmatrix / count_distribution else t(t(countmatrix) / count_distribution)
  prob[is.nan(prob) | is.infinite(prob)] <- 0

  # Compute entropy in nats
  entropy_terms <- ifelse(prob > 0, -prob * log(prob), 0)
  entropy_vals <- if (axis == 1) rowSums(entropy_terms) else colSums(entropy_terms)
  entropy <- entropy_vals
  entropy[is.nan(entropy) | is.infinite(entropy)] <- 0

  entropy
}


#' Jensen-Shannon Divergence
#'
#' Computes Jensen-Shannon divergence for genomic data with quality control filters.
#'
#' @param indata A data frame with hierarchical column structure (sampling_unit, feature)
#' @param debug Logical, if TRUE prints debug information
#' @param compute_MET Logical, if TRUE computes methylation level (MET)
#' @return A data frame with JSD values and associated statistics
#' @export


js_divergence <- function(indata, debug = FALSE, compute_MET = TRUE) {
  meta_cols <- c("chrom", "start", "end")
  feature_cols <- setdiff(colnames(indata), meta_cols)

  # Parse sample-feature names (split on last period)
  split_names <- strsplit(feature_cols, "\\.")
  sample_ids <- sapply(split_names, function(x) paste(head(x, -1), collapse = "."))
  features <- sapply(split_names, tail, 1)
  samples <- unique(sample_ids)
  feats <- unique(features)

  if (debug) {
    cat("Feature columns:", feature_cols, "\n")
    cat("Samples:", samples, "\n")
    cat("Features:", feats, "\n")
  }

  # Check for valid features
  if (length(feats) < 2 || !all(c("mC", "C") %in% feats)) {
    if (debug) cat("Insufficient or invalid features (need mC and C), returning empty data frame\n")
    return(data.frame())
  }

  # Reshape long to array: [positions × samples × features]
  df_long <- indata %>%
    pivot_longer(cols = all_of(feature_cols),
                 names_to = c("sample", "feature"),
                 names_sep = "\\.(?=[^\\.]+$)", # Split on last period
                 values_to = "count") %>%
    mutate(pos_id = paste(chrom, start, end, sep = "_"))

  pos_ids <- unique(df_long$pos_id)
  if (length(pos_ids) == 0) {
    if (debug) cat("No valid positions, returning empty data frame\n")
    return(data.frame())
  }

  count_array <- array(0,
                       dim = c(length(pos_ids), length(samples), length(feats)),
                       dimnames = list(pos_id = pos_ids, sample = samples, feature = feats))

  for (i in seq_len(nrow(df_long))) {
    count_array[df_long$pos_id[i], df_long$sample[i], df_long$feature[i]] <- df_long$count[i]
  }

  # Get mC + C per sample
  counts_per_sample <- apply(count_array, c(1, 2), sum, na.rm = TRUE)

  # QC filters
  min_count <- 3
  min_samples <- 2
  sample_size <- rowSums(counts_per_sample > 0)
  count_filter <- rowSums(counts_per_sample >= min_count) > 0
  size_filter <- sample_size >= min_samples
  keep <- which(count_filter & size_filter)

  if (debug) {
    cat("Sites passing filters:", length(keep), "\n")
    cat("Sample size summary:", summary(sample_size), "\n")
  }

  if (length(keep) == 0) return(data.frame())

  pos_valid <- pos_ids[keep]

  # Mixture entropy
  mix_counts <- apply(count_array[pos_valid, , , drop = FALSE], c(1, 3), sum, na.rm = TRUE)
  if (is.vector(mix_counts) || ncol(mix_counts) < 2) {
    if (debug) cat("Mix counts has insufficient features or is invalid, returning empty data frame\n")
    return(data.frame())
  }
  mix_entropy <- shannon_entropy(mix_counts, axis = 1)

  # Sample-wise entropies
  entropy_mat <- matrix(0, nrow = length(pos_valid), ncol = length(samples))
  colnames(entropy_mat) <- samples
  for (i in seq_along(samples)) {
    # Extract sample counts for current sample across all positions and features
    sample_counts <- count_array[pos_valid, i, , drop = FALSE]

    # Convert to 2D matrix: positions × features
    sample_counts_2d <- matrix(sample_counts, nrow = length(pos_valid), ncol = length(feats))
    colnames(sample_counts_2d) <- feats

    # Handle missing values
    sample_counts_2d[is.na(sample_counts_2d)] <- 0

    # Calculate entropy for this sample
    sample_entropy <- shannon_entropy(sample_counts_2d, axis = 1)

    # Ensure the entropy vector has the correct length
    if (length(sample_entropy) != length(pos_valid)) {
      if (debug) cat(sprintf("Warning: Sample %s entropy length mismatch. Expected %d, got %d\n",
                             samples[i], length(pos_valid), length(sample_entropy)))
      sample_entropy <- rep(0, length(pos_valid))
    }

    entropy_mat[, i] <- sample_entropy
  }

  # Sample-wise weights
  weights <- apply(count_array[pos_valid, , , drop = FALSE], c(1, 2), sum, na.rm = TRUE)
  total_coverage <- rowSums(weights, na.rm = TRUE)
  normalized_weights <- weights / total_coverage
  normalized_weights[is.na(normalized_weights) | is.infinite(normalized_weights)] <- 0

  if (debug) {
    cat("Mix counts head:", head(mix_counts), "\n")
    cat("Mix entropy summary:", summary(mix_entropy), "\n")
    cat("Sample entropy matrix head:", head(entropy_mat), "\n")
    cat("Weights summary:", summary(normalized_weights), "\n")
  }

  # Weighted average entropy ⟨H⟩
  avg_entropy <- rowSums(entropy_mat * normalized_weights, na.rm = TRUE)

  # Convert to bits
  LOG2E <- 1 / log(2)
  JSD_bit_ <- round(pmax(LOG2E * (mix_entropy - avg_entropy), 0), 3)
  HMIX_bit_ <- round(pmin(LOG2E * mix_entropy, 1), 3)
  samplesize <- sample_size[keep]

  # Compute mC, C, and MET
  mc_counts <- apply(count_array[pos_valid, , "mC", drop = FALSE], c(1), sum, na.rm = TRUE)
  c_counts <- apply(count_array[pos_valid, , "C", drop = FALSE], c(1), sum, na.rm = TRUE)
  total_counts <- mc_counts + c_counts
  met <- ifelse(total_counts > 0, mc_counts / total_counts, 0)
  met <- round(met, 3)

  # Output table
  out <- tibble::tibble(
    `#chrom` = sapply(strsplit(pos_valid, "_"), `[`, 1),
    start = as.integer(sapply(strsplit(pos_valid, "_"), `[`, 2)),
    end = as.integer(sapply(strsplit(pos_valid, "_"), `[`, 3)),
    JSD_bit_ = JSD_bit_,
    `sample size` = samplesize,
    HMIX_bit_ = HMIX_bit_,
    mC = mc_counts,
    C = c_counts
  )

  if (compute_MET) {
    out$MET <- met
  }

  if (debug) {
    cat("avg_entropy summary:", summary(avg_entropy), "\n")
    cat("JSD_bit_ summary:", summary(JSD_bit_), "\n")
    cat("HMIX_bit_ summary:", summary(HMIX_bit_), "\n")
    cat("mC summary:", summary(mc_counts), "\n")
    cat("C summary:", summary(c_counts), "\n")
    if (compute_MET) cat("MET summary:", summary(out$MET), "\n")
  }

  return(out)
}
