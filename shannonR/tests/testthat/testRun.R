# Test suite for the shannon R package
# Save this as tests/testthat/test-shannon.R

library(testthat)
library(shannon)  # Assuming your package is named 'shannon'

# Test data creation functions
create_test_data <- function() {
  # Create a small test dataset mimicking methylation count data
  # Make sure we have enough data points to pass quality filters
  test_data <- data.frame(
    sample1.CG = c(10, 5, 8, 15, 12),
    sample1.CHG = c(2, 8, 3, 1, 4),
    sample1.CHH = c(1, 6, 12, 4, 3),
    sample2.CG = c(12, 3, 9, 18, 11),
    sample2.CHG = c(1, 9, 4, 2, 5),
    sample2.CHH = c(3, 1, 10, 5, 2),
    sample3.CG = c(8, 7, 11, 12, 9),
    sample3.CHG = c(3, 6, 5, 2, 7),
    sample3.CHH = c(2, 2, 8, 6, 4)
  )
  rownames(test_data) <- paste0("chr1:", 1000 + 1:nrow(test_data) * 100, "-", 1000 + 1:nrow(test_data) * 100 + 99)
  return(test_data)
}

create_test_metadata <- function() {
  metadata <- data.frame(
    label = paste0("sample", 1:3),
    url = paste0("/path/to/sample", 1:3, "_CHG.bedGraph.gz"),
    treatment = c("control", "control", "treated"),
    stringsAsFactors = FALSE
  )
  return(metadata)
}

# Test Shannon entropy calculation
test_that("Shannon entropy calculation works correctly", {
  # Test with simple known case
  test_matrix <- matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2)
  entropy <- shannon_entropy(test_matrix, axis = 1)

  # For uniform distribution, entropy should be log(2)
  expected_entropy <- log(2)
  expect_equal(entropy, rep(expected_entropy, 2), tolerance = 1e-10)

  # Test edge case with zeros
  test_matrix_zeros <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  entropy_zeros <- shannon_entropy(test_matrix_zeros, axis = 1)
  expect_equal(entropy_zeros, c(0, 0), tolerance = 1e-10)
})

# Test JS divergence calculation
test_that("JS divergence calculation works", {
  test_data <- create_test_data()

  # Run JS divergence
  result <- js_divergence(test_data)

  # Check that result has expected structure
  expect_true(is.data.frame(result))
  expect_true("JSD_bit_" %in% colnames(result))
  expect_true("sample_size" %in% colnames(result))
  expect_true("HMIX_bit_" %in% colnames(result))

  # JSD should be non-negative
  expect_true(all(result$JSD_bit_ >= 0))

  # Sample size should be reasonable
  expect_true(all(result$sample_size > 0))
})

# Test with empty data
test_that("Functions handle empty data gracefully", {
  empty_data <- data.frame()
  result <- js_divergence(empty_data)
  expect_equal(nrow(result), 0)
})

# Test preprocessing functions
test_that("Preprocessing functions work", {
  # Test imputation
  imputed_val <- impute(data.frame(x = c(1, NA, 3)), method = "pseudocount")
  expect_equal(imputed_val, 1)

  # Test group naming
  new_name <- groupname(by = c("treatment", "tissue"),
                        name = c("control", "leaf"),
                        fname = "results.txt")
  expect_true(grepl("treatment_control_and_tissue_leaf", new_name))
})

# Integration test function (requires actual data files)
run_integration_test <- function(metadata_file = NULL, test_chrom = "chr1") {
  if (is.null(metadata_file)) {
    cat("Skipping integration test - no metadata file provided\n")
    return(invisible(NULL))
  }

  if (!file.exists(metadata_file)) {
    cat("Metadata file not found:", metadata_file, "\n")
    return(invisible(NULL))
  }

  # Read metadata
  sample_data <- read.csv(metadata_file, stringsAsFactors = FALSE)

  # Check required columns
  if (!all(c("url", "label") %in% colnames(sample_data))) {
    stop("Sample data must contain 'url' and 'label' columns")
  }

  # Test with a small subset of data
  sample_subset <- head(sample_data, 3)  # Use first 3 samples

  # Define data columns (assuming CHG methylation data)
  data_columns <- list(c(4, 5))  # Assuming count columns are 4 and 5

  # Create temporary output file
  temp_outfile <- tempfile(fileext = ".txt")

  tryCatch({
    # Run divergence calculation
    divergence(sample = sample_subset,
               chrom = test_chrom,
               data_columns = data_columns,
               outfile = temp_outfile,
               chunksize = 100)  # Small chunk for testing

    # Check if output file was created and has content
    if (file.exists(temp_outfile) && file.info(temp_outfile)$size > 0) {
      result <- read.table(temp_outfile, header = TRUE, sep = "\t")
      cat("Integration test successful! Results shape:", nrow(result), "x", ncol(result), "\n")
      cat("Column names:", paste(colnames(result), collapse = ", "), "\n")
      if (nrow(result) > 0) {
        cat("Sample JSD values:", head(result$JSD_bit_, 3), "\n")
      }
      return(result)
    } else {
      cat("No output generated - check if data files exist and are accessible\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("Integration test failed with error:", e$message, "\n")
    return(NULL)
  }, finally = {
    # Clean up
    if (file.exists(temp_outfile)) {
      unlink(temp_outfile)
    }
  })
}

# Manual test function for development
manual_test <- function() {
  cat("=== Manual Test Suite ===\n")

  # Test 1: Basic Shannon entropy
  cat("\n1. Testing Shannon entropy...\n")
  test_matrix <- matrix(c(2, 2, 1, 3), nrow = 2, ncol = 2)
  entropy_result <- shannon_entropy(test_matrix, axis = 1)
  cat("Entropy results:", entropy_result, "\n")

  # Test 2: JS divergence with synthetic data
  cat("\n2. Testing JS divergence...\n")
  test_data <- create_test_data()
  cat("Test data dimensions:", dim(test_data), "\n")
  cat("Test data preview:\n")
  print(head(test_data, 3))

  js_result <- js_divergence(test_data)
  cat("JS divergence results:\n")
  print(js_result)

  # Test 3: Constants
  cat("\n3. Testing constants...\n")
  cat("LOG2E value:", LOG2E, "\n")
  cat("Should be approximately:", log2(exp(1)), "\n")

  cat("\n=== Manual tests completed ===\n")
  return(invisible(list(entropy = entropy_result, js_div = js_result)))
}

# Usage examples:
# 1. Run unit tests: testthat::test_file("test-shannon.R")
# 2. Run manual tests: manual_test()
# 3. Run integration test: run_integration_test("path/to/metadata.csv", "chr1")
