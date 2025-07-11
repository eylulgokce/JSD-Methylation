library(testthat)
library(dplyr)
library(tidyr)
library(tibble)

# Test Constants
test_that("Constants are correctly defined", {
  expect_equal(LOG2E, log2(exp(1)))
  expect_true(is.numeric(LOG2E))
  expect_true(LOG2E > 0)
})

# Test Shannon Entropy Function
test_that("shannon_entropy works correctly", {
  # Test basic functionality
  test_matrix <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)

  # Test row-wise entropy (axis = 1)
  entropy_rows <- shannon_entropy(test_matrix, axis = 1)
  expect_length(entropy_rows, 2)
  expect_true(all(entropy_rows >= 0))

  # Test column-wise entropy (axis = 2)
  entropy_cols <- shannon_entropy(test_matrix, axis = 2)
  expect_length(entropy_cols, 3)
  expect_true(all(entropy_cols >= 0))

  # Test error cases
  expect_error(shannon_entropy(test_matrix, method = "invalid"))
  expect_error(shannon_entropy(NULL))

  # Test with zeros
  zero_matrix <- matrix(0, nrow = 2, ncol = 3)
  entropy_zero <- shannon_entropy(zero_matrix, axis = 1)
  expect_equal(entropy_zero, c(0, 0))

  # Test with NA values
  na_matrix <- matrix(c(1, 2, NA, 4, 5, 6), nrow = 2, ncol = 3)
  entropy_na <- shannon_entropy(na_matrix, axis = 1)
  expect_length(entropy_na, 2)
  expect_true(all(is.finite(entropy_na)))

  # Test 1D input
  vector_input <- c(1, 2, 3, 4, 5)
  entropy_1d <- shannon_entropy(vector_input)
  expect_equal(entropy_1d, rep(0, 5))

  # Test uniform distribution (maximum entropy)
  uniform_matrix <- matrix(rep(1, 6), nrow = 2, ncol = 3)
  entropy_uniform <- shannon_entropy(uniform_matrix, axis = 1)
  expect_true(all(entropy_uniform > 0))

  # Test single row/column edge cases
  single_row <- matrix(c(1, 2, 3), nrow = 1, ncol = 3)
  entropy_single <- shannon_entropy(single_row, axis = 1)
  expect_length(entropy_single, 1)
  expect_true(entropy_single >= 0)
})

# Test Jensen-Shannon Divergence Function
test_that("js_divergence works correctly", {
  # Create test data with proper structure
  test_data <- data.frame(
    chrom = rep("chr1", 5),
    start = c(100, 200, 300, 400, 500),
    end = c(150, 250, 350, 450, 550),
    sample1.mC = c(5, 10, 15, 20, 25),
    sample1.C = c(10, 20, 30, 40, 50),
    sample2.mC = c(3, 8, 12, 18, 22),
    sample2.C = c(12, 18, 28, 38, 48)
  )

  # Test basic functionality
  result <- js_divergence(test_data)
  expect_true(is.data.frame(result))
  expected_cols <- c("#chrom", "start", "end", "JSD_bit_", "sample size",
                     "HMIX_bit_", "mC", "C", "MET")
  expect_true(all(expected_cols %in% colnames(result)))

  # Test JSD values are non-negative
  if (nrow(result) > 0) {
    expect_true(all(result$JSD_bit_ >= 0))
    expect_true(all(result$HMIX_bit_ >= 0))
    expect_true(all(result$MET >= 0 & result$MET <= 1))
  }

  # Test with compute_MET = FALSE
  result_no_met <- js_divergence(test_data, compute_MET = FALSE)
  expect_false("MET" %in% colnames(result_no_met))

  # Test with insufficient features (missing C column)
  insufficient_data <- data.frame(
    chrom = rep("chr1", 3),
    start = c(100, 200, 300),
    end = c(150, 250, 350),
    sample1.mC = c(5, 10, 15)
  )
  result_insufficient <- js_divergence(insufficient_data)
  expect_equal(nrow(result_insufficient), 0)

  # Test with debug = TRUE
  expect_output(js_divergence(test_data, debug = TRUE))

  # Test with low coverage data (should be filtered out)
  low_coverage_data <- data.frame(
    chrom = rep("chr1", 3),
    start = c(100, 200, 300),
    end = c(150, 250, 350),
    sample1.mC = c(1, 1, 1),
    sample1.C = c(1, 1, 1)
  )
  result_low <- js_divergence(low_coverage_data)
  expect_equal(nrow(result_low), 0)

  # Test identical samples (should produce zero divergence)
  identical_data <- data.frame(
    chrom = rep("chr1", 3),
    start = c(100, 200, 300),
    end = c(150, 250, 350),
    sample1.mC = c(10, 20, 30),
    sample1.C = c(20, 40, 60),
    sample2.mC = c(10, 20, 30),
    sample2.C = c(20, 40, 60)
  )
  result_identical <- js_divergence(identical_data)
  if (nrow(result_identical) > 0) {
    expect_true(all(result_identical$JSD_bit_ == 0))
  }
})

# Test Parameter Validation Functions
test_that("parameter validation works", {
  # Test divergence function parameter validation
  sample_data <- data.frame(
    url = c("file1.bgz", "file2.bgz"),
    label = c("sample1", "sample2")
  )

  expect_error(divergence(sample_data, chrom = NULL))
  expect_error(divergence(sample_data, chrom = "chr1", outfile = NULL))
  expect_error(divergence(sample_data, chrom = "chr1", outfile = "test.txt",
                          data_columns = NULL))

  # Test invalid sample data frame
  invalid_sample <- data.frame(file = c("file1.bgz", "file2.bgz"))
  expect_error(divergence(invalid_sample, chrom = "chr1", outfile = "test.txt",
                          data_columns = list(c(4, 5))))
})

# Test Preprocessing Functions
test_that("preprocessing functions work correctly", {
  # Test impute function
  test_data <- data.frame(a = c(1, 2, NA), b = c(NA, 4, 5))
  result <- impute(test_data, method = "pseudocount")
  expect_equal(result, 1)

  result_default <- impute(test_data)
  expect_equal(result_default, 1)

  # Test groupname function
  result <- groupname(by = c("group", "condition"),
                      name = c("A", "control"),
                      fname = "output.txt")
  expect_true(is.character(result))
  expect_true(grepl("group_A_and_condition_control", result))

  # Test with NULL values
  result_null <- groupname(by = NULL, name = NULL, fname = "output.txt")
  expect_equal(result_null, "output.txt")

  # Test with no file extension
  result_no_ext <- groupname(by = "group", name = "A", fname = "output")
  expect_equal(result_no_ext, "output_group_A")
})

# Test Edge Cases
test_that("edge cases are handled correctly", {
  # Test shannon_entropy with empty matrix
  empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_equal(length(shannon_entropy(empty_matrix, axis = 1)), 0)

  # Test js_divergence with empty data
  empty_data <- data.frame()
  result_empty <- js_divergence(empty_data)
  expect_equal(nrow(result_empty), 0)

  # Test js_divergence with only metadata columns
  minimal_data <- data.frame(
    chrom = "chr1",
    start = 100,
    end = 150
  )
  result_minimal <- js_divergence(minimal_data)
  expect_equal(nrow(result_minimal), 0)

  # Test with single sample
  single_sample_data <- data.frame(
    chrom = rep("chr1", 3),
    start = c(100, 200, 300),
    end = c(150, 250, 350),
    sample1.mC = c(10, 20, 30),
    sample1.C = c(20, 40, 60)
  )
  result_single <- js_divergence(single_sample_data)
  expect_equal(nrow(result_single), 0)  # Should be filtered out due to min_samples = 2
})

# Test Data Type Handling
test_that("data type handling works correctly", {
  # Test shannon_entropy with character matrix (should convert to numeric)
  char_matrix <- matrix(c("1", "2", "3", "4"), nrow = 2, ncol = 2)
  result_char <- shannon_entropy(char_matrix, axis = 1)
  expect_true(is.numeric(result_char))
  expect_length(result_char, 2)

  # Test js_divergence with character counts (should handle conversion)
  char_data <- data.frame(
    chrom = rep("chr1", 3),
    start = c(100, 200, 300),
    end = c(150, 250, 350),
    sample1.mC = c("5", "10", "15"),
    sample1.C = c("10", "20", "30"),
    sample2.mC = c("3", "8", "12"),
    sample2.C = c("12", "18", "28"),
    stringsAsFactors = FALSE
  )
  result_char_js <- js_divergence(char_data)
  expect_true(is.data.frame(result_char_js))
})

# Test Mathematical Properties
test_that("mathematical properties are correct", {
  # Test that Shannon entropy is maximized for uniform distribution
  uniform_matrix <- matrix(rep(1, 20), nrow = 4, ncol = 5)
  entropy_uniform <- shannon_entropy(uniform_matrix, axis = 1)

  # Test that skewed distribution has lower entropy
  skewed_matrix <- matrix(c(10, 1, 1, 1, 1), nrow = 1, ncol = 5)
  entropy_skewed <- shannon_entropy(skewed_matrix, axis = 1)

  expect_true(entropy_uniform[1] > entropy_skewed[1])

  # Test that JSD is symmetric (same result regardless of sample order)
  test_data_order1 <- data.frame(
    chrom = rep("chr1", 3),
    start = c(100, 200, 300),
    end = c(150, 250, 350),
    sampleA.mC = c(10, 20, 30),
    sampleA.C = c(20, 40, 60),
    sampleB.mC = c(5, 15, 25),
    sampleB.C = c(15, 35, 55)
  )

  test_data_order2 <- data.frame(
    chrom = rep("chr1", 3),
    start = c(100, 200, 300),
    end = c(150, 250, 350),
    sampleB.mC = c(5, 15, 25),
    sampleB.C = c(15, 35, 55),
    sampleA.mC = c(10, 20, 30),
    sampleA.C = c(20, 40, 60)
  )

  result1 <- js_divergence(test_data_order1)
  result2 <- js_divergence(test_data_order2)

  if (nrow(result1) > 0 && nrow(result2) > 0) {
    expect_equal(result1$JSD_bit_, result2$JSD_bit_)
  }
})

# Performance test with realistic data
test_that("functions handle realistic datasets efficiently", {
  # Create a moderately sized dataset
  n_sites <- 50
  large_data <- data.frame(
    chrom = rep("chr1", n_sites),
    start = seq(1000, by = 100, length.out = n_sites),
    end = seq(1100, by = 100, length.out = n_sites),
    sample1.mC = sample(1:30, n_sites, replace = TRUE),
    sample1.C = sample(20:80, n_sites, replace = TRUE),
    sample2.mC = sample(1:30, n_sites, replace = TRUE),
    sample2.C = sample(20:80, n_sites, replace = TRUE),
    sample3.mC = sample(1:30, n_sites, replace = TRUE),
    sample3.C = sample(20:80, n_sites, replace = TRUE)
  )

  # Test that js_divergence completes in reasonable time
  start_time <- Sys.time()
  result_large <- js_divergence(large_data)
  end_time <- Sys.time()

  expect_true(difftime(end_time, start_time, units = "secs") < 5)  # Should complete in under 5 seconds
  expect_true(is.data.frame(result_large))

  if (nrow(result_large) > 0) {
    expect_true(all(result_large$JSD_bit_ >= 0))
    expect_true(all(result_large$HMIX_bit_ >= 0))
    expect_true(all(result_large$MET >= 0 & result_large$MET <= 1))
  }
})

cat("Tests completed successfully!\n")
