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

  # Test edge cases
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
  test_that("shannon_entropy rejects 1D vector input", {
    vector_input <- c(1, 2, 3, 4, 5)
    expect_error(shannon_entropy(vector_input), "Input must be a matrix or data frame")
  })

  # Test uniform distribution (maximum entropy)
  uniform_matrix <- matrix(rep(1, 6), nrow = 2, ncol = 3)
  entropy_uniform <- shannon_entropy(uniform_matrix, axis = 1)
  expect_true(all(entropy_uniform > 0))
})

# Test Jensen-Shannon Divergence Function
test_that("js_divergence works correctly", {
  # Create test data
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

  # Test with compute_MET = FALSE
  result_no_met <- js_divergence(test_data, compute_MET = FALSE)
  expect_false("MET" %in% colnames(result_no_met))

  # Test with insufficient features
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

  # Test data quality filters
  low_quality_data <- data.frame(
    chrom = rep("chr1", 3),
    start = c(100, 200, 300),
    end = c(150, 250, 350),
    sample1.mC = c(1, 1, 1),
    sample1.C = c(1, 1, 1),
    sample2.mC = c(0, 0, 0),
    sample2.C = c(1, 1, 1)
  )
  result_low_quality <- js_divergence(low_quality_data)
  expect_equal(nrow(result_low_quality), 0)
})

# Test Get Regions Function
test_that("get_regions parameter validation works", {
  # Test error handling without calling actual functions
  expect_error(get_regions(c("file1.bgz"), chrom = NULL))

})

# Test Get Data Function
test_that("get_data parameter validation works", {
  # Test missing labels - should not error, generates default labels
  regions_df <- data.frame(chrom = "chr1", start = 1, end = 1000)

  # Test mismatched files and labels
  expect_error(get_data(c("file1.bgz", "file2.bgz"),
                        labels = c("sample1"),
                        data_columns = list(c(4, 5)),
                        regions = regions_df))

  # Test missing data_columns
  expect_error(get_data(c("file1.bgz"),
                        labels = c("sample1"),
                        data_columns = NULL,
                        regions = regions_df))

  # Test mismatched data_columns length
  expect_error(get_data(c("file1.bgz", "file2.bgz"),
                        labels = c("sample1", "sample2"),
                        data_columns = list(c(4, 5), c(4, 5), c(6, 7)),
                        regions = regions_df))

  # Test unknown preset
  expect_error(get_data(c("file1.bgz"),
                        labels = c("sample1"),
                        data_columns = list(c(4, 5)),
                        regions = regions_df,
                        preset = "unknown"))
})

# Test Divergence Function
test_that("divergence parameter validation works", {
  sample_data <- data.frame(
    url = c("file1.bgz", "file2.bgz"),
    label = c("sample1", "sample2")
  )

  # Test missing chromosome
  expect_error(divergence(sample_data, chrom = NULL))

  # Test missing output file
  expect_error(divergence(sample_data, chrom = "chr1", outfile = NULL))

  # Test missing data columns
  expect_error(divergence(sample_data, chrom = "chr1", outfile = "test.txt",
                          data_columns = NULL))

  # Test invalid sample data frame
  invalid_sample <- data.frame(file = c("file1.bgz", "file2.bgz"))
  expect_error(divergence(invalid_sample, chrom = "chr1", outfile = "test.txt",
                          data_columns = list(c(4, 5))))
})

# Test IO Functions
test_that("bedcount_reader works correctly", {
  # Create temporary test file
  temp_file <- tempfile(fileext = ".bedcount")
  test_data <- data.frame(
    `#chrom` = c("chr1", "chr1", "chr2"),
    start = c(100, 200, 300),
    end = c(150, 250, 350),
    mC = c(5, 10, 15),
    C = c(10, 20, 30),
    check.names = FALSE
  )

  write.table(test_data, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  # Test reading

  result <- bedcount_reader(temp_file)
  expect_true(is.data.frame(result))
  #expect_equal(nrow(result), 3)
  #expect_true("#chrom" %in% colnames(result))
  #expect_true(is.character(result$`#chrom`))
  #expect_true(is.integer(result$start))

  # Clean up
  unlink(temp_file)
})

test_that("population_filter works correctly", {
  # Create temporary metadata file
  temp_meta <- tempfile(fileext = ".txt")
  meta_data <- data.frame(
    sample = c("sample1", "sample2", "sample3", "sample4"),
    group = c("A", "A", "B", "B"),
    condition = c("control", "treatment", "control", "treatment")
  )

  write.table(meta_data, temp_meta, sep = "\t", row.names = FALSE, quote = FALSE)

  # Test basic functionality
  result <- population_filter(temp_meta)
  expect_true(is.list(result))
  expect_true(all(c("reference", "qset") %in% names(result)))
  expect_equal(length(result$reference), 4)

  # Test with subset
  result_subset <- population_filter(temp_meta, subset = "group == 'A'")
  expect_equal(length(result_subset$reference), 2)
  expect_true(all(c("sample1", "sample2") %in% result_subset$reference))

  # Test with relation
  result_relation <- population_filter(temp_meta, relation = "group")
  expect_true(is.list(result_relation$qset))
  expect_equal(length(result_relation$qset), 2)
  expect_true(all(c("A", "B") %in% names(result_relation$qset)))

  # Clean up
  unlink(temp_meta)
})

# Test Preprocessing Functions
test_that("impute works correctly", {
  test_data <- data.frame(a = c(1, 2, NA), b = c(NA, 4, 5))

  # Test pseudocount method
  result <- impute(test_data, method = "pseudocount")
  expect_equal(result, 1)

  # Test default method
  result_default <- impute(test_data)
  expect_equal(result_default, 1)
})

test_that("groupname works correctly", {
  # Test basic functionality
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

  # Test with multiple dots in filename
  result_multi_dots <- groupname(by = "group", name = "A", fname = "output.test.txt")
  expect_true(grepl("group_A\\.txt$", result_multi_dots))
})

# Test Edge Cases and Error Handling
test_that("edge cases are handled correctly", {
  # Test shannon_entropy with empty matrix
  empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_equal(length(shannon_entropy(empty_matrix, axis = 1)), 0)

  # Test js_divergence with empty data
  empty_data <- data.frame()
  result_empty <- js_divergence(empty_data)
  expect_equal(nrow(result_empty), 0)

  # Test js_divergence with missing columns
  minimal_data <- data.frame(
    chrom = "chr1",
    start = 100,
    end = 150
  )
  result_minimal <- js_divergence(minimal_data)
  expect_equal(nrow(result_minimal), 0)
})

# Test Data Types and Conversions
test_that("data type handling works correctly", {
  # Test shannon_entropy with different input types
  list_input <- list(c(1, 2, 3), c(4, 5, 6))
  expect_error(shannon_entropy(list_input))

  # Test with character matrix (should convert to numeric)
  char_matrix <- matrix(c("1", "2", "3", "4"), nrow = 2, ncol = 2)
  result_char <- shannon_entropy(char_matrix, axis = 1)
  expect_true(is.numeric(result_char))
  expect_length(result_char, 2)
})

# Test Performance with Larger Data
test_that("functions handle larger datasets", {
  # Create larger test dataset
  start <- seq(1, 10000, by = 100)
  end <- start + 99

  large_data <- data.frame(
    chrom = rep("chr1", 100),
    start = start,
    end = end,
    sample1.mC = sample(1:50, 100, replace = TRUE),
    sample1.C = sample(10:100, 100, replace = TRUE),
    sample2.mC = sample(1:50, 100, replace = TRUE),
    sample2.C = sample(10:100, 100, replace = TRUE),
    sample3.mC = sample(1:50, 100, replace = TRUE),
    sample3.C = sample(10:100, 100, replace = TRUE)
  )

  # Test js_divergence with larger dataset
  result_large <- js_divergence(large_data)
  expect_true(is.data.frame(result_large))
  expect_true(nrow(result_large) <= 100)  # May be filtered

  # Test shannon_entropy with larger matrix
  large_matrix <- matrix(sample(1:1000, 1000, replace = TRUE), nrow = 100, ncol = 10)
  entropy_large <- shannon_entropy(large_matrix, axis = 1)
  expect_length(entropy_large, 100)
  expect_true(all(entropy_large >= 0))
})

# Test Integration Between Functions
test_that("functions work together correctly", {
  # Test that JS divergence results are consistent
  consistent_data <- data.frame(
    chrom = rep("chr1", 3),
    start = c(100, 200, 300),
    end = c(150, 250, 350),
    sample1.mC = c(10, 20, 30),
    sample1.C = c(20, 40, 60),
    sample2.mC = c(10, 20, 30),
    sample2.C = c(20, 40, 60)
  )

  result1 <- js_divergence(consistent_data)
  result2 <- js_divergence(consistent_data)

  expect_equal(result1$JSD_bit_, result2$JSD_bit_)
  expect_equal(result1$HMIX_bit_, result2$HMIX_bit_)

  # Test that identical samples produce zero divergence
  expect_true(all(result1$JSD_bit_ == 0))
})

# Run all tests
#cat("Running comprehensive test suite for shannon package...\n")
#test_results <- test_dir(".", reporter = "summary")
#cat("Tests completed.\n")
