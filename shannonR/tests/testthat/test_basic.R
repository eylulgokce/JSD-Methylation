test_that("Constants are correctly defined", {
  expect_equal(LOG2E, log2(exp(1)))
  expect_true(is.numeric(LOG2E))
  expect_true(LOG2E > 0)
})

test_that("Shannon entropy calculation works", {
  # Create a simple count matrix
  count_matrix <- matrix(c(10, 5, 3,
                           8, 7, 2,
                           12, 3, 4),
                         nrow = 3, byrow = TRUE)

  # Test entropy calculation
  entropy <- shannon_entropy(count_matrix, axis = 2, method = "plug-in")

  expect_true(is.numeric(entropy))
  expect_true(all(entropy >= 0))
  expect_equal(length(entropy), nrow(count_matrix))
})

test_that("Imputation function works", {
  test_data <- data.frame(a = c(1, 2, NA), b = c(NA, 4, 5))

  result <- impute(test_data, method = "pseudocount")
  expect_equal(result, 1)
})

test_that("Group name generation works", {
  result <- groupname(by = c("tissue", "condition"),
                      name = c("leaf", "control"),
                      fname = "results.txt")

  expect_true(is.character(result))
  expect_true(grepl("tissue_leaf_and_condition_control", result))
})

test_that("Population filter handles basic cases", {
  # Create temporary metadata file
  temp_file <- tempfile(fileext = ".txt")
  metadata <- data.frame(
    sample = c("s1", "s2", "s3"),
    condition = c("ctrl", "treat", "ctrl"),
    stringsAsFactors = FALSE
  )
  write.table(metadata, temp_file, sep = "\t", row.names = FALSE)

  # Test basic filtering
  result <- population_filter(temp_file)
  expect_equal(result$reference, c("s1", "s2", "s3"))

  # Cleanup
  unlink(temp_file)
})
