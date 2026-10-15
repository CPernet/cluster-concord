source(if (file.exists("test/helper-fixtures.R")) "test/helper-fixtures.R" else "helper-fixtures.R")

testthat::test_that("GAM imputation is deterministic, bounded, and leaves observed values unchanged", {
  set.seed(42)
  original <- data.frame(
    a = sin(seq_len(30) / 4),
    b = cos(seq_len(30) / 5),
    c = seq_len(30) / 10
  )
  incomplete <- original
  missing_a <- sample(seq_len(nrow(incomplete)), 2L)
  missing_b <- sample(setdiff(seq_len(nrow(incomplete)), missing_a), 2L)
  incomplete[missing_a, "a"] <- NA
  incomplete[missing_b, "b"] <- NA

  first <- gam_impute_chained(incomplete, names(incomplete), maxit = 1L, k = 3L, seed = 17L)
  second <- gam_impute_chained(incomplete, names(incomplete), maxit = 1L, k = 3L, seed = 17L)

  testthat::expect_equal(first$data, second$data)
  testthat::expect_false(anyNA(first$data))
  testthat::expect_equal(first$data$c, original$c)
  testthat::expect_true(all(first$data$a[missing_a] >= min(original$a)))
  testthat::expect_true(all(first$data$a[missing_a] <= max(original$a)))
  testthat::expect_true(all(first$data$b[missing_b] >= min(original$b)))
  testthat::expect_true(all(first$data$b[missing_b] <= max(original$b)))
  testthat::expect_equal(sum(first$audit$n_missing), 4L)
})

testthat::test_that("too few observed values use the audited observed-value fallback", {
  data <- data.frame(a = c(NA, 1:7), b = 11:18, c = 21:28)
  result <- gam_impute_chained(data, names(data), maxit = 1L, k = 3L, seed = 17L)

  testthat::expect_false(anyNA(result$data))
  testthat::expect_equal(result$audit$observed_resamples[result$audit$variable == "a"], 1L)
  testthat::expect_true(result$data$a[1] %in% 1:7)
})