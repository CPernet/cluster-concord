source(if (file.exists("test/helper-simulations.R")) {
  "test/helper-simulations.R"
} else {
  "helper-simulations.R"
})
source_module("generalization.R")

testthat::test_that("generalization reports calibration, coverage, and intervals", {
  development <- rep(c("1", "2", "3"), c(50, 30, 20))
  holdout <- c(rep("1", 10), rep("2", 6), rep("3", 4), NA, NA)
  result <- evaluate_generalization(
    development, holdout, method = "exact", seed = 161L
  )
  testthat::expect_equal(result$coverage, 20 / 22)
  testthat::expect_equal(sum(result$proportions$holdout_count), 20L)
  testthat::expect_equal(result$proportions$reference_proportion,
                         c(0.5, 0.3, 0.2))
  testthat::expect_identical(result$test$method,
                             "multinomial_exact_calibration")
  testthat::expect_equal(result$test$permutations, 0L)
  testthat::expect_true(result$test$p_value <= 1)
  testthat::expect_true(all(result$proportions$holdout_ci_lower >= 0))
  testthat::expect_true(all(result$proportions$holdout_ci_upper <= 1))
  testthat::expect_true(all(
    result$proportions$holdout_proportion >= result$proportions$holdout_ci_lower &
      result$proportions$holdout_proportion <= result$proportions$holdout_ci_upper
  ))
})

testthat::test_that("Monte Carlo calibration is explicit and reproducible", {
  development <- rep(c("1", "2", "3"), c(50, 30, 20))
  holdout <- c(rep("1", 5), rep("2", 2), rep("3", 1))
  first <- evaluate_generalization(
    development, holdout, method = "monte_carlo",
    permutations = 199L, seed = 162L
  )
  second <- evaluate_generalization(
    development, holdout, method = "monte_carlo",
    permutations = 199L, seed = 162L
  )
  testthat::expect_identical(first$test$method,
                             "multinomial_monte_carlo_calibration")
  testthat::expect_equal(first$test$p_value, second$test$p_value)
  testthat::expect_equal(first$test$p_resolution, 1 / 200)
  testthat::expect_length(first$null_statistics, 199L)
})

testthat::test_that("auto calibration uses exact enumeration when tractable", {
  result <- evaluate_generalization(
    rep(c("1", "2"), each = 20L),
    c(rep("1", 4L), rep("2", 2L)),
    method = "auto", exact_max_tables = 100L
  )
  testthat::expect_identical(result$test$method,
                             "multinomial_exact_calibration")
  testthat::expect_equal(result$test$possible_count_tables, 7)
})
