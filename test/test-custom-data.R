testthat::test_that("custom balancing supports arbitrary and absent covariates", {
  d <- data.frame(CIMBI = paste0("p", 1:30), value = seq_len(30),
                  site = rep(c("A", "B", NA), 10), group = rep(1:2, 15))
  a <- make_balanced_folds(d, sex_col = NULL, medication_col = NULL,
                           age_col = NULL, balance_cols = c("site", "group"))
  testthat::expect_equal(a$fold_sizes, rep(6L, 5))
  testthat::expect_true(all(c("site", "group") %in% names(a$balance_long)))
  b <- make_balanced_folds(d, sex_col = NULL, medication_col = NULL,
                           age_col = NULL, balance_cols = character())
  testthat::expect_equal(b$fold_sizes, rep(6L, 5))
  testthat::expect_error(make_balanced_folds(d, sex_col = NULL, medication_col = NULL,
                       age_col = NULL, balance_cols = "absent"), "Missing fold")
})

testthat::test_that("explicit legacy columns preserve fold allocation and age settings", {
  d <- data.frame(CIMBI = paste0("p", 1:30), value = seq_len(30),
                  sex = rep(c("F", "M"), 15), med = rep(c("Y", "N", "N"), 10), age = 18:47)
  a <- make_balanced_folds(d, sex_col = "sex", medication_col = "med", age_col = "age")
  b <- make_balanced_folds(d, sex_col = "sex", medication_col = "med", age_col = "age",
                           balance_cols = c("sex", "med"))
  testthat::expect_identical(a$fold_id, b$fold_id)
  c <- make_balanced_folds(d, sex_col = "sex", medication_col = "med", age_col = "age",
                           age_bin_width = 5, age_bin_origin = 18)
  testthat::expect_equal(diff(c$age_breaks), rep(5, length(c$age_breaks) - 1))
})

testthat::test_that("custom adapters run without behavioural columns and preserve training-only transforms", {
  root <- if (file.exists("src/custom_data_example.R")) "." else ".."
  source(file.path(root, "src/custom_data_example.R"))
  result <- run_custom_data_example()
  testthat::expect_length(result$cv$data_ids, 75L)
  testthat::expect_false(is.null(result$generalization))
  testthat::expect_null(result$inference)
  train <- data.frame(CIMBI=c("a", "b", "c"), x=c(1, NA, 3), y=c(3, 2, 8))
  test <- data.frame(CIMBI="d", x=1000, y=NA_real_)
  pipeline <- make_numeric_pipeline(c("x", "y"))
  fitted <- pipeline(train, test)
  testthat::expect_equal(unname(fitted$preprocessor$medians), c(2, 3))
  testthat::expect_equal(unname(fitted$preprocessor$center), c(2, 13/3))
  test$x <- -1000
  testthat::expect_identical(fitted$train, pipeline(train, test)$train)
})

testthat::test_that("custom loader retains its pipeline in the inference branch", {
  root <- if (file.exists("src/custom_data_example.R")) "." else ".."
  source(file.path(root, "src/custom_data_example.R"))
  set.seed(81)
  d <- data.frame(CIMBI=paste0("p", 1:40), x=rnorm(40), y=rnorm(40), z=rnorm(40))
  loader <- function(input_path, ...) list(screened_development_raw_data=d,
                                          screened_holdout_data=d[FALSE, ])
  # Exercise selection of the custom inference pipeline without running resampling.
  target_env <- environment(clustering_analysis)
  original_inference <- target_env$run_resampling_inference
  withr::defer(assign("run_resampling_inference", original_inference, envir=target_env))
  target_env$run_resampling_inference <- function(data, analysis_function, ...) {
    list(observed = analysis_function(data, 123L))
  }
  result <- clustering_analysis(input_path="unused", data_loader=loader,
    feature_pipeline=make_numeric_pipeline(c("x", "y", "z")),
    predictor_variables=c("x", "y", "z"), balance_cols=character(), age_col=NULL,
    methods=c("hca", "spherical_kmeans"), k_candidates=2L,
    task_blocks=list("x", "y", "z"), resamples=1L, write_outputs=FALSE)
  testthat::expect_false(is.null(result$inference$observed))
  testthat::expect_error(clustering_analysis(input_path="unused", data_loader=loader),
                         "requires feature_pipeline")
})
