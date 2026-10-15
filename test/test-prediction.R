source(if (file.exists("test/helper-simulations.R")) {
  "test/helper-simulations.R"
} else {
  "helper-simulations.R"
})

testthat::test_that("fitted feature preprocessing is unaffected by holdout sentinels", {
  task_spec <- data.frame(
    task = c("one", "two", "three"),
    perf_var = c("one_perf", "two_perf", "three_perf"),
    time_var = c("one_time", "two_time", "three_time"),
    stringsAsFactors = FALSE
  )
  train <- data.frame(
    CIMBI = paste0("T", 1:8),
    one_perf = 1:8, one_time = 8:1,
    two_perf = seq(2, 16, 2), two_time = seq(3, 24, 3),
    three_perf = c(NA, 2:8), three_time = seq(4, 32, 4)
  )
  holdout <- train[1:2, ]
  holdout$CIMBI <- c("H1", "H2")
  holdout[, setdiff(names(holdout), "CIMBI")] <- 1e9

  fit <- fit_cognitive_preprocessor(train, task_spec = task_spec)
  transformed_train <- apply_cognitive_preprocessor(fit, train)
  transformed_holdout <- apply_cognitive_preprocessor(fit, holdout)

  testthat::expect_equal(fit$center[["one_perf"]], mean(1:8))
  testthat::expect_equal(fit$impute[["three_perf"]], stats::median(2:8))
  testthat::expect_equal(ncol(transformed_train$features), choose(3, 2))
  testthat::expect_equal(nrow(transformed_holdout$features), 2L)
  testthat::expect_identical(rownames(transformed_holdout$features), c("H1", "H2"))
  testthat::expect_true(all(is.finite(transformed_holdout$features)))
  testthat::expect_equal(fit$center[["one_perf"]], mean(1:8))
})

testthat::test_that("loader exposes raw development and untouched holdout rows", {
  raw_data <- make_fake_raw_data(n = 40L, seed = 121L)
  workbook <- write_fake_workbook(raw_data, file.path(tempdir(), "preprocessor-holdout.xlsx"))
  result <- load_and_prepdata(
    workbook, holdout_fraction = 0.10,
    write_outputs = FALSE, return = "all", seed = 122L
  )
  testthat::expect_true("development_raw_data" %in% names(result))
  testthat::expect_false(any(result$holdout_data$CIMBI %in%
                               result$development_raw_data$CIMBI))
  testthat::expect_equal(
    nrow(result$development_raw_data) + nrow(result$holdout_data),
    result$preparation_audit$screened_total_n
  )
  testthat::expect_true(all(c(
    "screened_development_raw_data", "screened_holdout_data"
  ) %in% names(result)))
})

testthat::test_that("raw cognitive imputation is fitted on training rows only", {
  raw <- make_fake_raw_data(12L, seed = 612L)
  names(raw)[names(raw) == "CIMBI ID"] <- "CIMBI"
  # Use a source variable retained by derive_cognitive_measures(); some source
  # columns are intentionally renamed and therefore cannot be inspected after apply.
  column <- "Trail making A - # of sec"
  raw[[column]][1L] <- NA_real_
  train <- raw[1:10, , drop = FALSE]
  holdout <- raw[11:12, , drop = FALSE]
  holdout[[column]] <- 1e9

  fit <- fit_raw_cognitive_preprocessor(train)
  expected <- stats::median(train[[column]], na.rm = TRUE)
  transformed <- apply_raw_cognitive_preprocessor(fit, train)

  testthat::expect_equal(fit$imputation_values[[column]], expected)
  testthat::expect_equal(transformed[[column]][1L], expected)
  testthat::expect_false(any(fit$imputation_values == 1e9))
  testthat::expect_equal(nrow(apply_raw_cognitive_preprocessor(fit, holdout)), 2L)
})

testthat::test_that("held-out prediction stays in the fitted feature and consensus spaces", {
  source_module("base_clustering.R")
  source_module("consensus_clustering.R")
  source_module("prediction.R")
  fixture <- make_cluster_feature_result(n_per_cluster = 20L, seed = 131L)
  test_indices <- unlist(lapply(0:2, function(group) group * 20L + 18:20))
  train_indices <- setdiff(seq_len(nrow(fixture$features)), test_indices)
  subset_result <- function(indices) {
    features <- fixture$features[indices, , drop = FALSE]
    similarity <- cosine_similarity_between_subjects(features)
    dimnames(similarity) <- list(rownames(features), rownames(features))
    list(features = features, similarity = similarity,
         distance = stats::as.dist(1 - similarity),
         participant_ids = rownames(features))
  }
  train_result <- subset_result(train_indices)
  test_result <- subset_result(test_indices)
  base <- run_base_clustering(
    train_result,
    default_base_config(
      mode = "spec", methods = base_method_names(), k_candidates = 2:4,
      min_cluster_size = 3L, skmeans_restarts = 5L,
      ap_grid_size = 20L, hdbscan_min_pts = 5L,
      optional_failure = "error"
    ),
    seed = 132L
  )
  consensus <- run_consensus_clustering(
    base,
    default_consensus_config(mode = "spec", k_candidates = 2:4,
                             min_cluster_size = 3L)
  )
  fixed_labels <- consensus$labels
  base_prediction <- predict_base_clusters(base, test_result)
  prediction <- predict_consensus_cluster(
    consensus, base, base_prediction,
    min_methods = 2L, min_score = -Inf, min_margin = -Inf
  )

  testthat::expect_identical(names(base_prediction$assignments), base_method_names())
  testthat::expect_equal(nrow(base_prediction$assignments), length(test_indices))
  testthat::expect_equal(length(prediction$labels), length(test_indices))
  testthat::expect_identical(consensus$labels, fixed_labels)
  testthat::expect_equal(dim(prediction$cross_coassociation),
                         c(length(test_indices), length(train_indices)))
  testthat::expect_true(all(c("best_score", "second_score", "margin",
                              "contributing_methods", "rejected", "reason") %in%
                              names(prediction$diagnostics)))
})

testthat::test_that("consensus prediction uses pairwise available agreement", {
  source_module("base_clustering.R")
  source_module("consensus_clustering.R")
  source_module("prediction.R")
  ids <- c("A", "B", "C", "D")
  assignments <- data.frame(
    hca = c(1, 1, 2, 2),
    spherical_kmeans = c(1, 1, 2, 2),
    gmm = c(1, 2, 2, 2),
    row.names = ids
  )
  base <- structure(list(
    participant_ids = ids,
    assignments = assignments,
    methods = setNames(lapply(names(assignments), function(method) {
      list(status = "available", assignments = assignments[[method]])
    }), names(assignments))
  ), class = c("base_clustering_result", "list"))
  consensus <- structure(list(
    participant_ids = ids,
    active_methods = names(assignments),
    labels = setNames(c(1, 1, 2, 2), ids)
  ), class = c("consensus_clustering_result", "list"))
  predicted <- structure(list(
    participant_ids = "X",
    assignments = data.frame(hca = 1, spherical_kmeans = 1, gmm = NA,
                             row.names = "X"),
    diagnostics = list()
  ), class = c("base_prediction_result", "list"))
  result <- predict_consensus_cluster(
    consensus, base, predicted, min_methods = 2L,
    min_score = -Inf, min_margin = -Inf
  )
  testthat::expect_equal(unname(result$labels), 1L)
  testthat::expect_equal(unname(result$cross_coassociation[1, ]), c(1, 1, 0, 0))
})
