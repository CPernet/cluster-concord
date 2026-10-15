source(if (file.exists("test/helper-simulations.R")) {
  "test/helper-simulations.R"
} else {
  "helper-simulations.R"
})
source_module("clustering_analysis.R")

testthat::test_that("default analysis output is anchored to the repository root", {
  result_candidates <- c(
    file.path(getwd(), "results"),
    file.path(getwd(), "..", "results")
  )
  expected <- normalizePath(
    result_candidates[dir.exists(result_candidates)][1L],
    winslash = "/", mustWork = TRUE
  )
  original_directory <- getwd()
  on.exit(setwd(original_directory), add = TRUE)
  setwd(tempdir())

  testthat::expect_identical(default_clustering_output_dir(), expected)
  testthat::expect_identical(
    deparse(formals(clustering_analysis)$output_dir),
    "default_clustering_output_dir()"
  )
})

testthat::test_that("canonical entry point writes structured compatibility artifacts", {
  fixture <- make_cluster_feature_result(n_per_cluster = 10L, seed = 111L)
  data <- data.frame(
    CIMBI = fixture$participant_ids,
    Gender = rep(c("F", "M"), length.out = length(fixture$participant_ids)),
    Antidepressant = rep(c("No", "Yes"), length.out = length(fixture$participant_ids)),
    Age = seq(20, 68, length.out = length(fixture$participant_ids)),
    fixture$features,
    check.names = FALSE
  )
  feature_names <- colnames(fixture$features)
  builder <- function(frame) {
    features <- as.matrix(frame[, feature_names, drop = FALSE])
    rownames(features) <- frame$CIMBI
    similarity <- cosine_similarity_between_subjects(features)
    dimnames(similarity) <- list(frame$CIMBI, frame$CIMBI)
    list(features = features, similarity = similarity,
         distance = stats::as.dist(1 - similarity), participant_ids = frame$CIMBI)
  }
  output <- file.path(tempdir(), "consensus-entrypoint-test")
  result <- clustering_analysis(
    data = data,
    distance_builder = builder,
    predictor_variables = feature_names,
    sex_col = "Gender", medication_col = "Antidepressant", age_col = "Age",
    folds = 3L, mode = "compatibility", resamples = 0L, seed = 112L,
    output_dir = output,
    base_config = default_base_config(
      methods = base_method_names(), k_candidates = 2:4,
      min_cluster_size = 2L, skmeans_restarts = 3L,
      ap_grid_size = 15L, hdbscan_min_pts = 4L
    ),
    consensus_config = default_consensus_config(
      mode = "compatibility", k_candidates = 2:4, min_cluster_size = 2L
    )
  )
  testthat::expect_s3_class(result, "clustering_analysis_result")
  testthat::expect_true(file.exists(file.path(output, "clustering_analysis_result.rds")))
  testthat::expect_true(file.exists(file.path(output, "cv_cluster_counts.csv")))
  testthat::expect_true(file.exists(file.path(output, "cv_predictions.csv")))
  testthat::expect_true(file.exists(file.path(output, "cv_final_labels.csv")))
  final_labels <- utils::read.csv(file.path(output, "cv_final_labels.csv"))
  testthat::expect_true(all(c("CIMBI", "method", "fold_1", "fold_2",
                              "fold_3", "final_label", "stability",
                              "stable", "n_assignments") %in%
                              names(final_labels)))
  testthat::expect_false("stringsAsFactors" %in% names(final_labels))
  testthat::expect_true(all(c("method", "fold", "n_clusters") %in%
                              names(result$cv$summary$cluster_counts)))
})

testthat::test_that("SPEC entry point uses all methods and evaluates a fixed holdout", {
  fixture <- make_cluster_feature_result(n_per_cluster = 14L, seed = 171L)
  holdout_index <- unlist(lapply(0:2, function(group) group * 14L + 13:14))
  make_data <- function(indices) data.frame(
    CIMBI = fixture$participant_ids[indices],
    Gender = rep(c("F", "M"), length.out = length(indices)),
    Antidepressant = rep(c("No", "Yes"), length.out = length(indices)),
    Age = seq(20, 68, length.out = length(indices)),
    fixture$features[indices, , drop = FALSE],
    check.names = FALSE
  )
  development <- make_data(setdiff(seq_len(nrow(fixture$features)), holdout_index))
  holdout <- make_data(holdout_index)
  feature_names <- colnames(fixture$features)
  build_result <- function(frame) {
    features <- as.matrix(frame[, feature_names, drop = FALSE])
    rownames(features) <- frame$CIMBI
    similarity <- cosine_similarity_between_subjects(features)
    dimnames(similarity) <- list(frame$CIMBI, frame$CIMBI)
    list(features = features, similarity = similarity,
         distance = stats::as.dist(1 - similarity), participant_ids = frame$CIMBI)
  }
  pipeline <- function(train, test) {
    list(train = build_result(train), test = build_result(test), preprocessor = NULL)
  }
  output <- file.path(tempdir(), "consensus-spec-entrypoint-test")
  result <- clustering_analysis(
    data = development, generalization_data = holdout,
    distance_builder = build_result, feature_pipeline = pipeline,
    predictor_variables = feature_names,
    sex_col = "Gender", medication_col = "Antidepressant", age_col = "Age",
    folds = 3L, mode = "spec", resamples = 0L,
    generalization_permutations = 99L, seed = 172L, output_dir = output,
    base_config = default_base_config(
      mode = "spec", methods = base_method_names(), k_candidates = 2:4,
      min_cluster_size = 2L, skmeans_restarts = 3L,
      ap_grid_size = 12L, hdbscan_min_pts = 4L,
      optional_failure = "error"
    ),
    consensus_config = default_consensus_config(
      mode = "spec", k_candidates = 2:4, min_cluster_size = 2L
    )
  )
  testthat::expect_identical(result$final_fit$consensus$active_methods,
                             base_method_names())
  testthat::expect_equal(result$generalization$assigned_n +
                           result$generalization$excluded_n, nrow(holdout))
  testthat::expect_true(file.exists(file.path(output, "generalization_calibration.csv")))
  testthat::expect_true(all(c(
    "holdout_ci_lower", "holdout_ci_upper"
  ) %in% names(result$generalization$proportions)))
})
