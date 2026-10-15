source(if (file.exists("test/helper-simulations.R")) {
  "test/helper-simulations.R"
} else {
  "helper-simulations.R"
})
source_module("base_clustering.R")
source_module("consensus_clustering.R")
source_module("prediction.R")
source_module("cv_clustering.R")

testthat::test_that("stable labels require at least three of four assignments", {
  matrix <- rbind(
    four_of_four = c("A", "A", "A", "A", NA),
    three_of_four = c("B", "B", "B", "C", NA),
    two_of_four = c("A", "A", "B", "B", NA),
    insufficient = c("C", "C", NA, NA, NA)
  )
  result <- derive_stable_labels(matrix, threshold = 0.75,
                                 expected_assignments = 4L)
  testthat::expect_identical(result$stable_label,
                             c("A", "B", NA_character_, NA_character_))
  testthat::expect_equal(result$stability_proportion[1:3], c(1, 0.75, 0.5))
  testthat::expect_identical(result$stable[1:3], c(TRUE, TRUE, FALSE))
})

testthat::test_that("reference fold uses the modal cluster count deterministically", {
  testthat::expect_identical(select_reference_fold(c(4L, 3L, 3L, 5L, 3L)), 2L)
  testthat::expect_identical(select_reference_fold(c(2L, 3L, 2L, 3L)), 1L)
})

testthat::test_that("SPEC CV returns method-level assignments and metrics", {
  fixture <- make_cluster_feature_result(n_per_cluster = 15L, seed = 141L)
  data <- data.frame(
    CIMBI = fixture$participant_ids,
    sex = rep(c("F", "M"), length.out = nrow(fixture$features)),
    antidepressant = rep(c("No", "Yes", "No"), length.out = nrow(fixture$features)),
    age = seq(18, 70, length.out = nrow(fixture$features)),
    fixture$features,
    check.names = FALSE
  )
  feature_names <- colnames(fixture$features)
  build_result <- function(frame) {
    features <- as.matrix(frame[, feature_names, drop = FALSE])
    rownames(features) <- frame$CIMBI
    similarity <- cosine_similarity_between_subjects(features)
    dimnames(similarity) <- list(frame$CIMBI, frame$CIMBI)
    list(features = features, similarity = similarity,
         distance = stats::as.dist(1 - similarity), participant_ids = frame$CIMBI)
  }
  feature_pipeline <- function(train, test) {
    list(train = build_result(train), test = build_result(test), preprocessor = NULL)
  }
  result <- run_cross_validation(
    data = data,
    distance_builder = build_result,
    predictor_variables = feature_names,
    sex_col = "sex", medication_col = "antidepressant", age_col = "age",
    folds = 5L,
    base_config = default_base_config(
      mode = "spec", methods = base_method_names(), k_candidates = 2:4,
      min_cluster_size = 2L, skmeans_restarts = 3L,
      ap_grid_size = 12L, hdbscan_min_pts = 4L,
      optional_failure = "error"
    ),
    consensus_config = default_consensus_config(
      mode = "spec", k_candidates = 2:4, min_cluster_size = 2L
    ),
    prediction_mode = "feature_consensus",
    feature_pipeline = feature_pipeline,
    seed = 142L
  )
  methods <- c(base_method_names(), "consensus")
  testthat::expect_identical(names(result$summary$per_method), methods)
  testthat::expect_true(all(vapply(result$fold_fits, function(fold) {
    identical(names(fold$test_assignments), methods)
  }, logical(1))))
  testthat::expect_true(all(c("method", "mean_n_clusters", "mean_ari",
                              "accuracy", "coverage") %in%
                              names(result$summary$metrics)))
  testthat::expect_equal(nrow(result$summary$final_labels), nrow(data))
  testthat::expect_equal(nrow(result$summary$fold_labels),
                         nrow(data) * length(methods))
  testthat::expect_true(all(c("CIMBI", "method", "fold_1", "fold_2",
                              "fold_3", "fold_4", "fold_5", "final_label",
                              "stability", "stable", "n_assignments") %in%
                              names(result$summary$fold_labels)))
})

