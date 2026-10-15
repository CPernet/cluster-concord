source(if (file.exists("test/helper-simulations.R")) {
  "test/helper-simulations.R"
} else {
  "helper-simulations.R"
})
source_module("base_clustering.R")
source_module("consensus_clustering.R")
source_module("cv_clustering.R")

testthat::test_that("balanced folds are nearly equal and reproducible", {
  n <- 103L
  data <- data.frame(
    CIMBI = sprintf("P%03d", seq_len(n)),
    sex = rep(c("F", "M"), length.out = n),
    antidepressant = rep(c("No", "Yes", "No"), length.out = n),
    age = seq(19, 70, length.out = n)
  )
  first <- make_balanced_folds(
    data, folds = 5L, sex_col = "sex", medication_col = "antidepressant",
    age_col = "age", seed = 90L
  )
  second <- make_balanced_folds(
    data, folds = 5L, sex_col = "sex", medication_col = "antidepressant",
    age_col = "age", seed = 90L
  )
  testthat::expect_identical(first$fold_id, second$fold_id)
  testthat::expect_lte(diff(range(table(first$fold_id))), 1L)
  testthat::expect_equal(sort(unique(first$fold_id)), 1:5)
  testthat::expect_true(all(c("fold", "n", "sex", "medication", "age_bin") %in%
                              names(first$balance_long)))
})

testthat::test_that("legacy alignment recovers a pure label permutation", {
  reference <- c(A = 1L, B = 1L, C = 2L, D = 2L, E = 3L, F = 3L)
  candidate <- c(A = 3L, B = 3L, C = 1L, D = 1L, E = 2L, F = 2L)
  aligned <- align_fold_labels_legacy(reference, candidate)
  testthat::expect_identical(unname(aligned$labels), unname(as.character(reference)))
  testthat::expect_false(anyNA(aligned$labels))
})

testthat::test_that("one CV fold retains every base assignment", {
  fixture <- make_cluster_feature_result(n_per_cluster = 10L, seed = 92L)
  data <- data.frame(
    CIMBI = fixture$participant_ids,
    sex = rep(c("F", "M"), length.out = length(fixture$participant_ids)),
    antidepressant = rep(c("No", "Yes"), each = 3L, length.out = length(fixture$participant_ids)),
    age = seq(20, 65, length.out = length(fixture$participant_ids)),
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
  fold <- run_cv_fold(
    data = data,
    train_idx = 1:24,
    test_idx = 25:30,
    fold = 1L,
    distance_builder = builder,
    predictor_variables = feature_names,
    base_config = default_base_config(
      methods = base_method_names(), k_candidates = 2:4,
      min_cluster_size = 2L, skmeans_restarts = 3L,
      ap_grid_size = 15L, hdbscan_min_pts = 4L
    ),
    consensus_config = default_consensus_config(
      mode = "compatibility", k_candidates = 2:4, min_cluster_size = 2L
    ),
    prediction_mode = "legacy_raw_centroid",
    seed = 93L
  )
  testthat::expect_identical(names(fold$training_assignments), base_method_names())
  testthat::expect_equal(nrow(fold$training_assignments), 24L)
  testthat::expect_equal(length(fold$consensus_training), 24L)
  testthat::expect_equal(length(fold$consensus_prediction), 6L)
  testthat::expect_true(all(base_method_names() %in% names(fold$cluster_counts)))
})

