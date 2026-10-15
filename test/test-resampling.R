source(if (file.exists("test/helper-simulations.R")) {
  "test/helper-simulations.R"
} else {
  "helper-simulations.R"
})
source_module("clustering_analysis.R")

testthat::test_that("blocked permutation preserves within-task pairs", {
  data <- data.frame(
    CIMBI = paste0("P", 1:20),
    task1_score = 1:20, task1_time = 101:120,
    task2_score = 201:220, task2_time = 301:320
  )
  blocks <- list(c("task1_score", "task1_time"),
                 c("task2_score", "task2_time"))
  first <- permute_task_blocks(data, blocks, seed = 151L)
  second <- permute_task_blocks(data, blocks, seed = 151L)
  testthat::expect_identical(first, second)
  testthat::expect_equal(first$task1_time - first$task1_score, rep(100, 20))
  testthat::expect_equal(first$task2_time - first$task2_score, rep(100, 20))
  testthat::expect_false(identical(first$task1_score, data$task1_score))
})

testthat::test_that("participant bootstrap preserves shape and unique execution IDs", {
  data <- data.frame(CIMBI = paste0("P", 1:12), value = 1:12)
  sample <- bootstrap_participants(data, seed = 152L)
  testthat::expect_equal(nrow(sample$data), nrow(data))
  testthat::expect_equal(anyDuplicated(sample$data$CIMBI), 0L)
  testthat::expect_true(anyDuplicated(sample$source_indices) > 0)
})

testthat::test_that("empirical inference uses plus-one p-values and percentile CIs", {
  observed <- data.frame(
    method = c("hca", "consensus"),
    mean_n_clusters = c(3, 3), mean_ari = c(0.8, 0.9),
    accuracy = c(0.75, 0.85), coverage = c(1, 0.95)
  )
  null <- do.call(rbind, lapply(1:9, function(i) transform(
    observed, iteration = i, mean_ari = mean_ari - 0.4,
    accuracy = accuracy - 0.3
  )))
  bootstrap <- do.call(rbind, lapply(1:20, function(i) transform(
    observed, iteration = i,
    mean_ari = mean_ari + (i - 10) / 100,
    accuracy = accuracy + (i - 10) / 200
  )))
  inference <- summarize_resampling_inference(observed, null, bootstrap)
  testthat::expect_true(all(inference$p_values$p_value >= 0.1))
  testthat::expect_true(all(inference$p_values$p_value <= 1))
  testthat::expect_true(all(c("method", "metric", "lower", "upper") %in%
                              names(inference$confidence_intervals)))
  testthat::expect_true(all(inference$confidence_intervals$lower <=
                              inference$confidence_intervals$upper))
})

testthat::test_that("resampling seeds are deterministic and iteration-specific", {
  first <- resampling_iteration_seeds(100L, 1L)
  second <- resampling_iteration_seeds(100L, 1L)
  next_iteration <- resampling_iteration_seeds(100L, 2L)
  testthat::expect_identical(first, second)
  testthat::expect_false(any(first %in% next_iteration))
  testthat::expect_identical(
    names(first),
    c("permutation_data", "bootstrap_data", "permutation_analysis", "bootstrap_analysis")
  )
})

testthat::test_that("two-worker resampling matches sequential execution", {
  data <- data.frame(
    CIMBI = paste0("P", 1:18),
    task1_score = 1:18, task1_time = 101:118,
    task2_score = 201:218, task2_time = 301:318
  )
  blocks <- list(c("task1_score", "task1_time"),
                 c("task2_score", "task2_time"))
  analysis_function <- local({
    adjustment <- 0.25
    function(frame, seed) {
      set.seed(seed)
      data.frame(
        method = "test",
        score = mean(frame$task1_score) + stats::runif(1L) + adjustment
      )
    }
  })
  sequential <- run_resampling_inference(
    data, analysis_function, blocks, iterations = 3L, seed = 153L, workers = 1L
  )
  parallel <- run_resampling_inference(
    data, analysis_function, blocks, iterations = 3L, seed = 153L, workers = 2L,
    worker_source_files = unname(clustering_module_files(normalize = TRUE)),
    worker_required_functions = c("permute_task_blocks", "run_cross_validation"),
    worker_required_packages = "stats"
  )
  testthat::expect_identical(parallel$permutation, sequential$permutation)
  testthat::expect_identical(parallel$bootstrap, sequential$bootstrap)
  testthat::expect_identical(parallel$inference, sequential$inference)
})

testthat::test_that("parallel resampling reruns the cross-validation workflow", {
  fixture <- make_cluster_feature_result(
    n_per_cluster = 8L, dimensions = 4L, clusters = 2L, seed = 154L
  )
  data <- data.frame(
    CIMBI = fixture$participant_ids,
    sex = rep(c("F", "M"), length.out = nrow(fixture$features)),
    medication = rep(c("No", "Yes"), length.out = nrow(fixture$features)),
    age = seq(20, 60, length.out = nrow(fixture$features)),
    fixture$features,
    check.names = FALSE
  )
  feature_names <- colnames(fixture$features)
  build_result <- function(frame) {
    features <- as.matrix(frame[, feature_names, drop = FALSE])
    rownames(features) <- frame$CIMBI
    similarity <- cosine_similarity_between_subjects(features)
    dimnames(similarity) <- list(frame$CIMBI, frame$CIMBI)
    list(
      features = features, similarity = similarity,
      distance = stats::as.dist(1 - similarity), participant_ids = frame$CIMBI
    )
  }
  pipeline <- function(train, test) {
    list(train = build_result(train), test = build_result(test), preprocessor = NULL)
  }
  base_config <- default_base_config(
    mode = "spec", methods = c("hca", "spherical_kmeans"), k_candidates = 2L,
    min_cluster_size = 2L, skmeans_restarts = 2L, optional_failure = "error"
  )
  consensus_config <- default_consensus_config(
    mode = "spec", active_methods = c("hca", "spherical_kmeans"), k_candidates = 2L,
    min_cluster_size = 2L
  )
  analysis_function <- function(frame, seed) {
    fit <- run_cross_validation(
      data = frame, distance_builder = build_result,
      predictor_variables = feature_names,
      sex_col = "sex", medication_col = "medication", age_col = "age",
      folds = 2L, base_config = base_config,
      consensus_config = consensus_config,
      prediction_mode = "feature_consensus", seed = seed,
      feature_pipeline = pipeline
    )
    fit$summary$metrics
  }
  result <- run_resampling_inference(
    data, analysis_function,
    list(feature_names[1:2], feature_names[3:4]),
    iterations = 1L, seed = 155L, workers = 2L,
    worker_source_files = unname(clustering_module_files(normalize = TRUE)),
    worker_required_functions = c(
      "run_cross_validation", "run_base_clustering", "run_consensus_clustering"
    ),
    worker_required_packages = resampling_worker_packages(base_config$methods)
  )
  testthat::expect_equal(nrow(result$permutation), 3L)
  testthat::expect_equal(nrow(result$bootstrap), 3L)
  testthat::expect_true(all(
    c("hca", "spherical_kmeans", "consensus") %in% result$permutation$method
  ))
})
