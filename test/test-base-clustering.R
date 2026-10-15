source(if (file.exists("test/helper-simulations.R")) {
  "test/helper-simulations.R"
} else {
  "helper-simulations.R"
})
source_module("base_clustering.R")

testthat::test_that("base clustering exposes a deterministic structured contract", {
  fixture <- make_cluster_feature_result(seed = 12L)
  config <- default_base_config(
    mode = "compatibility",
    methods = c("hca", "spherical_kmeans", "gmm", "affinity_propagation", "hdbscan"),
    k_candidates = 2:4,
    min_cluster_size = 3L,
    skmeans_restarts = 5L,
    ap_grid_size = 20L,
    hdbscan_min_pts = 5L
  )
  result <- run_base_clustering(fixture, config = config, seed = 44L)

  testthat::expect_s3_class(result, "base_clustering_result")
  testthat::expect_identical(result$participant_ids, fixture$participant_ids)
  testthat::expect_identical(names(result$methods), config$methods)
  testthat::expect_equal(nrow(result$assignments), length(fixture$participant_ids))
  testthat::expect_identical(rownames(result$assignments), fixture$participant_ids)
  testthat::expect_true(all(config$methods %in% names(result$cluster_counts)))
  testthat::expect_identical(result$config$hca_linkage, "ward.D2")
  testthat::expect_identical(result$seed, 44L)
})

testthat::test_that("base clustering validates participant ordering and candidate k", {
  fixture <- make_cluster_feature_result(n_per_cluster = 3L, seed = 15L)
  broken <- fixture
  rownames(broken$features)[2] <- rownames(broken$features)[1]
  testthat::expect_error(
    run_base_clustering(broken, default_base_config(methods = "hca")),
    "unique"
  )
  testthat::expect_error(
    default_base_config(methods = "hca", k_candidates = c(1L, 2L)),
    "at least 2"
  )
})

testthat::test_that("HDBSCAN noise is unavailable rather than a regular cluster", {
  testthat::skip_if_not_installed("dbscan")
  fixture <- make_cluster_feature_result(n_per_cluster = 8L, seed = 18L)
  outlier <- matrix(-4, nrow = 1, ncol = ncol(fixture$features))
  rownames(outlier) <- "NOISE"
  colnames(outlier) <- colnames(fixture$features)
  features <- rbind(fixture$features, outlier)
  similarity <- cosine_similarity_between_subjects(features)
  dimnames(similarity) <- list(rownames(features), rownames(features))
  extended <- list(
    features = features,
    similarity = similarity,
    distance = stats::as.dist(1 - similarity),
    participant_ids = rownames(features)
  )
  result <- run_base_clustering(
    extended,
    default_base_config(methods = "hdbscan", hdbscan_min_pts = 6L),
    seed = 19L
  )
  raw <- result$methods$hdbscan$raw_assignments
  labels <- result$methods$hdbscan$assignments
  testthat::expect_true(all(is.na(labels[raw == 0L])))
})

testthat::test_that("SPEC mode uses UPGMA for the base HCA", {
  fixture <- make_cluster_feature_result(seed = 21L)
  config <- default_base_config(
    mode = "spec", methods = "hca", k_candidates = 2:4,
    min_cluster_size = 3L
  )
  result <- run_base_clustering(fixture, config, seed = 22L)
  testthat::expect_identical(config$hca_linkage, "average")
  testthat::expect_identical(result$methods$hca$prediction_metadata$linkage, "average")
  testthat::expect_identical(result$methods$hca$model$method, "average")
})
