source(if (file.exists("test/helper-simulations.R")) {
  "test/helper-simulations.R"
} else {
  "helper-simulations.R"
})
source_module("base_clustering.R")
source_module("consensus_clustering.R")

testthat::test_that("co-association uses pairwise available denominators", {
  assignments <- data.frame(
    hca = c(1L, 1L, 2L, 2L),
    spherical_kmeans = c(1L, 2L, 2L, NA),
    gmm = c(1L, 1L, NA, NA),
    row.names = paste0("P", 1:4)
  )
  result <- build_coassociation(assignments, no_evidence = "zero")
  testthat::expect_equal(result$available[1, 2], 3L)
  testthat::expect_equal(result$coassociation[1, 2], 2 / 3)
  testthat::expect_equal(result$available[3, 4], 1L)
  testthat::expect_equal(result$coassociation[3, 4], 1)
  testthat::expect_equal(result$available[2, 4], 1L)
  testthat::expect_equal(result$coassociation[2, 4], 0)
  testthat::expect_true(isSymmetric(result$coassociation))
})

testthat::test_that("no-evidence policy is explicit and diagnosed", {
  assignments <- data.frame(a = c(1L, NA), b = c(NA, 2L), row.names = c("A", "B"))
  preserved <- build_coassociation(assignments, no_evidence = "na")
  compatibility <- build_coassociation(assignments, no_evidence = "zero")
  testthat::expect_true(is.na(preserved$coassociation[1, 2]))
  testthat::expect_equal(compatibility$coassociation[1, 2], 0)
  testthat::expect_equal(preserved$no_evidence_pairs, 1L)
})

testthat::test_that("compatibility consensus defaults to the three legacy methods", {
  fixture <- make_cluster_feature_result(seed = 31L)
  base <- run_base_clustering(
    fixture,
    default_base_config(
      methods = base_method_names(), k_candidates = 2:4,
      min_cluster_size = 3L, skmeans_restarts = 5L,
      ap_grid_size = 20L, hdbscan_min_pts = 5L
    ),
    seed = 32L
  )
  consensus <- run_consensus_clustering(
    base,
    default_consensus_config(mode = "compatibility", k_candidates = 2:4,
                             min_cluster_size = 3L)
  )
  testthat::expect_s3_class(consensus, "consensus_clustering_result")
  testthat::expect_identical(
    consensus$active_methods,
    c("hca", "spherical_kmeans", "gmm")
  )
  testthat::expect_identical(consensus$config$linkage, "ward.D2")
  testthat::expect_equal(length(consensus$labels), nrow(base$assignments))
  testthat::expect_true(all(is.finite(consensus$distance)))
  testthat::expect_true(all(c("method_a", "method_b", "ari", "n_available") %in%
                              names(consensus$cross_method_ari)))
})

testthat::test_that("SPEC consensus includes all five available methods", {
  fixture <- make_cluster_feature_result(seed = 41L)
  base <- run_base_clustering(
    fixture,
    default_base_config(
      mode = "spec", methods = base_method_names(), k_candidates = 2:4,
      min_cluster_size = 3L, skmeans_restarts = 5L,
      ap_grid_size = 20L, hdbscan_min_pts = 5L,
      optional_failure = "error"
    ),
    seed = 42L
  )
  consensus <- run_consensus_clustering(
    base,
    default_consensus_config(
      mode = "spec", k_candidates = 2:4, min_cluster_size = 3L
    )
  )
  testthat::expect_identical(consensus$active_methods, base_method_names())
  testthat::expect_identical(consensus$config$linkage, "ward.D2")
  testthat::expect_identical(consensus$config$no_evidence, "na")
  testthat::expect_equal(unname(diag(consensus$coassociation)),
                         rep(1, nrow(base$assignments)))
})
