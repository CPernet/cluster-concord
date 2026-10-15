source_module("synthetic_validation.R")

testthat::test_that("full SPEC pipeline recovers known synthetic structure", {
  testthat::skip_if_not(identical(tolower(Sys.getenv("RUN_SLOW_TESTS")), "true"))
  validation <- run_synthetic_validation(
    seed = 2601L,
    n_per_cluster = 18L,
    folds = 5L,
    generalization_permutations = 99L,
    write_outputs = FALSE
  )

  testthat::expect_identical(
    validation$result$final_fit$consensus$active_methods,
    base_method_names()
  )
  testthat::expect_gte(validation$performance$final_consensus_ari, 0.80)
  testthat::expect_gte(validation$performance$stable_consensus_ari, 0.70)
  testthat::expect_gte(validation$performance$holdout_consensus_ari, 0.70)
  testthat::expect_gte(validation$performance$holdout_coverage, 0.80)
  testthat::expect_null(validation$result$inference)
  testthat::expect_equal(nrow(validation$subset_comparison), 26L)
  testthat::expect_true(all(c(
    "methods", "method_count", "ari", "selected_k", "is_all_five",
    "oracle_rank"
  ) %in% names(validation$subset_comparison)))
  testthat::expect_gte(
    validation$best_subset$ari,
    validation$subset_comparison$ari[validation$subset_comparison$is_all_five]
  )
  testthat::expect_identical(validation$subset_selection_scope,
                             "simulation_oracle_only")

  report_path <- file.path(tempdir(), "synthetic_validation_results.md")
  written <- write_synthetic_validation_report(validation, report_path)
  report <- paste(readLines(written, warn = FALSE), collapse = "\n")
  testthat::expect_true(file.exists(written))
  testthat::expect_match(report, "Best oracle consensus subset", fixed = TRUE)
  testthat::expect_match(report, "Participant-level recovered labels", fixed = TRUE)
  testthat::expect_match(report, validation$best_subset$methods, fixed = TRUE)
})
