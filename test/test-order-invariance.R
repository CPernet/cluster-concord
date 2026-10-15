
source(if (file.exists("test/helper-simulations.R")) "test/helper-simulations.R" else "helper-simulations.R")
source_module("clustering_analysis.R")

testthat::test_that("holdout and preparation survive row shuffle and ID replacement", {
  raw <- make_fake_raw_data(n = 40L, seed = 921L)
  a_path <- write_fake_workbook(raw, file.path(tempdir(), "invariance-original.xlsx"))
  shuffled <- raw[rev(seq_len(nrow(raw))), , drop = FALSE]
  original <- shuffled[["CIMBI ID"]]
  shuffled[["CIMBI ID"]] <- paste0("ANON", seq_len(nrow(shuffled)))
  inverse <- stats::setNames(original, shuffled[["CIMBI ID"]])
  b_path <- write_fake_workbook(shuffled, file.path(tempdir(), "invariance-anonymous.xlsx"))
  a <- load_and_prepdata(a_path, write_outputs = FALSE, return = "all")
  b <- load_and_prepdata(b_path, write_outputs = FALSE, return = "all")
  testthat::expect_identical(a$holdout_data$CIMBI, unname(inverse[b$holdout_data$CIMBI]))
  testthat::expect_identical(a$development_raw_data$CIMBI,
                            unname(inverse[b$development_raw_data$CIMBI]))
  testthat::expect_identical(unname(a$features), unname(b$features))
  zero <- load_and_prepdata(a_path, holdout_fraction = 0,
                           write_outputs = FALSE, return = "all")
  testthat::expect_equal(nrow(zero$development_raw_data), nrow(raw))
  testthat::expect_equal(nrow(zero$holdout_data), 0L)
  c <- load_and_prepdata(a_path, seed = 922L, write_outputs = FALSE, return = "all")
  testthat::expect_false(setequal(a$holdout_data$CIMBI, c$holdout_data$CIMBI))
})

testthat::test_that("fold membership is invariant across ten shuffles and renamed IDs", {
  fixture <- make_cluster_feature_result(n_per_cluster = 20L, seed = 923L)
  data <- data.frame(CIMBI = fixture$participant_ids,
                     sex = rep(c("F", "M"), length.out = 60L),
                     medication = rep(c("No", "Yes", NA), length.out = 60L),
                     age = seq(18, 66, length.out = 60L),
                     fixture$features, check.names = FALSE)
  plan <- function(d, seed = 2025L) make_balanced_folds(
    d, 5L, "sex", "medication", "age", seed
  )
  expected <- plan(data)$fold_id
  for (i in seq_len(10L)) {
    set.seed(100L + i)
    perm <- sample.int(nrow(data))
    changed <- data[perm, , drop = FALSE]
    changed$CIMBI <- paste0("replacement", seq_len(nrow(changed)))
    actual <- plan(changed)$fold_id
    testthat::expect_identical(actual[order(perm)], expected)
  }
  testthat::expect_false(identical(plan(data, 2026L)$fold_id, expected))
  testthat::expect_equal(sort(as.integer(table(expected))), c(12L, 12L, 12L, 12L, 12L),
                         ignore_attr = TRUE)
  # Age recoding must never change internal record ordering.
  rounded <- data
  rounded$age <- 18 + 3 * floor((rounded$age - 18) / 3) + 1.5
  testthat::expect_identical(canonical_record_order(data, age_col = "age"),
                            canonical_record_order(rounded, age_col = "age"))
})

testthat::test_that("actual SPEC model paths preserve partitions after ID replacement", {
  fixture <- make_cluster_feature_result(n_per_cluster = 12L, seed = 924L)
  data <- data.frame(CIMBI = fixture$participant_ids,
                     sex = rep(c("F", "M"), length.out = 36L),
                     medication = rep(c("No", "Yes"), length.out = 36L),
                     age = seq(18, 65, length.out = 36L),
                     fixture$features, check.names = FALSE)
  predictors <- colnames(fixture$features)
  builder <- function(d) {
    features <- as.matrix(d[, predictors, drop = FALSE])
    rownames(features) <- d$CIMBI
    similarity <- cosine_similarity_between_subjects(features)
    dimnames(similarity) <- list(d$CIMBI, d$CIMBI)
    list(features = features, similarity = similarity,
         distance = stats::as.dist(1 - similarity), participant_ids = d$CIMBI)
  }
  pipeline <- function(train, test) {
    list(train = builder(train), test = builder(test), preprocessor = NULL)
  }
  run <- function(d, holdout) clustering_analysis(
    data = d, generalization_data = holdout,
    distance_builder = builder, feature_pipeline = pipeline,
    predictor_variables = predictors, sex_col = "sex",
    medication_col = "medication", age_col = "age",
    folds = 3L, mode = "spec", resamples = 0L, seed = 2025L,
    write_outputs = FALSE, generalization_method = "exact",
    base_config = default_base_config(
      mode = "spec", methods = base_method_names(), k_candidates = 2:3,
      min_cluster_size = 2L, skmeans_restarts = 2L,
      ap_grid_size = 8L, hdbscan_min_pts = 4L, optional_failure = "error"),
    consensus_config = default_consensus_config(
      mode = "spec", k_candidates = 2:3, min_cluster_size = 2L))
  a <- run(data[4:36, , drop = FALSE], data[1:3, , drop = FALSE])
  changed <- data
  changed$CIMBI <- paste0("ANON", seq_len(nrow(changed)))
  inverse <- stats::setNames(data$CIMBI, changed$CIMBI)
  b <- run(changed[36:4, , drop = FALSE], changed[3:1, , drop = FALSE])
  # Compare membership by the private test mapping, not literal ID strings.
  testthat::expect_identical(a$cv$data_ids, unname(inverse[b$cv$data_ids]))
  testthat::expect_identical(a$cv$fold_id, b$cv$fold_id)
  for (i in seq_len(3L)) {
    testthat::expect_identical(a$cv$fold_fits[[i]]$base$cluster_counts,
                              b$cv$fold_fits[[i]]$base$cluster_counts)
    testthat::expect_identical(unname(as.matrix(a$cv$fold_fits[[i]]$training_assignments)),
                              unname(as.matrix(b$cv$fold_fits[[i]]$training_assignments)))
    testthat::expect_identical(unname(as.matrix(a$cv$fold_fits[[i]]$test_assignments)),
                              unname(as.matrix(b$cv$fold_fits[[i]]$test_assignments)))
  }
  testthat::expect_identical(unname(as.matrix(a$final_fit$base$assignments)),
                            unname(as.matrix(b$final_fit$base$assignments)))
  testthat::expect_identical(unname(a$final_fit$consensus$labels),
                            unname(b$final_fit$consensus$labels))
  testthat::expect_identical(a$cv$summary$metrics, b$cv$summary$metrics)
  testthat::expect_identical(a$generalization$proportions, b$generalization$proportions)
  testthat::expect_identical(a$generalization$test, b$generalization$test)
})

testthat::test_that("ambiguous records are reported without silently deleting rows", {
  d <- data.frame(CIMBI = c("a", "b"), age = c(20, 21), value = c(3, 3))
  testthat::expect_error(canonicalize_records(d, age_col = "age",
                                             require_distinct = TRUE),
                        "independent persistent analysis token")
  testthat::expect_equal(nrow(d), 2L)
})

testthat::test_that("fixed age bands cover arbitrary input ranges and tolerate nested recoding", {
  d <- data.frame(CIMBI = letters[1:9], value = seq_len(9L),
                  sex = "F", medication = "No",
                  age = c(17, 18, 23.9, 29.9, 30, 65, 75, NA, Inf))
  plan <- make_balanced_folds(d, 3L, "sex", "medication", "age", 2025L)
  testthat::expect_identical(plan$age_breaks, seq(6, 78, 12))
  testthat::expect_identical(plan$age_binning$width, 12)
  testthat::expect_identical(plan$age_binning$origin, 18)
  testthat::expect_true(all(c("[6,18)", "[18,30)", "[30,42)", "[54,66)",
                            "[66,78)", "Missing") %in% plan$age_band_counts$age_bin))
  testthat::expect_equal(sum(plan$age_band_counts$Freq), nrow(d))
  for (width in c(2, 3)) {
    recoded <- d
    finite <- is.finite(d$age)
    recoded$age[finite] <- 18 + width * floor((d$age[finite] - 18) / width) + width / 2
    altered <- make_balanced_folds(recoded, 3L, "sex", "medication", "age", 2025L)
    testthat::expect_identical(altered$fold_id, plan$fold_id)
  }
  d$age <- NA_real_
  missing <- make_balanced_folds(d, 3L, "sex", "medication", "age", 2025L)
  testthat::expect_identical(missing$age_breaks, numeric())
  testthat::expect_equal(missing$age_band_counts$Freq, 9L)
  testthat::expect_error(make_balanced_folds(d, 3L, "sex", "medication", "age",
                                            age_bin_width = 0),
                        "positive and finite")
  testthat::expect_true(all(plan$sparse_strata$n < 3L))
})

