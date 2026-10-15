source(if (file.exists("test/helper-fixtures.R")) "test/helper-fixtures.R" else "helper-fixtures.R")

testthat::test_that("two-dimensional task pairs use cosine similarity", {
  task_matrix <- rbind(c(1, 0), c(0, 1), c(1, 1))
  similarity <- task_similarity_matrix(task_matrix, matrix(TRUE, nrow = 3, ncol = 2))

  testthat::expect_equal(similarity[1, 2], 0)
  testthat::expect_equal(similarity[1, 3], 1 / sqrt(2))
  testthat::expect_equal(similarity[2, 3], 1 / sqrt(2))
})

testthat::test_that("structural one-measure pairs use Pearson-product features", {
  task_spec <- data.frame(
    task = c("A", "B", "C"),
    perf_var = c("a_perf", "b_perf", "c_perf"),
    time_var = c("a_time", "b_time", "c_time"),
    stringsAsFactors = FALSE
  )
  data <- data.frame(
    CIMBI = as.character(1:4),
    a_perf = c(1, 2, 3, 4),
    a_time = c(4, 3, 2, 1),
    b_perf = 0,
    b_time = c(1, 2, 3, 4),
    c_perf = c(2, 4, 1, 3),
    c_time = 0
  )

  result <- compute_subject_distance_matrix(data, task_spec = task_spec)
  expected_product <- z_safe(data$a_time) * z_safe(data$b_time)

  testthat::expect_equal(unname(result$features[, "t1_vs_t2"]), expected_product)
  testthat::expect_equal(
    result$feature_metadata$method,
    c("pearson_product", "pearson_product", "no_shared_dimension")
  )
  testthat::expect_equal(sum(expected_product) / (nrow(data) - 1), stats::cor(data$a_time, data$b_time))
  testthat::expect_equal(unname(result$features[, "t2_vs_t3"]), rep(0, nrow(data)))
  testthat::expect_s3_class(result$distance, "dist")
  testthat::expect_true(all(is.finite(result$distance)))
  testthat::expect_equal(unname(diag(as.matrix(result$distance))), rep(0, nrow(data)))
})

testthat::test_that("zero task vectors preserve self-similarity and zero other-similarity", {
  task_matrix <- rbind(c(0, 0), c(1, 1))
  similarity <- task_similarity_matrix(task_matrix, matrix(TRUE, nrow = 2, ncol = 2))

  testthat::expect_equal(similarity[1, 1], 1)
  testthat::expect_equal(similarity[1, 2], 0)
  testthat::expect_equal(similarity[2, 1], 0)
})