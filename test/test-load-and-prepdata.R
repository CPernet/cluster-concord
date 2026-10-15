source(if (file.exists("test/helper-fixtures.R")) "test/helper-fixtures.R" else "helper-fixtures.R")

testthat::test_that("loader creates a finite distance matrix and documented output artifacts", {
  raw_data <- make_fake_raw_data(n = 30L, seed = 11L)
  workbook_path <- write_fake_workbook(raw_data, file.path(tempdir(), "fake-cognitive-input.xlsx"))
  output_dir <- file.path(tempdir(), "fake-cognitive-results")

  result <- load_and_prepdata(
    workbook_path,
    holdout_fraction = 0.10,
    output_dir = output_dir,
    return = "all"
  )
  metadata <- jsonlite::read_json(result$output_paths$metadata_json, simplifyVector = TRUE)

  testthat::expect_s3_class(result$distance, "dist")
  testthat::expect_false(anyNA(result$distance))
  testthat::expect_true(all(is.finite(result$distance)))
  testthat::expect_setequal(c(result$prepared_data$CIMBI, result$holdout_data$CIMBI), as.character(raw_data[["CIMBI ID"]]))
  testthat::expect_equal(
    nrow(result$holdout_data),
    floor(result$preparation_audit$eligible_holdout_n * 0.10)
  )
  testthat::expect_equal(
    nrow(result$development_raw_data) + nrow(result$holdout_data),
    result$preparation_audit$screened_total_n
  )
  testthat::expect_identical(result$holdout_data,
                             result$screened_holdout_data)
  testthat::expect_equal(ncol(result$features), choose(nrow(result$task_spec), 2L))
  testthat::expect_true(file.exists(result$output_paths$distance_matrix_tsv))
  testthat::expect_true(file.exists(result$output_paths$metadata_json))
  testthat::expect_true(all(c("cosine", "pearson_product", "no_shared_dimension") %in% metadata$feature_metadata$method))
})


testthat::test_that("visible NA markers preserve missing-test screening and prepared data", {
  blank <- make_fake_raw_data(n = 30L, seed = 711L)
  cols <- cognitive_imputation_columns()
  for (column in cols) {
    blank[[column]] <- as.character(blank[[column]])
    blank[[column]][1L] <- NA_character_
  }
  marked <- blank
  for (column in cols) marked[[column]][1L] <- "NA"
  p1 <- write_fake_workbook(blank, file.path(tempdir(), "blank-marker-test.xlsx"))
  p2 <- write_fake_workbook(marked, file.path(tempdir(), "visible-marker-test.xlsx"))
  a <- load_and_prepdata(p1, write_outputs = FALSE, output_dir = tempdir(), return = "all")
  b <- load_and_prepdata(p2, write_outputs = FALSE, output_dir = tempdir(), return = "all")
  testthat::expect_false("ID00001" %in% b$screened_development_raw_data$CIMBI)
  testthat::expect_equal(b$preparation_audit$screened_total_n, 29L)
  testthat::expect_identical(a$screened_development_raw_data, b$screened_development_raw_data)
  testthat::expect_identical(a$screened_holdout_data, b$screened_holdout_data)
  testthat::expect_identical(a$prepared_data, b$prepared_data)
  testthat::expect_identical(a$features, b$features)
})


testthat::test_that("cleaned-input loader canonicalizes participants without removing IDs", {
  raw <- make_fake_raw_data(n = 30L, seed = 912L)
  raw <- raw[rev(seq_len(nrow(raw))), , drop = FALSE]
  path <- write_fake_workbook(raw, file.path(tempdir(), "cleaned-input-order.xlsx"))
  fit <- load_and_prepdata(path, write_outputs = FALSE, return = "all")
  expected <- as.character(raw[["CIMBI ID"]])
  canonical <- raw
  names(canonical)[names(canonical) == "CIMBI ID"] <- "CIMBI"
  canonical <- canonicalize_records(canonical)
  expected_dev <- canonical$CIMBI[!canonical$CIMBI %in% fit$holdout_data$CIMBI]
  testthat::expect_identical(fit$screened_development_raw_data$CIMBI, expected_dev)
  testthat::expect_setequal(c(fit$screened_development_raw_data$CIMBI, fit$screened_holdout_data$CIMBI), expected)
  testthat::expect_false("exclude_ids" %in% names(formals(load_and_prepdata)))
})


testthat::test_that("public ID header preserves analysis and cleaned input may omit source study", {
  raw <- make_fake_raw_data(n = 30L, seed = 1901L)
  old_path <- write_fake_workbook(raw, file.path(tempdir(), "original-id-header.xlsx"))
  names(raw)[names(raw) == "CIMBI ID"] <- "ID"
  new_path <- write_fake_workbook(raw, file.path(tempdir(), "public-id-header.xlsx"))
  old <- load_and_prepdata(old_path, write_outputs = FALSE, return = "all")
  renamed <- load_and_prepdata(new_path, write_outputs = FALSE, return = "all")
  testthat::expect_identical(old$prepared_data, renamed$prepared_data)
  testthat::expect_identical(old$features, renamed$features)
  testthat::expect_identical(old$holdout_data, renamed$holdout_data)
  raw[["Project ID (primary)"]] <- NULL
  clean_path <- write_fake_workbook(raw, file.path(tempdir(), "without-source-study.xlsx"))
  clean <- load_and_prepdata(clean_path, write_outputs = FALSE, return = "all")
  testthat::expect_false(clean$preparation_audit$project_filter_applied)
  testthat::expect_equal(clean$preparation_audit$initial_n, nrow(raw))
  testthat::expect_true(all(clean$holdout_data[["Overall comments"]] == "KOHORTE 2"))
  testthat::expect_true(all(is.finite(clean$features)))
  raw[["CIMBI ID"]] <- raw$ID
  ambiguous <- write_fake_workbook(raw, file.path(tempdir(), "ambiguous-id-columns.xlsx"))
  testthat::expect_error(load_and_prepdata(ambiguous, write_outputs = FALSE), "exactly one participant identifier")
})

