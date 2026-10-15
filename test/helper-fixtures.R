loader_path <- function() {
  candidates <- c(
    file.path(getwd(), "src", "R", "load_and_prepdata.R"),
    file.path(getwd(), "..", "src", "R", "load_and_prepdata.R")
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) {
    stop("Cannot find src/R/load_and_prepdata.R")
  }
  path
}

source(loader_path())

make_fake_raw_data <- function(n = 30L, seed = 1L) {
  set.seed(seed)
  columns <- cognitive_imputation_columns()
  values <- lapply(seq_along(columns), function(index) {
    stats::rnorm(n, mean = 5 + index / 10, sd = 1)
  })
  names(values) <- columns
  data <- as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)

  rate_columns <- grep("Hit rate|False alarm rate", columns, value = TRUE)
  for (column in rate_columns) {
    data[[column]] <- stats::runif(n, min = 0.15, max = 0.85)
  }
  time_columns <- grep("Mean rt|reaction latency", columns, value = TRUE)
  for (column in time_columns) {
    data[[column]] <- stats::runif(n, min = 0.4, max = 3)
  }
  intensity_columns <- grep("Mean intensity threshold", columns, value = TRUE)
  for (column in intensity_columns) {
    data[[column]] <- stats::runif(n, min = 1, max = 14)
  }
  rating_columns <- grep("Mean rating", columns, value = TRUE)
  for (column in rating_columns) {
    data[[column]] <- stats::runif(n, min = 1, max = 7)
  }

  data.frame(
    `CIMBI ID` = sprintf("ID%05d", seq_len(n)),
    `Project ID (primary)` = "Clinical cohort",
    `Age at neuropsych test` = seq(20, length.out = n),
    `Overall comments` = rep(c("KOHORTE 2", "Other"), length.out = n),
    data,
    check.names = FALSE
  )
}

write_fake_workbook <- function(data, path) {
  writexl::write_xlsx(data, path)
  path
}