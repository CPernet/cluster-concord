# Execute every test file independently and write a reproducible Markdown report.

project_root <- normalizePath(getwd(), winslash = "/")
if (!file.exists(file.path(project_root, "SPEC.md"))) {
  stop("Run test/generate_test_report.R from the repository root.")
}

old_slow_setting <- Sys.getenv("RUN_SLOW_TESTS", unset = NA_character_)
Sys.setenv(RUN_SLOW_TESTS = "true")

# test_file() does not automatically preload helper files. Loading both helpers
# here gives independently executed files the same support environment as
# test_dir().
sys.source(file.path(project_root, "test", "helper-fixtures.R"), envir = .GlobalEnv)
sys.source(file.path(project_root, "test", "helper-simulations.R"), envir = .GlobalEnv)

test_files <- sort(list.files(
  file.path(project_root, "test"), pattern = "^test-.*\\.R$", full.names = TRUE
))
if (!length(test_files)) stop("No test files were found.")

execute_file <- function(path) {
  started <- proc.time()[["elapsed"]]
  result <- testthat::test_file(path, reporter = "silent")
  elapsed <- proc.time()[["elapsed"]] - started
  table <- as.data.frame(result, stringsAsFactors = FALSE)
  table$file <- basename(path)
  list(path = path, table = table, elapsed = elapsed)
}

executions <- lapply(test_files, execute_file)
case_results <- do.call(rbind, lapply(executions, `[[`, "table"))

case_status <- function(row) {
  if (isTRUE(row[["error"]]) || row[["failed"]] > 0L) return("FAIL")
  if (isTRUE(row[["skipped"]])) return("SKIP")
  if (row[["warning"]] > 0L) return("PASS WITH WARNING")
  "PASS"
}
case_results$status <- vapply(seq_len(nrow(case_results)), function(index) {
  case_status(case_results[index, , drop = FALSE])
}, character(1))

file_results <- do.call(rbind, lapply(executions, function(execution) {
  table <- execution$table
  failed <- sum(table$failed) + sum(table$error)
  data.frame(
    file = basename(execution$path),
    cases = nrow(table),
    assertions = sum(table$nb),
    passed = sum(table$passed),
    failed = failed,
    skipped = sum(table$skipped),
    warnings = sum(table$warning),
    seconds = round(execution$elapsed, 2L),
    status = if (failed > 0L) "FAIL" else if (all(table$skipped)) "SKIP" else "PASS",
    stringsAsFactors = FALSE
  )
}))

escape_markdown <- function(value) {
  value <- ifelse(is.na(value), "", as.character(value))
  gsub("|", "\\|", value, fixed = TRUE)
}

markdown_table <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  formatted <- data
  formatted[] <- lapply(formatted, escape_markdown)
  rows <- apply(formatted, 1L, function(row) {
    paste0("| ", paste(row, collapse = " | "), " |")
  })
  c(
    paste0("| ", paste(names(formatted), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(formatted)), collapse = " | "), " |"),
    rows
  )
}

git_commit <- tryCatch(
  system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)[1L],
  error = function(error) "unavailable"
)
working_tree_changes <- tryCatch(
  length(system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE)) > 0L,
  error = function(error) NA
)
total_failed <- sum(file_results$failed)
overall <- if (total_failed == 0L) "PASS" else "FAIL"

summary_table <- data.frame(
  result = overall,
  files = nrow(file_results),
  cases = sum(file_results$cases),
  assertions = sum(file_results$assertions),
  passed_expectations = sum(file_results$passed),
  failed_expectations_or_cases = total_failed,
  skipped_cases = sum(file_results$skipped),
  warnings = sum(file_results$warnings),
  elapsed_seconds = round(sum(file_results$seconds), 2L),
  stringsAsFactors = FALSE
)

report_lines <- c(
  "# Consensus-clustering test report",
  "",
  "## Execution metadata",
  "",
  sprintf("- Generated: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("- Git base commit: `%s`", git_commit),
  sprintf("- Working-tree changes included in this run: `%s`", working_tree_changes),
  sprintf("- R: `%s`", R.version.string),
  sprintf("- testthat: `%s`", as.character(utils::packageVersion("testthat"))),
  "- Command: `Rscript test/generate_test_report.R`",
  "- `RUN_SLOW_TESTS=true`: every test file, including the known-truth full pipeline, was executed.",
  "- Inputs: deterministic generated fixtures and temporary files; no production workbook was read or modified.",
  "",
  "## Overall result",
  "",
  markdown_table(summary_table),
  "",
  "## Result by test file",
  "",
  markdown_table(file_results),
  "",
  "## Validation cases",
  ""
)

for (execution in executions) {
  table <- execution$table
  table$status <- vapply(seq_len(nrow(table)), function(index) {
    case_status(table[index, , drop = FALSE])
  }, character(1))
  detail <- data.frame(
    validation = table$test,
    assertions = table$nb,
    passed = table$passed,
    failed = table$failed,
    warnings = table$warning,
    seconds = round(table$real, 2L),
    status = table$status,
    stringsAsFactors = FALSE
  )
  report_lines <- c(
    report_lines,
    sprintf("### `%s`", basename(execution$path)),
    "",
    markdown_table(detail),
    ""
  )
}

report_lines <- c(
  report_lines,
  "## Interpretation",
  "",
  if (overall == "PASS") {
    "Every validation case passed, including the all-five-method synthetic pipeline and held-out consensus prediction."
  } else {
    "At least one validation failed. The failing file and case are identified above; this report was still written before returning a nonzero status."
  },
  "",
  "Detailed test intent and coverage are documented in `test/README.md`."
)

report_path <- file.path(project_root, "test", "test_report.md")
writeLines(report_lines, report_path, useBytes = TRUE)

if (is.na(old_slow_setting)) {
  Sys.unsetenv("RUN_SLOW_TESTS")
} else {
  Sys.setenv(RUN_SLOW_TESTS = old_slow_setting)
}

cat("Wrote", normalizePath(report_path, winslash = "/"), "\n")
if (total_failed > 0L) quit(status = 1L)
