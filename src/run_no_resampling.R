
# Run from the project root; optionally supply an input workbook as first argument.
source("src/R/clustering_analysis.R")
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[[1L]] else "data/data.xlsx"
result <- clustering_analysis(
  input_path = input, resamples = 0L,
  generalization_method = "exact", generalization_permutations = 0L,
  seed = 2025L, workers = 1L,
  age_bin_width = 12, age_bin_origin = 18
)

