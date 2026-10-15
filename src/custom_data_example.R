# Run from the repository root. This example uses synthetic numeric data only.
source(if (file.exists("src/R/clustering_analysis.R")) {
  "src/R/clustering_analysis.R"
} else {
  "../src/R/clustering_analysis.R"
})

# Choose a representation appropriate to your question. This example deliberately
# uses cosine similarity to agree with the existing held-out prediction adapters.
make_numeric_pipeline <- function(columns) {
  force(columns)
  function(train, test) {
    x <- as.matrix(train[, columns, drop = FALSE])
    if (!is.numeric(x)) stop("Example features must be numeric.")
    medians <- apply(x, 2, median, na.rm = TRUE)
    if (any(!is.finite(medians))) stop("Training features cannot be entirely missing.")
    impute <- function(frame) {
      z <- as.matrix(frame[, columns, drop = FALSE])
      for (j in seq_along(columns)) z[is.na(z[, j]), j] <- medians[j]
      if (any(!is.finite(z))) stop("Features must be finite after imputation.")
      z
    }
    x <- impute(train)
    center <- colMeans(x)
    scale <- apply(x, 2, sd)
    scale[!is.finite(scale) | scale == 0] <- 1
    transform <- function(frame) {
      features <- sweep(sweep(impute(frame), 2, center), 2, scale, "/")
      rownames(features) <- as.character(frame$CIMBI)
      norms <- sqrt(rowSums(features^2))
      if (any(norms == 0)) stop("Cosine representation requires nonzero feature rows.")
      unit <- features / norms
      similarity <- tcrossprod(unit)
      similarity[similarity > 1] <- 1
      similarity[similarity < -1] <- -1
      diag(similarity) <- 1
      list(features = features, similarity = similarity,
           distance = as.dist(1 - similarity), participant_ids = rownames(features))
    }
    list(train = transform(train), test = transform(test),
         preprocessor = list(medians = medians, center = center, scale = scale))
  }
}

run_custom_data_example <- function() {
  set.seed(51)
  example <- data.frame(CIMBI = paste0("example_", 1:90),
                        site = rep(c("A", "B", "C"), 30),
                        x = rnorm(90), y = rnorm(90), z = rnorm(90))
  clustering_analysis(
    data = example[1:75, ], generalization_data = example[76:90, ],
    predictor_variables = c("x", "y", "z"),
    feature_pipeline = make_numeric_pipeline(c("x", "y", "z")),
    balance_cols = "site", age_col = NULL, holdout_cohort = NULL,
    methods = c("hca", "spherical_kmeans"), k_candidates = 2:3,
    resamples = 0L, write_outputs = FALSE
  )
}
