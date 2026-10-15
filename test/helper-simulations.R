source(if (file.exists("src/R/load_and_prepdata.R")) {
  "src/R/load_and_prepdata.R"
} else {
  "../src/R/load_and_prepdata.R"
})

source_module <- function(filename) {
  candidates <- c(
    file.path(getwd(), "src", "R", filename),
    file.path(getwd(), "..", "src", "R", filename)
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) stop("Cannot find module: ", filename)
  source(path)
}

make_cluster_feature_result <- function(
  n_per_cluster = 18L,
  dimensions = 8L,
  clusters = 3L,
  noise_sd = 0.08,
  seed = 101L
) {
  set.seed(seed)
  centers <- matrix(0.05, nrow = clusters, ncol = dimensions)
  for (group in seq_len(clusters)) {
    start <- ((group - 1L) * 2L) %% dimensions + 1L
    centers[group, start] <- 1
    centers[group, (start %% dimensions) + 1L] <- 0.7
  }
  truth <- rep(seq_len(clusters), each = n_per_cluster)
  features <- centers[truth, , drop = FALSE] +
    matrix(stats::rnorm(length(truth) * dimensions, sd = noise_sd),
           nrow = length(truth))
  ids <- sprintf("SIM%04d", seq_along(truth))
  rownames(features) <- ids
  colnames(features) <- sprintf("feature_%02d", seq_len(dimensions))
  similarity <- cosine_similarity_between_subjects(features)
  dimnames(similarity) <- list(ids, ids)
  distance <- stats::as.dist(1 - similarity)
  attr(distance, "Labels") <- ids
  list(
    distance = distance,
    similarity = similarity,
    features = features,
    participant_ids = ids,
    truth = truth
  )
}

