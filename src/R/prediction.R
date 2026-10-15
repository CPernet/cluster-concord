# Out-of-sample prediction adapters in the fitted task-similarity feature space.

# Compute cosine similarity from held-out 45-feature rows to fixed training rows,
# enforcing the training feature order and preserving unavailable similarities.
cross_cosine_similarity <- function(new_features, training_features, eps = 1e-12) {
  new_features <- as.matrix(new_features)
  training_features <- as.matrix(training_features)
  if (!identical(colnames(new_features), colnames(training_features))) {
    missing <- setdiff(colnames(training_features), colnames(new_features))
    if (length(missing)) stop("New features are missing training columns: ", paste(missing, collapse = ", "))
    new_features <- new_features[, colnames(training_features), drop = FALSE]
  }
  new_norm <- sqrt(rowSums(new_features^2))
  train_norm <- sqrt(rowSums(training_features^2))
  similarity <- (new_features %*% t(training_features)) /
    (new_norm %o% train_norm)
  similarity[!is.finite(similarity)] <- NA_real_
  similarity[similarity > 1] <- 1
  similarity[similarity < -1] <- -1
  similarity
}

# Extend UPGMA/HCA out of sample by assigning each participant to the cluster
# with smallest average cosine dissimilarity to its training members.
predict_hca_clusters <- function(method_result, new_features) {
  metadata <- method_result$prediction_metadata
  training <- metadata$training_features
  labels <- as.character(metadata$training_assignments)
  similarity <- cross_cosine_similarity(new_features, training)
  levels <- sort(unique(labels[!is.na(labels)]))
  # Compute mean held-out-to-member cosine dissimilarity for one HCA cluster.
  mean_distance <- sapply(levels, function(label) {
    rowMeans(1 - similarity[, labels == label, drop = FALSE], na.rm = TRUE)
  })
  if (is.null(dim(mean_distance))) mean_distance <- matrix(mean_distance, ncol = length(levels))
  prediction <- levels[max.col(-mean_distance, ties.method = "first")]
  list(assignments = prediction, diagnostics = list(mean_distance = mean_distance))
}

# Assign held-out participants to the fitted normalized spherical-k-means center
# with greatest cosine similarity.
predict_spherical_kmeans_clusters <- function(method_result, new_features) {
  centers <- method_result$prediction_metadata$centers
  similarity <- cross_cosine_similarity(new_features, centers)
  prediction <- rownames(centers)[max.col(similarity, ties.method = "first")]
  list(assignments = prediction, diagnostics = list(center_similarity = similarity))
}

# Apply the training feature selection, scaling, and PCA rotation, then use the
# fitted Mclust posterior prediction without refitting on held-out data.
predict_gmm_clusters <- function(method_result, new_features) {
  metadata <- method_result$prediction_metadata
  required <- metadata$kept_features
  missing <- setdiff(required, colnames(new_features))
  if (length(missing)) stop("New GMM features are missing: ", paste(missing, collapse = ", "))
  matrix <- as.matrix(new_features[, required, drop = FALSE])
  scaled <- scale(matrix, center = metadata$center, scale = metadata$scale)
  scores <- scaled %*% metadata$rotation
  prediction <- stats::predict(method_result$model, newdata = scores)
  list(
    assignments = as.character(prediction$classification),
    diagnostics = list(
      posterior = prediction$z,
      maximum_posterior = apply(prediction$z, 1L, max)
    )
  )
}

# Assign each held-out participant to the affinity-propagation exemplar with
# greatest cosine similarity and inherit that exemplar's training label.
predict_affinity_propagation_clusters <- function(method_result, new_features) {
  metadata <- method_result$prediction_metadata
  exemplar_features <- metadata$exemplar_features
  similarity <- cross_cosine_similarity(new_features, exemplar_features)
  closest <- max.col(similarity, ties.method = "first")
  exemplar_ids <- rownames(exemplar_features)[closest]
  training_labels <- method_result$assignments
  prediction <- as.character(training_labels[exemplar_ids])
  list(assignments = prediction,
       diagnostics = list(exemplar_similarity = similarity,
                          exemplar_id = exemplar_ids))
}

# Use HDBSCAN's supported fitted-model prediction and retain predicted noise as
# NA so it is unavailable, rather than an ordinary consensus cluster.
predict_hdbscan_clusters <- function(method_result, new_features) {
  training <- method_result$prediction_metadata$training_features
  prediction <- stats::predict(method_result$model, newdata = as.matrix(new_features),
                               data = training)
  raw <- as.integer(prediction)
  labels <- raw
  labels[labels == 0L] <- NA_integer_
  list(assignments = labels,
       diagnostics = list(raw_assignments = raw, predicted_noise = raw == 0L))
}

# Dispatch every fitted method-specific adapter in the common training-derived
# feature space and return aligned assignments plus uncertainty diagnostics.
predict_base_clusters <- function(base_result, new_distance_result,
                                  methods = names(base_result$methods)) {
  if (!inherits(base_result, "base_clustering_result")) {
    stop("base_result must be a base_clustering_result.")
  }
  new_features <- as.matrix(new_distance_result$features)
  ids <- new_distance_result$participant_ids
  if (is.null(ids)) ids <- rownames(new_features)
  ids <- as.character(ids)
  if (length(ids) != nrow(new_features) || anyNA(ids) || anyDuplicated(ids)) {
    stop("New participant IDs must uniquely identify feature rows.")
  }
  rownames(new_features) <- ids
  missing_methods <- setdiff(methods, names(base_result$methods))
  if (length(missing_methods)) stop("Base result lacks methods: ", paste(missing_methods, collapse = ", "))
  adapters <- list(
    hca = predict_hca_clusters,
    spherical_kmeans = predict_spherical_kmeans_clusters,
    gmm = predict_gmm_clusters,
    affinity_propagation = predict_affinity_propagation_clusters,
    hdbscan = predict_hdbscan_clusters
  )
  results <- vector("list", length(methods))
  names(results) <- methods
  for (method in methods) {
    fitted <- base_result$methods[[method]]
    if (!identical(fitted$status, "available")) {
      results[[method]] <- list(
        assignments = rep(NA_character_, length(ids)),
        diagnostics = list(reason = "training method unavailable")
      )
    } else {
      results[[method]] <- tryCatch(
        adapters[[method]](fitted, new_features),
        # Convert one adapter failure into unavailable predictions with an
        # explicit reason while preserving other method predictions.
        error = function(error) list(
          assignments = rep(NA_character_, length(ids)),
          diagnostics = list(reason = conditionMessage(error))
        )
      )
    }
  }
  # Align one method's predicted assignments to the validated held-out ID order.
  assignments <- as.data.frame(
    lapply(results, function(result) result$assignments),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  rownames(assignments) <- ids
  structure(list(
    participant_ids = ids,
    assignments = assignments,
    diagnostics = lapply(results, `[[`, "diagnostics"),
    methods = methods
  ), class = c("base_prediction_result", "list"))
}

# Extend method agreement from each held-out participant to the fixed training
# consensus partition, assign by greatest mean co-association, and retain
# rejection/uncertainty diagnostics without reclustering the combined sample.
predict_consensus_cluster <- function(consensus_result, base_result,
                                      base_prediction, min_methods = 2L,
                                      min_score = 0, min_margin = 0) {
  if (!inherits(consensus_result, "consensus_clustering_result")) {
    stop("consensus_result must be a consensus_clustering_result.")
  }
  if (!inherits(base_prediction, "base_prediction_result")) {
    stop("base_prediction must be a base_prediction_result.")
  }
  methods <- consensus_result$active_methods
  training <- validate_assignments(base_result$assignments, methods)
  predicted <- as.data.frame(base_prediction$assignments,
                             stringsAsFactors = FALSE, check.names = FALSE)
  missing_methods <- setdiff(methods, names(predicted))
  if (length(missing_methods)) {
    stop("Predictions lack consensus methods: ", paste(missing_methods, collapse = ", "))
  }
  predicted <- predicted[, methods, drop = FALSE]
  # Normalize one predicted method column to character labels and canonical NA.
  predicted[] <- lapply(predicted, function(values) {
    values <- as.character(values)
    values[values == ""] <- NA_character_
    values
  })
  if (is.null(rownames(predicted)) || anyNA(rownames(predicted)) ||
      anyDuplicated(rownames(predicted))) {
    stop("Prediction rows require unique, non-missing participant IDs.")
  }
  test_ids <- rownames(predicted)
  training_ids <- rownames(training)
  agreement <- matrix(NA_real_, nrow(predicted), nrow(training),
                      dimnames = list(test_ids, training_ids))
  available_count <- matrix(0L, nrow(predicted), nrow(training),
                            dimnames = dimnames(agreement))
  co_count <- available_count
  for (method in methods) {
    test_labels <- predicted[[method]]
    train_labels <- training[[method]]
    both <- outer(!is.na(test_labels), !is.na(train_labels), "&")
    same <- outer(test_labels, train_labels, "==")
    same[!both | is.na(same)] <- FALSE
    available_count <- available_count + both
    co_count <- co_count + same
  }
  observed <- available_count > 0L
  agreement[observed] <- co_count[observed] / available_count[observed]
  consensus_labels <- as.character(consensus_result$labels[training_ids])
  cluster_levels <- sort(unique(consensus_labels[!is.na(consensus_labels)]))
  # Score one fixed training consensus cluster by mean held-out co-association.
  scores <- sapply(cluster_levels, function(label) {
    rowMeans(agreement[, consensus_labels == label, drop = FALSE], na.rm = TRUE)
  })
  if (is.null(dim(scores))) scores <- matrix(scores, ncol = length(cluster_levels))
  colnames(scores) <- cluster_levels
  rownames(scores) <- test_ids
  best_index <- max.col(scores, ties.method = "first")
  best_score <- scores[cbind(seq_len(nrow(scores)), best_index)]
  # Extract the runner-up cluster score for one held-out participant so the
  # assignment margin can support rejection diagnostics.
  second_score <- vapply(seq_len(nrow(scores)), function(index) {
    ordered <- sort(scores[index, ], decreasing = TRUE)
    if (length(ordered) >= 2L) ordered[2L] else NA_real_
  }, numeric(1))
  margin <- best_score - second_score
  contributing <- rowSums(!is.na(predicted))
  rejected <- contributing < min_methods | !is.finite(best_score) |
    best_score < min_score | (!is.na(margin) & margin < min_margin)
  labels <- as.integer(factor(cluster_levels[best_index], levels = cluster_levels))
  labels[rejected] <- NA_integer_
  names(labels) <- test_ids
  reason <- rep(NA_character_, length(test_ids))
  reason[contributing < min_methods] <- "insufficient_methods"
  reason[!is.finite(best_score)] <- "no_consensus_evidence"
  reason[is.na(reason) & best_score < min_score] <- "low_consensus_score"
  reason[is.na(reason) & !is.na(margin) & margin < min_margin] <- "low_margin"
  diagnostics <- data.frame(
    CIMBI = test_ids,
    best_score = best_score,
    second_score = second_score,
    margin = margin,
    contributing_methods = contributing,
    hdbscan_noise = if ("hdbscan" %in% names(predicted)) is.na(predicted$hdbscan) else NA,
    gmm_maximum_posterior = if (!is.null(base_prediction$diagnostics$gmm$maximum_posterior)) {
      base_prediction$diagnostics$gmm$maximum_posterior
    } else NA_real_,
    rejected = rejected,
    reason = reason,
    stringsAsFactors = FALSE
  )
  structure(list(
    participant_ids = test_ids,
    labels = labels,
    native_label = stats::setNames(cluster_levels[best_index], test_ids),
    cross_coassociation = agreement,
    available_methods = available_count,
    cluster_scores = scores,
    diagnostics = diagnostics,
    thresholds = list(min_methods = min_methods, min_score = min_score,
                      min_margin = min_margin)
  ), class = c("consensus_prediction_result", "list"))
}
