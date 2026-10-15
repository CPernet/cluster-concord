# Function-based base clustering extracted from the legacy analysis scripts.

# Return the five canonical base-method identifiers used consistently in
# configuration, assignments, consensus inputs, and output tables.
base_method_names <- function() {
  c("hca", "spherical_kmeans", "gmm", "affinity_propagation", "hdbscan")
}

# Build and validate explicit base-clustering defaults. SPEC mode selects UPGMA
# for base HCA; all algorithm tuning values are persisted for auditability.
default_base_config <- function(
  mode = c("compatibility", "spec"),
  methods = base_method_names(),
  k_candidates = 2:6,
  min_cluster_size = 5L,
  skmeans_restarts = 50L,
  ap_grid_size = 300L,
  hdbscan_min_pts = 10L,
  gmm_variance = 0.95,
  gmm_pca_cap = 20L,
  optional_failure = c("unavailable", "error")
) {
  mode <- match.arg(mode)
  optional_failure <- match.arg(optional_failure)
  methods <- unique(as.character(methods))
  unknown <- setdiff(methods, base_method_names())
  if (length(unknown)) stop("Unknown clustering methods: ", paste(unknown, collapse = ", "))
  k_candidates <- sort(unique(as.integer(k_candidates)))
  if (!length(k_candidates) || anyNA(k_candidates) || any(k_candidates < 2L)) {
    stop("k_candidates must contain integers of at least 2.", call. = FALSE)
  }
  if (min_cluster_size < 1L) stop("min_cluster_size must be positive.", call. = FALSE)
  list(
    mode = mode,
    methods = methods,
    k_candidates = k_candidates,
    min_cluster_size = as.integer(min_cluster_size),
    skmeans_restarts = as.integer(skmeans_restarts),
    ap_grid_size = as.integer(ap_grid_size),
    hdbscan_min_pts = as.integer(hdbscan_min_pts),
    gmm_variance = as.numeric(gmm_variance),
    gmm_pca_cap = as.integer(gmm_pca_cap),
    hca_linkage = if (identical(mode, "spec")) "average" else "ward.D2",
    optional_failure = optional_failure
  )
}

# Compute mean silhouette safely for candidate-partition selection, returning
# NA when a valid multi-cluster silhouette cannot be formed.
safe_silhouette_mean <- function(labels, distance) {
  if (length(unique(labels[!is.na(labels)])) < 2L) return(NA_real_)
  result <- tryCatch(
    cluster::silhouette(as.integer(labels), distance),
    # Convert numerical or structural silhouette failures into an unavailable
    # score so other candidate partitions remain evaluable.
    error = function(e) NULL
  )
  if (is.null(result)) return(NA_real_)
  mean(result[, "sil_width"], na.rm = TRUE)
}

# Validate the shared 45-feature/similarity/distance contract and participant
# ordering before any base algorithm is fitted.
validate_distance_result <- function(distance_result) {
  required <- c("features", "similarity", "distance")
  missing <- setdiff(required, names(distance_result))
  if (length(missing)) stop("distance_result is missing: ", paste(missing, collapse = ", "))
  features <- as.matrix(distance_result$features)
  similarity <- as.matrix(distance_result$similarity)
  distance <- distance_result$distance
  if (!inherits(distance, "dist")) stop("distance_result$distance must inherit from dist.")
  ids <- distance_result$participant_ids
  if (is.null(ids)) ids <- rownames(features)
  ids <- as.character(ids)
  if (length(ids) != nrow(features) || anyNA(ids) || anyDuplicated(ids)) {
    stop("Participant IDs must be unique, non-missing, and match feature rows.")
  }
  if (is.null(rownames(features)) || anyNA(rownames(features)) ||
      anyDuplicated(rownames(features))) {
    stop("Feature row names must be unique and non-missing.")
  }
  if (is.null(rownames(features)) || !identical(rownames(features), ids)) {
    stop("Feature rows must be ordered by participant IDs.")
  }
  if (!identical(dim(similarity), c(length(ids), length(ids)))) {
    stop("Similarity dimensions must match participant IDs.")
  }
  if (!is.null(rownames(similarity)) && !identical(rownames(similarity), ids)) {
    stop("Similarity ordering must match participant IDs.")
  }
  if (attr(distance, "Size") != length(ids)) stop("Distance size must match participant IDs.")
  if (anyNA(features) || any(!is.finite(features))) stop("Features must be finite and complete.")
  if (anyNA(similarity) || any(!is.finite(similarity))) stop("Similarity must be finite and complete.")
  list(features = features, similarity = similarity, distance = distance, ids = ids)
}

# Construct a typed unavailable-method result so one optional algorithm failure
# does not corrupt assignments or diagnostics from the other methods.
unavailable_method_result <- function(method, ids, reason) {
  structure(list(
    method = method,
    status = "unavailable",
    assignments = stats::setNames(rep(NA_integer_, length(ids)), ids),
    selected_k = NA_integer_,
    model = NULL,
    prediction_metadata = NULL,
    diagnostics = list(),
    warnings = as.character(reason)
  ), class = c("base_method_result", "list"))
}

# Construct the common successful base-method result, including fitted model,
# prediction metadata, selected k, diagnostics, and optional raw assignments.
available_method_result <- function(method, assignments, selected_k, model,
                                    prediction_metadata = NULL,
                                    diagnostics = list(), warnings = character(),
                                    raw_assignments = NULL) {
  out <- list(
    method = method,
    status = "available",
    assignments = assignments,
    selected_k = as.integer(selected_k),
    model = model,
    prediction_metadata = prediction_metadata,
    diagnostics = diagnostics,
    warnings = warnings
  )
  if (!is.null(raw_assignments)) out$raw_assignments <- raw_assignments
  structure(out, class = c("base_method_result", "list"))
}

# Select a feasible dendrogram cut by mean silhouette, first enforcing the
# configured minimum cluster size and then using a documented relaxed fallback.
select_partition_k <- function(tree, distance, candidates, min_cluster_size) {
  candidates <- candidates[candidates < attr(distance, "Size")]
  if (!length(candidates)) stop("No feasible candidate k for this sample size.")
  # Score every feasible k under either strict or relaxed size constraints.
  evaluate <- function(enforce_size) {
    # Cut and score one candidate partition without changing the fitted tree.
    scores <- vapply(candidates, function(k) {
      labels <- stats::cutree(tree, k = k)
      if (enforce_size && min(table(labels)) < min_cluster_size) return(NA_real_)
      safe_silhouette_mean(labels, distance)
    }, numeric(1))
    scores
  }
  scores <- evaluate(TRUE)
  relaxed <- FALSE
  if (all(!is.finite(scores))) {
    scores <- evaluate(FALSE)
    relaxed <- TRUE
  }
  if (all(!is.finite(scores))) {
    selected <- min(candidates)
  } else {
    selected <- candidates[which.max(scores)]
  }
  list(k = selected, diagnostics = data.frame(
    k = candidates,
    mean_silhouette = scores,
    selected = candidates == selected,
    size_constraint_relaxed = relaxed
  ))
}

# Fit SPEC base hierarchical clustering in cosine-distance space. UPGMA is used
# in SPEC mode; the fitted tree and training members support held-out prediction.
cluster_hca <- function(distance_result, config, seed = 1L) {
  validated <- validate_distance_result(distance_result)
  tree <- stats::hclust(validated$distance, method = config$hca_linkage)
  selection <- select_partition_k(
    tree, validated$distance, config$k_candidates, config$min_cluster_size
  )
  assignments <- stats::setNames(stats::cutree(tree, k = selection$k), validated$ids)
  available_method_result(
    "hca", assignments, selection$k, tree,
    prediction_metadata = list(
      linkage = config$hca_linkage,
      training_features = validated$features,
      training_assignments = assignments,
      distance = "cosine"
    ),
    diagnostics = list(k_selection = selection$diagnostics)
  )
}

# Run deterministic repeated spherical-k-means starts for one k and retain the
# successful fit with the strongest reported objective.
fit_skmeans_repeated <- function(features, k, restarts, seed) {
  fits <- vector("list", restarts)
  objective <- rep(NA_real_, restarts)
  for (index in seq_len(restarts)) {
    set.seed(seed + index)
    fits[[index]] <- tryCatch(
      skmeans::skmeans(features, k = k, method = "genetic"),
      # Mark a failed stochastic start unavailable while preserving other starts.
      error = function(e) NULL
    )
    if (!is.null(fits[[index]])) {
      objective[index] <- tryCatch(
        tail(fits[[index]]$value, 1L),
        # Some package versions omit a usable objective; silhouette selection
        # remains available even when this restart tie-break is missing.
        error = function(e) NA_real_
      )
    }
  }
  valid <- which(!vapply(fits, is.null, logical(1)))
  if (!length(valid)) return(NULL)
  if (all(!is.finite(objective[valid]))) return(fits[[valid[1L]]])
  fits[[which.max(objective)]]
}

# Fit spherical k-means to normalized 45-feature rows, choose k by cosine-space
# silhouette, and retain normalized centers for held-out assignment.
cluster_spherical_kmeans <- function(distance_result, config, seed = 1L) {
  validated <- validate_distance_result(distance_result)
  norms <- sqrt(rowSums(validated$features^2))
  if (any(!is.finite(norms) | norms == 0)) stop("Spherical k-means requires non-zero feature rows.")
  features <- validated$features / norms
  cosine_similarity <- features %*% t(features)
  cosine_similarity[cosine_similarity > 1] <- 1
  cosine_similarity[cosine_similarity < -1] <- -1
  cosine_distance <- stats::as.dist(1 - cosine_similarity)
  candidates <- config$k_candidates[config$k_candidates < nrow(features)]
  if (!length(candidates)) stop("No feasible candidate k for spherical k-means.")
  fits <- vector("list", length(candidates))
  scores <- rep(NA_real_, length(candidates))
  for (index in seq_along(candidates)) {
    fit <- fit_skmeans_repeated(features, candidates[index], config$skmeans_restarts,
                                seed + 1000L * index)
    fits[[index]] <- fit
    if (is.null(fit) || min(table(fit$cluster)) < config$min_cluster_size) next
    scores[index] <- safe_silhouette_mean(fit$cluster, cosine_distance)
  }
  if (all(!is.finite(scores))) {
    for (index in seq_along(candidates)) {
      if (!is.null(fits[[index]])) scores[index] <- safe_silhouette_mean(
        fits[[index]]$cluster, cosine_distance
      )
    }
  }
  if (all(!is.finite(scores))) stop("All spherical k-means candidates failed.")
  selected_index <- which.max(scores)
  fit <- fits[[selected_index]]
  assignments <- stats::setNames(as.integer(fit$cluster), validated$ids)
  # Compute and renormalize the fitted cosine centroid for one selected cluster.
  centers <- do.call(rbind, lapply(sort(unique(assignments)), function(label) {
    center <- colMeans(features[assignments == label, , drop = FALSE])
    center / sqrt(sum(center^2))
  }))
  rownames(centers) <- sort(unique(assignments))
  available_method_result(
    "spherical_kmeans", assignments, candidates[selected_index], fit,
    prediction_metadata = list(centers = centers, distance = "cosine"),
    diagnostics = list(k_selection = data.frame(
      k = candidates, mean_silhouette = scores,
      selected = seq_along(candidates) == selected_index
    ))
  )
}

# Fit training-only scaling, PCA, and mclust model selection, retaining every
# transform required for posterior prediction of held-out participants.
cluster_gmm <- function(distance_result, config, seed = 1L) {
  validated <- validate_distance_result(distance_result)
  # mclust 6.1.3 under R 4.6 resolves several internal helpers from the
  # attached package environment rather than only from its namespace.
  suppressPackageStartupMessages(
    require("mclust", character.only = TRUE, quietly = TRUE)
  )
  features <- validated$features
  keep <- apply(features, 2L, stats::sd) > 0
  if (!any(keep)) stop("GMM has no non-constant features.")
  features <- features[, keep, drop = FALSE]
  center <- colMeans(features)
  scale_values <- apply(features, 2L, stats::sd)
  scaled <- scale(features, center = center, scale = scale_values)
  pca <- stats::prcomp(scaled, center = FALSE, scale. = FALSE)
  cumulative <- summary(pca)$importance[3L, ]
  components <- which(cumulative >= config$gmm_variance)[1L]
  if (is.na(components)) components <- ncol(pca$x)
  components <- min(max(1L, components), config$gmm_pca_cap)
  scores <- pca$x[, seq_len(components), drop = FALSE]
  candidates <- config$k_candidates[config$k_candidates < nrow(scores)]
  if (!length(candidates)) stop("No feasible component count for GMM.")
  set.seed(seed)
  fit <- mclust::Mclust(
    scores,
    G = candidates,
    modelNames = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE"),
    verbose = FALSE
  )
  assignments <- stats::setNames(as.integer(fit$classification), validated$ids)
  available_method_result(
    "gmm", assignments, fit$G, fit,
    prediction_metadata = list(
      kept_features = colnames(features), center = center, scale = scale_values,
      rotation = pca$rotation[, seq_len(components), drop = FALSE],
      max_posterior = apply(fit$z, 1L, max)
    ),
    diagnostics = list(model_name = fit$modelName, components = components,
                       uncertainty = fit$uncertainty)
  )
}

# Fit affinity propagation across an explicit preference grid, select a bounded
# cluster count, and retain exemplars for greatest-similarity prediction.
cluster_affinity_propagation <- function(distance_result, config, seed = 1L) {
  validated <- validate_distance_result(distance_result)
  similarity <- (validated$similarity + t(validated$similarity)) / 2
  diag(similarity) <- 1
  values <- similarity[upper.tri(similarity)]
  span <- diff(range(values))
  if (!is.finite(span) || span < .Machine$double.eps) span <- 1
  low_n <- max(5L, ceiling(config$ap_grid_size * 2 / 3))
  mid_n <- max(5L, config$ap_grid_size - low_n)
  preferences <- sort(unique(c(
    seq(min(values) - 20 * span, min(values), length.out = low_n),
    stats::quantile(values, probs = seq(0.001, 0.5, length.out = mid_n), names = FALSE)
  )))
  models <- vector("list", length(preferences))
  counts <- rep(NA_integer_, length(preferences))
  set.seed(seed)
  for (index in seq_along(preferences)) {
    models[[index]] <- tryCatch(
      suppressWarnings(apcluster::apcluster(
        s = similarity, p = preferences[index], maxits = 3000,
        convits = 200, lam = 0.9
      )),
      # Preserve the preference grid when one AP candidate fails to converge.
      error = function(e) NULL
    )
    if (!is.null(models[[index]])) counts[index] <- length(models[[index]]@clusters)
  }
  valid <- which(is.finite(counts))
  if (!length(valid)) stop("All affinity-propagation candidates failed.")
  in_range <- valid[counts[valid] %in% config$k_candidates]
  choices <- if (length(in_range)) in_range else valid
  target <- stats::median(config$k_candidates)
  selected <- choices[which.min(abs(counts[choices] - target))]
  fit <- models[[selected]]
  assignments <- rep(NA_integer_, length(validated$ids))
  for (cluster_index in seq_along(fit@clusters)) {
    assignments[fit@clusters[[cluster_index]]] <- cluster_index
  }
  assignments <- stats::setNames(assignments, validated$ids)
  exemplars <- as.integer(fit@exemplars)
  available_method_result(
    "affinity_propagation", assignments, counts[selected], fit,
    prediction_metadata = list(
      exemplar_indices = exemplars,
      exemplar_ids = validated$ids[exemplars],
      exemplar_features = validated$features[exemplars, , drop = FALSE],
      preference = preferences[selected], distance = "cosine"
    ),
    diagnostics = list(preferences = preferences, cluster_counts = counts)
  )
}

# Fit density-based HDBSCAN in the common feature space. Noise label zero is
# converted to NA so consensus denominators treat it as unavailable evidence.
cluster_hdbscan <- function(distance_result, config, seed = 1L) {
  validated <- validate_distance_result(distance_result)
  set.seed(seed)
  fit <- dbscan::hdbscan(validated$features, minPts = config$hdbscan_min_pts)
  raw <- stats::setNames(as.integer(fit$cluster), validated$ids)
  assignments <- raw
  assignments[assignments == 0L] <- NA_integer_
  available_method_result(
    "hdbscan", assignments, length(unique(stats::na.omit(assignments))), fit,
    prediction_metadata = list(
      training_features = validated$features,
      membership_probability = fit$membership_prob,
      min_pts = config$hdbscan_min_pts
    ),
    diagnostics = list(noise_fraction = mean(raw == 0L)),
    raw_assignments = raw
  )
}

# Dispatch all requested base algorithms with deterministic method-specific
# seeds and return one aligned, serializable result for consensus and prediction.
run_base_clustering <- function(distance_result, config = default_base_config(), seed = 1L) {
  validated <- validate_distance_result(distance_result)
  seed <- as.integer(seed)
  dispatch <- list(
    hca = cluster_hca,
    spherical_kmeans = cluster_spherical_kmeans,
    gmm = cluster_gmm,
    affinity_propagation = cluster_affinity_propagation,
    hdbscan = cluster_hdbscan
  )
  packages <- c(
    hca = "cluster", spherical_kmeans = "skmeans", gmm = "mclust",
    affinity_propagation = "apcluster", hdbscan = "dbscan"
  )
  results <- vector("list", length(config$methods))
  names(results) <- config$methods
  for (index in seq_along(config$methods)) {
    method <- config$methods[index]
    if (!requireNamespace(packages[[method]], quietly = TRUE)) {
      reason <- paste0("Required package is not installed: ", packages[[method]])
      if (identical(config$optional_failure, "error")) stop(reason)
      results[[method]] <- unavailable_method_result(method, validated$ids, reason)
      next
    }
    results[[method]] <- tryCatch(
      dispatch[[method]](distance_result, config, seed + 10000L * index),
      # Either propagate required-method failure or convert an optional failure
      # into a typed unavailable result without stopping unrelated algorithms.
      error = function(error) {
        if (identical(config$optional_failure, "error")) stop(error)
        unavailable_method_result(method, validated$ids, conditionMessage(error))
      }
    )
  }
  # Align each method's assignment vector to the validated participant order.
  assignments <- as.data.frame(
    lapply(results, function(result) unname(result$assignments)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  rownames(assignments) <- validated$ids
  # Extract selected cluster counts for fold and resampling summaries.
  cluster_counts <- vapply(results, function(result) result$selected_k, integer(1))
  warnings <- unlist(lapply(results, `[[`, "warnings"), use.names = FALSE)
  structure(list(
    participant_ids = validated$ids,
    features = validated$features,
    similarity = validated$similarity,
    distance = validated$distance,
    methods = results,
    assignments = assignments,
    cluster_counts = cluster_counts,
    diagnostics = lapply(results, `[[`, "diagnostics"),
    warnings = warnings,
    config = config,
    seed = seed
  ), class = c("base_clustering_result", "list"))
}
