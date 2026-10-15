# Label-invariant consensus clustering for named base-method assignments.

# Build explicit consensus defaults. SPEC mode requests all five methods,
# preserves no-evidence pairs as NA, and uses Ward D2 only on distance 1-C.
default_consensus_config <- function(
  mode = c("compatibility", "spec"),
  active_methods = NULL,
  k_candidates = 2:6,
  min_cluster_size = 5L,
  no_evidence = NULL,
  score_weights = c(
    silhouette = 0.30,
    separation = 0.30,
    within = 0.20,
    balance = 0.10,
    pac = 0.10
  )
) {
  mode <- match.arg(mode)
  if (is.null(active_methods)) {
    active_methods <- if (identical(mode, "compatibility")) {
      c("hca", "spherical_kmeans", "gmm")
    } else {
      base_method_names()
    }
  }
  active_methods <- unique(as.character(active_methods))
  unknown <- setdiff(active_methods, base_method_names())
  if (length(unknown)) stop("Unknown consensus methods: ", paste(unknown, collapse = ", "))
  k_candidates <- sort(unique(as.integer(k_candidates)))
  if (!length(k_candidates) || anyNA(k_candidates) || any(k_candidates < 2L)) {
    stop("Consensus k candidates must be integers of at least 2.")
  }
  if (is.null(no_evidence)) no_evidence <- if (identical(mode, "compatibility")) "zero" else "na"
  no_evidence <- match.arg(no_evidence, c("zero", "na", "error"))
  score_weights <- as.numeric(score_weights)
  names(score_weights) <- c("silhouette", "separation", "within", "balance", "pac")
  if (any(!is.finite(score_weights)) || any(score_weights < 0) || sum(score_weights) <= 0) {
    stop("Consensus score weights must be finite non-negative values with positive sum.")
  }
  score_weights <- score_weights / sum(score_weights)
  list(
    mode = mode,
    active_methods = active_methods,
    k_candidates = k_candidates,
    min_cluster_size = as.integer(min_cluster_size),
    no_evidence = no_evidence,
    score_weights = score_weights,
    linkage = "ward.D2"
  )
}

# Validate and normalize an ID-indexed method-assignment table while preserving
# NA values used for HDBSCAN noise and unavailable method assignments.
validate_assignments <- function(assignments, methods = names(assignments)) {
  assignments <- as.data.frame(assignments, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(assignments) < 2L) stop("At least two participants are required.")
  if (is.null(rownames(assignments)) || anyNA(rownames(assignments)) ||
      anyDuplicated(rownames(assignments))) {
    stop("Assignment rows require unique, non-missing participant IDs.")
  }
  missing <- setdiff(methods, names(assignments))
  if (length(missing)) stop("Missing assignment methods: ", paste(missing, collapse = ", "))
  assignments <- assignments[, methods, drop = FALSE]
  # Normalize one method column to character labels and canonical missing values.
  assignments[] <- lapply(assignments, function(values) {
    values <- as.character(values)
    values[values == ""] <- NA_character_
    values
  })
  assignments
}

# Compute label-invariant pairwise ARI between base methods using only
# participants assigned by both methods.
compute_method_ari <- function(assignments) {
  assignments <- validate_assignments(assignments)
  methods <- names(assignments)
  pairs <- expand.grid(method_a = methods, method_b = methods, stringsAsFactors = FALSE)
  pairs$n_available <- 0L
  pairs$ari <- NA_real_
  for (index in seq_len(nrow(pairs))) {
    first <- assignments[[pairs$method_a[index]]]
    second <- assignments[[pairs$method_b[index]]]
    available <- !is.na(first) & !is.na(second)
    pairs$n_available[index] <- sum(available)
    if (sum(available) >= 2L) {
      pairs$ari[index] <- mclust::adjustedRandIndex(first[available], second[available])
    }
  }
  pairs
}

# Build the evidence-accumulation co-association matrix. Each participant pair's
# denominator counts only methods with two available assignments.
build_coassociation <- function(assignments, no_evidence = c("na", "zero", "error")) {
  no_evidence <- match.arg(no_evidence)
  assignments <- validate_assignments(assignments)
  count <- matrix(0L, nrow(assignments), nrow(assignments),
                  dimnames = list(rownames(assignments), rownames(assignments)))
  available <- count
  for (method in names(assignments)) {
    labels <- assignments[[method]]
    present <- !is.na(labels)
    both <- outer(present, present, "&")
    same <- outer(labels, labels, "==")
    same[!both | is.na(same)] <- FALSE
    count <- count + same
    available <- available + both
  }
  coassociation <- matrix(NA_real_, nrow(count), ncol(count), dimnames = dimnames(count))
  observed <- available > 0L
  coassociation[observed] <- count[observed] / available[observed]
  diag(coassociation) <- 1
  unavailable <- !observed
  diag(unavailable) <- FALSE
  no_evidence_pairs <- sum(unavailable[upper.tri(unavailable)])
  if (no_evidence_pairs && identical(no_evidence, "error")) {
    stop(no_evidence_pairs, " participant pairs have no available method agreement.")
  }
  if (identical(no_evidence, "zero")) coassociation[unavailable] <- 0
  list(
    coassociation = coassociation,
    co_clustered = count,
    available = available,
    no_evidence_pairs = as.integer(no_evidence_pairs),
    no_evidence_policy = no_evidence
  )
}

# Rescale one consensus selection metric to [0,1], with deterministic handling
# of constant or completely unavailable metrics.
scale_consensus_metric <- function(values) {
  if (all(!is.finite(values))) return(rep(NA_real_, length(values)))
  limits <- range(values[is.finite(values)])
  if (diff(limits) < 1e-12) return(ifelse(is.finite(values), 0.5, NA_real_))
  (values - limits[1L]) / diff(limits)
}

# Calculate the prespecified k-selection metrics for each cut of the fixed
# Ward-D2 consensus tree: silhouette, agreement, separation, balance, and PAC.
consensus_partition_metrics <- function(tree, distance, coassociation, candidates,
                                        min_cluster_size) {
  # Cut and score one candidate consensus cluster count.
  do.call(rbind, lapply(candidates, function(k) {
    labels <- stats::cutree(tree, k = k)
    sizes <- table(labels)
    within <- numeric()
    between <- numeric()
    for (first in seq_len(length(labels) - 1L)) {
      for (second in (first + 1L):length(labels)) {
        if (labels[first] == labels[second]) {
          within <- c(within, coassociation[first, second])
        } else {
          between <- c(between, coassociation[first, second])
        }
      }
    }
    data.frame(
      k = k,
      min_size = min(sizes),
      max_size = max(sizes),
      valid_size = min(sizes) >= min_cluster_size,
      silhouette = safe_silhouette_mean(labels, distance),
      within = mean(within),
      between = mean(between),
      separation = mean(within) - mean(between),
      balance = min(sizes) / max(sizes),
      pac = mean(within > 0.1 & within < 0.9),
      stringsAsFactors = FALSE
    )
  }))
}

# Select consensus k through the configured weighted score, enforcing minimum
# cluster size first and recording any necessary relaxed fallback.
select_consensus_k <- function(tree, distance, coassociation, config) {
  candidates <- config$k_candidates[config$k_candidates < nrow(coassociation)]
  if (!length(candidates)) stop("No feasible consensus k for this sample size.")
  metrics <- consensus_partition_metrics(
    tree, distance, coassociation, candidates, config$min_cluster_size
  )
  eligible <- metrics$valid_size & apply(
    metrics[, c("silhouette", "within", "between", "separation", "balance", "pac")],
    # Require every primary selection metric to be finite for one candidate row.
    1L, function(row) all(is.finite(row))
  )
  relaxed <- FALSE
  if (!any(eligible)) {
    eligible <- apply(
      metrics[, c("silhouette", "within", "between", "separation", "balance", "pac")],
      # Under relaxed size constraints, retain only rows with complete metrics.
      1L, function(row) all(is.finite(row))
    )
    relaxed <- TRUE
  }
  if (!any(eligible)) stop("No consensus candidate has finite selection metrics.")
  scored <- metrics[eligible, , drop = FALSE]
  scored$silhouette_scaled <- scale_consensus_metric(scored$silhouette)
  scored$separation_scaled <- scale_consensus_metric(scored$separation)
  scored$within_scaled <- scale_consensus_metric(scored$within)
  scored$balance_scaled <- scale_consensus_metric(scored$balance)
  scored$pac_scaled <- scale_consensus_metric(-scored$pac)
  weights <- config$score_weights
  scored$score <-
    weights[["silhouette"]] * scored$silhouette_scaled +
    weights[["separation"]] * scored$separation_scaled +
    weights[["within"]] * scored$within_scaled +
    weights[["balance"]] * scored$balance_scaled +
    weights[["pac"]] * scored$pac_scaled
  selected <- scored$k[which.max(scored$score)]
  metrics$selected <- metrics$k == selected
  list(k = as.integer(selected), metrics = metrics, scored = scored,
       size_constraint_relaxed = relaxed)
}

# Convert a finite complete co-association matrix to distance 1-C, fit the fixed
# Ward-D2 tree, select k, and return the training consensus partition.
fit_consensus <- function(coassociation, config) {
  coassociation <- as.matrix(coassociation)
  if (nrow(coassociation) < 3L || nrow(coassociation) != ncol(coassociation)) {
    stop("Consensus requires a square matrix for at least three participants.")
  }
  if (anyNA(coassociation) || any(!is.finite(coassociation))) {
    stop("A finite complete consensus distance cannot be formed; inspect no-evidence pairs.")
  }
  if (!isSymmetric(coassociation, tol = 1e-10)) stop("Co-association must be symmetric.")
  if (any(coassociation < 0 | coassociation > 1)) stop("Co-association must lie in [0, 1].")
  distance <- stats::as.dist(1 - coassociation)
  tree <- stats::hclust(distance, method = config$linkage)
  selection <- select_consensus_k(tree, distance, coassociation, config)
  labels <- stats::setNames(stats::cutree(tree, k = selection$k), rownames(coassociation))
  list(distance = distance, tree = tree, labels = labels,
       selected_k = selection$k, selection = selection)
}

# Combine available requested base methods into the structured consensus result,
# retaining assignments, denominators, diagnostics, selected k, and omissions.
run_consensus_clustering <- function(base_result,
                                     config = default_consensus_config()) {
  if (!inherits(base_result, "base_clustering_result")) {
    stop("base_result must be a base_clustering_result.")
  }
  requested <- config$active_methods
  # Keep one requested method only when its fit succeeded and at least one
  # participant has a non-noise assignment.
  available <- requested[vapply(requested, function(method) {
    method %in% names(base_result$methods) &&
      identical(base_result$methods[[method]]$status, "available") &&
      any(!is.na(base_result$methods[[method]]$assignments))
  }, logical(1))]
  omitted <- setdiff(requested, available)
  if (length(available) < 2L) stop("At least two available base methods are required.")
  assignments <- validate_assignments(base_result$assignments, available)
  coassociation <- build_coassociation(assignments, config$no_evidence)
  fitted <- fit_consensus(coassociation$coassociation, config)
  warnings <- if (length(omitted)) {
    paste("Unavailable consensus methods omitted:", paste(omitted, collapse = ", "))
  } else character()
  structure(list(
    participant_ids = rownames(assignments),
    active_methods = available,
    omitted_methods = omitted,
    assignments = assignments,
    co_clustered = coassociation$co_clustered,
    available = coassociation$available,
    coassociation = coassociation$coassociation,
    no_evidence_pairs = coassociation$no_evidence_pairs,
    distance = fitted$distance,
    tree = fitted$tree,
    labels = fitted$labels,
    selected_k = fitted$selected_k,
    selection_diagnostics = fitted$selection,
    cross_method_ari = compute_method_ari(assignments),
    warnings = warnings,
    config = config
  ), class = c("consensus_clustering_result", "list"))
}
