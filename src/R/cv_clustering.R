# Cross-validation orchestration using the modular base and consensus APIs.

# Create deterministic, nearly equal folds using joint categorical strata and
# optional age bins. NULL balance_cols retains the behavioural defaults.
make_balanced_folds <- function(data, folds = 5L, sex_col,
                                medication_col, age_col, seed = 1L,
                                age_bin_width = 12, age_bin_origin = 18,
                                balance_cols = NULL) {
  folds <- as.integer(folds)
  if (folds < 2L || folds > nrow(data)) stop("folds must be between 2 and sample size.")
  categorical_cols <- if (is.null(balance_cols)) c(sex_col, medication_col) else balance_cols
  if (!is.character(categorical_cols) || anyNA(categorical_cols) ||
      anyDuplicated(categorical_cols)) stop("balance_cols must be unique column names.")
  required <- c(categorical_cols, age_col)
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Missing fold-balancing columns: ", paste(missing, collapse = ", "))
  # Normalize one categorical balancing variable and retain missingness as an
  # explicit stratum rather than dropping participants.
  clean_category <- function(values) {
    values <- as.character(values)
    values[is.na(values) | values == ""] <- "Missing"
    values
  }
  categories <- lapply(categorical_cols, function(column) clean_category(data[[column]]))
  names(categories) <- if (is.null(balance_cols)) c("sex", "medication") else categorical_cols
  if (!is.null(balance_cols) && any(names(categories) %in% c("fold", "n", "age_bin", "all_rows", "can_populate_every_fold"))) {
    stop("Rename balancing columns that conflict with diagnostic fields.")
  }
  age <- if (is.null(age_col)) rep(NA_real_, nrow(data)) else suppressWarnings(as.numeric(data[[age_col]]))
  if (length(age_bin_width) != 1L || !is.finite(age_bin_width) ||
      age_bin_width <= 0) stop("age_bin_width must be positive and finite.")
  if (length(age_bin_origin) != 1L || !is.finite(age_bin_origin)) {
    stop("age_bin_origin must be finite.")
  }
  finite <- is.finite(age)
  lower <- rep(NA_real_, length(age))
  lower[finite] <- age_bin_origin +
    floor((age[finite] - age_bin_origin) / age_bin_width) * age_bin_width
  age_bin <- rep("Missing", length(age))
  age_bin[finite] <- paste0("[", lower[finite], ",",
                            lower[finite] + age_bin_width, ")")
  breaks <- if (any(finite)) {
    seq(min(lower[finite]), max(lower[finite]) + age_bin_width,
        by = age_bin_width)
  } else numeric()
  if (!is.null(age_col)) categories$age_bin <- age_bin
  if (!length(categories)) categories$all_rows <- rep("All", nrow(data))
  strata <- do.call(interaction, c(categories, list(drop = TRUE, lex.order = TRUE)))
  n <- nrow(data)
  target <- rep(n %/% folds, folds)
  if (n %% folds) target[seq_len(n %% folds)] <- target[seq_len(n %% folds)] + 1L
  fold_id <- rep(NA_integer_, n)
  fold_sizes <- integer(folds)
  canonical_positions <- canonical_record_order(data, age_col = age_col)
  fold_seed <- analysis_stage_seed(seed, 200000L)
  set.seed(fold_seed)
  stratum_levels <- levels(strata)
  for (stratum in stratum_levels[sample.int(length(stratum_levels))]) {
    stratum_indices <- canonical_positions[strata[canonical_positions] == stratum]
    indices <- stratum_indices[sample.int(length(stratum_indices))]
    for (index in indices) {
      eligible <- which(fold_sizes < target)
      smallest <- min(fold_sizes[eligible])
      choices <- eligible[fold_sizes[eligible] == smallest]
      chosen <- choices[sample.int(length(choices), 1L)]
      fold_id[index] <- chosen
      fold_sizes[chosen] <- fold_sizes[chosen] + 1L
    }
  }
  balance_rows <- data.frame(fold = fold_id, categories, check.names = FALSE,
                             stringsAsFactors = FALSE)
  balance_long <- stats::aggregate(
    rep(1L, nrow(balance_rows)), by = balance_rows, FUN = sum
  )
  names(balance_long)[names(balance_long) == "x"] <- "n"
  balance_long <- balance_long[, c("fold", "n", names(categories)), drop = FALSE]
  stratum_counts <- stats::aggregate(
    rep(1L, n), by = categories,
    FUN = sum
  )
  names(stratum_counts)[names(stratum_counts) == "x"] <- "n"
  stratum_counts$can_populate_every_fold <- stratum_counts$n >= folds
  list(
    fold_id = fold_id,
    target_sizes = target,
    fold_sizes = fold_sizes,
    balance_long = balance_long,
    age_breaks = breaks,
    age_binning = list(strategy = if (is.null(age_col)) "disabled" else "fixed_width", width = age_bin_width,
                       origin = age_bin_origin, intervals = "left_closed_right_open"),
    age_band_counts = as.data.frame(table(age_bin), stringsAsFactors = FALSE),
    stratum_counts = stratum_counts,
    sparse_strata = stratum_counts[stratum_counts$n < folds, , drop = FALSE],
    fold_seed = fold_seed,
    seed = as.integer(seed)
  )
}

# Align arbitrary candidate labels to a reference space through padded
# rectangular Hungarian matching, preserving unmatched cluster labels.
align_fold_labels_legacy <- function(reference, candidate) {
  reference <- stats::setNames(as.character(reference), names(reference))
  candidate <- stats::setNames(as.character(candidate), names(candidate))
  overlap <- intersect(names(reference), names(candidate))
  ref_levels <- sort(unique(reference[!is.na(reference)]))
  candidate_levels <- sort(unique(candidate[!is.na(candidate)]))
  usable <- overlap[!is.na(reference[overlap]) & !is.na(candidate[overlap])]
  if (!length(usable) || !length(ref_levels) || !length(candidate_levels)) {
    return(list(labels = rep(NA_character_, length(candidate)), mapping = character(),
                overlap_n = length(usable)))
  }
  contingency <- matrix(0L, length(ref_levels), length(candidate_levels),
                        dimnames = list(ref_levels, candidate_levels))
  for (id in usable) {
    contingency[reference[id], candidate[id]] <- contingency[reference[id], candidate[id]] + 1L
  }
  size <- max(dim(contingency))
  padded <- matrix(0, size, size)
  padded[seq_len(nrow(contingency)), seq_len(ncol(contingency))] <- contingency
  assignment <- as.integer(clue::solve_LSAP(max(padded) - padded))
  mapping <- stats::setNames(rep(NA_character_, length(candidate_levels)), candidate_levels)
  for (ref_index in seq_len(length(ref_levels))) {
    candidate_index <- assignment[ref_index]
    if (candidate_index <= length(candidate_levels)) {
      mapping[candidate_index] <- ref_levels[ref_index]
    }
  }
  unmatched <- which(is.na(mapping))
  if (length(unmatched)) {
    mapping[unmatched] <- paste0("unmatched_", candidate_levels[unmatched])
  }
  labels <- unname(mapping[candidate])
  names(labels) <- names(candidate)
  list(labels = labels, mapping = mapping, overlap_n = length(usable),
       contingency = contingency)
}

# Historical reproduction helper: fit raw-predictor scaling on train and apply
# it unchanged to test for the compatibility centroid rule.
legacy_prepare_predictors <- function(train, test, variables) {
  train_matrix <- as.matrix(train[, variables, drop = FALSE])
  test_matrix <- as.matrix(test[, variables, drop = FALSE])
  storage.mode(train_matrix) <- "double"
  storage.mode(test_matrix) <- "double"
  medians <- apply(train_matrix, 2L, stats::median, na.rm = TRUE)
  medians[!is.finite(medians)] <- 0
  for (column in seq_along(variables)) {
    train_matrix[is.na(train_matrix[, column]), column] <- medians[column]
    test_matrix[is.na(test_matrix[, column]), column] <- medians[column]
  }
  center <- colMeans(train_matrix)
  scale_values <- apply(train_matrix, 2L, stats::sd)
  scale_values[!is.finite(scale_values) | scale_values == 0] <- 1
  list(
    train = scale(train_matrix, center = center, scale = scale_values),
    test = scale(test_matrix, center = center, scale = scale_values),
    center = center, scale = scale_values, medians = medians
  )
}

# Historical reproduction helper: assign test rows to nearest Euclidean
# training centroids and return distance margins.
legacy_centroid_prediction <- function(train_matrix, test_matrix, labels) {
  labels <- as.character(labels)
  valid <- !is.na(labels)
  levels <- sort(unique(labels[valid]))
  # Compute one available training cluster's raw-variable centroid.
  centers <- do.call(rbind, lapply(levels, function(label) {
    colMeans(train_matrix[valid & labels == label, , drop = FALSE])
  }))
  rownames(centers) <- levels
  # Compute distances from one held-out row to every training centroid.
  distances <- t(apply(test_matrix, 1L, function(row) {
    # Score one candidate centroid for the current held-out row.
    apply(centers, 1L, function(center) sum((row - center)^2))
  }))
  if (is.null(dim(distances))) distances <- matrix(distances, nrow = nrow(test_matrix))
  predictions <- levels[max.col(-distances, ties.method = "first")]
  ordered <- t(apply(distances, 1L, sort))
  margin <- if (ncol(distances) >= 2L) ordered[, 2L] - ordered[, 1L] else rep(NA_real_, nrow(test_matrix))
  list(prediction = predictions, margin = margin, centers = centers, distances = distances)
}

# Fit one outer fold with training-only preprocessing, all requested base
# methods, consensus, held-out adapters, and fixed-partition consensus extension.
run_cv_fold <- function(data, train_idx, test_idx, fold, distance_builder,
                        predictor_variables,
                        base_config = default_base_config(),
                        consensus_config = default_consensus_config(),
                        prediction_mode = c("legacy_raw_centroid", "feature_consensus"),
                        seed = 1L, feature_pipeline = NULL,
                        age_col = "Age at neuropsych test") {
  prediction_mode <- match.arg(prediction_mode)
  train_idx <- train_idx[canonical_record_order(data[train_idx, , drop = FALSE],
                                              age_col = age_col)]
  test_idx <- test_idx[canonical_record_order(data[test_idx, , drop = FALSE],
                                            age_col = age_col)]
  train <- data[train_idx, , drop = FALSE]
  test <- data[test_idx, , drop = FALSE]
  feature_fit <- NULL
  if (identical(prediction_mode, "feature_consensus")) {
    if (is.null(feature_pipeline)) {
      stop("feature_consensus prediction requires feature_pipeline(train, test).")
    }
    feature_fit <- feature_pipeline(train, test)
    if (!all(c("train", "test") %in% names(feature_fit))) {
      stop("feature_pipeline must return named train and test distance results.")
    }
    distance_result <- feature_fit$train
  } else {
    distance_result <- distance_builder(train)
  }
  base <- run_base_clustering(distance_result, base_config, seed = seed)
  consensus <- run_consensus_clustering(base, consensus_config)
  if (identical(prediction_mode, "legacy_raw_centroid")) {
    prepared <- legacy_prepare_predictors(train, test, predictor_variables)
    prediction <- legacy_centroid_prediction(
      prepared$train, prepared$test, consensus$labels
    )
    consensus_prediction <- stats::setNames(prediction$prediction, test$CIMBI)
    test_assignments <- data.frame(
      consensus = prediction$prediction,
      row.names = as.character(test$CIMBI), stringsAsFactors = FALSE
    )
    prediction_diagnostics <- list(
      margin = stats::setNames(prediction$margin, test$CIMBI),
      preprocessing = prepared,
      centers = prediction$centers
    )
  } else {
    base_prediction <- predict_base_clusters(base, feature_fit$test)
    consensus_prediction_result <- predict_consensus_cluster(
      consensus, base, base_prediction,
      min_methods = max(2L, ceiling(length(consensus$active_methods) / 2L)),
      min_score = 0, min_margin = 0
    )
    test_assignments <- base_prediction$assignments
    test_assignments$consensus <- unname(consensus_prediction_result$native_label)
    consensus_prediction <- stats::setNames(
      as.character(consensus_prediction_result$native_label), test$CIMBI
    )
    prediction_diagnostics <- list(
      base = base_prediction$diagnostics,
      consensus = consensus_prediction_result$diagnostics,
      cross_coassociation = consensus_prediction_result$cross_coassociation,
      preprocessor = feature_fit$preprocessor
    )
  }
  training_assignments <- base$assignments
  rownames(training_assignments) <- as.character(train$CIMBI)
  consensus_training <- stats::setNames(as.character(consensus$labels), train$CIMBI)
  all_training_assignments <- training_assignments
  all_training_assignments$consensus <- unname(consensus_training)
  list(
    fold = as.integer(fold),
    train_idx = train_idx,
    test_idx = test_idx,
    train_ids = as.character(train$CIMBI),
    test_ids = as.character(test$CIMBI),
    base = base,
    consensus = consensus,
    training_assignments = training_assignments,
    all_training_assignments = all_training_assignments,
    test_assignments = test_assignments,
    consensus_training = consensus_training,
    consensus_prediction = consensus_prediction,
    prediction_diagnostics = prediction_diagnostics,
    cluster_counts = c(base$cluster_counts, consensus = consensus$selected_k),
    config = list(base = base_config, consensus = consensus_config,
                  prediction_mode = prediction_mode),
    seed = as.integer(seed)
  )
}

# Select the common label-space reference fold by modal cluster count and a
# deterministic lowest-index tie-break.
select_reference_fold <- function(cluster_counts) {
  cluster_counts <- as.integer(cluster_counts)
  valid <- which(is.finite(cluster_counts))
  if (!length(valid)) return(NA_integer_)
  frequency <- table(cluster_counts[valid])
  modal_frequency <- max(frequency)
  modal_counts <- as.integer(names(frequency)[frequency == modal_frequency])
  selected_count <- min(modal_counts)
  valid[which(cluster_counts[valid] == selected_count)[1L]]
}

# Derive final subject labels directly from four training-fold assignments,
# requiring at least 3/4 agreement and returning NA for ties or instability.
derive_stable_labels <- function(label_matrix, threshold = 0.75,
                                 expected_assignments = 4L) {
  label_matrix <- as.matrix(label_matrix)
  # Summarize modal label, agreement fraction, and stability for one participant.
  rows <- lapply(seq_len(nrow(label_matrix)), function(index) {
    labels <- label_matrix[index, ]
    labels <- labels[!is.na(labels) & labels != ""]
    counts <- sort(table(labels), decreasing = TRUE)
    modal <- if (length(counts)) names(counts)[1L] else NA_character_
    proportion <- if (length(counts)) unname(counts[1L]) / length(labels) else NA_real_
    stable <- length(labels) >= expected_assignments && is.finite(proportion) &&
      proportion >= threshold
    c(label = modal, proportion = proportion, stable = stable,
      n_assignments = length(labels))
  })
  rows <- do.call(rbind, rows)
  stable <- as.logical(rows[, "stable"])
  stable_label <- rows[, "label"]
  stable_label[!stable] <- NA_character_
  list(
    modal_label = rows[, "label"],
    stable_label = unname(stable_label),
    stability_proportion = as.numeric(rows[, "proportion"]),
    stable = unname(stable),
    n_assignments = as.integer(rows[, "n_assignments"])
  )
}

# Align and summarize one method across folds, returning overlap ARI, held-out
# accuracy, coverage, confusion, balanced accuracy, and stable labels.
summarize_one_method_spec <- function(fold_fits, participant_ids, method,
                                      stability_threshold = 0.75) {
  # Extract this method's selected cluster count from one fold.
  counts <- vapply(fold_fits, function(fit) unname(fit$cluster_counts[[method]]), integer(1))
  reference_index <- select_reference_fold(counts)
  reference <- stats::setNames(
    as.character(fold_fits[[reference_index]]$all_training_assignments[[method]]),
    fold_fits[[reference_index]]$train_ids
  )
  aligned_training <- vector("list", length(fold_fits))
  aligned_prediction <- vector("list", length(fold_fits))
  mappings <- vector("list", length(fold_fits))
  reference_levels <- sort(unique(reference[!is.na(reference)]))
  for (index in seq_along(fold_fits)) {
    native_training <- stats::setNames(
      as.character(fold_fits[[index]]$all_training_assignments[[method]]),
      fold_fits[[index]]$train_ids
    )
    native_prediction <- stats::setNames(
      as.character(fold_fits[[index]]$test_assignments[[method]]),
      fold_fits[[index]]$test_ids
    )
    if (index == reference_index) {
      mapping <- stats::setNames(reference_levels, reference_levels)
      aligned_training[[index]] <- native_training
      mappings[[index]] <- list(mapping = mapping, labels = native_training)
    } else {
      mappings[[index]] <- align_fold_labels_legacy(reference, native_training)
      mapping <- mappings[[index]]$mapping
      aligned_training[[index]] <- mappings[[index]]$labels
    }
    aligned_prediction[[index]] <- stats::setNames(
      unname(mapping[native_prediction]), names(native_prediction)
    )
  }
  label_matrix <- matrix(NA_character_, length(participant_ids), length(fold_fits),
                         dimnames = list(participant_ids, paste0("fold_", seq_along(fold_fits))))
  for (index in seq_along(aligned_training)) {
    label_matrix[names(aligned_training[[index]]), index] <- aligned_training[[index]]
  }
  stability <- derive_stable_labels(
    label_matrix, threshold = stability_threshold,
    expected_assignments = length(fold_fits) - 1L
  )
  # Compare one fold's aligned prediction with the subject's stable reference.
  prediction_rows <- do.call(rbind, lapply(seq_along(fold_fits), function(index) {
    prediction <- aligned_prediction[[index]]
    participant_index <- match(names(prediction), participant_ids)
    target <- stability$stable_label[participant_index]
    data.frame(
      method = method,
      fold = fold_fits[[index]]$fold,
      CIMBI = names(prediction),
      prediction = unname(prediction),
      stable_target = target,
      scored = !is.na(prediction) & !is.na(target),
      correct = ifelse(!is.na(prediction) & !is.na(target),
                       unname(prediction) == target, NA),
      stringsAsFactors = FALSE
    )
  }))
  pairwise <- list()
  counter <- 0L
  for (first in seq_len(length(fold_fits) - 1L)) {
    for (second in (first + 1L):length(fold_fits)) {
      overlap <- intersect(names(aligned_training[[first]]), names(aligned_training[[second]]))
      usable <- overlap[!is.na(aligned_training[[first]][overlap]) &
                          !is.na(aligned_training[[second]][overlap])]
      counter <- counter + 1L
      pairwise[[counter]] <- data.frame(
        method = method, fold_a = fold_fits[[first]]$fold,
        fold_b = fold_fits[[second]]$fold, n_overlap = length(usable),
        ari = if (length(usable) >= 2L) mclust::adjustedRandIndex(
          aligned_training[[first]][usable], aligned_training[[second]][usable]
        ) else NA_real_
      )
    }
  }
  overlap_ari <- do.call(rbind, pairwise)
  scored <- prediction_rows[prediction_rows$scored, , drop = FALSE]
  recalls <- if (nrow(scored)) tapply(scored$correct, scored$stable_target, mean) else numeric()
  list(
    method = method,
    reference_fold = fold_fits[[reference_index]]$fold,
    mappings = mappings,
    aligned_training = aligned_training,
    aligned_prediction = aligned_prediction,
    label_matrix = label_matrix,
    stability = stability,
    predictions = prediction_rows,
    overlap_ari = overlap_ari,
    confusion = table(scored$stable_target, scored$prediction),
    metrics = data.frame(
      method = method,
      mean_n_clusters = mean(counts, na.rm = TRUE),
      mean_ari = mean(overlap_ari$ari, na.rm = TRUE),
      accuracy = if (nrow(scored)) mean(scored$correct) else NA_real_,
      balanced_accuracy = if (length(recalls)) mean(recalls) else NA_real_,
      coverage = mean(!is.na(prediction_rows$prediction)),
      rejected_n = sum(is.na(prediction_rows$prediction)),
      stable_n = sum(stability$stable),
      stringsAsFactors = FALSE
    )
  )
}

# Build the complete SPEC summary across every base method and consensus,
# including final stable labels and relational metric tables.
summarize_cv_spec <- function(fold_fits, participant_ids,
                              stability_threshold = 0.75) {
  methods <- names(fold_fits[[1L]]$test_assignments)
  # Apply the identical alignment and metric contract to one method.
  per_method <- stats::setNames(lapply(methods, function(method) {
    summarize_one_method_spec(fold_fits, participant_ids, method,
                              stability_threshold)
  }), methods)
  final_labels <- data.frame(CIMBI = participant_ids, stringsAsFactors = FALSE)
  for (method in methods) {
    final_labels[[paste0(method, "_label")]] <- per_method[[method]]$stability$stable_label
    final_labels[[paste0(method, "_stability")]] <-
      per_method[[method]]$stability$stability_proportion
  }
  fold_labels <- do.call(rbind, lapply(methods, function(method) {
    method_result <- per_method[[method]]
    training_labels <- as.data.frame(method_result$label_matrix,
                                     stringsAsFactors = FALSE)
    data.frame(
      CIMBI = participant_ids,
      method = method,
      training_labels,
      final_label = method_result$stability$stable_label,
      stability = method_result$stability$stability_proportion,
      stable = method_result$stability$stable,
      n_assignments = method_result$stability$n_assignments,
      stringsAsFactors = FALSE
    )
  }))
  list(
    per_method = per_method,
    metrics = do.call(rbind, lapply(per_method, `[[`, "metrics")),
    final_labels = final_labels,
    fold_labels = fold_labels,
    predictions = do.call(rbind, lapply(per_method, `[[`, "predictions")),
    overlap_ari = do.call(rbind, lapply(per_method, `[[`, "overlap_ari")),
    # Convert one fold's named cluster counts to relational long form.
    cluster_counts = do.call(rbind, lapply(fold_fits, function(fit) {
      data.frame(fold = fit$fold, method = names(fit$cluster_counts),
                 n_clusters = unname(fit$cluster_counts))
    }))
  )
}

# Reproduce the historical consensus-only CV summary for compatibility mode;
# new SPEC analysis uses summarize_cv_spec().
summarize_cv_legacy <- function(fold_fits, participant_ids) {
  valid <- !vapply(fold_fits, is.null, logical(1))
  fold_fits <- fold_fits[valid]
  if (length(fold_fits) < 2L) stop("At least two successful folds are required.")
  # Convert one compatibility fold's consensus count to a summary row.
  cluster_counts <- do.call(rbind, lapply(fold_fits, function(fit) {
    data.frame(fold = fit$fold, method = names(fit$cluster_counts),
               n_clusters = unname(fit$cluster_counts))
  }))
  reference <- fold_fits[[1L]]$consensus_training
  aligned_training <- vector("list", length(fold_fits))
  aligned_prediction <- vector("list", length(fold_fits))
  aligned_training[[1L]] <- reference
  aligned_prediction[[1L]] <- fold_fits[[1L]]$consensus_prediction
  for (index in seq.int(2L, length(fold_fits))) {
    mapping <- align_fold_labels_legacy(reference, fold_fits[[index]]$consensus_training)
    aligned_training[[index]] <- mapping$labels
    native_prediction <- as.character(fold_fits[[index]]$consensus_prediction)
    aligned_prediction[[index]] <- stats::setNames(
      unname(mapping$mapping[native_prediction]), names(fold_fits[[index]]$consensus_prediction)
    )
  }
  label_matrix <- matrix(NA_character_, length(participant_ids), length(fold_fits),
                         dimnames = list(participant_ids, paste0("fold_", seq_along(fold_fits))))
  for (index in seq_along(aligned_training)) {
    label_matrix[names(aligned_training[[index]]), index] <- aligned_training[[index]]
  }
  # Derive the historical fold-modal label for one participant.
  subject_consensus <- do.call(rbind, lapply(seq_along(participant_ids), function(index) {
    labels <- stats::na.omit(label_matrix[index, ])
    counts <- sort(table(labels), decreasing = TRUE)
    data.frame(
      CIMBI = participant_ids[index],
      n_training = length(labels),
      consensus_label = if (length(counts)) names(counts)[1L] else NA_character_,
      consensus_proportion = if (length(counts)) unname(counts[1L]) / length(labels) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  fold_labels <- data.frame(
    CIMBI = participant_ids,
    method = "consensus",
    as.data.frame(label_matrix, stringsAsFactors = FALSE),
    final_label = subject_consensus$consensus_label,
    stability = subject_consensus$consensus_proportion,
    stable = NA,
    n_assignments = subject_consensus$n_training,
    stringsAsFactors = FALSE
  )
  # Score one compatibility fold's raw-centroid held-out predictions.
  prediction_rows <- do.call(rbind, lapply(seq_along(fold_fits), function(index) {
    prediction <- aligned_prediction[[index]]
    target <- subject_consensus$consensus_label[match(names(prediction), subject_consensus$CIMBI)]
    data.frame(
      fold = fold_fits[[index]]$fold, CIMBI = names(prediction),
      prediction = unname(prediction), target = target,
      correct = unname(prediction) == target, stringsAsFactors = FALSE
    )
  }))
  pairwise <- list()
  counter <- 0L
  for (first in seq_len(length(aligned_training) - 1L)) {
    for (second in (first + 1L):length(aligned_training)) {
      overlap <- intersect(names(aligned_training[[first]]), names(aligned_training[[second]]))
      usable <- overlap[!is.na(aligned_training[[first]][overlap]) &
                          !is.na(aligned_training[[second]][overlap])]
      counter <- counter + 1L
      pairwise[[counter]] <- data.frame(
        fold_a = fold_fits[[first]]$fold, fold_b = fold_fits[[second]]$fold,
        n_overlap = length(usable),
        ari = if (length(usable) >= 2L) mclust::adjustedRandIndex(
          aligned_training[[first]][usable], aligned_training[[second]][usable]
        ) else NA_real_
      )
    }
  }
  list(
    cluster_counts = cluster_counts,
    mean_cluster_counts = stats::aggregate(n_clusters ~ method, cluster_counts, mean),
    aligned_training = aligned_training,
    label_matrix = label_matrix,
    subject_consensus = subject_consensus,
    fold_labels = fold_labels,
    predictions = prediction_rows,
    accuracy = mean(prediction_rows$correct, na.rm = TRUE),
    overlap_ari = do.call(rbind, pairwise)
  )
}

# Orchestrate every outer fold, retaining fold fits and applying the SPEC or
# compatibility summarizer without global state.
run_cross_validation <- function(data, distance_builder, predictor_variables,
                                 sex_col, medication_col, age_col, folds = 5L,
                                 base_config = default_base_config(),
                                 consensus_config = default_consensus_config(),
                                 prediction_mode = "legacy_raw_centroid",
                                 seed = 1L, feature_pipeline = NULL,
                                 age_bin_width = 12, age_bin_origin = 18,
                                 balance_cols = NULL) {
  data <- canonicalize_records(data, age_col = age_col)
  fold_plan <- make_balanced_folds(data, folds, sex_col, medication_col, age_col,
                                   seed, age_bin_width, age_bin_origin, balance_cols)
  # Fit one deterministic fold using its balanced held-out membership.
  fits <- lapply(seq_len(folds), function(fold) {
    run_cv_fold(
      data, which(fold_plan$fold_id != fold), which(fold_plan$fold_id == fold), fold,
      distance_builder, predictor_variables, base_config, consensus_config,
      prediction_mode, seed + fold * 10000L, feature_pipeline, age_col
    )
  })
  summary <- if (identical(prediction_mode, "feature_consensus")) {
    summarize_cv_spec(fits, as.character(data$CIMBI))
  } else {
    summarize_cv_legacy(fits, as.character(data$CIMBI))
  }
  structure(list(
    data_ids = as.character(data$CIMBI), fold_id = fold_plan$fold_id,
    fold_plan = fold_plan, fold_fits = fits, summary = summary,
    config = list(folds = folds, base = base_config, consensus = consensus_config,
                  prediction_mode = prediction_mode,
                  age_bin_width = age_bin_width,
                  age_bin_origin = age_bin_origin), seed = as.integer(seed)
  ), class = c("cv_clustering_result", "list"))
}

# Historical compatibility wrapper for blocked permutations; complete SPEC
# inference lives in resampling.R.
run_blocked_permutations_legacy <- function(data, task_blocks, iterations,
                                            analysis_function, seed = 1L) {
  # Permute and analyze one legacy Monte Carlo iteration.
  lapply(seq_len(iterations), function(iteration) {
    set.seed(seed + iteration)
    permuted <- data
    for (block in task_blocks) {
      donors <- sample.int(nrow(data))
      permuted[, block] <- data[donors, block, drop = FALSE]
    }
    analysis_function(permuted, seed + iteration)
  })
}

# Historical compatibility wrapper for participant bootstraps with unique
# execution IDs.
run_bootstraps_legacy <- function(data, iterations, analysis_function, seed = 1L) {
  # Draw and analyze one with-replacement participant sample.
  lapply(seq_len(iterations), function(iteration) {
    set.seed(seed + iteration)
    indices <- sample.int(nrow(data), nrow(data), replace = TRUE)
    sample <- data[indices, , drop = FALSE]
    if ("CIMBI" %in% names(sample)) {
      sample$CIMBI <- make.unique(paste0(as.character(sample$CIMBI), "_bootstrap"))
    }
    analysis_function(sample, seed + iteration)
  })
}
