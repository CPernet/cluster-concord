# Deterministic known-structure validation for the canonical SPEC workflow.

# Capture this module's source path for stable project/module discovery.
.synthetic_validation_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

# Locate the repository root so generated validation results are written under
# the canonical results directory regardless of the caller's working directory.
synthetic_project_root <- function() {
  candidates <- unique(c(
    if (!is.null(.synthetic_validation_file)) {
      dirname(dirname(dirname(normalizePath(.synthetic_validation_file))))
    },
    getwd(), dirname(getwd()), dirname(dirname(getwd()))
  ))
  root <- candidates[file.exists(file.path(candidates, "SPEC.md"))][1L]
  if (is.na(root)) stop("Cannot locate project root containing SPEC.md.")
  normalizePath(root, winslash = "/")
}

# Load the canonical analysis entry point and its dependencies when validation is
# invoked standalone, without re-sourcing an already loaded workflow.
ensure_synthetic_modules <- function() {
  if (exists("clustering_analysis", mode = "function", inherits = TRUE)) {
    return(invisible(TRUE))
  }
  candidates <- unique(c(
    if (!is.null(.synthetic_validation_file)) dirname(normalizePath(.synthetic_validation_file)),
    file.path(getwd(), "src", "R"),
    file.path(getwd(), "..", "src", "R"),
    getwd()
  ))
  directory <- candidates[file.exists(file.path(candidates, "clustering_analysis.R"))][1L]
  path <- file.path(directory, "clustering_analysis.R")
  if (!file.exists(path)) stop("Cannot find clustering_analysis.R: ", path)
  # These canonical modules are shared by several helper functions below, so
  # load them into the caller session rather than a transient function frame.
  sys.source(path, envir = .GlobalEnv)
  invisible(TRUE)
}

# Generate independent known cluster templates over ten paired cognitive task
# measures; truth labels are retained only for post-fit recovery evaluation.
generate_known_structure_data <- function(n_per_cluster = 24L, clusters = 3L,
                                          noise_sd = 0.08, seed = 2601L) {
  ensure_synthetic_modules()
  n_per_cluster <- as.integer(n_per_cluster)
  clusters <- as.integer(clusters)
  if (clusters < 2L || n_per_cluster < 10L) {
    stop("Use at least two clusters and ten participants per cluster.")
  }
  set.seed(as.integer(seed))
  task_spec <- cognitive_task_spec()
  variables <- unique(c(task_spec$perf_var, task_spec$time_var))
  truth <- rep(seq_len(clusters), each = n_per_cluster)

  # Independent cluster templates in the paired task-measure space. Applying
  # fit_cognitive_preprocessor() converts these into the intended 45 pairwise
  # within-participant task-similarity features.
  templates <- matrix(stats::rnorm(clusters * length(variables), sd = 1.8),
                      nrow = clusters)
  values <- templates[truth, , drop = FALSE] +
    matrix(stats::rnorm(length(truth) * length(variables), sd = noise_sd),
           nrow = length(truth))
  colnames(values) <- variables
  ids <- sprintf("SYN%04d", seq_along(truth))
  data <- data.frame(
    CIMBI = ids,
    Gender = rep(c("F", "M"), length.out = length(truth)),
    Antidepressant = rep(c("No", "Yes", "No"), length.out = length(truth)),
    Age = 20 + ((seq_along(truth) * 7L) %% 45L),
    truth = truth,
    values,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  list(data = data, truth = stats::setNames(truth, ids), task_spec = task_spec,
       seed = as.integer(seed), noise_sd = noise_sd)
}

# Compute ARI on available labels with explicit perfect handling for identical
# degenerate partitions where the standard ARI denominator is zero.
safe_partition_ari <- function(truth, prediction) {
  usable <- !is.na(truth) & !is.na(prediction)
  if (sum(usable) < 2L) return(NA_real_)
  truth <- as.character(truth[usable])
  prediction <- as.character(prediction[usable])
  value <- mclust::adjustedRandIndex(truth, prediction)
  if (is.finite(value)) return(value)
  # ARI has a zero denominator for some identical degenerate partitions (for
  # example, every observation is a singleton). Their pairwise equivalence
  # relations still provide an unambiguous perfect-partition result.
  if (identical(outer(truth, truth, "=="), outer(prediction, prediction, "=="))) {
    return(1)
  }
  NA_real_
}

# Exhaustively fit every 2-to-5-method consensus subset against the same base
# outputs and rank recovery using simulation truth as an oracle diagnostic only.
compare_consensus_method_subsets <- function(base_result, truth,
                                             k_candidates = 2:4,
                                             min_cluster_size = 3L) {
  methods <- base_method_names()
  truth <- truth[base_result$participant_ids]
  # Score one base method's ARI, coverage, and recovered cluster count.
  method_performance <- do.call(rbind, lapply(methods, function(method) {
    labels <- base_result$assignments[[method]]
    data.frame(
      method = method,
      ari = safe_partition_ari(truth, labels),
      coverage = mean(!is.na(labels)),
      n_clusters = length(unique(labels[!is.na(labels)])),
      stringsAsFactors = FALSE
    )
  }))
  # Enumerate all method combinations of one requested subset size.
  subsets <- unlist(lapply(2:length(methods), function(size) {
    utils::combn(methods, size, simplify = FALSE)
  }), recursive = FALSE)
  # Fit and score one fixed method subset through the standard consensus API.
  rows <- lapply(subsets, function(subset) {
    consensus <- tryCatch(
      run_consensus_clustering(
        base_result,
        default_consensus_config(
          mode = "spec", active_methods = subset,
          k_candidates = k_candidates, min_cluster_size = min_cluster_size
        )
      ),
      error = identity
    )
    if (inherits(consensus, "error")) {
      return(data.frame(
        methods = paste(subset, collapse = "+"), method_count = length(subset),
        ari = NA_real_, selected_k = NA_integer_, coverage = NA_real_,
        error = conditionMessage(consensus), stringsAsFactors = FALSE
      ))
    }
    data.frame(
      methods = paste(subset, collapse = "+"), method_count = length(subset),
      ari = safe_partition_ari(truth, consensus$labels),
      selected_k = consensus$selected_k,
      coverage = mean(!is.na(consensus$labels)), error = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  comparison <- do.call(rbind, rows)
  comparison$is_all_five <- comparison$method_count == length(methods)
  ordering <- order(-comparison$ari, comparison$method_count, comparison$methods,
                    na.last = TRUE)
  comparison$oracle_rank <- match(seq_len(nrow(comparison)), ordering)
  comparison <- comparison[order(comparison$oracle_rank), , drop = FALSE]
  rownames(comparison) <- NULL
  best_methods <- strsplit(comparison$methods[1L], "+", fixed = TRUE)[[1L]]
  best_consensus <- run_consensus_clustering(
    base_result,
    default_consensus_config(
      mode = "spec", active_methods = best_methods,
      k_candidates = k_candidates, min_cluster_size = min_cluster_size
    )
  )
  list(
    subset_comparison = comparison,
    base_method_performance = method_performance,
    best_subset = comparison[1L, , drop = FALSE],
    best_subset_labels = best_consensus$labels,
    scope = "simulation_oracle_only"
  )
}

# Convert a data frame to a compact GitHub-flavored Markdown table for the
# generated known-truth report.
markdown_table <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  # Format one cell deterministically and escape Markdown pipe characters.
  format_cell <- function(value) {
    if (length(value) == 0L || is.na(value)) return("")
    if (is.numeric(value)) value <- format(round(value, 4L), trim = TRUE)
    gsub("|", "\\|", as.character(value), fixed = TRUE)
  }
  # Render one formatted data row between Markdown table delimiters.
  rows <- apply(data, 1L, function(row) {
    paste0("| ", paste(vapply(row, format_cell, character(1)), collapse = " | "), " |")
  })
  c(
    paste0("| ", paste(names(data), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(data)), collapse = " | "), " |"),
    rows
  )
}

# Write the complete known-truth validation narrative, method/subset rankings,
# and participant-level recovered labels to results Markdown.
write_synthetic_validation_report <- function(validation, path = NULL) {
  if (!inherits(validation, "synthetic_validation_result")) {
    stop("validation must be a synthetic_validation_result.")
  }
  if (is.null(path)) {
    path <- file.path(synthetic_project_root(), "results", "synthetic_validation_results.md")
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  best <- validation$best_subset[1L, , drop = FALSE]
  all_five <- validation$subset_comparison[
    validation$subset_comparison$is_all_five, , drop = FALSE
  ][1L, , drop = FALSE]
  improvement <- best$ari - all_five$ari
  improvement_statement <- if (!is.finite(improvement)) {
    "The ARI difference could not be estimated."
  } else if (improvement > 1e-12) {
    sprintf("The best subset improves ARI over all five methods by %.4f.", improvement)
  } else if (improvement < -1e-12) {
    sprintf("The all-five consensus exceeds the ranked subset by %.4f.", -improvement)
  } else {
    "The best-ranked subset and the all-five consensus have equal ARI; the subset wins only the deterministic fewer-method tie-break."
  }

  truth <- validation$simulation$truth
  ids <- names(truth)
  final <- stats::setNames(
    as.character(validation$result$final_fit$consensus$labels),
    validation$result$final_fit$base$participant_ids
  )
  best_labels <- stats::setNames(
    as.character(validation$best_subset_labels),
    validation$result$final_fit$base$participant_ids
  )
  stable <- validation$result$cv$summary$final_labels
  stable_labels <- stats::setNames(as.character(stable$consensus_label), stable$CIMBI)
  heldout <- validation$result$generalization$consensus_prediction$native_label
  labels <- data.frame(
    CIMBI = ids,
    population = ifelse(ids %in% validation$holdout_ids, "holdout", "development"),
    truth = unname(truth),
    all_five_final = unname(final[ids]),
    best_subset_final = unname(best_labels[ids]),
    stable_cv = unname(stable_labels[ids]),
    holdout_prediction = unname(heldout[ids]),
    stringsAsFactors = FALSE
  )
  performance <- data.frame(
    metric = names(validation$performance),
    value = unlist(validation$performance, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  lines <- c(
    "# Synthetic consensus-clustering validation results",
    "",
    "## Synthetic data and validation design",
    "",
    sprintf("- Seed: `%d`.", validation$simulation$seed),
    sprintf("- Participants: %d development and %d held out (floor 10%% overall, allocated nearly equally across truth clusters).",
            validation$performance$development_n, validation$performance$holdout_n),
    "- Known structure: three clusters generated from independent templates over ten paired performance/reaction-time task measures, with Gaussian within-cluster noise.",
    "- Representation: training-fitted preprocessing transforms the 20 task measures into 45 pairwise task-similarity features.",
    "- Validation: five-fold CV, 3-of-4 stable labels, all five clustering algorithms, fixed-partition held-out prediction, and every 2-to-5-method consensus subset.",
    "- Bootstrap and permutation loops are excluded because this report evaluates recovery against known truth; those mechanisms have dedicated tests.",
    "",
    "## Recovery summary",
    "",
    markdown_table(performance),
    "",
    "## Individual clustering methods",
    "",
    markdown_table(validation$base_method_performance),
    "",
    "## Best oracle consensus subset",
    "",
    sprintf("The best-ranked method combination is **%s** with ARI %.4f and selected k = %s.",
            best$methods, best$ari, best$selected_k),
    sprintf("The all-five consensus has ARI %.4f and selected k = %s. %s",
            all_five$ari, all_five$selected_k, improvement_statement),
    "",
    "This is a simulation oracle result: known truth is used to rank subsets. It must not be used to select methods on the empirical dataset without a prespecified or nested validation rule.",
    "",
    "## Consensus subset ranking",
    "",
    markdown_table(validation$subset_comparison),
    "",
    "## Participant-level recovered labels",
    "",
    sprintf("`best_subset_final` uses: %s.", best$methods),
    "Blank cells are structurally unavailable: development participants have no holdout prediction, and held-out participants are not part of development fits.",
    "",
    markdown_table(labels)
  )
  writeLines(lines, path, useBytes = TRUE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

# Run the deterministic full SPEC validation with five-fold CV, all five base
# methods, fixed 10% holdout prediction, and oracle subset comparison. Resampling
# is intentionally excluded because known-truth recovery is the target here.
run_synthetic_validation <- function(seed = 2601L, n_per_cluster = 24L,
                                     folds = 5L,
                                     generalization_permutations = 999L,
                                     output_dir = "results/synthetic_validation",
                                     write_outputs = TRUE) {
  ensure_synthetic_modules()
  simulation <- generate_known_structure_data(
    n_per_cluster = n_per_cluster, seed = seed
  )
  all_data <- simulation$data
  cluster_levels <- sort(unique(all_data$truth))
  set.seed(as.integer(seed + 1L))
  holdout_n <- max(length(cluster_levels), floor(nrow(all_data) * 0.10))
  allocation <- rep(holdout_n %/% length(cluster_levels), length(cluster_levels))
  remainder <- holdout_n %% length(cluster_levels)
  if (remainder) allocation[seq_len(remainder)] <- allocation[seq_len(remainder)] + 1L
  # Draw one cluster's allocated share of the stratified floor-10% holdout.
  holdout_indices <- unlist(lapply(seq_along(cluster_levels), function(index) {
    group <- cluster_levels[index]
    candidates <- which(all_data$truth == group)
    sample(candidates, allocation[index])
  }), use.names = FALSE)
  development <- all_data[-holdout_indices, , drop = FALSE]
  holdout <- all_data[holdout_indices, , drop = FALSE]
  task_spec <- simulation$task_spec
  # Resolve one simulated task to its paired performance/reaction-time columns.
  task_blocks <- lapply(seq_len(nrow(task_spec)), function(index) {
    c(task_spec$perf_var[index], task_spec$time_var[index])
  })

  result <- clustering_analysis(
    data = development,
    generalization_data = holdout,
    predictor_variables = unique(unlist(task_blocks, use.names = FALSE)),
    sex_col = "Gender", medication_col = "Antidepressant", age_col = "Age",
    folds = folds, resamples = 0L, seed = seed + 2L, workers = 1L,
    mode = "spec", task_blocks = task_blocks,
    generalization_permutations = generalization_permutations,
    output_dir = output_dir, write_outputs = write_outputs,
    base_config = default_base_config(
      mode = "spec", methods = base_method_names(), k_candidates = 2:4,
      min_cluster_size = 3L, skmeans_restarts = 5L, ap_grid_size = 15L,
      hdbscan_min_pts = 4L, optional_failure = "error"
    ),
    consensus_config = default_consensus_config(
      mode = "spec", active_methods = base_method_names(),
      k_candidates = 2:4, min_cluster_size = 3L
    )
  )

  development_truth <- simulation$truth[result$final_fit$base$participant_ids]
  final_labels <- result$final_fit$consensus$labels
  stable_table <- result$cv$summary$final_labels
  stable_truth <- simulation$truth[stable_table$CIMBI]
  holdout_prediction <- result$generalization$consensus_prediction$native_label
  holdout_truth <- simulation$truth[names(holdout_prediction)]
  subset_evaluation <- compare_consensus_method_subsets(
    result$final_fit$base, simulation$truth,
    k_candidates = 2:4, min_cluster_size = 3L
  )
  performance <- list(
    final_consensus_ari = safe_partition_ari(development_truth, final_labels),
    stable_consensus_ari = safe_partition_ari(
      stable_truth, stable_table$consensus_label
    ),
    stable_coverage = mean(!is.na(stable_table$consensus_label)),
    holdout_consensus_ari = safe_partition_ari(holdout_truth, holdout_prediction),
    holdout_coverage = result$generalization$coverage,
    development_n = nrow(development),
    holdout_n = nrow(holdout),
    feature_count = ncol(result$final_fit$base$features)
  )
  validation <- structure(list(
    result = result,
    performance = performance,
    subset_comparison = subset_evaluation$subset_comparison,
    base_method_performance = subset_evaluation$base_method_performance,
    best_subset = subset_evaluation$best_subset,
    best_subset_labels = subset_evaluation$best_subset_labels,
    subset_selection_scope = subset_evaluation$scope,
    simulation = simulation,
    development_ids = development$CIMBI,
    holdout_ids = holdout$CIMBI
  ), class = c("synthetic_validation_result", "list"))
  if (isTRUE(write_outputs)) {
    validation$report_path <- write_synthetic_validation_report(validation)
  }
  validation
}
