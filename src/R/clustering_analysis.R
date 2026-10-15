# CONSENSUS-CLUSTERING ANALYSIS ENTRY POINT
# =========================================
#
# This file provides clustering_analysis(), the canonical implementation of the
# workflow described in SPEC.md. Source this file, then call the function with
# either a raw Excel workbook (`input_path`) or an in-memory data frame (`data`).
# Do not provide both unless `data` is intended to take precedence.
#
# Production call using a raw workbook and ten parallel Monte Carlo workers:
#
#   source("src/R/clustering_analysis.R")
#   result <- clustering_analysis(
#     input_path = "data/data.xlsx",
#     workers = 10L
#   )
#
# `workers` controls parallel execution of the permutation/bootstrap loop. Ten
# workers are shown because the default 1,000 resamples are computationally
# expensive; choose a value supported by the server's CPU and memory capacity.
#
# A shorter deterministic development run (no permutation/bootstrap inference):
#
#   result <- clustering_analysis(
#     input_path = "data/data.xlsx",
#     resamples = 0L,
#     write_outputs = FALSE
#   )
#
# The principal defaults implement the SPEC workflow:
#
# - `mode = "spec"`: leakage-free preprocessing, feature-space prediction,
#   five-method consensus, final development fit, inference, and generalization.
# - `methods = base_method_names()`: HCA/UPGMA, spherical k-means, GMM,
#   affinity propagation, and HDBSCAN are all requested.
# - `k_candidates = 2:6`: candidate cluster counts considered by bounded methods
#   and consensus selection.
# - `folds = 5L`: five nearly equal folds balanced by sex, antidepressant use,
#   and age distribution.
# - `resamples = 1000L`: 1,000 blocked permutations and 1,000 participant
#   bootstraps. Set to zero to skip inference during development or testing.
# - `seed = 2025L`: master deterministic seed; distinct offsets are used for
#   final fitting, generalization, and resampling.
# - `workers = 1L`: sequential, Windows-safe execution. Parallel execution is
#   opt-in and uses deterministic per-iteration seeds.
# - `output_dir = default_clustering_output_dir()` and `write_outputs = TRUE`:
#   save the complete RDS result and CSV summaries in the repository-root
#   `results/` directory, independently of the caller's working directory.
# - `generalization_method = "auto"`: enumerate the exact multinomial null
#   when tractable and otherwise use `generalization_permutations = 9999L`
#   Monte Carlo replicates for the holdout calibration test.
#
# Input and filtering arguments:
#
# - `input_path`: optional raw Excel workbook. The loader applies cohort rules,
#   screens missingness, creates the seeded eligible-cohort holdout, and returns
#   raw development data that passed the same screen.
# - `data`: optional in-memory development data. It must contain unique `CIMBI`
#   IDs, the balancing covariates, and either the configured cognitive variables
#   or columns accepted by a custom feature pipeline.
# - `generalization_data`: optional fixed holdout data when `data` is supplied.
#   With `input_path` in SPEC mode, the loader's screened holdout is used.
# - `holdout_fraction` / `holdout_cohort`: define the post-screen sampling
#   fraction and eligible population. A NULL cohort uses all screened rows.
# - `study_filter`: optional allowed values of `Project ID (primary)`.
# - `sex_col`, `medication_col`, and `age_col`: names used to balance CV folds.
# - age_bin_width = 12 / age_bin_origin = 18: fixed-width age groups
#   extending automatically over the input range; missing age is its own group.
#
# Optional extension arguments:
#
# - `data_loader`: optional loader implementing the contract in src/custom_data.md.
#   Requires a custom feature_pipeline; otherwise the behavioural loader is used.
# - `balance_cols`: NULL preserves sex_col/medication_col; a character vector
#   selects categorical balancing columns (character(0) selects none).
#   Set age_col = NULL to disable numeric age balancing.
# - `distance_builder`: compatibility-mode prepared-frame distance builder.
#   SPEC mode instead obtains distance results from feature_pipeline.
# - `predictor_variables`: explicit variables required by the analysis. Defaults
#   to the performance and reaction-time columns in cognitive_task_spec().
# - `base_config` / `consensus_config`: complete configuration overrides. When
#   NULL, they are created from mode, methods, and k_candidates.
# - `feature_pipeline`: function(train, test) returning `train`, `test`, and
#   `preprocessor`. The SPEC default fits every learned transform on train only.
# - `task_blocks`: list of ten performance/reaction-time column pairs used by
#   blocked permutations. The cognitive task specification is the default.
# - `mode = "compatibility"`: historical reproduction only; it retains the old
#   raw-centroid prediction behavior and omits SPEC-only final inference stages.
#
# Return value:
#
# A `clustering_analysis_result` list containing cross-validation fits and
# summaries (`cv`), the fixed full-development base/consensus fit (`final_fit`),
# loader audit state (`loaded`), permutation/bootstrap inference (`inference`),
# holdout calibration (`generalization`), run metadata, and written paths.
# HDBSCAN noise and rejected/unstable assignments remain `NA`, not extra clusters.

# Capture the path used to source this file so dependent modules can be located
# without relying on a caller's working directory.
.clustering_analysis_file <- tryCatch({
  source_file <- NULL
  for (frame in rev(sys.frames())) {
    candidate <- frame$ofile
    if (!is.null(candidate) &&
        identical(basename(candidate), "clustering_analysis.R")) {
      source_file <- candidate
      break
    }
  }
  if (is.null(source_file)) NULL else normalizePath(
    source_file, winslash = "/", mustWork = TRUE
  )
}, error = function(e) NULL)

# Locate the canonical src/R module directory when this file is sourced from the
# repository root, the test directory, or src/R itself.
clustering_module_directory <- function() {
  current_file <- .clustering_analysis_file
  if (!is.null(current_file)) return(dirname(normalizePath(current_file)))
  candidates <- c(
    file.path(getwd(), "src", "R"),
    file.path(getwd(), "..", "src", "R"),
    getwd()
  )
  candidates[file.exists(file.path(candidates, "load_and_prepdata.R"))][1L]
}

# Locate the repository root from the canonical src/R module directory so
# default artifacts never depend on the caller's current working directory.
clustering_project_root <- function() {
  module_directory <- clustering_module_directory()
  if (is.na(module_directory) || !nzchar(module_directory)) {
    stop("Cannot locate the clustering module directory.")
  }
  root <- normalizePath(
    file.path(module_directory, "..", ".."),
    winslash = "/", mustWork = TRUE
  )
  if (!file.exists(file.path(root, "SPEC.md"))) {
    stop("Cannot locate the repository root containing SPEC.md.")
  }
  root
}

# Return the canonical repository-root results directory used when callers do
# not explicitly override output_dir.
default_clustering_output_dir <- function() {
  file.path(clustering_project_root(), "results")
}

# Source every module required by the end-to-end SPEC workflow. The function
# fails early with the missing path instead of allowing a partial analysis.
clustering_module_files <- function(normalize = FALSE) {
  directory <- clustering_module_directory()
  modules <- c(
    "load_and_prepdata.R", "base_clustering.R",
    "consensus_clustering.R", "prediction.R", "cv_clustering.R",
    "resampling.R", "generalization.R"
  )
  paths <- file.path(directory, modules)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Cannot find workflow module: ", paste(missing, collapse = ", "))
  }
  if (isTRUE(normalize)) {
    paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  }
  stats::setNames(paths, modules)
}

source_clustering_modules <- function() {
  module_files <- clustering_module_files()
  for (module in names(module_files)) {
    if (!exists(sub("\\.R$", "", module), mode = "function", inherits = TRUE) ||
        identical(module, "load_and_prepdata.R")) {
      source(module_files[[module]], local = parent.frame())
    }
  }
  invisible(dirname(module_files[[1L]]))
}

source_clustering_modules()

# Return the default twenty cognitive inputs: one performance and one reaction-
# time variable for each of the ten tasks defined by cognitive_task_spec().
default_cognitive_predictors <- function(task_spec = cognitive_task_spec()) {
  unique(c(task_spec$perf_var, task_spec$time_var))
}

# Return the package namespaces needed by the selected base methods plus the
# packages used by cross-validation and consensus diagnostics for every method.
resampling_worker_packages <- function(methods) {
  method_packages <- c(
    spherical_kmeans = "skmeans",
    gmm = "mclust",
    affinity_propagation = "apcluster",
    hdbscan = "dbscan"
  )
  selected <- unname(method_packages[intersect(methods, names(method_packages))])
  unique(c("cluster", "clue", "mclust", selected))
}

# Persist the structured analysis result and relational CSV summaries required
# by the SPEC: fold cluster counts, predictions, overlap ARI, inference, and
# generalization proportions. Returns the paths of every artifact written.
write_clustering_artifacts <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list(
    result_rds = file.path(output_dir, "clustering_analysis_result.rds"),
    cv_cluster_counts = file.path(output_dir, "cv_cluster_counts.csv"),
    cv_predictions = file.path(output_dir, "cv_predictions.csv"),
    cv_overlap_ari = file.path(output_dir, "cv_overlap_ari.csv"),
    cv_final_labels = file.path(output_dir, "cv_final_labels.csv"),
    cv_fold_balance = file.path(output_dir, "cv_fold_balance.csv"),
    cv_stratum_counts = file.path(output_dir, "cv_stratum_counts.csv")
  )
  saveRDS(result, paths$result_rds)
  utils::write.csv(result$cv$summary$cluster_counts, paths$cv_cluster_counts,
                   row.names = FALSE)
  utils::write.csv(result$cv$summary$predictions, paths$cv_predictions,
                   row.names = FALSE)
  utils::write.csv(result$cv$summary$overlap_ari, paths$cv_overlap_ari,
                   row.names = FALSE)
  utils::write.csv(result$cv$summary$fold_labels, paths$cv_final_labels,
                   row.names = FALSE)
  utils::write.csv(result$cv$fold_plan$balance_long, paths$cv_fold_balance,
                   row.names = FALSE)
  utils::write.csv(result$cv$fold_plan$stratum_counts, paths$cv_stratum_counts,
                   row.names = FALSE)
  if (!is.null(result$inference)) {
    paths$permutation_p_values <- file.path(output_dir, "permutation_p_values.csv")
    paths$bootstrap_confidence_intervals <- file.path(output_dir, "bootstrap_confidence_intervals.csv")
    utils::write.csv(result$inference$inference$p_values,
                     paths$permutation_p_values, row.names = FALSE)
    utils::write.csv(result$inference$inference$confidence_intervals,
                     paths$bootstrap_confidence_intervals, row.names = FALSE)
  }
  if (!is.null(result$generalization)) {
    paths$generalization_calibration <- file.path(output_dir, "generalization_calibration.csv")
    utils::write.csv(result$generalization$proportions,
                     paths$generalization_calibration, row.names = FALSE)
  }
  paths
}

# Run the complete consensus-clustering workflow described in SPEC.md.
#
# The function validates/loads input, constructs balanced folds, fits and
# predicts all requested base methods, creates pairwise-available consensus,
# derives 3-of-4 stable labels, fits a fixed development partition, evaluates
# the untouched holdout, optionally runs blocked permutations and participant
# bootstraps, and writes structured results.
clustering_analysis <- function(
  input_path = NULL,
  data = NULL,
  generalization_data = NULL,
  distance_builder = NULL,
  predictor_variables = NULL,
  study_filter = NULL,
  holdout_fraction = 0.10,
  holdout_cohort = "KOHORTE 2",
  methods = base_method_names(),
  k_candidates = 2:6,
  folds = 5L,
  resamples = 1000L,
  seed = 2025L,
  workers = 1L,
  output_dir = default_clustering_output_dir(),
  mode = c("spec", "compatibility"),
  sex_col = "Gender",
  medication_col = "BDDin Q9 - Antidepressants",
  age_col = "Age at neuropsych test",
  base_config = NULL,
  consensus_config = NULL,
  feature_pipeline = NULL,
  task_blocks = NULL,
  generalization_method = c("auto", "exact", "monte_carlo"),
  generalization_permutations = 9999L,
  generalization_confidence = 0.95,
  generalization_exact_max_tables = 1000000L,
  write_outputs = TRUE,
  age_bin_width = 12,
  age_bin_origin = 18,
  balance_cols = NULL,
  data_loader = NULL
) {
  # Normalize the public mode and reject invalid execution settings before any
  # data are read or output directories are touched.
  mode <- match.arg(mode)
  generalization_method <- match.arg(generalization_method)
  if (workers < 1L) stop("workers must be positive.")
  if (resamples < 0L) stop("resamples cannot be negative.")
  loaded <- NULL
  raw_input <- FALSE
  custom_loader <- !is.null(data_loader)
  custom_feature_pipeline <- !is.null(feature_pipeline)
  if (custom_loader && !is.function(data_loader)) stop("data_loader must be a function.")
  if (custom_loader && !is.null(data)) stop("Use data_loader with input_path, or supply data directly.")
  if (custom_loader && !custom_feature_pipeline) stop("A custom data_loader requires feature_pipeline.")

  # When a workbook is supplied, screen rows and create the holdout before any
  # learned preprocessing. SPEC mode keeps screened raw rows so each CV training
  # fold can own imputation and transformation parameters.
  if (is.null(data)) {
    if (is.null(input_path)) stop("Provide either input_path or data.")
    loader <- if (custom_loader) data_loader else load_and_prepdata
    loaded <- loader(
      input_path = input_path,
      holdout_fraction = holdout_fraction,
      holdout_cohort = holdout_cohort,
      seed = seed,
      output_dir = output_dir,
      write_outputs = FALSE,
      return = "all"
    )
    if (custom_loader && identical(mode, "spec")) {
      required <- c("screened_development_raw_data", "screened_holdout_data")
      if (!is.list(loaded) || !all(required %in% names(loaded))) {
        stop("data_loader must return screened_development_raw_data and screened_holdout_data.")
      }
      if (!all(vapply(loaded[required], is.data.frame, logical(1)))) {
        stop("Custom loader development and holdout outputs must be data frames.")
      }
    }
    if (identical(mode, "spec")) {
      data <- loaded$screened_development_raw_data
      raw_input <- TRUE
      if (is.null(generalization_data) && nrow(loaded$screened_holdout_data)) {
        generalization_data <- loaded$screened_holdout_data
      }
    } else {
      data <- loaded$prepared_data
    }
  }

  # Enforce the participant-key contract shared by every downstream module.
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (!"CIMBI" %in% names(data)) stop("Input data must contain CIMBI.")
  data$CIMBI <- as.character(data$CIMBI)
  if (anyNA(data$CIMBI) || any(!nzchar(data$CIMBI)) || anyDuplicated(data$CIMBI)) {
    stop("Input participant IDs must be unique and non-missing.")
  }
  if (!is.null(study_filter)) {
    if (!"Project ID (primary)" %in% names(data)) stop("Study filter requires Project ID (primary).")
    data <- data[data[["Project ID (primary)"]] %in% study_filter, , drop = FALSE]
  }

  # Canonical ordering precedes every learned transform and stochastic fit.
  data <- canonicalize_records(data, age_col = age_col, require_distinct = TRUE)
  if (!is.null(generalization_data) && nrow(generalization_data)) {
    generalization_data <- canonicalize_records(
      as.data.frame(generalization_data), age_col = age_col, require_distinct = TRUE
    )
  }

  # Default distance builder: convert prepared cognitive measures into the 45
  # task-pair features and cosine subject distance required by the SPEC.
  if (is.null(distance_builder)) {
    # Convert one prepared frame into the structured distance result expected by
    # run_base_clustering().
    distance_builder <- function(frame) compute_subject_distance_matrix(frame, id_col = "CIMBI")
  }

  # Validate that direct/prepared inputs expose the configured cognitive
  # predictors. Raw workbook input is validated later by its fitted pipeline.
  if (is.null(predictor_variables)) predictor_variables <- default_cognitive_predictors()
  missing_predictors <- setdiff(predictor_variables, names(data))
  if (length(missing_predictors) && !(identical(mode, "spec") && raw_input)) {
    stop("Missing predictor variables: ", paste(missing_predictors, collapse = ", "))
  }

  # Build explicit base-method defaults unless the caller supplies a complete
  # override. SPEC defaults request all five algorithms and UPGMA base HCA.
  if (is.null(base_config)) {
    base_config <- default_base_config(
      mode = mode, methods = methods, k_candidates = k_candidates
    )
  }

  # Build explicit consensus defaults unless overridden. SPEC consensus uses
  # all requested methods, pairwise availability, and Ward D2 on 1-C.
  if (is.null(consensus_config)) {
    consensus_config <- default_consensus_config(
      mode = mode, active_methods = methods, k_candidates = k_candidates
    )
  }
  prediction_mode <- if (identical(mode, "compatibility")) {
    "legacy_raw_centroid"
  } else {
    "feature_consensus"
  }

  # Default raw-data feature pipeline: fit raw median imputation and all derived
  # feature transformations on the current training fold, then apply them
  # unchanged to that fold's validation participants.
  if (identical(mode, "spec") && is.null(feature_pipeline)) {
    if (raw_input) {
      # Fit raw and feature preprocessors on `train`; apply both fitted objects
      # unchanged to `test` and return auditable fitted state.
      feature_pipeline <- function(train, test) {
        raw_preprocessor <- fit_raw_cognitive_preprocessor(train, id_col = "CIMBI")
        derived_train <- apply_raw_cognitive_preprocessor(raw_preprocessor, train)
        derived_test <- apply_raw_cognitive_preprocessor(raw_preprocessor, test)
        feature_preprocessor <- fit_cognitive_preprocessor(derived_train, id_col = "CIMBI")
        list(
          train = apply_cognitive_preprocessor(feature_preprocessor, derived_train),
          test = apply_cognitive_preprocessor(feature_preprocessor, derived_test),
          preprocessor = list(raw = raw_preprocessor, features = feature_preprocessor)
        )
      }
    } else {
      # Default prepared-data feature pipeline: fit imputation, centering,
      # scaling, availability, and 45-feature construction on training only.
      feature_pipeline <- function(train, test) {
        preprocessor <- fit_cognitive_preprocessor(train, id_col = "CIMBI")
        list(
          train = apply_cognitive_preprocessor(preprocessor, train),
          test = apply_cognitive_preprocessor(preprocessor, test),
          preprocessor = preprocessor
        )
      }
    }
  }

  # Phase 1 of the SPEC analysis: five balanced outer folds. Each fold retains
  # every base assignment, fixed-partition consensus prediction, overlap ARI,
  # accuracy, coverage, and aligned labels used for the 3-of-4 stability rule.
  cv <- run_cross_validation(
    data = data,
    distance_builder = distance_builder,
    predictor_variables = predictor_variables,
    sex_col = sex_col,
    medication_col = medication_col,
    age_col = age_col,
    folds = folds,
    base_config = base_config,
    consensus_config = consensus_config,
    prediction_mode = prediction_mode,
    seed = seed,
    feature_pipeline = feature_pipeline,
    age_bin_width = age_bin_width, age_bin_origin = age_bin_origin,
    balance_cols = balance_cols
  )
  final_fit <- NULL
  generalization <- NULL
  inference <- NULL
  if (identical(mode, "spec")) {
    # Phase 2: fit preprocessing and all clustering methods once on the entire
    # development population. The placeholder lets the same train/apply API run
    # when no external holdout was requested; it never changes the fitted model.
    placeholder <- if (!is.null(generalization_data) && nrow(generalization_data)) {
      generalization_data
    } else {
      data[1L, , drop = FALSE]
    }
    full_features <- feature_pipeline(data, placeholder)
    final_base <- run_base_clustering(full_features$train, base_config, seed + 700000L)
    final_consensus <- run_consensus_clustering(final_base, consensus_config)
    final_fit <- list(
      base = final_base, consensus = final_consensus,
      preprocessor = full_features$preprocessor
    )

    # Phase 3: transform the untouched screened holdout with development
    # parameters, predict each base method in the 45-feature space, extend
    # co-association to the fixed consensus partition, and diagnose marginal
    # assignment calibration. When the loader sampled a configured cohort, use
    # that same cohort as the development reference population.
    if (!is.null(generalization_data) && nrow(generalization_data)) {
      reference_labels <- final_consensus$labels
      reference_population <- "all development participants"
      if (raw_input && !custom_loader && !is.null(holdout_cohort) &&
          "Overall comments" %in% names(data)) {
        reference_rows <- match(final_consensus$participant_ids, data$CIMBI)
        in_reference <- data[["Overall comments"]][reference_rows] %in% holdout_cohort
        reference_labels <- reference_labels[in_reference]
        reference_population <- paste0(
          "screened development participants with Overall comments in: ",
          paste(holdout_cohort, collapse = ", ")
        )
      }
      generalization <- run_generalization_prediction(
        final_base, final_consensus, full_features$test,
        development_labels = reference_labels,
        method = generalization_method,
        permutations = generalization_permutations,
        seed = seed + 800000L,
        min_methods = max(2L, ceiling(length(final_consensus$active_methods) / 2L)),
        confidence = generalization_confidence,
        exact_max_tables = generalization_exact_max_tables,
        reference_population = reference_population
      )
    }

    # Phase 4: if requested, run blocked permutation tests and participant-level
    # bootstraps for every base method and consensus metric.
    if (resamples > 0L) {
      inference_data <- data
      inference_feature_pipeline <- feature_pipeline
      if (raw_input && !custom_feature_pipeline) {
        # Resampling operates on one consistently derived development table;
        # each resampled CV still refits feature scaling within its training fold.
        inference_raw_fit <- fit_raw_cognitive_preprocessor(data, id_col = "CIMBI")
        inference_data <- apply_raw_cognitive_preprocessor(inference_raw_fit, data)
        # Within each resampled CV fold, refit the prepared-feature transform on
        # that fold's training participants and apply it to its test participants.
        inference_feature_pipeline <- function(train, test) {
          preprocessor <- fit_cognitive_preprocessor(train, id_col = "CIMBI")
          list(
            train = apply_cognitive_preprocessor(preprocessor, train),
            test = apply_cognitive_preprocessor(preprocessor, test),
            preprocessor = preprocessor
          )
        }
      }
      if (is.null(task_blocks)) {
        spec <- cognitive_task_spec()
        # Create the ten shared-donor performance/reaction-time permutation
        # blocks prescribed by the SPEC.
        proposed <- lapply(seq_len(nrow(spec)), function(index) {
          c(spec$perf_var[index], spec$time_var[index])
        })
        task_blocks <- proposed[vapply(proposed, function(block) all(block %in% names(inference_data)), logical(1))]
      }
      if (!length(task_blocks)) {
        stop("SPEC resampling requires task_blocks present in the analysis data.")
      }

      # Re-run the complete CV pipeline for one permuted or bootstrapped data
      # set and return its method-level primary metrics to the inference engine.
      analysis_function <- function(resampled_data, resample_seed) {
        fit <- run_cross_validation(
          data = resampled_data,
          distance_builder = distance_builder,
          predictor_variables = predictor_variables,
          sex_col = sex_col, medication_col = medication_col, age_col = age_col,
          folds = folds, base_config = base_config,
          consensus_config = consensus_config,
          prediction_mode = "feature_consensus",
          seed = resample_seed, feature_pipeline = inference_feature_pipeline,
          age_bin_width = age_bin_width, age_bin_origin = age_bin_origin,
          balance_cols = balance_cols
        )
        fit$summary$metrics
      }
      inference <- run_resampling_inference(
        inference_data, analysis_function, task_blocks,
        iterations = resamples, seed = seed + 900000L, workers = workers,
        worker_source_files = unname(clustering_module_files(normalize = TRUE)),
        worker_required_functions = c(
          "run_cross_validation", "run_base_clustering",
          "run_consensus_clustering", "predict_base_clusters",
          "fit_cognitive_preprocessor"
        ),
        worker_required_packages = resampling_worker_packages(base_config$methods)
      )
    }
  }

  # Assemble a serializable result with explicit run metadata. Compatibility
  # mode leaves SPEC-only final_fit, inference, and generalization entries NULL.
  result <- structure(list(
    cv = cv,
    final_fit = final_fit,
    loaded = loaded,
    inference = inference,
    generalization = generalization,
    metadata = list(
      mode = mode, methods = methods, k_candidates = k_candidates,
      folds = folds, requested_resamples = as.integer(resamples),
      age_bin_width = age_bin_width, age_bin_origin = age_bin_origin,
      holdout_fraction = holdout_fraction, holdout_cohort = holdout_cohort,
      generalization_method = generalization_method,
      generalization_permutations = as.integer(generalization_permutations),
      generalization_confidence = generalization_confidence,
      seed = as.integer(seed), workers = as.integer(workers),
      fit_population = "development",
      row_order_policy = "non-ID non-age record contents; canonical_v1",
      created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    output_paths = NULL
  ), class = c("clustering_analysis_result", "list"))

  # Write artifacts only after successful completion so failed analyses do not
  # leave a misleading partial result set.
  if (isTRUE(write_outputs)) {
    result$output_paths <- write_clustering_artifacts(result, output_dir)
    saveRDS(result, result$output_paths$result_rds)
  }
  result
}
