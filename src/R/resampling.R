# Deterministic blocked permutation and participant-bootstrap inference.

# Shuffle each cognitive performance/reaction-time pair with one shared donor
# permutation, disrupting between-task profiles while preserving each pair.
permute_task_blocks <- function(data, task_blocks, seed) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  variables <- unique(unlist(task_blocks, use.names = FALSE))
  missing <- setdiff(variables, names(data))
  if (length(missing)) stop("Permutation blocks contain missing variables: ", paste(missing, collapse = ", "))
  set.seed(as.integer(seed))
  permuted <- data
  for (block in task_blocks) {
    donors <- sample.int(nrow(data), nrow(data), replace = FALSE)
    permuted[, block] <- data[donors, block, drop = FALSE]
  }
  permuted
}

# Draw one participant-level bootstrap sample with replacement and create unique
# execution IDs while retaining source-participant audit information.
bootstrap_participants <- function(data, seed, id_col = "CIMBI") {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (!id_col %in% names(data)) stop("Bootstrap data lacks ID column: ", id_col)
  set.seed(as.integer(seed))
  indices <- sample.int(nrow(data), nrow(data), replace = TRUE)
  sample <- data[indices, , drop = FALSE]
  source_ids <- as.character(sample[[id_col]])
  sample[[id_col]] <- make.unique(paste0(source_ids, "_bootstrap"))
  rownames(sample) <- NULL
  list(data = sample, source_indices = indices, source_ids = source_ids,
       seed = as.integer(seed))
}

# Identify numeric primary metrics shared by observed, permutation, and
# bootstrap tables while excluding the Monte Carlo iteration key.
metric_columns <- function(data) {
  setdiff(names(data)[vapply(data, is.numeric, logical(1))], "iteration")
}

# Convert null and bootstrap distributions into plus-one empirical p-values and
# 2.5%/97.5% percentile intervals for every method and metric.
summarize_resampling_inference <- function(observed, permutation, bootstrap) {
  required <- "method"
  if (!all(required %in% names(observed)) || !all(required %in% names(permutation)) ||
      !all(required %in% names(bootstrap))) {
    stop("Observed and resampling metrics must contain method.")
  }
  metrics <- Reduce(intersect, list(metric_columns(observed), metric_columns(permutation),
                                    metric_columns(bootstrap)))
  methods <- intersect(observed$method, intersect(permutation$method, bootstrap$method))
  p_rows <- list()
  ci_rows <- list()
  counter <- 0L
  for (method in methods) {
    observed_row <- observed[observed$method == method, , drop = FALSE][1L, ]
    for (metric in metrics) {
      observed_value <- observed_row[[metric]]
      null <- permutation[permutation$method == method, metric]
      null <- null[is.finite(null)]
      boot <- bootstrap[bootstrap$method == method, metric]
      boot <- boot[is.finite(boot)]
      if (!length(null) || !is.finite(observed_value)) {
        p_value <- NA_real_
        tail <- if (identical(metric, "mean_n_clusters")) "two_sided_null_departure" else "upper"
      } else if (identical(metric, "mean_n_clusters")) {
        center <- stats::median(null)
        p_value <- (1 + sum(abs(null - center) >= abs(observed_value - center))) /
          (length(null) + 1)
        tail <- "two_sided_null_departure"
      } else {
        p_value <- (1 + sum(null >= observed_value)) / (length(null) + 1)
        tail <- "upper"
      }
      counter <- counter + 1L
      p_rows[[counter]] <- data.frame(
        method = method, metric = metric, observed = observed_value,
        null_n = length(null), tail = tail, p_value = p_value,
        stringsAsFactors = FALSE
      )
      interval <- if (length(boot)) {
        stats::quantile(boot, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      } else c(NA_real_, NA_real_)
      ci_rows[[counter]] <- data.frame(
        method = method, metric = metric, bootstrap_n = length(boot),
        lower = interval[1L], upper = interval[2L],
        stringsAsFactors = FALSE
      )
    }
  }
  list(p_values = do.call(rbind, p_rows),
       confidence_intervals = do.call(rbind, ci_rows))
}

# Run the observed CV metrics plus deterministic blocked permutations and
# participant bootstraps, sequentially by default or on optional workers.
# `analysis_function(data, seed)` is the caller-supplied complete CV callback.
initialize_resampling_workers <- function(cluster, source_files = character(),
                                          library_paths = .libPaths()) {
  source_files <- as.character(source_files)
  missing_files <- source_files[!file.exists(source_files)]
  if (length(missing_files)) {
    stop("Cannot initialize resampling workers; missing source files: ",
         paste(missing_files, collapse = ", "))
  }
  parallel::clusterCall(
    cluster,
    function(files, libraries) {
      .libPaths(unique(c(libraries, .libPaths())))
      for (path in files) sys.source(path, envir = .GlobalEnv)
      TRUE
    },
    source_files, as.character(library_paths)
  )
}

preflight_resampling_workers <- function(cluster, required_functions = character(),
                                         required_packages = character()) {
  checks <- parallel::clusterCall(
    cluster,
    function(function_names, package_names) {
      missing_functions <- function_names[!vapply(
        function_names, exists, logical(1), mode = "function", inherits = TRUE
      )]
      missing_packages <- package_names[!vapply(
        package_names, requireNamespace, logical(1), quietly = TRUE
      )]
      list(
        missing_functions = missing_functions,
        missing_packages = missing_packages,
        library_paths = .libPaths()
      )
    },
    unique(as.character(required_functions)),
    unique(as.character(required_packages))
  )
  failures <- vapply(checks, function(check) {
    length(check$missing_functions) || length(check$missing_packages)
  }, logical(1))
  if (any(failures)) {
    details <- vapply(which(failures), function(index) {
      check <- checks[[index]]
      paste0(
        "worker ", index,
        if (length(check$missing_functions)) paste0(
          " missing functions [", paste(check$missing_functions, collapse = ", "), "]"
        ) else "",
        if (length(check$missing_packages)) paste0(
          " missing packages [", paste(check$missing_packages, collapse = ", "), "]"
        ) else "",
        "; .libPaths()=", paste(check$library_paths, collapse = ";")
      )
    }, character(1))
    stop("Parallel resampling worker preflight failed: ", paste(details, collapse = " | "))
  }
  invisible(checks)
}

resampling_iteration_seeds <- function(seed, index) {
  seed <- as.double(seed)
  index <- as.double(index)
  if (length(seed) != 1L || !is.finite(seed) || seed != trunc(seed)) {
    stop("seed must be one finite integer.")
  }
  if (length(index) != 1L || !is.finite(index) || index < 1 || index != trunc(index)) {
    stop("resampling iteration index must be a positive integer.")
  }
  seeds <- seed + index * 1000 + seq_len(4L)
  if (any(seeds > .Machine$integer.max) || any(seeds < -.Machine$integer.max)) {
    stop("Derived resampling seed exceeds R's integer range.")
  }
  stats::setNames(
    as.integer(seeds),
    c("permutation_data", "bootstrap_data", "permutation_analysis", "bootstrap_analysis")
  )
}

run_resampling_inference <- function(data, analysis_function, task_blocks,
                                     iterations = 1000L, seed = 1L,
                                     workers = 1L,
                                     worker_source_files = character(),
                                     worker_required_functions = character(),
                                     worker_required_packages = character()) {
  iterations <- as.integer(iterations)
  workers <- as.integer(workers)
  if (iterations < 1L) stop("iterations must be positive.")
  if (workers < 1L) stop("workers must be positive.")
  if (!is.function(analysis_function)) stop("analysis_function must be a function.")
  if (!is.list(task_blocks) || !length(task_blocks)) {
    stop("task_blocks must be a non-empty list.")
  }
  seed <- as.integer(seed)
  if (length(seed) != 1L || is.na(seed)) stop("seed must be one finite integer.")
  observed <- analysis_function(data, seed)
  # Keep task-specific state and the caller callback inside the serialized
  # closure. Workflow functions referenced globally by the callback are loaded
  # separately during PSOCK worker initialization.
  permute_function <- permute_task_blocks
  bootstrap_function <- bootstrap_participants
  seed_function <- resampling_iteration_seeds
  # Execute one paired permutation/bootstrap iteration with non-overlapping,
  # reproducible seeds and attach its iteration key to every metric row.
  run_iteration <- function(index) {
    iteration_seeds <- seed_function(seed, index)
    permuted <- permute_function(data, task_blocks, iteration_seeds[["permutation_data"]])
    bootstrap <- bootstrap_function(data, iteration_seeds[["bootstrap_data"]])$data
    permutation_metrics <- analysis_function(
      permuted, iteration_seeds[["permutation_analysis"]]
    )
    bootstrap_metrics <- analysis_function(
      bootstrap, iteration_seeds[["bootstrap_analysis"]]
    )
    permutation_metrics$iteration <- index
    bootstrap_metrics$iteration <- index
    list(permutation = permutation_metrics, bootstrap = bootstrap_metrics)
  }
  if (workers == 1L) {
    results <- lapply(seq_len(iterations), run_iteration)
  } else {
    cluster <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    initialize_resampling_workers(cluster, worker_source_files, .libPaths())
    preflight_resampling_workers(
      cluster, worker_required_functions, worker_required_packages
    )
    results <- parallel::parLapplyLB(cluster, seq_len(iterations), run_iteration)
  }
  permutation <- do.call(rbind, lapply(results, `[[`, "permutation"))
  bootstrap <- do.call(rbind, lapply(results, `[[`, "bootstrap"))
  list(
    observed = observed,
    permutation = permutation,
    bootstrap = bootstrap,
    inference = summarize_resampling_inference(observed, permutation, bootstrap),
    iterations = iterations,
    seed = as.integer(seed),
    workers = workers
  )
}
