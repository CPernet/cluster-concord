# Fixed-development-partition generalization calibration.

# Return conservative simultaneous Wilson intervals for multinomial cell
# proportions. Bonferroni adjustment gives at least the requested family-wise
# coverage while retaining sensible limits for zero-count cells.
simultaneous_wilson_intervals <- function(counts, total, confidence = 0.95) {
  counts <- as.numeric(counts)
  cells <- length(counts)
  if (!length(total) || total <= 0L || cells < 1L) {
    return(data.frame(lower = rep(NA_real_, cells), upper = rep(NA_real_, cells)))
  }
  alpha <- 1 - confidence
  z <- stats::qnorm(1 - alpha / (2 * cells))
  estimate <- counts / total
  denominator <- 1 + z^2 / total
  center <- (estimate + z^2 / (2 * total)) / denominator
  half_width <- z * sqrt(
    estimate * (1 - estimate) / total + z^2 / (4 * total^2)
  ) / denominator
  data.frame(
    lower = pmax(0, center - half_width),
    upper = pmin(1, center + half_width)
  )
}

# Enumerate every allocation of `total` observations among `cells` categories.
# This is used only when the number of weak compositions is below a configured
# limit, so the returned matrix remains bounded in size.
enumerate_multinomial_counts <- function(total, cells) {
  total <- as.integer(total)
  cells <- as.integer(cells)
  if (cells < 1L || total < 0L) stop("Invalid multinomial dimensions.")
  if (cells == 1L) return(matrix(total, nrow = 1L))
  bars <- utils::combn(total + cells - 1L, cells - 1L)
  apply(bars, 2L, function(value) {
    diff(c(0L, value, total + cells)) - 1L
  })
}

# Compare assigned holdout proportions with a prespecified development
# reference population. This is a calibration/distribution-shift diagnostic;
# it does not measure subject-level assignment correctness.
evaluate_generalization <- function(
  development_labels,
  holdout_labels,
  method = c("auto", "exact", "monte_carlo"),
  permutations = 9999L,
  seed = 1L,
  confidence = 0.95,
  exact_max_tables = 1000000L,
  reference_population = "development"
) {
  method <- match.arg(method)
  permutations <- as.integer(permutations)
  exact_max_tables <- as.numeric(exact_max_tables)
  if (permutations < 1L && identical(method, "monte_carlo")) {
    stop("Monte Carlo calibration requires at least one replicate.")
  }
  if (!is.numeric(confidence) || length(confidence) != 1L ||
      !is.finite(confidence) || confidence <= 0 || confidence >= 1) {
    stop("confidence must be a single number strictly between 0 and 1.")
  }

  development_labels <- as.character(development_labels)
  holdout_labels <- as.character(holdout_labels)
  development_labels <- development_labels[!is.na(development_labels)]
  levels <- sort(unique(development_labels))
  if (length(levels) < 2L) stop("Development labels require at least two clusters.")
  development_counts <- table(factor(development_labels, levels = levels))
  development_proportions <- as.numeric(development_counts) / sum(development_counts)
  assigned <- holdout_labels[!is.na(holdout_labels)]
  holdout_counts <- table(factor(assigned, levels = levels))
  expected <- length(assigned) * development_proportions
  statistic <- if (length(assigned) && all(expected > 0)) {
    sum((as.numeric(holdout_counts) - expected)^2 / expected)
  } else NA_real_

  possible_tables <- if (length(assigned)) {
    choose(length(assigned) + length(levels) - 1L, length(levels) - 1L)
  } else 0
  resolved_method <- method
  if (identical(method, "auto")) {
    resolved_method <- if (is.finite(possible_tables) &&
                            possible_tables <= exact_max_tables) {
      "exact"
    } else {
      "monte_carlo"
    }
  }

  null_statistics <- numeric()
  exceedances <- NA_integer_
  monte_carlo_se <- NA_real_
  p_resolution <- NA_real_
  p_value <- NA_real_
  if (is.finite(statistic) && identical(resolved_method, "exact")) {
    count_tables <- enumerate_multinomial_counts(length(assigned), length(levels))
    probabilities <- apply(
      count_tables, 2L, stats::dmultinom, prob = development_proportions
    )
    exact_statistics <- colSums((count_tables - expected)^2 / expected)
    # The tolerance prevents floating-point ties at the observed statistic from
    # being omitted from the upper-tail probability.
    p_value <- sum(probabilities[exact_statistics >= statistic - 1e-12])
  } else if (is.finite(statistic) && identical(resolved_method, "monte_carlo")) {
    set.seed(as.integer(seed))
    null_statistics <- replicate(permutations, {
      simulated <- as.numeric(stats::rmultinom(
        1L, length(assigned), development_proportions
      ))
      sum((simulated - expected)^2 / expected)
    })
    exceedances <- sum(null_statistics >= statistic)
    p_value <- (exceedances + 1) / (permutations + 1)
    monte_carlo_se <- sqrt(p_value * (1 - p_value) / (permutations + 1))
    p_resolution <- 1 / (permutations + 1)
  }

  intervals <- simultaneous_wilson_intervals(
    holdout_counts, length(assigned), confidence = confidence
  )
  list(
    proportions = data.frame(
      cluster = levels,
      reference_count = as.integer(development_counts),
      reference_proportion = development_proportions,
      holdout_count = as.integer(holdout_counts),
      holdout_proportion = if (length(assigned)) {
        as.numeric(holdout_counts) / length(assigned)
      } else {
        NA_real_
      },
      holdout_ci_lower = intervals$lower,
      holdout_ci_upper = intervals$upper,
      stringsAsFactors = FALSE
    ),
    coverage = length(assigned) / length(holdout_labels),
    assigned_n = length(assigned),
    excluded_n = sum(is.na(holdout_labels)),
    reference_n = length(development_labels),
    reference_population = as.character(reference_population),
    confidence_level = confidence,
    test = list(
      method = paste0("multinomial_", resolved_method, "_calibration"),
      requested_method = method,
      statistic = statistic,
      p_value = p_value,
      possible_count_tables = possible_tables,
      permutations = if (identical(resolved_method, "monte_carlo")) permutations else 0L,
      exceedances = exceedances,
      monte_carlo_se = monte_carlo_se,
      p_resolution = p_resolution,
      seed = if (identical(resolved_method, "monte_carlo")) as.integer(seed) else NA_integer_,
      interpretation = paste(
        "Tests marginal assignment proportions against the fixed development",
        "reference; it does not test subject-level assignment correctness or",
        "independent recovery of the cluster structure."
      )
    ),
    null_statistics = null_statistics
  )
}

# Predict every base method and the fixed consensus partition for the holdout,
# retain rejection diagnostics, and evaluate marginal assignment calibration.
run_generalization_prediction <- function(
  base_result,
  consensus_result,
  holdout_distance_result,
  development_labels = consensus_result$labels,
  method = c("auto", "exact", "monte_carlo"),
  permutations = 9999L,
  seed = 1L,
  min_methods = 2L,
  confidence = 0.95,
  exact_max_tables = 1000000L,
  reference_population = "development"
) {
  method <- match.arg(method)
  base_prediction <- predict_base_clusters(base_result, holdout_distance_result)
  consensus_prediction <- predict_consensus_cluster(
    consensus_result, base_result, base_prediction,
    min_methods = min_methods, min_score = 0, min_margin = 0
  )
  predicted_labels <- consensus_prediction$native_label
  predicted_labels[consensus_prediction$diagnostics$rejected] <- NA_character_
  evaluation <- evaluate_generalization(
    development_labels = development_labels,
    holdout_labels = predicted_labels,
    method = method,
    permutations = permutations,
    seed = seed,
    confidence = confidence,
    exact_max_tables = exact_max_tables,
    reference_population = reference_population
  )
  evaluation$base_prediction <- base_prediction
  evaluation$consensus_prediction <- consensus_prediction
  evaluation
}
