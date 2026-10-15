
# Stable internal ordering uses non-ID record contents, never spreadsheet rows.
# Age is excluded so future age recoding does not alter the ordering key.
# Keys are transient, reversible encodings, not anonymization tokens or outputs.
canonical_record_keys <- function(data, id_col = "CIMBI",
                                  age_col = "Age at neuropsych test") {
  columns <- sort(setdiff(names(data), c(id_col, "ID", "CIMBI ID", age_col)),
                  method = "radix")
  if (!length(columns)) stop("Canonical ordering requires non-ID, non-age data.")
  values <- lapply(data[, columns, drop = FALSE], function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (inherits(x, "Date")) x <- as.character(x)
    if (inherits(x, "POSIXt")) x <- format(x, tz = "UTC", usetz = TRUE)
    if (is.numeric(x)) {
      x <- as.double(x)
      x[is.nan(x)] <- NA_real_
      x[!is.na(x) & x == 0] <- 0
    }
    if (is.character(x)) x <- enc2utf8(x)
    unname(x)
  })
  vapply(seq_len(nrow(data)), function(i) {
    row <- lapply(values, function(x) x[[i]])
    paste(sprintf("%02x", as.integer(serialize(row, NULL, version = 2))),
          collapse = "")
  }, character(1))
}

canonical_record_order <- function(data, id_col = "CIMBI",
                                   age_col = "Age at neuropsych test",
                                   require_distinct = FALSE) {
  keys <- canonical_record_keys(data, id_col, age_col)
  if (require_distinct && anyDuplicated(keys)) {
    stop("Records indistinguishable without ID and age require an independent persistent analysis token; no rows were removed.")
  }
  order(keys, method = "radix")
}

canonicalize_records <- function(data, id_col = "CIMBI",
                                 age_col = "Age at neuropsych test",
                                 require_distinct = FALSE) {
  out <- data[canonical_record_order(data, id_col, age_col, require_distinct),
              , drop = FALSE]
  rownames(out) <- NULL
  out
}

analysis_stage_seed <- function(seed, offset) {
  if (length(seed) != 1L || !is.finite(seed)) stop("seed must be finite.")
  as.integer((as.double(seed) + as.double(offset)) %% .Machine$integer.max)
}

# Load and prepare cognitive data from the raw neuropsychology workbook.
# Example:
# prepared <- load_and_prepdata(
#   "data/data.xlsx",
# )

# Return the complete raw workbook columns eligible for cognitive imputation and
# downstream derivation under the study protocol.
cognitive_imputation_columns <- function() {
  c(
    "SRT - Total mean reaction latency",
    "Trail making B - # of sec",
    "Trail making A - # of sec",
    "SDMT total score",
    "Letter-number sequencing total score",
    "D-KEFS (STROOP) - Color-word interference condition 4 raw score",
    "RAVLT - correct A1-A5",
    "RAVLT - correct A6 (immediate recall)",
    "RAVLT - correct A7 (delayed recall)",
    "Emoticom-ER - Hit rate for Eyes-Happy",
    "Emoticom-ER - Hit rate for Eyes-Sad",
    "Emoticom-ER - False alarm rate for Eyes-Happy",
    "Emoticom-ER - False alarm rate for Eyes-Sad",
    "Emoticom-ER - Hit rate for Eyes-Anger",
    "Emoticom-ER - Hit rate for Eyes-Fear",
    "Emoticom-ER - False alarm rate for Eyes-Anger",
    "Emoticom-ER - False alarm rate for Eyes-Fear",
    "Emoticom-ER - Mean rt (s) for Eyes-Happy hits",
    "Emoticom-ER - Mean rt (s) for Eyes-Sad hits",
    "Emoticom-ER - Mean rt (s) for Eyes-Anger hits",
    "Emoticom-ER - Mean rt (s) for Eyes-Fear hits",
    "Emoticom-ER - Mean rt (s) for Eyes-Anger false alarms",
    "Emoticom-ER - Mean rt (s) for Eyes-Fear false alarms",
    "Emoticom-ER - Mean rt (s) for Eyes-Happy false alarms",
    "Emoticom-ER - Mean rt (s) for Eyes-Sad false alarms",
    "Mean intensity threshold - Sad (ascending)",
    "Mean intensity threshold - Happy (ascending)",
    "Mean intensity threshold - Sad (descending)",
    "Mean intensity threshold - Happy (descending)",
    "Mean intensity threshold - Anger (ascending)",
    "Mean intensity threshold - Anger (descending)",
    "Mean intensity threshold - Fear (ascending)",
    "Mean intensity threshold - Fear (descending)",
    "Mean intensity threshold - Disgust (ascending)",
    "Mean intensity threshold - Disgust (descending)",
    "Emoticom-MJ - Mean rating (1-7), Guilty, Victim Intentional",
    "Emoticom-MJ - Mean rating (1-7), Ashamed, Victim Intentional",
    "Emoticom-MJ - Mean rating (1-7), Guilty, Agent Intentional",
    "Emoticom-MJ - Mean rating (1-7), Ashamed, Agent Intentional",
    "Emoticom-MJ - Mean rating (1-7), Guilty, Victim Unintentional",
    "Emoticom-MJ - Mean rating (1-7), Ashamed, Victim Unintentional",
    "Emoticom-MJ - Mean rating (1-7), Guilty, Agent Unintentional",
    "Emoticom-MJ - Mean rating (1-7), Ashamed, Agent Unintentional",
    "Emoticom-MJ - Mean rt (sec), Guilty, Victim Intentional",
    "Emoticom-MJ - Mean rt (sec), Ashamed, Victim Intentional",
    "Emoticom-MJ - Mean rt (sec), Guilty, Agent Intentional",
    "Emoticom-MJ - Mean rt (sec), Ashamed, Agent Intentional",
    "Emoticom-MJ - Mean rt (sec), Guilty, Victim Unintentional",
    "Emoticom-MJ - Mean rt (sec), Ashamed, Victim Unintentional",
    "Emoticom-MJ - Mean rt (sec), Guilty, Agent Unintentional",
    "Emoticom-MJ - Mean rt (sec), Ashamed, Agent Unintentional"
  )
}

# Group raw variables into the ten cognitive tests used by the participant-level
# missing-test screening rule.
cognitive_missingness_groups <- function() {
  list(
    SRT = "SRT - Total mean reaction latency",
    TMT = c("Trail making B - # of sec", "Trail making A - # of sec"),
    SDMT = "SDMT total score",
    LNS = "Letter-number sequencing total score",
    Stroop = "D-KEFS (STROOP) - Color-word interference condition 4 raw score",
    RAVLT = c(
      "RAVLT - correct A1-A5",
      "RAVLT - correct A6 (immediate recall)",
      "RAVLT - correct A7 (delayed recall)"
    ),
    ERT = c(
      "Emoticom-ER - Hit rate for Eyes-Happy",
      "Emoticom-ER - Hit rate for Eyes-Sad",
      "Emoticom-ER - False alarm rate for Eyes-Happy",
      "Emoticom-ER - False alarm rate for Eyes-Sad",
      "Emoticom-ER - Hit rate for Eyes-Anger",
      "Emoticom-ER - Hit rate for Eyes-Fear",
      "Emoticom-ER - False alarm rate for Eyes-Anger",
      "Emoticom-ER - False alarm rate for Eyes-Fear",
      "Emoticom-ER - Mean rt (s) for Eyes-Happy hits",
      "Emoticom-ER - Mean rt (s) for Eyes-Sad hits",
      "Emoticom-ER - Mean rt (s) for Eyes-Anger hits",
      "Emoticom-ER - Mean rt (s) for Eyes-Fear hits",
      "Emoticom-ER - Mean rt (s) for Eyes-Anger false alarms",
      "Emoticom-ER - Mean rt (s) for Eyes-Fear false alarms",
      "Emoticom-ER - Mean rt (s) for Eyes-Happy false alarms",
      "Emoticom-ER - Mean rt (s) for Eyes-Sad false alarms"
    ),
    IMT = c(
      "Mean intensity threshold - Sad (ascending)",
      "Mean intensity threshold - Happy (ascending)",
      "Mean intensity threshold - Sad (descending)",
      "Mean intensity threshold - Happy (descending)",
      "Mean intensity threshold - Anger (ascending)",
      "Mean intensity threshold - Anger (descending)",
      "Mean intensity threshold - Fear (ascending)",
      "Mean intensity threshold - Fear (descending)",
      "Mean intensity threshold - Disgust (ascending)",
      "Mean intensity threshold - Disgust (descending)"
    ),
    MET = c(
      "Emoticom-MJ - Mean rating (1-7), Guilty, Victim Intentional",
      "Emoticom-MJ - Mean rating (1-7), Ashamed, Victim Intentional",
      "Emoticom-MJ - Mean rating (1-7), Guilty, Agent Intentional",
      "Emoticom-MJ - Mean rating (1-7), Ashamed, Agent Intentional",
      "Emoticom-MJ - Mean rating (1-7), Guilty, Victim Unintentional",
      "Emoticom-MJ - Mean rating (1-7), Ashamed, Victim Unintentional",
      "Emoticom-MJ - Mean rating (1-7), Guilty, Agent Unintentional",
      "Emoticom-MJ - Mean rating (1-7), Ashamed, Agent Unintentional",
      "Emoticom-MJ - Mean rt (sec), Guilty, Victim Intentional",
      "Emoticom-MJ - Mean rt (sec), Ashamed, Victim Intentional",
      "Emoticom-MJ - Mean rt (sec), Guilty, Agent Intentional",
      "Emoticom-MJ - Mean rt (sec), Ashamed, Agent Intentional",
      "Emoticom-MJ - Mean rt (sec), Guilty, Victim Unintentional",
      "Emoticom-MJ - Mean rt (sec), Ashamed, Victim Unintentional",
      "Emoticom-MJ - Mean rt (sec), Guilty, Agent Unintentional",
      "Emoticom-MJ - Mean rt (sec), Ashamed, Agent Unintentional"
    )
  )
}

# Fail early with a complete list of unavailable runtime packages required by
# loading, preparation, or output serialization.
require_cognitive_packages <- function(packages) {
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Required package(s) are not installed: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }
}

# Enforce an explicit required-column contract and report all missing columns in
# the caller-supplied context.
assert_required_columns <- function(data, columns, context) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      context,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

# Exclude participants missing more than the configured number of entire tests,
# before learned imputation or feature transformation.
filter_missing_tests <- function(data, groups, max_missing_tests = 3L) {
  missing_by_group <- vapply(
    groups,
    # Mark one cognitive test missing only when all of its source variables are NA.
    function(columns) rowSums(!is.na(data[, columns, drop = FALSE])) == 0,
    logical(nrow(data))
  )
  if (is.null(dim(missing_by_group))) {
    missing_by_group <- matrix(missing_by_group, ncol = 1L)
  }
  keep <- rowSums(missing_by_group) <= max_missing_tests
  data[keep, , drop = FALSE]
}

# Perform the audited historical chained-GAM imputation used for compatibility
# preparation, with deterministic seeds and bounded stochastic draws.
gam_impute_chained <- function(data, columns, maxit = 3L, k = 10L, seed = 2025L) {
  set.seed(seed)
  data <- as.data.frame(data)

  if (!is.numeric(maxit) || length(maxit) != 1L || maxit < 1L) {
    stop("maxit must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || k < 3L) {
    stop("k must be an integer of at least 3.", call. = FALSE)
  }

  for (column in columns) {
    data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
  }

  # Preserve the original missingness mask for one imputation variable.
  was_missing <- lapply(columns, function(column) is.na(data[[column]]))
  names(was_missing) <- columns
  audit <- data.frame(
    variable = columns,
    n_missing = vapply(was_missing, sum, integer(1)),
    gam_fits = integer(length(columns)),
    gam_fit_failures = integer(length(columns)),
    predictive_draws = integer(length(columns)),
    residual_bootstraps = integer(length(columns)),
    observed_resamples = integer(length(columns)),
    stringsAsFactors = FALSE
  )

  for (column in columns) {
    missing_index <- which(was_missing[[column]])
    observed <- data[[column]][!was_missing[[column]]]
    if (length(missing_index) > 0 && length(observed) == 0) {
      stop("Cannot impute an entirely missing variable: ", column, call. = FALSE)
    }
    if (length(missing_index) > 0 && length(observed) > 0) {
      data[[column]][missing_index] <- sample(observed, length(missing_index), replace = TRUE)
    }
  }

  # Quote one workbook column safely for dynamically constructed GAM formulas.
  quote_name <- function(x) paste0("`", x, "`")
  for (iteration in seq_len(maxit)) {
    for (target in columns) {
      missing_index <- which(was_missing[[target]])
      fit_index <- which(!was_missing[[target]])
      audit_row <- match(target, audit$variable)
      if (length(missing_index) == 0) {
        next
      }

      observed_target <- data[[target]][fit_index]
      # Fall back to observed target resampling when a valid GAM prediction
      # cannot be fitted; increment the corresponding audit counter.
      assign_observed_resample <- function() {
        data[missing_index, target] <<- sample(observed_target, length(missing_index), replace = TRUE)
        audit$observed_resamples[audit_row] <<- audit$observed_resamples[audit_row] + length(missing_index)
      }

      if (length(fit_index) < 10L) {
        assign_observed_resample()
        next
      }

      predictors <- setdiff(columns, target)
      formula_rhs <- paste0("s(", quote_name(predictors), ", k=", k, ")", collapse = " + ")
      model_formula <- stats::as.formula(paste(quote_name(target), "~", formula_rhs))
      fit <- try(
        mgcv::gam(
          model_formula,
          data = data[fit_index, , drop = FALSE],
          family = stats::gaussian(),
          select = TRUE
        ),
        silent = TRUE
      )
      if (inherits(fit, "try-error")) {
        audit$gam_fit_failures[audit_row] <- audit$gam_fit_failures[audit_row] + 1L
        assign_observed_resample()
        next
      }
      audit$gam_fits[audit_row] <- audit$gam_fits[audit_row] + 1L

      prediction <- try(
        stats::predict(
          fit,
          newdata = data[missing_index, , drop = FALSE],
          type = "response",
          se.fit = TRUE
        ),
        silent = TRUE
      )

      prediction_is_usable <- !inherits(prediction, "try-error") &&
        is.list(prediction) &&
        length(prediction$fit) == length(missing_index) &&
        length(prediction$se.fit) == length(missing_index) &&
        all(is.finite(prediction$fit)) &&
        all(is.finite(prediction$se.fit))
      lower <- min(observed_target, na.rm = TRUE)
      upper <- max(observed_target, na.rm = TRUE)

      if (prediction_is_usable) {
        sigma2 <- try(summary(fit)$scale, silent = TRUE)
        if (inherits(sigma2, "try-error") || !is.finite(sigma2)) {
          sigma2 <- stats::var(stats::residuals(fit, type = "response"), na.rm = TRUE)
        }
        sigma2 <- if (is.finite(sigma2) && sigma2 > 0) sigma2 else 0
        sd_prediction <- sqrt(pmax(0, sigma2 + as.numeric(prediction$se.fit)^2))
        sd_cap <- stats::quantile(sd_prediction, 0.95, na.rm = TRUE)
        if (is.finite(sd_cap)) {
          sd_prediction <- pmin(sd_prediction, sd_cap)
        }
        draw <- stats::rnorm(length(prediction$fit), 0, sd_prediction)
        imputed <- as.numeric(prediction$fit) + draw
        imputed <- pmin(pmax(imputed, lower), upper)
        data[missing_index, target] <- imputed
        audit$predictive_draws[audit_row] <- audit$predictive_draws[audit_row] + length(missing_index)
      } else {
        predicted_mean <- try(stats::predict(
          fit,
          newdata = data[missing_index, , drop = FALSE],
          type = "response"
        ), silent = TRUE)
        residuals <- stats::residuals(fit, type = "response")
        residuals <- residuals[is.finite(residuals)]
        if (!inherits(predicted_mean, "try-error") &&
            length(predicted_mean) == length(missing_index) &&
            all(is.finite(predicted_mean)) &&
            length(residuals) > 0) {
          imputed <- as.numeric(predicted_mean) + sample(residuals, length(predicted_mean), replace = TRUE)
          data[missing_index, target] <- pmin(pmax(imputed, lower), upper)
          audit$residual_bootstraps[audit_row] <- audit$residual_bootstraps[audit_row] + length(missing_index)
        } else {
          assign_observed_resample()
        }
      }
    }
  }

  list(data = data, was_missing = was_missing, audit = audit)
}

# A deterministic train/apply imputer for outer-CV and generalization splits.
# The historical GAM routine remains available for compatibility preparation,
# but cannot validly be refitted on validation observations.
fit_raw_cognitive_preprocessor <- function(data, id_col = "CIMBI",
                                            columns = cognitive_imputation_columns()) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  assert_required_columns(data, c(id_col, columns), "Raw cognitive training data")
  # Learn one deterministic training median per raw variable; held-out values
  # never contribute to these parameters.
  values <- vapply(columns, function(column) {
    numeric_values <- suppressWarnings(as.numeric(data[[column]]))
    value <- stats::median(numeric_values, na.rm = TRUE)
    if (!is.finite(value)) {
      stop("Cannot fit raw preprocessing: entirely missing variable: ", column,
           call. = FALSE)
    }
    value
  }, numeric(1))
  structure(list(
    id_col = id_col,
    columns = columns,
    imputation_values = values,
    fit_ids = as.character(data[[id_col]]),
    imputation_method = "training_median"
  ), class = c("raw_cognitive_preprocessor", "list"))
}

# Apply fixed training-derived raw imputation values, then derive the configured
# cognitive measures without learning from validation/generalization rows.
apply_raw_cognitive_preprocessor <- function(preprocessor, data) {
  if (!inherits(preprocessor, "raw_cognitive_preprocessor")) {
    stop("preprocessor must be a raw_cognitive_preprocessor.", call. = FALSE)
  }
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  assert_required_columns(
    data, c(preprocessor$id_col, preprocessor$columns), "Raw cognitive apply data"
  )
  for (column in preprocessor$columns) {
    values <- suppressWarnings(as.numeric(data[[column]]))
    values[!is.finite(values)] <- preprocessor$imputation_values[[column]]
    data[[column]] <- values
  }
  derive_cognitive_measures(data, id_col = preprocessor$id_col)
}

# Convert one emotion's hit/false-alarm rates to counts and compute signal-
# detection indices used by the ERT performance measure.
compute_sdt_emotion <- function(data, emotion_label, emotion_name, id_col = "CIMBI") {
  hit_column <- paste0("Emoticom-ER - Hit rate for Eyes-", emotion_name)
  false_alarm_column <- paste0("Emoticom-ER - False alarm rate for Eyes-", emotion_name)
  hit_count_column <- paste0("ERT_n_hit_", tolower(emotion_label))
  false_alarm_count_column <- paste0("ERT_n_false_", tolower(emotion_label))
  miss_count_column <- paste0("ERT_n_miss_", tolower(emotion_label))
  correct_rejection_column <- paste0("ERT_n_corr_", tolower(emotion_label))

  data[[hit_column]] <- as.numeric(data[[hit_column]])
  data[[false_alarm_column]] <- as.numeric(data[[false_alarm_column]])
  data[[hit_count_column]] <- as.integer(round(data[[hit_column]] * 20L))
  data[[false_alarm_count_column]] <- as.integer(round(data[[false_alarm_column]] * 60L))
  data[[miss_count_column]] <- 20L - data[[hit_count_column]]
  data[[correct_rejection_column]] <- 60L - data[[false_alarm_count_column]]

  count_columns <- c(
    hit_count_column,
    false_alarm_count_column,
    miss_count_column,
    correct_rejection_column
  )
  keep <- stats::complete.cases(data[, count_columns, drop = FALSE])
  indices <- psycho::dprime(
    n_hit = data[[hit_count_column]][keep],
    n_fa = data[[false_alarm_count_column]][keep],
    n_miss = data[[miss_count_column]][keep],
    n_cr = data[[correct_rejection_column]][keep]
  )
  indices <- as.data.frame(indices, stringsAsFactors = FALSE)
  names(indices)[names(indices) == "dprime"] <- paste0("SDT_", emotion_label, "_dprime")
  names(indices)[names(indices) == "beta"] <- paste0("SDT_", emotion_label, "_beta")
  names(indices)[names(indices) == "aprime"] <- paste0("SDT_", emotion_label, "_aprime")
  names(indices)[names(indices) == "bppd"] <- paste0("SDT_", emotion_label, "_bppd")
  names(indices)[names(indices) == "c"] <- paste0("SDT_", emotion_label, "_c")

  list(
    data = data,
    indices = cbind(data[keep, id_col, drop = FALSE], indices)
  )
}

# Compute and ID-join signal-detection indices for all four configured emotions.
add_sdt_indices <- function(data, id_col = "CIMBI") {
  emotion_map <- c(Happy = "Happy", Sad = "Sad", Angry = "Anger", Fear = "Fear")
  indices <- list()
  for (emotion_label in names(emotion_map)) {
    result <- compute_sdt_emotion(data, emotion_label, emotion_map[[emotion_label]], id_col)
    data <- result$data
    indices[[emotion_label]] <- result$indices
  }

  for (emotion_label in names(indices)) {
    index_data <- indices[[emotion_label]]
    new_columns <- setdiff(names(index_data), id_col)
    matching_rows <- match(data[[id_col]], index_data[[id_col]])
    for (column in new_columns) {
      data[[column]] <- index_data[[column]][matching_rows]
    }
  }
  data
}

# Rename one required raw source column to its canonical derived-measure name,
# failing clearly if the source is absent.
rename_column <- function(data, source, target) {
  if (!(source %in% names(data))) {
    stop("Cannot derive ", target, ": source column is missing: ", source, call. = FALSE)
  }
  names(data)[names(data) == source] <- target
  data
}

# Derive the ten cognitive task performance/reaction-time measures from imputed
# raw workbook columns according to the fixed study definitions.
derive_cognitive_measures <- function(data, id_col = "CIMBI") {
  data <- rename_column(data, "SRT - Total mean reaction latency", "SRT")
  data$log_a <- log(data$SRT)
  data$SRT_scale <- 0
  data$TMT_Cor_B <- data[["Trail making B - # of sec"]] - data[["Trail making A - # of sec"]]
  data$TMT_scale <- 0
  data <- rename_column(data, "SDMT total score", "SDMT")
  data$SDMT_time <- 0
  data <- rename_column(data, "Letter-number sequencing total score", "LNS")
  data$LNS_time <- 0
  data$Stroop_time <- data[["D-KEFS (STROOP) - Color-word interference condition 4 raw score"]]
  data$Stroop_scale <- 0
  data$RAVLT <- (
    data[["RAVLT - correct A1-A5"]] +
      data[["RAVLT - correct A6 (immediate recall)"]] +
      data[["RAVLT - correct A7 (delayed recall)"]]
  ) / 7
  data$RAVLT_time <- 0

  data$ERT_hitrate_HS <- data[["Emoticom-ER - Hit rate for Eyes-Happy"]] -
    data[["Emoticom-ER - Hit rate for Eyes-Sad"]]
  data$ERT_falsealarm_HS <- data[["Emoticom-ER - False alarm rate for Eyes-Happy"]] -
    data[["Emoticom-ER - False alarm rate for Eyes-Sad"]]
  data$ERT_hitrate_HS_time <- (
    data[["Emoticom-ER - Mean rt (s) for Eyes-Happy hits"]] +
      data[["Emoticom-ER - Mean rt (s) for Eyes-Sad hits"]]
  ) / 2
  data$ERT_falsealarm_HS_time <- (
    data[["Emoticom-ER - Mean rt (s) for Eyes-Happy false alarms"]] +
      data[["Emoticom-ER - Mean rt (s) for Eyes-Sad false alarms"]]
  ) / 2
  data$ERT_time_avg_HS <- (data$ERT_hitrate_HS_time + data$ERT_falsealarm_HS_time) / 2
  data <- add_sdt_indices(data, id_col)
  data$ERT_hitrate_HSAF_time <- (
    data[["Emoticom-ER - Mean rt (s) for Eyes-Happy hits"]] +
      data[["Emoticom-ER - Mean rt (s) for Eyes-Sad hits"]] +
      data[["Emoticom-ER - Mean rt (s) for Eyes-Anger hits"]] +
      data[["Emoticom-ER - Mean rt (s) for Eyes-Fear hits"]]
  ) / 4
  data$ERT_falsealarm_HSAF_time <- (
    data[["Emoticom-ER - Mean rt (s) for Eyes-Happy false alarms"]] +
      data[["Emoticom-ER - Mean rt (s) for Eyes-Sad false alarms"]] +
      data[["Emoticom-ER - Mean rt (s) for Eyes-Anger false alarms"]] +
      data[["Emoticom-ER - Mean rt (s) for Eyes-Fear false alarms"]]
  ) / 4
  data$ERT_time_avg_HSAF <- (data$ERT_hitrate_HSAF_time + data$ERT_falsealarm_HSAF_time) / 2

  data$IMT_asc_SH_num <- data[["Mean intensity threshold - Sad (ascending)"]] -
    data[["Mean intensity threshold - Happy (ascending)"]]
  data$IMT_desc_SH_num <- data[["Mean intensity threshold - Sad (descending)"]] -
    data[["Mean intensity threshold - Happy (descending)"]]
  data$IMT_asc_desc_avg_SH_num <- (data$IMT_asc_SH_num + data$IMT_desc_SH_num) / 2
  data$IMT_asc_SH_PP <- (
    ((data[["Mean intensity threshold - Sad (ascending)"]] - 1) / 14) * 100
  ) - (((data[["Mean intensity threshold - Happy (ascending)"]] - 1) / 14) * 100)
  data$IMT_desc_SH_PP <- (
    ((data[["Mean intensity threshold - Sad (descending)"]] - 1) / 14) * 100
  ) - (((data[["Mean intensity threshold - Happy (descending)"]] - 1) / 14) * 100)
  data$IMT_asc_desc_avg_SH_PP <- (data$IMT_asc_SH_PP + data$IMT_desc_SH_PP) / 2
  data$IMT_asc_time <- 0
  data$IMT_desc_time <- 0
  data$IMT_time <- 0

  data$IMT_H <- (data[["Mean intensity threshold - Happy (ascending)"]] +
    data[["Mean intensity threshold - Happy (descending)"]]) / 2
  data$IMT_S <- (data[["Mean intensity threshold - Sad (ascending)"]] +
    data[["Mean intensity threshold - Sad (descending)"]]) / 2
  data$IMT_A <- (data[["Mean intensity threshold - Anger (ascending)"]] +
    data[["Mean intensity threshold - Anger (descending)"]]) / 2
  data$IMT_F <- (data[["Mean intensity threshold - Fear (ascending)"]] +
    data[["Mean intensity threshold - Fear (descending)"]]) / 2
  data$IMT_D <- (data[["Mean intensity threshold - Disgust (ascending)"]] +
    data[["Mean intensity threshold - Disgust (descending)"]]) / 2
  data$IMT_H_bias <- abs(data$IMT_H - 7)
  data$IMT_S_bias <- abs(data$IMT_S - 7)
  data$IMT_A_bias <- abs(data$IMT_A - 7)
  data$IMT_D_bias <- abs(data$IMT_D - 7)
  data$IMT_F_bias <- abs(data$IMT_F - 7)
  data$IMT_HSADF_Bias_old <- (
    data$IMT_H_bias + data$IMT_S_bias + data$IMT_A_bias + data$IMT_D_bias + data$IMT_F_bias
  ) / 5
  data$IMT_HSADF_Bias <- ((data$IMT_S + data$IMT_A + data$IMT_D + data$IMT_F) / 4) - data$IMT_H

  met_names <- c(
    MET_GVI = "Emoticom-MJ - Mean rating (1-7), Guilty, Victim Intentional",
    MET_SVI = "Emoticom-MJ - Mean rating (1-7), Ashamed, Victim Intentional",
    MET_GAI = "Emoticom-MJ - Mean rating (1-7), Guilty, Agent Intentional",
    MET_SAI = "Emoticom-MJ - Mean rating (1-7), Ashamed, Agent Intentional",
    MET_GVU = "Emoticom-MJ - Mean rating (1-7), Guilty, Victim Unintentional",
    MET_SVU = "Emoticom-MJ - Mean rating (1-7), Ashamed, Victim Unintentional",
    MET_GAU = "Emoticom-MJ - Mean rating (1-7), Guilty, Agent Unintentional",
    MET_SAU = "Emoticom-MJ - Mean rating (1-7), Ashamed, Agent Unintentional",
    MET_GVI_time = "Emoticom-MJ - Mean rt (sec), Guilty, Victim Intentional",
    MET_SVI_time = "Emoticom-MJ - Mean rt (sec), Ashamed, Victim Intentional",
    MET_GAI_time = "Emoticom-MJ - Mean rt (sec), Guilty, Agent Intentional",
    MET_SAI_time = "Emoticom-MJ - Mean rt (sec), Ashamed, Agent Intentional",
    MET_GVU_time = "Emoticom-MJ - Mean rt (sec), Guilty, Victim Unintentional",
    MET_SVU_time = "Emoticom-MJ - Mean rt (sec), Ashamed, Victim Unintentional",
    MET_GAU_time = "Emoticom-MJ - Mean rt (sec), Guilty, Agent Unintentional",
    MET_SAU_time = "Emoticom-MJ - Mean rt (sec), Ashamed, Agent Unintentional"
  )
  for (target in names(met_names)) {
    data <- rename_column(data, met_names[[target]], target)
  }
  data$MET_mean_guilt <- (data$MET_GVI + data$MET_GAI + data$MET_GVU + data$MET_GAU) / 4
  data$MET_mean_ashamed <- (data$MET_SVI + data$MET_SAI + data$MET_SVU + data$MET_SAU) / 4
  data$MET_mean_guilt_time <- (
    data$MET_GVI_time + data$MET_GAI_time + data$MET_GVU_time + data$MET_GAU_time
  ) / 4
  data$MET_mean_ashamed_time <- (
    data$MET_SVI_time + data$MET_SAI_time + data$MET_SVU_time + data$MET_SAU_time
  ) / 4

  data$ERT_dprime_average_HS <- (data$SDT_Happy_dprime + data$SDT_Sad_dprime) / 2
  data$ERT_dprime_average_HSAF <- (
    data$SDT_Happy_dprime + data$SDT_Sad_dprime + data$SDT_Angry_dprime + data$SDT_Fear_dprime
  ) / 4
  data$ERT_dprime_bias <- data$SDT_Happy_dprime - rowMeans(
    data[, c("SDT_Sad_dprime", "SDT_Angry_dprime", "SDT_Fear_dprime"), drop = FALSE],
    na.rm = TRUE
  )

  z_score_columns <- c(
    SRT_Z_Score = "SRT",
    TMT_Corr_B_Z_Score = "TMT_Cor_B",
    SDMT_Z_Score = "SDMT",
    LNS_Z_Score = "LNS",
    Stroop_Z_Score = "Stroop_time",
    RAVLT_Z_Score = "RAVLT",
    ERT_hitrate_HS_Z_Score = "ERT_hitrate_HS",
    ERT_falsealarm_HS_Z_Score = "ERT_falsealarm_HS",
    SDT_Happy_dprime_Z_Score = "SDT_Happy_dprime",
    SDT_Sad_dprime_Z_Score = "SDT_Sad_dprime",
    IMT_asc_SH_PP_Z_Score = "IMT_asc_SH_PP",
    IMT_desc_SH_PP_Z_Score = "IMT_desc_SH_PP",
    IMT_asc_desc_avg_SH_PP_Z_Score = "IMT_asc_desc_avg_SH_PP",
    MET_mean_guilt_Z_score = "MET_mean_guilt",
    MET_mean_ashamed_Z_score = "MET_mean_ashamed"
  )
  for (target in names(z_score_columns)) {
    source <- z_score_columns[[target]]
    data[[target]] <- (data[[source]] - mean(data[[source]], na.rm = TRUE)) /
      stats::sd(data[[source]], na.rm = TRUE)
  }

  data
}

# Return the fixed ten-task mapping from canonical performance/reaction-time
# measures to structural availability used to construct 45 pairwise features.
cognitive_task_spec <- function() {
  data.frame(
    task = c(
      "SRT", "TMT", "SDMT", "LNS", "Stroop", "RAVLT",
      "ERT", "IMT", "MET_guilt", "MET_ashamed"
    ),
    perf_var = c(
      "SRT_scale", "TMT_scale", "SDMT", "LNS", "Stroop_scale", "RAVLT",
      "ERT_dprime_bias", "IMT_HSADF_Bias", "MET_mean_guilt", "MET_mean_ashamed"
    ),
    time_var = c(
      "SRT", "TMT_Cor_B", "SDMT_time", "LNS_time", "Stroop_time", "RAVLT_time",
      "ERT_time_avg_HSAF", "IMT_time", "MET_mean_guilt_time", "MET_mean_ashamed_time"
    ),
    perf_available = c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    time_available = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
}

# Standardize one task measure safely; structurally constant dimensions become
# zeros so availability can be handled explicitly rather than producing NaN.
z_safe <- function(x, eps = 1e-12) {
  standard_deviation <- stats::sd(x, na.rm = TRUE)
  mean_value <- mean(x, na.rm = TRUE)
  if (!is.finite(standard_deviation) || standard_deviation < eps) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mean_value) / standard_deviation)
}

# Generate stable names for every lower-triangle task pair.
make_pair_names <- function(n_tasks) {
  # Format one two-task index pair as a relational feature identifier.
  utils::combn(seq_len(n_tasks), 2L, FUN = function(pair) {
    sprintf("t%d_vs_t%d", pair[1], pair[2])
  })
}

# Describe every task-pair feature, its two source tasks, shared dimensions, and
# cosine/Pearson-product/no-shared-dimension construction rule.
task_pair_metadata <- function(task_spec, availability) {
  pair_indices <- utils::combn(seq_len(nrow(task_spec)), 2L)
  entries <- vector("list", ncol(pair_indices))
  for (pair_number in seq_len(ncol(pair_indices))) {
    pair <- pair_indices[, pair_number]
    shared_dimensions <- c("perf", "time")[availability[pair[1], ] & availability[pair[2], ]]
    method <- if (length(shared_dimensions) == 2L) {
      "cosine"
    } else if (length(shared_dimensions) == 1L) {
      "pearson_product"
    } else {
      "no_shared_dimension"
    }
    entries[[pair_number]] <- data.frame(
      feature = sprintf("t%d_vs_t%d", pair[1], pair[2]),
      task_a = task_spec$task[pair[1]],
      task_b = task_spec$task[pair[2]],
      shared_dimensions = paste(shared_dimensions, collapse = ","),
      method = method,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, entries)
}

# Compute one participant's ten-by-ten task similarity matrix: cosine with two
# shared dimensions, standardized product with one, and zero with none.
task_similarity_matrix <- function(task_matrix, availability, eps = 1e-12) {
  task_matrix <- as.matrix(task_matrix)
  storage.mode(task_matrix) <- "double"
  n_tasks <- nrow(task_matrix)
  similarity <- diag(1, n_tasks)

  for (first_task in seq_len(n_tasks - 1L)) {
    for (second_task in (first_task + 1L):n_tasks) {
      shared <- availability[first_task, ] & availability[second_task, ]
      if (!any(shared)) {
        value <- 0
      } else {
        first_vector <- task_matrix[first_task, shared, drop = TRUE]
        second_vector <- task_matrix[second_task, shared, drop = TRUE]
        if (any(!is.finite(first_vector)) || any(!is.finite(second_vector))) {
          value <- 0
        } else if (sum(shared) == 2L) {
          first_norm <- sqrt(sum(first_vector^2))
          second_norm <- sqrt(sum(second_vector^2))
          if (first_norm <= eps || second_norm <= eps) {
            value <- 0
          } else {
            value <- sum(first_vector * second_vector) / (first_norm * second_norm)
            value <- max(-1, min(1, value))
          }
        } else {
          value <- sum(first_vector * second_vector)
        }
      }
      similarity[first_task, second_task] <- value
      similarity[second_task, first_task] <- value
    }
  }
  similarity
}

# Convert one participant's paired task values into the 45 lower-triangle
# task-similarity features used by every clustering and prediction method.
subject_features <- function(row_data, task_pairs, availability, eps = 1e-12) {
  task_matrix <- do.call(
    rbind,
    # Extract the performance/reaction-time vector for one task.
    lapply(task_pairs, function(columns) as.numeric(row_data[1, columns]))
  )
  similarity <- task_similarity_matrix(task_matrix, availability, eps)
  similarity[lower.tri(similarity)]
}

# Fit the learned portion of task-feature preprocessing on development/training
# participants only. The object can then be applied unchanged to validation or
# generalization participants.
fit_cognitive_preprocessor <- function(data, id_col = "CIMBI",
                                       task_spec = cognitive_task_spec(),
                                       eps = 1e-12) {
  required <- unique(c(id_col, task_spec$perf_var, task_spec$time_var))
  assert_required_columns(data, required, "Preprocessor training data")
  ids <- as.character(data[[id_col]])
  if (anyNA(ids) || anyDuplicated(ids)) {
    stop("Preprocessor training IDs must be unique and non-missing.", call. = FALSE)
  }
  variables <- unique(c(task_spec$perf_var, task_spec$time_var))
  numeric_data <- data[, variables, drop = FALSE]
  numeric_data[] <- lapply(numeric_data, as.numeric)
  # Learn the training median for one canonical cognitive input.
  impute <- vapply(numeric_data, function(values) {
    value <- stats::median(values, na.rm = TRUE)
    if (is.finite(value)) value else 0
  }, numeric(1))
  completed <- numeric_data
  for (variable in variables) {
    missing <- !is.finite(completed[[variable]])
    completed[[variable]][missing] <- impute[[variable]]
  }
  center <- vapply(completed, mean, numeric(1))
  scale <- vapply(completed, stats::sd, numeric(1))
  scale[!is.finite(scale) | scale < eps] <- 1
  availability <- cbind(
    # Mark one task's performance dimension available from training variation.
    perf = vapply(task_spec$perf_var, function(variable) {
      stats::sd(completed[[variable]]) >= eps
    }, logical(1)),
    # Mark one task's reaction-time dimension available from training variation.
    time = vapply(task_spec$time_var, function(variable) {
      stats::sd(completed[[variable]]) >= eps
    }, logical(1))
  )
  rownames(availability) <- task_spec$task
  structure(list(
    id_col = id_col,
    task_spec = task_spec,
    variables = variables,
    impute = impute,
    center = center,
    scale = scale,
    availability = availability,
    eps = eps,
    training_ids = ids,
    imputation_policy = "training_median",
    fitted_n = nrow(data)
  ), class = c("cognitive_preprocessor", "list"))
}

# Apply fixed training medians, centers, scales, and availability decisions to
# new rows and return their aligned 45-feature/similarity/distance representation.
apply_cognitive_preprocessor <- function(preprocessor, data) {
  if (!inherits(preprocessor, "cognitive_preprocessor")) {
    stop("preprocessor must be a cognitive_preprocessor.", call. = FALSE)
  }
  required <- c(preprocessor$id_col, preprocessor$variables)
  assert_required_columns(data, required, "Preprocessor application data")
  ids <- as.character(data[[preprocessor$id_col]])
  if (anyNA(ids) || anyDuplicated(ids)) {
    stop("Application IDs must be unique and non-missing.", call. = FALSE)
  }
  completed <- data[, preprocessor$variables, drop = FALSE]
  completed[] <- lapply(completed, as.numeric)
  for (variable in preprocessor$variables) {
    missing <- !is.finite(completed[[variable]])
    completed[[variable]][missing] <- preprocessor$impute[[variable]]
  }
  standardized <- sweep(as.matrix(completed), 2L, preprocessor$center, "-")
  standardized <- sweep(standardized, 2L, preprocessor$scale, "/")
  standardized_data <- data.frame(CIMBI = ids, stringsAsFactors = FALSE)
  task_spec <- preprocessor$task_spec
  for (task_number in seq_len(nrow(task_spec))) {
    standardized_data[[paste0(task_spec$task[task_number], "_perf")]] <-
      standardized[, task_spec$perf_var[task_number]]
    standardized_data[[paste0(task_spec$task[task_number], "_time")]] <-
      standardized[, task_spec$time_var[task_number]]
  }
  # Resolve one task to its standardized performance/reaction-time columns.
  task_pairs <- lapply(task_spec$task, function(task) {
    c(paste0(task, "_perf"), paste0(task, "_time"))
  })
  n_features <- choose(nrow(task_spec), 2L)
  features <- t(vapply(
    seq_len(nrow(standardized_data)),
    # Construct the fixed 45-feature vector for one participant.
    function(index) subject_features(
      standardized_data[index, , drop = FALSE], task_pairs,
      preprocessor$availability, preprocessor$eps
    ),
    numeric(n_features)
  ))
  features <- as.matrix(features)
  rownames(features) <- ids
  colnames(features) <- make_pair_names(nrow(task_spec))
  similarity <- cosine_similarity_between_subjects(features, preprocessor$eps)
  dimnames(similarity) <- list(ids, ids)
  distance <- stats::as.dist(1 - similarity)
  attr(distance, "Labels") <- ids
  list(
    distance = distance,
    similarity = similarity,
    features = features,
    standardized_data = standardized_data,
    task_spec = task_spec,
    task_availability = preprocessor$availability,
    participant_ids = ids,
    preprocessor = preprocessor
  )
}

# Compute the finite cosine similarity matrix between participant feature rows,
# preserving self-similarity and explicit zero-vector behavior.
cosine_similarity_between_subjects <- function(features, eps = 1e-12) {
  features <- as.matrix(features)
  storage.mode(features) <- "double"
  norms <- sqrt(rowSums(features^2))
  valid_rows <- is.finite(norms) & norms > eps
  similarity <- matrix(0, nrow(features), nrow(features))
  if (any(valid_rows)) {
    valid_features <- features[valid_rows, , drop = FALSE]
    similarity[valid_rows, valid_rows] <-
      (valid_features %*% t(valid_features)) / (norms[valid_rows] %o% norms[valid_rows])
  }
  diag(similarity) <- 1
  similarity[similarity > 1] <- 1
  similarity[similarity < -1] <- -1
  similarity
}

# Build the complete cluster-ready representation from prepared cognitive data:
# standardized task values, 45 task-pair features, cosine similarity, and 1-S distance.
compute_subject_distance_matrix <- function(data, id_col = "CIMBI", task_spec = cognitive_task_spec(), eps = 1e-12) {
  required_columns <- unique(c(id_col, task_spec$perf_var, task_spec$time_var))
  assert_required_columns(data, required_columns, "Prepared data")
  identifiers <- as.character(data[[id_col]])
  if (anyNA(identifiers) || anyDuplicated(identifiers)) {
    stop("Prepared data must have unique, non-missing participant IDs.", call. = FALSE)
  }

  selected_data <- data[, required_columns, drop = FALSE]
  raw_summary <- data.frame(
    variable = required_columns,
    # Count missing values for one input column before standardization.
    n_missing = vapply(selected_data, function(x) sum(is.na(x)), integer(1)),
    # Count non-finite numeric values for one input column.
    n_nonfinite = vapply(selected_data, function(x) {
      if (is.numeric(x)) sum(!is.finite(x)) else 0L
    }, integer(1)),
    # Record one input column's observed standard deviation for auditability.
    sd_value = vapply(selected_data, function(x) {
      if (is.numeric(x)) stats::sd(x, na.rm = TRUE) else NA_real_
    }, numeric(1)),
    stringsAsFactors = FALSE
  )

  standardized_data <- data.frame(CIMBI = identifiers, stringsAsFactors = FALSE)
  for (task_number in seq_len(nrow(task_spec))) {
    standardized_data[[paste0(task_spec$task[task_number], "_perf")]] <-
      z_safe(selected_data[[task_spec$perf_var[task_number]]], eps)
    standardized_data[[paste0(task_spec$task[task_number], "_time")]] <-
      z_safe(selected_data[[task_spec$time_var[task_number]]], eps)
  }
  standardized_columns <- setdiff(names(standardized_data), "CIMBI")
  if (anyNA(standardized_data[, standardized_columns, drop = FALSE]) ||
      any(!is.finite(as.matrix(standardized_data[, standardized_columns, drop = FALSE])))) {
    stop("Standardized task variables contain missing or non-finite values.", call. = FALSE)
  }

  availability <- cbind(
    # Detect observed performance variation for one standardized task.
    perf = vapply(task_spec$task, function(task) {
      any(abs(standardized_data[[paste0(task, "_perf")]]) > eps)
    }, logical(1)),
    # Detect observed reaction-time variation for one standardized task.
    time = vapply(task_spec$task, function(task) {
      any(abs(standardized_data[[paste0(task, "_time")]]) > eps)
    }, logical(1))
  )
  rownames(availability) <- task_spec$task
  # Resolve one task to its standardized performance/reaction-time pair.
  task_pairs <- lapply(task_spec$task, function(task) c(paste0(task, "_perf"), paste0(task, "_time")))
  n_features <- choose(nrow(task_spec), 2L)
  feature_matrix <- t(vapply(
    seq_len(nrow(standardized_data)),
    # Construct one participant's 45 task-similarity features.
    function(index) subject_features(standardized_data[index, , drop = FALSE], task_pairs, availability, eps),
    numeric(n_features)
  ))
  feature_matrix <- as.matrix(feature_matrix)
  colnames(feature_matrix) <- make_pair_names(nrow(task_spec))
  rownames(feature_matrix) <- identifiers
  if (anyNA(feature_matrix) || any(!is.finite(feature_matrix))) {
    stop("Participant feature matrix contains missing or non-finite values.", call. = FALSE)
  }

  subject_similarity <- cosine_similarity_between_subjects(feature_matrix, eps)
  rownames(subject_similarity) <- identifiers
  colnames(subject_similarity) <- identifiers
  subject_distance <- stats::as.dist(1 - subject_similarity)
  attr(subject_distance, "Labels") <- identifiers

  list(
    distance = subject_distance,
    similarity = subject_similarity,
    features = feature_matrix,
    standardized_data = standardized_data,
    task_spec = task_spec,
    task_availability = availability,
    feature_metadata = task_pair_metadata(task_spec, availability),
    raw_summary = raw_summary
  )
}

# Write the auditable subject-distance TSV and JSON metadata describing inputs,
# task availability, feature construction, seeds, and package versions.
write_distance_outputs <- function(distance_result, preparation_audit, input_path, output_dir, seed) {
  require_cognitive_packages("jsonlite")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  matrix_path <- file.path(output_dir, "subject_distance_matrix.tsv")
  metadata_path <- file.path(output_dir, "subject_distance_metadata.json")
  distance_matrix <- as.matrix(distance_result$distance)
  utils::write.table(
    data.frame(CIMBI = rownames(distance_matrix), distance_matrix, check.names = FALSE),
    matrix_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  metadata <- list(
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    input_path = normalizePath(input_path, winslash = "/", mustWork = TRUE),
    seed = seed,
    r_version = R.version.string,
    package_versions = list(
      mgcv = as.character(utils::packageVersion("mgcv")),
      readxl = as.character(utils::packageVersion("readxl")),
      psycho = as.character(utils::packageVersion("psycho")),
      jsonlite = as.character(utils::packageVersion("jsonlite"))
    ),
    preparation_audit = preparation_audit,
    task_specification = distance_result$task_spec,
    task_availability = data.frame(task = rownames(distance_result$task_availability), distance_result$task_availability),
    feature_metadata = distance_result$feature_metadata,
    participant_ids = rownames(distance_result$features),
    output_files = list(distance_matrix_tsv = normalizePath(matrix_path, winslash = "/", mustWork = FALSE))
  )
  jsonlite::write_json(metadata, metadata_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  list(distance_matrix_tsv = matrix_path, metadata_json = metadata_path)
}

#' Load, screen, impute, derive, and cluster-ready distance-encode cognitive data.
#'
#' The default return value is the subject-level `dist` object. Use `return = "all"`
#' to obtain the prepared data, feature matrix, metadata, and holdout records.
#' The seeded 10% generalization holdout is removed before learned preparation,
#' as required by the SPEC for this cleaned cohort.
load_and_prepdata <- function(
  input_path,
  holdout_fraction = 0.10,
  holdout_cohort = "KOHORTE 2",
  seed = 2025L,
  max_missing_tests = 3L,
  output_dir = "results",
  write_outputs = TRUE,
  return = c("distance", "all")
) {
  require_cognitive_packages(c("readxl", "mgcv", "psycho"))
  return <- match.arg(return)
  if (!file.exists(input_path)) {
    stop("Input workbook does not exist: ", input_path, call. = FALSE)
  }
  if (!is.numeric(holdout_fraction) || length(holdout_fraction) != 1L ||
      is.na(holdout_fraction) || holdout_fraction < 0 || holdout_fraction > 1) {
    stop("holdout_fraction must be a single number from 0 to 1.", call. = FALSE)
  }

  data <- as.data.frame(readxl::read_excel(input_path, na = c("", "NA")), stringsAsFactors = FALSE)
  id_columns <- intersect(c("ID", "CIMBI ID", "CIMBI"), names(data))
  if (length(id_columns) != 1L) {
    stop("Input workbook must contain exactly one participant identifier column: ID, CIMBI ID or CIMBI.")
  }
  names(data)[names(data) == id_columns] <- "CIMBI"
  required_columns <- c(
    "CIMBI",
    "Age at neuropsych test",
    "Overall comments",
    cognitive_imputation_columns()
  )
  assert_required_columns(data, required_columns, "Input workbook")
  data$CIMBI <- as.character(data$CIMBI)
  if (anyNA(data$CIMBI) || any(!nzchar(data$CIMBI)) || anyDuplicated(data$CIMBI)) {
    stop("Input participant IDs must be unique and non-missing; no rows were removed.")
  }
  data <- canonicalize_records(data, require_distinct = TRUE)
  initial_n <- nrow(data)

  # Already-cleaned inputs may omit the source study column.
  # If supplied, retain the original study exclusion rule.
  project_filter_applied <- "Project ID (primary)" %in% names(data)
  if (project_filter_applied) {
  data <- data[data[["Project ID (primary)"]] != "NeuroPharm NP1 (Depression)", , drop = FALSE]
  }
  after_project_filter_n <- nrow(data)

  # Apply the analysis missingness screen before sampling the holdout. This
  # makes the requested fraction describe the participants who can actually
  # enter the analysis rather than a larger pre-screen population.
  groups <- cognitive_missingness_groups()
  screened_data <- filter_missing_tests(data, groups, max_missing_tests)
  screened_total_n <- nrow(screened_data)

  # By default the protocol targets KOHORTE 2. A NULL cohort makes every
  # screened participant eligible. The downstream calibration reference is
  # restricted to the same population when a cohort is configured.
  eligible <- if (is.null(holdout_cohort)) {
    seq_len(screened_total_n)
  } else {
    which(screened_data[["Overall comments"]] %in% holdout_cohort)
  }
  holdout_n <- floor(length(eligible) * holdout_fraction)
  holdout_seed <- analysis_stage_seed(seed, 100000L)
  set.seed(holdout_seed)
  holdout_rows <- if (holdout_n > 0L) {
    eligible[sample.int(length(eligible), holdout_n, replace = FALSE)]
  } else integer()
  holdout_data <- canonicalize_records(screened_data[holdout_rows, , drop = FALSE])
  development_rows <- setdiff(seq_len(screened_total_n), holdout_rows)
  data <- screened_data[development_rows, , drop = FALSE]
  development_raw_data <- data
  screened_development_raw_data <- data
  screened_holdout_data <- holdout_data
  after_missingness_filter_n <- nrow(data)
  imputation_columns <- cognitive_imputation_columns()
  imputed <- gam_impute_chained(data, imputation_columns, seed = seed)
  prepared_data <- derive_cognitive_measures(imputed$data)
  preparation_audit <- list(
    initial_n = initial_n,
    after_project_filter_n = after_project_filter_n,
    project_filter_applied = project_filter_applied,
    screened_total_n = screened_total_n,
    eligible_holdout_n = length(eligible),
    holdout_n = nrow(holdout_data),
    requested_holdout_fraction = holdout_fraction,
    realized_eligible_holdout_fraction = if (length(eligible)) {
      nrow(holdout_data) / length(eligible)
    } else {
      NA_real_
    },
    holdout_cohort = holdout_cohort,
    holdout_ids = holdout_data$CIMBI,
    holdout_seed = holdout_seed,
    row_order_policy = "non-ID non-age record contents; canonical_v1",
    after_missingness_filter_n = after_missingness_filter_n,
    imputation_missingness = vapply(imputed$was_missing, sum, integer(1)),
    imputation_audit = imputed$audit
  )
  distance_result <- compute_subject_distance_matrix(prepared_data)
  output_paths <- if (isTRUE(write_outputs)) {
    write_distance_outputs(distance_result, preparation_audit, input_path, output_dir, seed)
  } else {
    NULL
  }
  result <- c(
    list(
      distance = distance_result$distance,
      prepared_data = prepared_data,
      development_raw_data = development_raw_data,
      holdout_data = holdout_data,
      screened_development_raw_data = screened_development_raw_data,
      screened_holdout_data = screened_holdout_data,
      preparation_audit = preparation_audit,
      output_paths = output_paths
    ),
    distance_result[names(distance_result) != "distance"]
  )
  if (identical(return, "distance")) result$distance else result
}
