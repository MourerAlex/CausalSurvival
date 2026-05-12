# ARCHIVED — `validate_person_time()`
#
# Ported from CausalCompetingRisks (CCR) on 2026-05-06 (Phase 1).
# Removed from the live codebase on 2026-05-12 because spec §3.0.9
# (LOCKED 2026-05-08) drops the BYO ("bring your own pt_data") path:
#
#   "Users **cannot** supply already-discretized person-time data;
#   the function always discretizes from subject-level input."
#
# With BYO gone, this validator never runs in the canonical flow
# (`to_person_time()` → `causal_survival()`). Kept here in case
# §3.0.9 is reversed in v1.x and BYO is reintroduced.
#
# Last refresh before archival:
#   - Schema updated to y_event / dep_cens / indep_cens.
#   - k integer >= 1 enforced (spec §3.0.2).
#   - Left-truncation check shifted to k = 1.
#   - NA in flag columns hard-errored (spec §3.0.3).
#
# Live validators retained: validate_input_shape(),
# validate_subject_level(), check_covariate_quality().

#' Validate Person-Time Input
#'
#' Checks that user-supplied person-time data (one row per subject-interval)
#' has the structure required by the package's hazard pipeline. Only
#' called when the input does NOT already come from [to_person_time()]
#' (i.e., no `"person_time"` class).
#'
#' @param pt_data A data.frame in person-time format.
#' @param id,treatment Character column names.
#' @param covariates Character vector.
#'
#' @return Invisibly returns TRUE if all checks pass.
#'
#' @details
#' ## Hard errors
#' - NULL column names
#' - Required columns missing (id, treatment, k, y_event, dep_cens,
#'   indep_cens)
#' - NAs in id or treatment
#' - NAs in any of `{y_event, dep_cens, indep_cens}` (spec §3.0.3)
#' - Flag columns contain values other than 0 or 1
#' - Mutual-exclusivity invariant violated (more than one of
#'   `{y_event, dep_cens, indep_cens}` equals 1 in the same row)
#' - Duplicate (id, k) pairs
#' - Treatment not coded as {0, 1}
#' - `k` not integer, negative, or zero (must be `>= 1`)
#' - Left-truncated subjects (no row at k = 1)
#' - Covariate NAs or unsupported types (via [check_covariate_quality()])
#'
#' ## Warnings
#' - Constant or high-cardinality covariates (via [check_covariate_quality()])
#'
#' @keywords internal
validate_person_time <- function(pt_data,
                                 id,
                                 treatment,
                                 covariates) {

  # --- Column name arguments must not be NULL ---
  if (is.null(id) || is.null(treatment)) {
    stop(
      "id and treatment must be column names, not NULL.",
      call. = FALSE
    )
  }

  # --- Required person-time columns ---
  flag_cols <- c("y_event", "dep_cens", "indep_cens")
  required  <- c(id, treatment, "k", flag_cols, covariates)
  missing_cols <- setdiff(required, names(pt_data))
  if (length(missing_cols) > 0) {
    stop(
      "Person-time data is missing required column(s): ",
      paste(missing_cols, collapse = ", "), ". ",
      "Use to_person_time() to prepare your data.",
      call. = FALSE
    )
  }

  # --- Critical columns must have no NAs ---
  for (col_name in c(id, treatment)) {
    n_na <- sum(is.na(pt_data[[col_name]]))
    if (n_na > 0) {
      stop(
        "Column '", col_name, "' contains ", n_na, " NA value(s). ",
        "Critical columns must have no missing values.",
        call. = FALSE
      )
    }
  }

  # --- Flag columns must contain only {0, 1} (spec §3.0.3 — no NA) ---
  for (flag_col in flag_cols) {
    vals <- pt_data[[flag_col]]
    if (any(is.na(vals))) {
      n_na <- sum(is.na(vals))
      stop(
        "Column '", flag_col, "' contains ", n_na, " NA value(s). ",
        "Flag columns must be 0 or 1 (spec \u00a73.0.3).",
        call. = FALSE
      )
    }
    unique_vals <- unique(vals)
    if (!all(unique_vals %in% c(0, 1))) {
      stop(
        "Column '", flag_col, "' must contain only 0 or 1. ",
        "Found: ", paste(unique_vals, collapse = ", "),
        call. = FALSE
      )
    }
  }

  # --- Mutual-exclusivity: at most one of {y_event, dep_cens,
  #     indep_cens} can be 1 in any given row (spec §3.0.3). NAs are
  #     already rejected above; rowSums runs fail-loud. ---
  flag_sum <- rowSums(
    pt_data[, flag_cols, drop = FALSE] == 1L
  )
  if (any(flag_sum > 1)) {
    n_bad <- sum(flag_sum > 1)
    stop(
      "Mutual-exclusivity invariant violated in ", n_bad, " row(s): ",
      "more than one of {y_event, dep_cens, indep_cens} equals 1. ",
      "Each row must encode at most one terminal event.",
      call. = FALSE
    )
  }

  # --- Treatment must be {0, 1} ---
  trt_vals <- unique(pt_data[[treatment]])
  if (!setequal(trt_vals, c(0, 1))) {
    stop(
      "Treatment column '", treatment, "' must be coded as {0, 1}. ",
      "Found: ", paste(trt_vals, collapse = ", "), ". ",
      "Recode before calling, or use to_person_time() which standardizes ",
      "to {0, 1}.",
      call. = FALSE
    )
  }

  # --- k must be integer-valued and >= 1 (spec §3.0.2) ---
  k_vals <- pt_data$k
  if (any(is.na(k_vals))) {
    stop("Column 'k' contains NA values.", call. = FALSE)
  }
  if (!is.numeric(k_vals) || any(k_vals != as.integer(k_vals))) {
    stop("Column 'k' must contain integer values.", call. = FALSE)
  }
  if (any(k_vals < 1)) {
    stop(
      "Column 'k' must be >= 1 (interval index starts at 1 under the ",
      "(0, t_1], ..., (t_{K_max-1}, T_max] convention). ",
      "Found min k = ", min(k_vals), ".",
      call. = FALSE
    )
  }

  # --- Duplicate (id, k) check ---
  dupes <- duplicated(pt_data[, c(id, "k")])
  if (any(dupes)) {
    stop(
      "Duplicate (", id, ", k) pairs detected in person-time data. ",
      "Each subject must have at most one row per interval.",
      call. = FALSE
    )
  }

  # --- Left-truncation check: every subject must have a k = 1 row ---
  # Left-truncated data (subjects first observed at k > 1) is not
  # supported in v1 and will not be supported in future versions. This
  # is a structural assumption of the discrete-time pooled logistic
  # framework.
  min_k_per_subject <- tapply(pt_data$k, pt_data[[id]], min)
  if (any(min_k_per_subject > 1)) {
    n_affected <- sum(min_k_per_subject > 1)
    stop(
      "Left-truncated data is not supported. ",
      n_affected, " subject(s) have no row at k = 1. ",
      "Every subject must be observed from the first interval onward.",
      call. = FALSE
    )
  }

  # --- Covariate checks (errors + warnings) ---
  check_covariate_quality(pt_data, covariates)

  invisible(TRUE)
}
