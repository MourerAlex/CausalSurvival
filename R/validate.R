#' Validate Basic Input Data Frame Shape
#'
#' Checks that an input object is a non-empty data.frame.
#'
#' @param x The object to check.
#' @param name Character. Name of the argument (for error messages).
#'
#' @return Invisibly TRUE if all checks pass; otherwise throws a helpful
#'   error.
#' @keywords internal
validate_input_shape <- function(x, name) {
  if (is.null(x)) {
    stop(name, " is NULL.", call. = FALSE)
  }
  if (!is.data.frame(x)) {
    stop(
      name, " must be a data.frame. Got class: ",
      paste(class(x), collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(x) == 0) {
    stop(name, " has 0 rows.", call. = FALSE)
  }
  invisible(TRUE)
}


#' Shared Covariate Checks
#'
#' Hard errors on (a) NAs, (b) unsupported types (anything other than
#' numeric / integer / factor / character / logical — these break glm
#' downstream). Warnings on (c) constant columns, (d) high-cardinality
#' factors.
#'
#' @param data A data.frame.
#' @param covariates Character vector of covariate column names.
#' @return Invisibly NULL.
#' @keywords internal
check_covariate_quality <- function(data, covariates) {

  errors        <- character()
  warnings_msgs <- character()

  for (cov in covariates) {
    col <- data[[cov]]
    col_class <- class(col)[1]

    # --- Type (ERROR) ---
    if (!col_class %in% c("numeric", "integer", "factor",
                          "character", "logical")) {
      errors <- c(errors,
        paste0("  - '", cov, "': unsupported type '", col_class,
               "' (expected numeric, integer, factor, character, or logical)")
      )
      next
    }

    # --- NAs (ERROR) ---
    if (any(is.na(col))) {
      errors <- c(errors,
        paste0("  - '", cov, "': contains ", sum(is.na(col)), " NA value(s)")
      )
    }

    # --- Constant (WARNING) ---
    n_unique_nona <- length(unique(col[!is.na(col)]))
    if (n_unique_nona < 2) {
      warnings_msgs <- c(warnings_msgs,
        paste0("  - '", cov, "': fewer than 2 unique values")
      )
    }

    # --- Cardinality, categorical only (WARNING) ---
    if (col_class %in% c("factor", "character")) {
      n_levels <- if (col_class == "factor") length(levels(col))
                  else length(unique(col))
      if (n_levels > 20) {
        warnings_msgs <- c(warnings_msgs,
          paste0("  - '", cov, "': ", n_levels, " levels ",
                 "(high cardinality - consider grouping)")
        )
      }
    }
  }

  if (length(errors) > 0) {
    stop(
      "Covariate validation failed:\n",
      paste(errors, collapse = "\n"),
      call. = FALSE
    )
  }
  if (length(warnings_msgs) > 0) {
    warning(
      "Covariate quality issues:\n",
      paste(warnings_msgs, collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(NULL)
}


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
#' - Required columns missing (id, treatment, k, y_flag, c_flag)
#' - NAs in id or treatment
#' - Flag columns contain values other than 0, 1, or NA
#' - Duplicate (id, k) pairs
#' - Treatment not coded as {0, 1}
#' - Left-truncated subjects (no row at k = 0)
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
  required <- c(id, treatment, "k", "y_flag", "c_flag", covariates)
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

  # --- Flag columns must contain only {0, 1, NA} ---
  for (flag_col in c("y_flag", "c_flag")) {
    vals <- unique(pt_data[[flag_col]])
    vals_nona <- vals[!is.na(vals)]
    if (!all(vals_nona %in% c(0L, 1L, 0, 1))) {
      stop(
        "Column '", flag_col, "' must contain only 0, 1, or NA. ",
        "Found: ", paste(vals, collapse = ", "),
        call. = FALSE
      )
    }
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

  # --- Duplicate (id, k) check ---
  dupes <- duplicated(pt_data[, c(id, "k")])
  if (any(dupes)) {
    stop(
      "Duplicate (", id, ", k) pairs detected in person-time data. ",
      "Each subject must have at most one row per interval.",
      call. = FALSE
    )
  }

  # --- Left-truncation check: every subject must have a k = 0 row ---
  # Left-truncated data (subjects first observed at k > 0) is not supported
  # in v1 and will not be supported in future versions. This is a structural
  # assumption of the discrete-time pooled logistic framework.
  min_k_per_subject <- tapply(pt_data$k, pt_data[[id]], min)
  if (any(min_k_per_subject > 0)) {
    n_affected <- sum(min_k_per_subject > 0)
    stop(
      "Left-truncated data is not supported. ",
      n_affected, " subject(s) have no row at k = 0. ",
      "Every subject must be observed from the first interval onward.",
      call. = FALSE
    )
  }

  # --- Covariate checks (errors + warnings) ---
  check_covariate_quality(pt_data, covariates)

  invisible(TRUE)
}
