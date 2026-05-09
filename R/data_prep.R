#' Standardize Treatment to `{0, 1}`
#'
#' Coerces a 2-level treatment column to integer `{0, 1}` and returns
#' both the new vector and the original level mapping. The mapping is
#' stashed on the person-time output so display layers (`print`,
#' `summary`, `plot`) can relabel arms back to the user's original
#' codes.
#'
#' Mapping rule by input type:
#' - factor: `levels(.)[1]` -> 0, `levels(.)[2]` -> 1
#' - character: sorted unique values; lower -> 0, higher -> 1
#' - logical: `FALSE` -> 0, `TRUE` -> 1
#' - numeric: must already be `{0, 1}`; passed through unchanged
#'
#' @param x A 2-level treatment vector.
#' @return A list with two elements:
#'   - `values`: integer vector in `{0, 1}` of the same length as `x`.
#'   - `levels`: length-2 character vector. `levels[1]` is the original
#'     label for `A = 0`; `levels[2]` for `A = 1`.
#' @keywords internal
standardize_treatment <- function(x) {

  if (is.factor(x)) {
    if (length(levels(x)) != 2) {
      stop("Treatment factor must have exactly 2 levels. Got: ",
           length(levels(x)), call. = FALSE)
    }
    list(
      values = as.integer(x) - 1L,
      levels = levels(x)
    )

  } else if (is.character(x)) {
    lv <- sort(unique(x))
    if (length(lv) != 2) {
      stop("Treatment character vector must have exactly 2 unique ",
           "values. Got: ", length(lv), call. = FALSE)
    }
    list(
      values = as.integer(match(x, lv)) - 1L,
      levels = lv
    )

  } else if (is.logical(x)) {
    list(
      values = as.integer(x),
      levels = c("FALSE", "TRUE")
    )

  } else if (is.numeric(x)) {
    u <- unique(x)
    if (!setequal(u, c(0, 1))) {
      stop("Numeric treatment must be coded as {0, 1}. Found: ",
           paste(u, collapse = ", "), call. = FALSE)
    }
    list(
      values = as.integer(x),
      levels = c("0", "1")
    )

  } else {
    stop("Unsupported treatment type: ", class(x)[1],
         ". Use factor, character, numeric, or logical.",
         call. = FALSE)
  }
}


#' Prepare Subject-Level Data for Causal Survival Analysis
#'
#' Converts subject-level data (one row per subject) into discrete-time
#' person-time format using [survival::survSplit()]. Encodes a three-way
#' event/censoring split (`y_event`, `dep_cens`, `indep_cens`) per spec
#' §3.0, standardizes treatment to `{0, 1}`, and attaches metadata for
#' downstream estimation.
#'
#' @param data A subject-level data.frame (one row per subject).
#' @param id Character. Subject identifier column.
#' @param time Character. Event/censoring time column.
#' @param status Character. Binary event indicator column
#'   (`1` = event, `0` = censored).
#' @param ipcw Logical. Scalar (`TRUE` / `FALSE`) or length-`nrow(data)`
#'   logical vector classifying each subject's natural censoring
#'   mechanism. `TRUE` = dependent (subject's censoring contributes to
#'   the c-hazard fit and gets weighted via IPCW); `FALSE` = independent
#'   (weight 1, treated as cause-specific competitor). Default `TRUE`
#'   ("all natural censoring is dependent"). Per-subject vector form is
#'   the user's a-priori labeling, never inferred from data. Inert under
#'   `method = "gformula"` (`causal_survival()` warns).
#' @param T_max Numeric. Administrative truncation horizon. `NULL`
#'   (default) → `max(data[[time]])`.
#' @param treatment Character. 2-level treatment column (factor,
#'   character, logical, or numeric `{0, 1}`).
#' @param covariates Character vector. Baseline covariate columns.
#' @param cut_points Time grid specification:
#'   - `NULL` (default): 12 equi-spaced intervals over `(0, T_max]`.
#'   - Single positive integer: that many equi-spaced intervals.
#'   - Numeric vector of length >= 2: explicit interior cut points,
#'     strictly within `(0, T_max)`.
#' @param time_varying Reserved for future use. Must be `NULL` in v1.
#'
#' @return A data.frame of class `c("person_time", "data.frame")` with
#'   columns `id, k, treatment, <covariates>, y_event, dep_cens,
#'   indep_cens` and attributes `cut_times`, `T_max`, `K_max`,
#'   `treatment_levels`, `id_col`, `treatment_col`, `covariates`.
#'
#' @details
#' ## Time grid
#' Intervals are left-open right-closed:
#' `(0, t_1], (t_1, t_2], ..., (t_{K_max-1}, T_max]`. `k` is the
#' integer interval index (`1..K_max`); a subject's event time `t`
#' lands in interval `k` iff `t_{k-1} < t <= t_k`. Events at
#' `time = 0` are not supported (hard error). Subjects with
#' `time > T_max` are administratively censored: at-risk rows up to
#' `k = K_max`, no exit row materialized.
#'
#' ## Row encoding
#' Per row, at most one of `{y_event, dep_cens, indep_cens}` is `1`
#' (exit row); otherwise all three are `0` (at-risk row, including the
#' final row of admin-censored subjects). Encoding rules per spec
#' §3.0.4.
#'
#' ## Treatment standardization
#' The treatment column is coerced to integer `{0, 1}` via
#' [standardize_treatment()]. The original level mapping is stashed as
#' `attr(<output>, "treatment_levels")` for later display by Phase 3
#' accessors (`print`, `summary`, `plot`).
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   id = 1:100, time = rexp(100), status = rbinom(100, 1, 0.6),
#'   A = sample(c("ctrl", "trt"), 100, TRUE),
#'   age = rnorm(100, 60, 10)
#' )
#' pt <- to_person_time(
#'   df, id = "id", time = "time", status = "status",
#'   treatment = "A", covariates = "age", cut_points = 10
#' )
#' attr(pt, "treatment_levels")  # c("ctrl", "trt")
#' }
#'
#' @seealso [validate_subject_level()] for the input contract,
#'   [validate_person_time()] for the output schema (used when users
#'   bypass this builder).
#' @export
to_person_time <- function(data,
                           id           = "id",
                           time         = "time",
                           status       = "status",
                           ipcw         = TRUE,
                           T_max        = NULL,
                           treatment    = "A",
                           covariates   = character(),
                           cut_points   = NULL,
                           time_varying = NULL) {

  # --- Input shape + subject-level validation ---
  validate_input_shape(data, "data")
  validate_subject_level(
    data, id = id, time = time, status = status, treatment = treatment,
    covariates = covariates, cut_points = cut_points,
    time_varying = time_varying
  )
  # TODO: extend validate_subject_level for `ipcw` shape + T_max range.

  # Inline ipcw + T_max validation pending the validator extension.
  if (!is.logical(ipcw) || any(is.na(ipcw))) {
    stop("`ipcw` must be a non-NA logical scalar or vector.",
         call. = FALSE)
  }
  if (length(ipcw) != 1L && length(ipcw) != nrow(data)) {
    stop("`ipcw` must have length 1 or nrow(data).", call. = FALSE)
  }

  # --- Resolve T_max ---
  if (is.null(T_max)) {
    T_max <- max(data[[time]])
  } else if (T_max > max(data[[time]])) {
    stop("`T_max` must be <= max(data[[time]]).", call. = FALSE)
  }

  # --- Resolve cut_points -> cut_times (length K_max, last = T_max) ---
  if (is.null(cut_points)) {
    cut_times <- seq(0, T_max, length.out = 12L + 1L)[-1L]
  } else if (length(cut_points) == 1L) {
    cut_times <- seq(0, T_max, length.out = cut_points + 1L)[-1L]
  } else {
    cut_times <- sort(cut_points)
    if (cut_times[length(cut_times)] != T_max) {
      cut_times <- c(cut_times, T_max)
    }
  }
  K_max <- length(cut_times)

  # --- Resolve ipcw to length-n logical vector ---
  if (length(ipcw) == 1L) ipcw <- rep(ipcw, nrow(data))

  # --- Hard error: events at time = 0 (no home interval under (a, b]) ---
  events_at_zero <- which(data[[time]] == 0 & data[[status]] == 1)
  if (length(events_at_zero)) {
    stop("Event(s) at `time = 0` are not supported: intervals are ",
         "left-open ((0, t_1] excludes 0). Affected subject id(s): ",
         paste(data[[id]][events_at_zero], collapse = ", "), ".",
         call. = FALSE)
  }

  # --- Standardize treatment ---
  trt <- standardize_treatment(data[[treatment]])

  # --- Prepare for survSplit (intervals are (tstart, tstop]) ---
  # Subjects with `time > T_max` are admin-censored: tstop capped at
  # T_max and event_indicator = 0, so survSplit produces at-risk rows
  # up to k = K_max with no exit row. Events at exactly time = T_max
  # fire normally at k = K_max (T_max is inside the last interval
  # under the (a, b] convention).
  df <- data
  df[[treatment]]    <- trt$values
  df$.ipcw_label     <- ipcw
  df$tstart          <- 0
  df$tstop           <- pmin(df[[time]], T_max)
  df$event_indicator <- as.integer(df[[time]] <= T_max)

  # --- survSplit ---
  pt_data <- survival::survSplit(
    data  = df,
    cut   = cut_times,
    start = "tstart",
    end   = "tstop",
    event = "event_indicator"
  )

  # --- Integer interval index k = 1..K_max ---
  cut_breaks <- c(0, cut_times[-K_max])
  pt_data$k  <- match(pt_data$tstart, cut_breaks)

  # --- Derive y_event / dep_cens / indep_cens ---
  # Terminal rows carry event_indicator = 1 (subject exited before T_max).
  # Combine with status and per-subject ipcw label to split the three
  # mutually-exclusive indicators per spec §3.0.4.
  status_col <- pt_data[[status]]
  ipcw_col   <- pt_data$.ipcw_label
  pt_data$y_event    <- as.integer(
    pt_data$event_indicator == 1L & status_col == 1
  )
  pt_data$dep_cens   <- as.integer(
    pt_data$event_indicator == 1L & status_col == 0 & ipcw_col
  )
  pt_data$indep_cens <- as.integer(
    pt_data$event_indicator == 1L & status_col == 0 & !ipcw_col
  )

  # --- Drop scaffolding + original time/status ---
  pt_data$tstart          <- NULL
  pt_data$tstop           <- NULL
  pt_data$event_indicator <- NULL
  pt_data$.ipcw_label     <- NULL
  pt_data[[time]]         <- NULL
  pt_data[[status]]       <- NULL

  # --- Reorder columns ---
  ordered_cols <- c(id, "k", treatment, covariates,
                    "y_event", "dep_cens", "indep_cens")
  pt_data <- pt_data[, ordered_cols, drop = FALSE]

  # --- Attach metadata + class ---
  attr(pt_data, "cut_times")        <- cut_times
  attr(pt_data, "T_max")            <- T_max
  attr(pt_data, "K_max")            <- K_max
  attr(pt_data, "treatment_levels") <- trt$levels
  attr(pt_data, "id_col")           <- id
  attr(pt_data, "treatment_col")    <- treatment
  attr(pt_data, "covariates")       <- covariates
  class(pt_data) <- c("person_time", class(pt_data))

  pt_data
}
