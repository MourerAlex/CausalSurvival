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
#' person-time format using [survival::survSplit()]. Derives event-flag
#' columns (`y_flag`, `c_flag`) on terminal rows, standardizes treatment
#' to `{0, 1}`, and attaches metadata for downstream estimation.
#'
#' @param data A subject-level data.frame (one row per subject).
#' @param id Character. Subject identifier column.
#' @param time Character. Event/censoring time column.
#' @param status Character. Binary event indicator column
#'   (`1` = event, `0` = censored).
#' @param treatment Character. 2-level treatment column (factor,
#'   character, logical, or numeric `{0, 1}`).
#' @param covariates Character vector. Baseline covariate columns.
#' @param cut_points Time grid specification:
#'   - `NULL` (default): 12 equi-spaced intervals over `[0, max(time)]`.
#'   - Single positive integer: that many equi-spaced intervals.
#'   - Numeric vector of length >= 2: explicit interior cut points,
#'     strictly within `(0, max(time))`.
#' @param time_varying Reserved for future use. Must be `NULL` in v1.
#'
#' @return A data.frame of class `c("person_time", "data.frame")` with
#'   columns `id, k, treatment, <covariates>, y_flag, c_flag` and
#'   attributes `cut_times`, `treatment_levels`, `id_col`,
#'   `treatment_col`, `covariates`.
#'
#' @details
#' ## Time grid
#' Event times are shifted by `+1` (`tstop = time + 1`) so that `k = 0`
#' represents the first at-risk interval. Subjects with `time = 0` land
#' in `k = 0`.
#'
#' For example, with `cut_points = 4` and `max(time) = 3` (yearly events
#' encoded as `time` in `{0, 1, 2, 3}`), the +1 shift produces 4
#' unit-width intervals labeled by `k`:
#' - `k = 0` -> `[0, 1)`   "year 1"   (subjects with `time = 0` land here)
#' - `k = 1` -> `[1, 2)`   "year 2"
#' - `k = 2` -> `[2, 3)`   "year 3"
#' - `k = 3` -> `[3, 4)`   "year 4"
#'
#' ## survSplit convention
#' `event_indicator = 1L` is set for ALL subjects (including censored)
#' before survSplit. This flags the terminal row of every subject so
#' that `y_flag` and `c_flag` can be derived by combining the indicator
#' with the original `status` column.
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
                           id = "id",
                           time = "time",
                           status = "status",
                           treatment = "A",
                           covariates = character(),
                           cut_points = NULL,
                           time_varying = NULL) {

  # --- Input shape + subject-level validation ---
  validate_input_shape(data, "data")
  validate_subject_level(
    data, id = id, time = time, status = status, treatment = treatment,
    covariates = covariates, cut_points = cut_points,
    time_varying = time_varying
  )

  # --- Resolve cut_points -> cut_times ---
  max_time <- max(data[[time]])
  if (is.null(cut_points)) {
    cut_times <- seq(0, max_time, length.out = 12L + 1L)[-1L]
  } else if (length(cut_points) == 1L) {
    cut_times <- seq(0, max_time, length.out = cut_points + 1L)[-1L]
  } else {
    cut_times <- sort(cut_points)
  }

  # --- Standardize treatment ---
  trt <- standardize_treatment(data[[treatment]])

  # --- Prepare for survSplit ---
  df <- data
  df[[treatment]] <- trt$values
  df$event_indicator <- 1L
  df$tstart <- 0
  df$tstop <- df[[time]] + 1

  # --- survSplit ---
  pt_data <- survival::survSplit(
    data  = df,
    cut   = cut_times,
    start = "tstart",
    end   = "tstop",
    event = "event_indicator"
  )

  # --- k = 0-based interval index ---
  pt_data$k <- pt_data$tstart

  # --- Derive y_flag, c_flag from terminal rows ---
  # event_indicator == 1 marks the row where the subject leaves the
  # study; combine with original status to split into y_flag vs c_flag.
  pt_data$y_flag <- as.integer(
    pt_data$event_indicator == 1L & pt_data[[status]] == 1
  )
  pt_data$c_flag <- as.integer(
    pt_data$event_indicator == 1L & pt_data[[status]] == 0
  )

  # --- Drop survSplit scaffolding + original time/status ---
  # Original `time` and `status` are dropped: they are redundant with
  # `k + y_flag + c_flag` once discretized. See dev/TODO.md
  # "Optional pass-through of original time / status..." — when we
  # build `simulate_causal_survival_data()`, add a flag to keep them
  # broadcast on each row for DGP cross-checks against continuous truth.
  pt_data$tstart          <- NULL
  pt_data$tstop           <- NULL
  pt_data$event_indicator <- NULL
  pt_data[[time]]         <- NULL
  pt_data[[status]]       <- NULL

  # --- Reorder columns: id, k, treatment, covariates, y_flag, c_flag ---
  ordered_cols <- c(id, "k", treatment, covariates, "y_flag", "c_flag")
  pt_data <- pt_data[, ordered_cols, drop = FALSE]

  # --- Attach metadata + class ---
  attr(pt_data, "cut_times")        <- cut_times
  attr(pt_data, "treatment_levels") <- trt$levels
  attr(pt_data, "id_col")           <- id
  attr(pt_data, "treatment_col")    <- treatment
  attr(pt_data, "covariates")       <- covariates
  class(pt_data) <- c("person_time", class(pt_data))

  pt_data
}
