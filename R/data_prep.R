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
#' person-time format using [survival::survSplit()]. Derives the
#' three-way censoring split (`y_event`, `dep_cens`, `indep_cens`) on
#' terminal rows, standardizes treatment to `{0, 1}`, materializes the
#' admin-truncation convention, and attaches metadata for downstream
#' estimation.
#'
#' Schema and structural ordering are LOCKED — see
#' `dev/CAUSAL_SURVIVAL_SPEC.md` §3.0.
#'
#' @param data A subject-level data.frame (one row per subject).
#' @param id Character. Subject identifier column.
#' @param time Character. Event/censoring time column. Must be strictly
#'   positive — `time = 0` events are a hard error (no home interval
#'   under the `(0, t_1]` convention).
#' @param status Character. Binary event indicator column
#'   (`1` = event, `0` = censored).
#' @param ipcw Logical, scalar OR length-`nrow(data)` vector.
#'   `TRUE` = the subject's censoring is dependent (LTFU /
#'   treatment-switch / etc.) — contributes to the c-hazard fit and
#'   gets IPCW-weighted. `FALSE` = independent (admin-style) —
#'   excluded from the c-hazard fit and assigned weight 1. This is the
#'   user's a-priori labeling of each subject's censoring mechanism;
#'   never inferred from data. Honored only when `status = 0` AND
#'   `time <= T_max`; otherwise ignored. Default `TRUE`.
#' @param T_max Numeric scalar, end of the analyzable time grid. The
#'   reporting grid is `(0, T_max]` partitioned into `K_max` intervals
#'   (see `cut_points`). Subjects with `time > T_max` are
#'   administratively truncated (no exit row, contribute at-risk rows
#'   up to `k = K_max`). Defaults to `max(data[[time]])`. Hard error
#'   if `T_max > max(data[[time]])` (would induce empty trailing
#'   intervals).
#' @param treatment Character. 2-level treatment column (factor,
#'   character, logical, or numeric `{0, 1}`).
#' @param covariates Character vector. Baseline covariate columns.
#' @param cut_points Time grid specification over `(0, T_max]`:
#'   - `NULL` (default): 12 equi-spaced intervals.
#'   - Single positive integer: that many equi-spaced intervals.
#'   - Numeric vector of length >= 2: explicit interior cut points,
#'     strictly within `(0, T_max)`.
#' @param time_varying Reserved for future use. Must be `NULL` in v1.
#'
#' @return A data.frame of class `c("person_time", "data.frame")` with
#'   columns `id, k, A, <covariates>, y_event, dep_cens, indep_cens`
#'   and attributes `cut_times`, `T_max`, `K_max`, `treatment_levels`,
#'   `id_col`, `treatment_col`, `covariates`. Per row, at most one of
#'   `{y_event, dep_cens, indep_cens}` is `1` (exit row); otherwise all
#'   three are `0` (at-risk row, including the final row of admin-
#'   truncated subjects).
#'
#' @details
#' ## Time grid (LOCKED §3.0.2)
#' Continuous time over `(0, T_max]` is partitioned into `K_max`
#' analyzable intervals. Intervals are left-open right-closed: a
#' subject's event time `t` lands in interval `k` iff
#' `t_{k-1} < t <= t_k`. `k = 0` is pre-baseline and not in the data.
#' Estimates are reported at `k = 1, ..., K_max`.
#'
#' Boundary handling at `t = T_max`: `time = T_max` is inside the last
#' interval `(t_{K_max-1}, T_max]` (events fire normally at
#' `k = K_max`). `time > T_max` is outside all intervals — those
#' subjects are admin-truncated.
#'
#' ## Structural ordering (LOCKED §3.0.1)
#' Within each interval: `C_admin -> C_dep -> Y`. Admin censoring is
#' placed first within the interval — subjects with `time > T_max`
#' have no exit row regardless of `status`.
#'
#' ## Row encoding (LOCKED §3.0.4)
#'
#' | input                                              | indicator        |
#' |----------------------------------------------------|------------------|
#' | `status = 1, time <= T_max`                        | `y_event = 1`    |
#' | `status = 1, time > T_max`                         | (no exit row)    |
#' | `status = 0, time <= T_max, ipcw[i] = TRUE`        | `dep_cens = 1`   |
#' | `status = 0, time <= T_max, ipcw[i] = FALSE`       | `indep_cens = 1` |
#' | `status = 0, time > T_max`                         | (no exit row)    |
#' | at-risk (incl. admin-truncated subjects)           | all `0`          |
#'
#' ## Diagnostics
#' - **Hard error** if `T_max > max(data[[time]])` — empty trailing
#'   intervals.
#' - **Hard error** if any `status = 1` row has `time = 0` — no home
#'   interval. Lists affected subject ids.
#' - **Warning** if `mean(admin_reach) < 0.5` after encoding — "hazard
#'   at K_max from thin risk set; CIF at K_max unreliable". Fit
#'   proceeds.
#'
#' ## Treatment standardization
#' The treatment column is coerced to integer `{0, 1}` via
#' [standardize_treatment()]. The original level mapping is stashed as
#' `attr(<output>, "treatment_levels")` for later display by accessors
#' (`print`, `summary`, `plot`).
#'
#' ## Pre-split mode
#' Dropped from v1 (see spec §3.0.9). Users with already-classified
#' censoring must run their classification through `status` + `ipcw`
#' (per-subject vector). The function always discretizes from
#' subject-level input.
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
#'
#' # Per-subject censoring classification:
#' ipcw_vec <- df$reason %in% c("ltfu", "switch")  # TRUE = dep_cens
#' pt2 <- to_person_time(df, id = "id", time = "time", status = "status",
#'                       treatment = "A", ipcw = ipcw_vec)
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
                           ipcw = TRUE,
                           T_max = NULL,
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
  # NOTE: ipcw / T_max validation lives here for now; will move into
  # validate_subject_level() in Phase 3 step 4 of the refactor.

  # --- Resolve T_max ---
  max_time <- max(data[[time]])
  if (is.null(T_max)) {
    T_max <- max_time
  } else if (T_max > max_time) {
    stop("T_max (", T_max, ") exceeds max(data[[time]]) (", max_time,
         "). Would create empty trailing intervals.", call. = FALSE)
  }

  # --- Hard error: any time = 0 (no home interval under (0, t_1]) ---
  # Subjects with time = 0 have no follow-up window, so they cannot
  # be placed in any analyzable interval. Errored regardless of
  # status to avoid a silent drop by survSplit (tstop = 0 produces
  # no rows). NA-safe: validate_subject_level already rejects NA in
  # time/status, but check defensively in case validator is bypassed.
  zero_time <- which(data[[time]] == 0)
  if (length(zero_time) > 0) {
    ids <- data[[id]][zero_time]
    stop("time = 0 not supported (no home interval under (0, t_1] ",
         "convention; would be silently dropped by survSplit). ",
         "Affected ", id, ": ",
         paste(ids, collapse = ", "), ". ",
         "Drop or recode these subjects upstream.", call. = FALSE)
  }

  # --- Resolve ipcw shape ---
  if (length(ipcw) == 1L) {
    ipcw_vec <- rep(as.logical(ipcw), nrow(data))
  } else if (length(ipcw) == nrow(data)) {
    ipcw_vec <- as.logical(ipcw)
  } else {
    stop("ipcw must be a scalar logical or length-", nrow(data),
         " logical vector. Got length ", length(ipcw), ".",
         call. = FALSE)
  }
  if (any(is.na(ipcw_vec))) {
    stop("ipcw contains NA values. Must be TRUE/FALSE per subject.",
         call. = FALSE)
  }

  # --- Resolve cut_points -> cut_times over (0, T_max] ---
  # cut_times holds the right endpoints t_1 < ... < t_{K_max} = T_max.
  if (is.null(cut_points)) {
    cut_times <- seq(0, T_max, length.out = 12L + 1L)[-1L]
  } else if (length(cut_points) == 1L) {
    cut_times <- seq(0, T_max, length.out = cut_points + 1L)[-1L]
  } else {
    cut_times <- sort(cut_points)
    if (cut_times[1] <= 0 || cut_times[length(cut_times)] >= T_max) {
      stop("Explicit cut_points must lie strictly within (0, T_max). ",
           "T_max = ", T_max, ".", call. = FALSE)
    }
    cut_times <- c(cut_times, T_max)
  }
  K_max <- length(cut_times)

  # --- Standardize treatment ---
  trt <- standardize_treatment(data[[treatment]])

  # --- Identify admin-truncated subjects (time > T_max).  These get
  #     no exit row; their final at-risk row is at k = K_max. ---
  admin_trunc <- data[[time]] > T_max

  # --- Prepare for survSplit ---
  # Cap times at T_max so survSplit produces rows up to k = K_max for
  # admin-truncated subjects. event_indicator = 1 will be suppressed
  # on those rows later.
  df <- data
  df[[treatment]] <- trt$values
  df[[time]]      <- pmin(df[[time]], T_max)
  df$event_indicator <- 1L
  df$tstart <- 0
  df$tstop  <- df[[time]]

  # Broadcast subject-level metadata for per-row flag assembly.
  df$.ipcw_subject        <- ipcw_vec
  df$.admin_trunc_subject <- admin_trunc
  df$.status_subject      <- df[[status]]

  # --- survSplit ---
  pt_data <- survival::survSplit(
    data  = df,
    cut   = cut_times,
    start = "tstart",
    end   = "tstop",
    event = "event_indicator"
  )

  # --- Integer k under (0, t_1], ..., (t_{K_max-1}, T_max] ---
  pt_data$k <- findInterval(pt_data$tstop, c(0, cut_times),
                            left.open = TRUE)

  # --- Three-way exit-flag assembly ---
  # event_indicator == 1 marks each subject's terminal row. For admin-
  # truncated subjects, the terminal row is forced back to at-risk
  # (all three flags 0). Otherwise the row splits by status + ipcw.
  is_terminal <- pt_data$event_indicator == 1L &
                 !pt_data$.admin_trunc_subject
  is_event    <- is_terminal & pt_data$.status_subject == 1
  is_dep      <- is_terminal & pt_data$.status_subject == 0 &
                   pt_data$.ipcw_subject
  is_indep    <- is_terminal & pt_data$.status_subject == 0 &
                   !pt_data$.ipcw_subject

  pt_data$y_event    <- as.integer(is_event)
  pt_data$dep_cens   <- as.integer(is_dep)
  pt_data$indep_cens <- as.integer(is_indep)

  # --- Drop survSplit scaffolding + original time/status ---
  pt_data$tstart          <- NULL
  pt_data$tstop           <- NULL
  pt_data$event_indicator <- NULL
  pt_data[[time]]         <- NULL
  pt_data[[status]]       <- NULL
  pt_data$.ipcw_subject        <- NULL
  pt_data$.admin_trunc_subject <- NULL
  pt_data$.status_subject      <- NULL

  # --- Reorder columns ---
  ordered_cols <- c(id, "k", treatment, covariates,
                    "y_event", "dep_cens", "indep_cens")
  pt_data <- pt_data[, ordered_cols, drop = FALSE]

  # --- admin_reach diagnostic ---
  n_subjects   <- length(unique(data[[id]]))
  n_reach_kmax <- length(unique(pt_data[[id]][pt_data$k == K_max]))
  admin_reach_frac <- n_reach_kmax / n_subjects
  if (admin_reach_frac < 0.5) {
    warning("mean(admin_reach) = ",
            formatC(admin_reach_frac, digits = 3, format = "f"),
            " (< 0.5): hazard at K_max from thin risk set; ",
            "CIF at K_max unreliable.", call. = FALSE)
  }

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
