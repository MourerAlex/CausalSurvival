#' Null-Coalescing Operator
#'
#' Returns `x` if not NULL, otherwise `y`.
#'
#' @param x,y Values to coalesce.
#' @return `x` if not NULL, else `y`.
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


#' Clone Baseline Across Interval Indices
#'
#' Broadcast each subject's baseline row across all `K_max` interval
#' indices, setting the treatment column to a fixed value `a`. Used by
#' the g-formula and IPW-MSM workers to predict counterfactual hazards
#' at every (subject, k) regardless of the subject's observed
#' event/censoring time.
#'
#' The `k` column on the clone holds the integer interval index
#' (`1..K_max`), matching the fit-time encoding of the Y-hazard / Y-MSM
#' model. The time grid `cut_times` is preserved on the calling side
#' (e.g. via `attr(pt_data, "cut_times")`) for report-time alignment.
#'
#' @param baseline data.frame. One row per subject (typically
#'   `pt_data[pt_data$k == 1, ]`).
#' @param cut_times Numeric vector of interval-end times `t_1, ..., T_max`.
#'   Used here only for its length `K_max`.
#' @param treatment_col Character. Treatment column name.
#' @param a Numeric (0 or 1). Counterfactual treatment value.
#' @return data.frame with `nrow(baseline) * K_max` rows.
#' @keywords internal
make_clone <- function(baseline, cut_times, treatment_col, a) {
  n <- nrow(baseline)
  K <- length(cut_times)
  clone <- baseline[rep(seq_len(n), each = K), , drop = FALSE]
  clone$k                <- rep(seq_len(K), times = n)
  clone[[treatment_col]] <- a
  rownames(clone) <- NULL
  clone
}
