#' Risk-table counts by observed treatment arm
#'
#' Returns counts at each cut time, grouped by the observed treatment
#' level (e.g., A = 0 and A = 1). Complements the cumulative
#' incidence curves for reporting (adjustedCurves-style risk table).
#'
#' The user picks **one** count type per call. Counts are computed
#' on observed subjects — there are no counterfactual rows under
#' v0.1.0's binary point-treatment scope.
#'
#' @param fit A `"causal_survival_fit"` object from
#'   [causal_survival()].
#' @param count Character. One of:
#'   \describe{
#'     \item{`"at_risk"`}{Number of subjects with a row at `k` (still
#'       under observation through interval `k`).}
#'     \item{`"events_y"`}{Number of `Y` events in interval `k`
#'       (`y_event == 1`).}
#'     \item{`"censored"`}{Number of subjects censored in interval
#'       `k` (`dep_cens == 1` or `indep_cens == 1`).}
#'   }
#'
#' @return A data.frame with column `k` and one column per observed
#'   treatment value, named e.g. `A_0` and `A_1`.
#'
#' @seealso [causal_survival()], [plot.causal_survival_risk()]
#' @family accessors
#' @export
causal_risk_table <- function(fit,
                              count = c("at_risk", "events_y",
                                        "censored")) {
  stopifnot(inherits(fit, "causal_survival_fit"))
  count <- match.arg(count)
  if (is.null(fit$pt_data)) {
    stop("`fit$pt_data` is NULL; refit with `keep_data = TRUE` to ",
         "compute a risk table.", call. = FALSE)
  }
  risk_table_internal(
    pt_data   = fit$pt_data,
    id_col    = fit$id_col,
    trt_col   = fit$treatment_col,
    cut_times = fit$cut_times,
    count     = count
  )
}


#' Internal risk-table worker (no class check)
#'
#' Used by both [causal_risk_table()] and the plot method's optional
#' table panel. The split lets the plot method compute the table from
#' a `causal_survival_risk` object's stored references without
#' needing a back-reference to the full fit.
#'
#' @param pt_data Person-time data.frame (must include `k`, the
#'   treatment column, `y_event`, `dep_cens`, `indep_cens`).
#' @param id_col,trt_col Character column names.
#' @param cut_times Numeric vector of cut times. The function looks
#'   up rows by the integer interval index `1..length(cut_times)`,
#'   then reports counts aligned with the cut-time values.
#' @param count One of `"at_risk"`, `"events_y"`, `"censored"`.
#' @return A data.frame with `k` (= cut time value) plus one
#'   per-arm count column.
#' @family internal
#' @keywords internal
risk_table_internal <- function(pt_data, id_col, trt_col,
                                cut_times, count) {
  trt_vals <- sort(unique(pt_data[[trt_col]]))
  K_max    <- length(cut_times)

  result <- data.frame(k = cut_times)
  for (a in trt_vals) {
    col_name <- paste0(trt_col, "_", a)
    vals <- vapply(seq_len(K_max), function(k_int) {
      rows_a_k <- pt_data[pt_data[[trt_col]] == a &
                            pt_data$k == k_int, , drop = FALSE]
      if (count == "at_risk") {
        length(unique(rows_a_k[[id_col]]))
      } else if (count == "events_y") {
        as.integer(sum(rows_a_k$y_event, na.rm = TRUE))
      } else {
        as.integer(sum(rows_a_k$dep_cens, na.rm = TRUE) +
                     sum(rows_a_k$indep_cens, na.rm = TRUE))
      }
    }, integer(1))
    result[[col_name]] <- vals
  }
  result
}
