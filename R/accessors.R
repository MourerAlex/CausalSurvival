#' Cumulative incidence / survival under each arm
#'
#' Extract the counterfactual cumulative incidence
#' \eqn{F^a(k) = E[Y_k^{a, c = 0}]} (or the survival
#' \eqn{S^a(k) = 1 - F^a(k)}) under each level `a` of the binary
#' treatment, indexed by the reporting grid `fit$cut_times`. Optionally
#' folds bootstrap confidence bands into the same table.
#'
#' @param fit A `"causal_survival_fit"` object from
#'   [causal_survival()].
#' @param scale One of `"incidence"` (default, returns
#'   \eqn{F^a(k)}) or `"survival"` (returns \eqn{S^a(k)}).
#' @param ci Optional. A `"causal_survival_bootstrap"` object from
#'   [bootstrap()]. When provided, the `lower` / `upper` columns of
#'   `$risk` are populated; otherwise they are `NA_real_`.
#'
#' @return An S3 object of class `"causal_survival_risk"` with:
#'   \describe{
#'     \item{risk}{Long-format data.frame with columns `method`,
#'       `treatment`, `k`, `time`, `value`, `lower`, `upper`. One row
#'       per `(method, treatment, k)` triple.}
#'     \item{scale}{The selected scale.}
#'     \item{replicates, alpha}{Carried from `ci` for per-contrast
#'       bands in [plot.causal_survival_risk()] (or `NULL`).}
#'     \item{pt_data, id_col, treatment_col, cut_times}{References to
#'       the fit, used by the plot method's optional risk-table panel.}
#'   }
#'
#' @seealso [causal_contrast()], [plot.causal_survival_risk()]
#' @family accessors
#' @export
causal_risk <- function(fit, scale = c("incidence", "survival"),
                        ci = NULL) {
  stopifnot(inherits(fit, "causal_survival_fit"))
  scale <- match.arg(scale)
  if (!is.null(ci)) {
    stopifnot(inherits(ci, "causal_survival_bootstrap"))
  }
  structure(
    list(
      risk          = build_risk_long(fit, ci, scale),
      scale         = scale,
      replicates    = if (!is.null(ci)) ci$replicates else NULL,
      alpha         = if (!is.null(ci)) ci$alpha      else NULL,
      pt_data       = fit$pt_data,
      id_col        = fit$id_col,
      treatment_col = fit$treatment_col,
      cut_times     = fit$cut_times
    ),
    class = "causal_survival_risk"
  )
}


#' Build the long-format `$risk` data.frame
#'
#' Pivot `fit$cumulative_incidence` (per-method long table with
#' columns `treatment`, `k`, `time`, `surv`, `inc`) and optionally
#' `ci$ci_curves` (per-method long bands) into the canonical accessor
#' shape used by [causal_risk()].
#'
#' @param fit A `"causal_survival_fit"` object.
#' @param ci A `"causal_survival_bootstrap"` object or `NULL`.
#' @param scale One of `"incidence"` or `"survival"`.
#' @return Long-format data.frame with columns `method`, `treatment`,
#'   `k`, `time`, `value`, `lower`, `upper`. `lower` / `upper` are
#'   `NA_real_` when `ci` is `NULL`.
#' @family internal
#' @keywords internal
build_risk_long <- function(fit, ci, scale) {
  cum_inc_list <- fit$cumulative_incidence
  value_col    <- if (scale == "incidence") "inc" else "surv"

  rows <- list()
  for (m in names(cum_inc_list)) {
    est <- cum_inc_list[[m]]
    if (is.null(est)) next  # method not fit this run

    out <- data.frame(
      method    = m,
      treatment = est$treatment,
      k         = est$k,
      time      = est$time,
      value     = est[[value_col]],
      lower     = NA_real_,
      upper     = NA_real_,
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    # Join in bootstrap bands when present. Per spec §4.3 the
    # bootstrap object stores `ci_lower` / `ci_upper` as
    # data.frames with columns (treatment, k, lower / upper) for the
    # method that produced the fit; we match by (treatment, k).
    if (!is.null(ci)) {
      key_out <- paste(out$treatment, out$k, sep = "|")
      key_lo  <- paste(ci$ci_lower$treatment, ci$ci_lower$k, sep = "|")
      key_up  <- paste(ci$ci_upper$treatment, ci$ci_upper$k, sep = "|")
      out$lower <- ci$ci_lower$lower[match(key_out, key_lo)]
      out$upper <- ci$ci_upper$upper[match(key_out, key_up)]
      # When `scale = "survival"`, the cumulative incidence bands
      # map to survival bands via `S = 1 - F` — flipped order.
      if (scale == "survival") {
        flipped_lower <- 1 - out$upper
        flipped_upper <- 1 - out$lower
        out$lower <- flipped_lower
        out$upper <- flipped_upper
      }
    }
    rows[[length(rows) + 1L]] <- out
  }
  if (length(rows) == 0L) {
    return(data.frame(
      method = character(), treatment = integer(),
      k = integer(), time = numeric(),
      value = numeric(), lower = numeric(), upper = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}


#' Counterfactual contrasts on the reporting grid
#'
#' Compute pairwise contrasts of the counterfactual cumulative
#' incidence between treatment arms on the reporting grid. The risk
#' difference `RD(k) = F^{a}(k) - F^{a'}(k)` (operator `"-"`) and the
#' risk ratio `RR(k) = F^{a}(k) / F^{a'}(k)` (operator `"/"`) are the
#' two default operators.
#'
#' Under the v0.1.0 binary point-treatment scope there is a single
#' arm pair `(max, min)` of `fit$treatment_levels`. The `reference` /
#' `contrasts` arguments are carried in the signature for forward
#' compatibility with multi-arm (v2) and separable-effects (downstream
#' package) scopes; in v0.1.0 supplying a non-`NULL` `contrasts`
#' raises an error.
#'
#' @param fit A `"causal_survival_fit"` object from
#'   [causal_survival()].
#' @param reference Optional reference arm (one of the levels in
#'   `fit$treatment_levels`). `NULL` defaults to the lowest level.
#' @param contrasts Optional named list of custom contrasts (deferred
#'   to v2 multi-arm; see `dev/unused_code/resolve_contrast_pairs_v2.R`).
#'   Must be `NULL` in v0.1.0.
#' @param ci Optional. A `"causal_survival_bootstrap"` object from
#'   [bootstrap()]. When `NULL`, point estimates only are returned
#'   with `lower` / `upper` set to `NA_real_`. A loud `warning()` is
#'   emitted in that case: contrasts without confidence intervals are
#'   not meaningful and their interpretation is strongly discouraged.
#' @param time `NULL` (default, resolves to the final cut time per
#'   spec §3.4) or a numeric scalar resolved on the reporting grid
#'   via [snap_time()]. Out-of-bounds values raise an error.
#'
#' @return An S3 object of class `"causal_survival_contrast"` with:
#'   \describe{
#'     \item{contrasts}{Long-format data.frame with columns `method`,
#'       `name`, `treatment_a`, `treatment_b`, `op`, `k`, `time`,
#'       `estimate`, `lower`, `upper`.}
#'     \item{method}{The fit's estimation method.}
#'     \item{alpha}{Bootstrap significance level (or `NULL`).}
#'     \item{time}{The resolved cut time.}
#'   }
#'
#' @seealso [causal_risk()], [plot.causal_survival_contrast()],
#'   [bootstrap()]
#' @family accessors
#' @export
causal_contrast <- function(fit, reference = NULL, contrasts = NULL,
                            ci = NULL, time = NULL) {
  stopifnot(inherits(fit, "causal_survival_fit"))

  if (is.null(ci)) {
    warning(
      "`ci = NULL`: causal_contrast() is returning point estimates ",
      "without confidence intervals. Contrasts without uncertainty ",
      "are not meaningful — their interpretation is strongly ",
      "discouraged. Compute a bootstrap first:\n  ",
      "boot <- bootstrap(fit, n_boot = 500)\n  ",
      "Then: causal_contrast(fit, ci = boot)",
      call. = FALSE
    )
  } else {
    stopifnot(inherits(ci, "causal_survival_bootstrap"))
  }

  if (!is.null(contrasts)) {
    stop(
      "`contrasts` (custom list) is not supported in v0.1.0 ",
      "(binary point treatment only). The drafted helper is parked ",
      "at dev/unused_code/resolve_contrast_pairs_v2.R; see TODO ",
      "'Multi-arm treatment'.",
      call. = FALSE
    )
  }

  levels_vec <- fit$treatment_levels
  if (length(levels_vec) != 2L) {
    stop("v0.1.0 supports binary treatment only; ",
         "fit$treatment_levels has length ", length(levels_vec), ".",
         call. = FALSE)
  }
  if (is.null(reference)) reference <- min(levels_vec)
  if (!reference %in% levels_vec) {
    stop("`reference` must be one of the fit's treatment levels: ",
         paste(levels_vec, collapse = ", "), ".", call. = FALSE)
  }

  # Binary inline expansion: one pair, both operators.
  comparator <- setdiff(levels_vec, reference)
  pairs <- list(
    name        = rep(paste0(comparator, "_vs_", reference), 2L),
    treatment_a = rep(comparator, 2L),
    treatment_b = rep(reference, 2L),
    op          = c("-", "/")
  )

  out_df <- compute_contrast_table(fit, pairs, ci)

  # Spec §3.4: `time = NULL` resolves to the final cut time.
  t_at <- snap_time(time, fit$cut_times)
  out_df <- out_df[out_df$time == t_at, , drop = FALSE]

  structure(
    list(
      contrasts = out_df,
      method    = fit$method,
      alpha     = if (!is.null(ci)) ci$alpha else NULL,
      time      = t_at
    ),
    class = "causal_survival_contrast"
  )
}


#' Build the long-format `$contrasts` data.frame
#'
#' Iterate over the methods present in `fit$cumulative_incidence` and
#' the arm pairs in `pairs`, computing the per-`k` contrast under the
#' requested operator. Bootstrap bands (when supplied) are obtained
#' by computing the contrast per replicate and taking the requested
#' lower/upper quantile.
#'
#' @param fit A `"causal_survival_fit"` object.
#' @param pairs A list with parallel vectors `name`, `treatment_a`,
#'   `treatment_b`, `op` — one element per (pair x operator) cell.
#' @param ci A `"causal_survival_bootstrap"` object or `NULL`.
#' @return Long-format data.frame with columns `method`, `name`,
#'   `treatment_a`, `treatment_b`, `op`, `k`, `time`, `estimate`,
#'   `lower`, `upper`. `lower` / `upper` are `NA_real_` when `ci` is
#'   `NULL`.
#' @family internal
#' @keywords internal
compute_contrast_table <- function(fit, pairs, ci) {
  alpha <- if (!is.null(ci)) ci$alpha else NA_real_
  reps  <- if (!is.null(ci)) ci$replicates else NULL

  out <- list()
  for (m in names(fit$cumulative_incidence)) {
    est <- fit$cumulative_incidence[[m]]
    if (is.null(est)) next  # method not fit this run

    for (j in seq_along(pairs$name)) {
      a   <- pairs$treatment_a[j]
      b   <- pairs$treatment_b[j]
      op  <- pairs$op[j]
      nm  <- pairs$name[j]

      F_a <- est$inc[est$treatment == a]
      F_b <- est$inc[est$treatment == b]

      estimate <- contrast_op(F_a, F_b, op)

      if (!is.null(reps)) {
        # `reps` is long-format `data.frame(boot_id, treatment, k,
        # value)` per spec §4.3. Reshape to a `[k, boot_id]` matrix
        # per arm, compute the contrast per (k, boot_id), then take
        # quantiles along the boot_id axis.
        wide_a <- tapply(
          reps$value[reps$treatment == a],
          list(reps$k[reps$treatment == a],
               reps$boot_id[reps$treatment == a]),
          identity
        )
        wide_b <- tapply(
          reps$value[reps$treatment == b],
          list(reps$k[reps$treatment == b],
               reps$boot_id[reps$treatment == b]),
          identity
        )
        per_rep <- contrast_op(wide_a, wide_b, op)
        lower   <- apply(per_rep, 1L, stats::quantile,
                         probs = alpha / 2,       na.rm = TRUE)
        upper   <- apply(per_rep, 1L, stats::quantile,
                         probs = 1 - alpha / 2,   na.rm = TRUE)
      } else {
        lower <- rep(NA_real_, length(estimate))
        upper <- rep(NA_real_, length(estimate))
      }

      out[[length(out) + 1L]] <- data.frame(
        method      = m,
        name        = nm,
        treatment_a = a,
        treatment_b = b,
        op          = op,
        k           = est$k[est$treatment == a],
        time        = est$time[est$treatment == a],
        estimate    = estimate,
        lower       = lower,
        upper       = upper,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }
  if (length(out) == 0L) {
    return(data.frame(
      method = character(), name = character(),
      treatment_a = numeric(), treatment_b = numeric(),
      op = character(), k = integer(), time = numeric(),
      estimate = numeric(), lower = numeric(), upper = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, out)
}


#' Tiny operator dispatcher used by [compute_contrast_table()]
#'
#' Applies the contrast operator (`"-"` or `"/"`) to either a pair of
#' numeric vectors (point estimate) or a pair of matrices indexed by
#' `(k, boot_id)` (per-replicate bands). Kept as a single function so
#' the two shapes share one implementation.
#'
#' @param x,y Numeric vector or matrix.
#' @param op One of `"-"` or `"/"`.
#' @return Numeric of the same shape as `x` / `y`.
#' @family internal
#' @keywords internal
contrast_op <- function(x, y, op) {
  switch(op,
    "-" = x - y,
    "/" = x / y,
    stop("Unknown contrast operator: ", op, call. = FALSE)
  )
}


#' Resolve a user-supplied time to the reporting grid
#'
#' Used by [causal_contrast()] and [summary.causal_survival_fit()] to
#' interpret the optional `time` argument. `NULL` resolves to the
#' final cut time. A numeric value outside `[min(cut_times),
#' max(cut_times)]` is rejected with a hard error (spec §3.5). An
#' in-bounds value is snapped to the nearest entry of `cut_times`;
#' a `message()` is emitted if snapping changes it.
#'
#' @param time `NULL` or a numeric scalar.
#' @param cut_times Numeric vector of available cut times.
#' @return A single numeric value drawn from `cut_times`.
#' @family internal
#' @keywords internal
snap_time <- function(time, cut_times) {
  if (is.null(time)) return(max(cut_times))
  if (!is.numeric(time) || length(time) != 1L || is.na(time)) {
    stop("`time` must be NULL or a single non-missing numeric.",
         call. = FALSE)
  }
  if (time > max(cut_times) || time < min(cut_times)) {
    stop(sprintf(
      "`time = %g` is outside the reporting grid [%g, %g].",
      time, min(cut_times), max(cut_times)
    ), call. = FALSE)
  }
  idx <- which.min(abs(cut_times - time))
  k_at <- cut_times[idx]
  if (!isTRUE(all.equal(k_at, time))) {
    message(sprintf("`time = %g` snapped to nearest cut time: %g.",
                    time, k_at))
  }
  k_at
}
