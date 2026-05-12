#' Print a causal_survival_fit object
#'
#' One-screen summary: estimand framing, method (+ IPW engine when
#' relevant), cohort size, time grid, per-arm cumulative incidence at
#' the final cut time, and the warning / model-check tallies.
#'
#' @param x A `"causal_survival_fit"` object.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_survival_fit <- function(x, ...) {
  cat("Counterfactual survival fit (causal_survival_fit)\n")
  cat("-------------------------------------------------\n")
  cat("Estimand: E[Y_k^{a, c = 0}] under exchangeability ",
      "conditional on baseline L_0\n", sep = "")
  cat("See `causal_assumptions(fit)` for the full identification block.\n")

  engine_str <- if (x$method == "ipw")
    paste0(" (engine: ", x$ipw_engine, ")") else ""
  cat("Method: ", x$method, engine_str, "\n", sep = "")
  cat("N subjects: ", length(unique(x$pt_data[[x$id_col]])), "\n", sep = "")
  cat("Cut times: ", length(x$cut_times),
      " (T_max = ", max(x$cut_times), ")\n", sep = "")

  # Per-arm cumulative incidence at final cut time
  est <- x$cumulative_incidence[[x$method]]
  if (!is.null(est)) {
    K_max <- max(est$k)
    last  <- est[est$k == K_max, , drop = FALSE]
    cat("\nCumulative incidence at k = ", K_max,
        " (t = ", max(x$cut_times), "):\n", sep = "")
    for (i in seq_len(nrow(last))) {
      cat(sprintf("  a = %s : F^a = %.4f\n",
                  format(last$treatment[i]), last$inc[i]))
    }
  }

  # Model-checks tally
  if (!is.null(x$model_checks)) {
    issues <- 0L
    for (chk in x$model_checks) {
      if (is.null(chk)) next
      if (!isTRUE(chk$converged)) issues <- issues + 1L
      if (length(chk$glm_warnings) > 0L) issues <- issues + 1L
    }
    if (issues > 0L) {
      cat("\nModel checks: ", issues,
          " issue(s) - use `fit$model_checks` to inspect.\n", sep = "")
    }
  }

  # Warnings - spec §3.5 line 537
  if (length(x$warnings) > 0L) {
    cat("\nFit completed with ", length(x$warnings),
        " warning(s) (see fit$warnings).\n", sep = "")
  }

  cat("\nUse causal_risk(), causal_contrast(), bootstrap() to extract components.\n")
  invisible(x)
}


#' Summary of a causal_survival_fit object
#'
#' Per-arm cumulative incidence at the selected cut time, the
#' model-check tally, and (when a bootstrap is supplied) the RD/RR
#' contrast at the same time. The identification block is delegated
#' to [causal_assumptions()] - not inlined here.
#'
#' @param object A `"causal_survival_fit"` object.
#' @param ci Optional. A `"causal_survival_bootstrap"` object from
#'   [bootstrap()]. When supplied, the contrast row at the selected
#'   `time` is printed beneath the per-arm risk block.
#' @param time `NULL` (default, resolves to the final cut time per
#'   spec §3.4) or a numeric scalar resolved on the reporting grid
#'   via [snap_time()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns the cumulative incidence rows at the
#'   selected time.
#' @export
summary.causal_survival_fit <- function(object, ci = NULL,
                                        time = NULL, ...) {
  if (!is.null(ci)) {
    stopifnot(inherits(ci, "causal_survival_bootstrap"))
  }
  t_at <- snap_time(time, object$cut_times)
  k_at <- which(object$cut_times == t_at)

  cat("Counterfactual survival - summary\n")
  cat("=================================\n\n")
  engine_str <- if (object$method == "ipw")
    paste0(" (engine: ", object$ipw_engine, ")") else ""
  cat("Method: ", object$method, engine_str,
      " | N: ", length(unique(object$pt_data[[object$id_col]])),
      "\n", sep = "")
  cat("Cut times: ", length(object$cut_times),
      " (T_max = ", max(object$cut_times), ")\n\n", sep = "")

  est <- object$cumulative_incidence[[object$method]]
  row <- est[est$k == k_at, , drop = FALSE]
  cat(sprintf("Cumulative incidence at k = %d (t = %g):\n", k_at, t_at))
  for (i in seq_len(nrow(row))) {
    cat(sprintf("  a = %s : F^a = %.4f\n",
                format(row$treatment[i]), row$inc[i]))
  }

  # Contrast block (only when bootstrap attached)
  if (!is.null(ci)) {
    ctr <- causal_contrast(object, ci = ci, time = t_at)$contrasts
    cat(sprintf(
      "\nContrasts at k = %d (%.0f%% CIs):\n", k_at, (1 - ci$alpha) * 100
    ))
    for (i in seq_len(nrow(ctr))) {
      r <- ctr[i, ]
      label <- if (r$op == "-") "RD" else "RR"
      cat(sprintf(
        "  %s  %s = %6.3f  [%6.3f, %6.3f]\n",
        r$name, label, r$estimate, r$lower, r$upper
      ))
    }
  }

  # Model-check tally
  if (!is.null(object$model_checks)) {
    non_converged <- 0L
    for (chk in object$model_checks) {
      if (is.null(chk)) next
      if (!isTRUE(chk$converged)) non_converged <- non_converged + 1L
    }
    if (non_converged > 0L) {
      cat(sprintf(
        "\nModel checks: %d non-converged model(s). See `fit$model_checks`.\n",
        non_converged
      ))
    }
  }

  # Pointers
  cat("\nIdentification block: causal_assumptions(fit).\n")
  if (is.null(ci)) {
    cat("For contrasts with CIs: boot <- bootstrap(fit, n_boot = 500); ",
        "summary(fit, ci = boot)\n", sep = "")
  }
  invisible(row)
}


#' Print a causal_survival_risk object
#'
#' Per-method per-arm value at the final cut time on the requested
#' scale; indicates whether bootstrap bands are attached.
#'
#' @param x A `"causal_survival_risk"` object from [causal_risk()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_survival_risk <- function(x, ...) {
  cat("Counterfactual ", x$scale,
      " curves (causal_survival_risk)\n", sep = "")
  cat("------------------------------------------------\n")

  methods_avail <- unique(x$risk$method)
  cat("Methods: ", paste(methods_avail, collapse = ", "), "\n", sep = "")
  cat("Bootstrap CIs: ",
      if (any(!is.na(x$risk$lower))) "yes" else "no", "\n\n", sep = "")

  value_label <- if (x$scale == "incidence") "F^a" else "S^a"

  for (m in methods_avail) {
    sub <- x$risk[x$risk$method == m, , drop = FALSE]
    K_max <- max(sub$k)
    last  <- sub[sub$k == K_max, , drop = FALSE]
    cat(sprintf("[%s] at final time (k = %d, t = %g):\n",
                m, K_max, last$time[1]))
    for (i in seq_len(nrow(last))) {
      band <- if (is.na(last$lower[i])) "" else
        sprintf("  [%.4f, %.4f]", last$lower[i], last$upper[i])
      cat(sprintf("  a = %s : %s = %.4f%s\n",
                  format(last$treatment[i]), value_label,
                  last$value[i], band))
    }
    cat("\n")
  }
  cat("Use plot(causal_risk(fit)) to visualize.\n")
  invisible(x)
}


#' Print a causal_survival_contrast object
#'
#' Renders the contrast table at the selected time, with the method
#' and significance level.
#'
#' @param x A `"causal_survival_contrast"` object from
#'   [causal_contrast()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_survival_contrast <- function(x, ...) {
  cat("Counterfactual contrasts (causal_survival_contrast)\n")
  cat("---------------------------------------------------\n")
  cat("Method: ", x$method, "\n", sep = "")
  if (!is.null(x$alpha)) {
    cat(sprintf("Significance level: %g (%.0f%% CIs)\n\n",
                x$alpha, (1 - x$alpha) * 100))
  } else {
    cat("Significance level: - (no bootstrap supplied)\n\n")
  }

  cat(sprintf("At t = %g:\n", x$time))
  out <- x$contrasts[, c("name", "op", "estimate", "lower", "upper"),
                     drop = FALSE]
  print(out, row.names = FALSE)
  invisible(x)
}


#' Print a causal_survival_bootstrap object
#'
#' Replicate count (requested vs effective), significance level,
#' failed-replicate count, and a pointer to the accessor that pairs
#' the bands with a fit.
#'
#' @param x A `"causal_survival_bootstrap"` object from
#'   [bootstrap()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_survival_bootstrap <- function(x, ...) {
  cat("Bootstrap confidence bands (causal_survival_bootstrap)\n")
  cat("------------------------------------------------------\n")
  cat("Replicates requested: ", x$n_boot_requested, "\n", sep = "")
  cat("Replicates effective: ", x$n_boot_effective, "\n", sep = "")
  if (length(x$failed_reps) > 0L) {
    cat("Failed replicates: ", length(x$failed_reps),
        " (see $failed_reps for indices)\n", sep = "")
  }
  if (!is.na(x$warnings_count) && x$warnings_count > 0L) {
    cat("Warnings inside replicates: ", x$warnings_count, "\n", sep = "")
  }
  cat(sprintf("Significance level: %g (%.0f%% CIs)\n",
              x$alpha, (1 - x$alpha) * 100))
  cat("\nUse `causal_contrast(fit, ci = <this>)` for contrast bands,\n")
  cat("or `plot(causal_risk(fit, ci = <this>))` for curve bands.\n")
  invisible(x)
}


#' Confidence intervals for causal_survival_fit (intentionally not provided)
#'
#' This method exists to redirect callers to the supported pattern.
#' Confidence intervals are not stored on the fit object; they live
#' on a separate `"causal_survival_bootstrap"` object that pairs with
#' the fit at accessor time.
#'
#' @param object A `"causal_survival_fit"` object.
#' @param parm,level,... Unused.
#' @return Always errors - use the [bootstrap()] + [causal_contrast()]
#'   pattern instead.
#' @export
confint.causal_survival_fit <- function(object, parm = NULL,
                                        level = 0.95, ...) {
  stop(
    "Confidence intervals are not stored on `fit` in this package. ",
    "Pair the fit with a bootstrap object:\n  ",
    "boot <- bootstrap(fit, n_boot = 500)\n  ",
    "causal_contrast(fit, ci = boot)",
    call. = FALSE
  )
}
