#' Check Fitted Probabilities (Positivity Signal)
#'
#' Inspects predicted probabilities from a fitted glm and returns a
#' diagnostics list with the continuous signal (min and max fitted). Emits
#' warnings for: non-convergence, NA fitted values, and predicted
#' probabilities near 0 or 1 (the user interprets the continuous extremes;
#' no binary positivity-violation flag is reported because any specific
#' threshold would be arbitrary).
#'
#' @param model A fitted glm object.
#' @param model_label Character label for the model (e.g., "C-hazard").
#' @param warn_eps Threshold below 0 (and above 1-warn_eps) for emitting an
#'   "extreme probabilities predicted" warning. Default 1e-6. Only controls
#'   whether a warning fires — the continuous min/max are always returned.
#' @param glm_warnings Character vector of glm warnings captured at fit time.
#'
#' @return A list with:
#'   \describe{
#'     \item{label}{Model label.}
#'     \item{converged}{Logical.}
#'     \item{min_fitted, max_fitted}{Extremes of fitted probabilities (or NA
#'       if fitted values themselves contain NAs). Always reported so users
#'       can judge positivity for themselves.}
#'     \item{glm_warnings}{Character vector passed through from the fit.}
#'   }
#' @keywords internal
check_fitted_positivity <- function(model, model_label,
                                    warn_eps = 1e-6,
                                    glm_warnings = character()) {

  converged <- isTRUE(model$converged)
  if (!converged) {
    warning(
      model_label, " model did not converge.",
      call. = FALSE
    )
  }

  probs <- stats::fitted(model)

  if (any(is.na(probs))) {
    warning(
      model_label, " model has NA fitted values - fit may have failed.",
      call. = FALSE
    )
    return(list(
      label = model_label,
      converged = converged,
      min_fitted = NA_real_,
      max_fitted = NA_real_,
      glm_warnings = glm_warnings
    ))
  }

  min_p <- min(probs)
  max_p <- max(probs)

  # Warn about extreme probabilities (user interprets what "extreme" means).
  # No binary positivity-violation flag is stored — the min/max are the
  # continuous signal, and any cutoff would be arbitrary.
  if (min_p < warn_eps || max_p > 1 - warn_eps) {
    warning(
      model_label,
      " model predicted extreme probabilities ",
      "(min=", signif(min_p, 3), ", max=", signif(max_p, 3), "). ",
      "Inspect diagnostic(fit) to assess positivity.",
      call. = FALSE
    )
  }

  list(
    label = model_label,
    converged = converged,
    min_fitted = min_p,
    max_fitted = max_p,
    glm_warnings = glm_warnings
  )
}


#' Fit a Logistic GLM with Warning Capture and Positivity Check
#'
#' Internal helper. Fits a binomial logistic glm, captures any `glm()`
#' warnings via `withCallingHandlers()`, re-emits them so the calling
#' fitter's collector grabs them, and builds a diagnostics list via
#' [check_fitted_positivity()].
#'
#' Used by [fit_hazard_models()] (Y, C hazards on person-time), by
#' [fit_propensity()] (treatment model on baseline rows), and by the IPW
#' worker for the weighted Y-MSM fit.
#'
#' @param formula A fitted model formula.
#' @param data Data frame passed to `glm()`.
#' @param label Human-readable label for warnings and diagnostics.
#' @param weights Optional numeric vector of `glm` case weights, length
#'   `nrow(data)`. NULL (default) fits an unweighted glm. When supplied,
#'   the vector is stashed on `data` under an internal column name
#'   (`.fit_logistic_w_`) and passed to `glm()` by that name; this
#'   side-steps NSE on the `weights` formal.
#'
#' @return A list with `model` (the glm object) and `check` (diagnostics).
#' @keywords internal
fit_logistic <- function(formula, data, label, weights = NULL) {
  glm_warnings <- character()

  handler <- function(w) {
    glm_warnings <<- c(glm_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }

  if (is.null(weights)) {
    model <- withCallingHandlers(
      stats::glm(formula, data = data,
                 family = stats::binomial(link = "logit")),
      warning = handler
    )
  } else {
    if (!is.numeric(weights) || length(weights) != nrow(data)) {
      stop("`weights` must be a numeric vector of length nrow(data) (",
           nrow(data), "); got length ", length(weights), ".",
           call. = FALSE)
    }
    data$.fit_logistic_w_ <- weights
    model <- withCallingHandlers(
      stats::glm(formula, data = data,
                 family = stats::binomial(link = "logit"),
                 weights = .fit_logistic_w_),
      warning = handler
    )
  }

  # Re-emit glm warnings so they are collected by the calling fitter's handler
  for (w in glm_warnings) warning(w, call. = FALSE)

  check <- check_fitted_positivity(model, label, glm_warnings = glm_warnings)

  list(model = model, check = check)
}


#' Predict with `withCallingHandlers` Warning Capture
#'
#' Internal helper wrapping `stats::predict()` so that any warnings (e.g.,
#' rank-deficient fit) are captured and re-emitted with a model label
#' prefix. This keeps the surfaces visible to upstream warning collectors.
#'
#' @param model A fitted glm object.
#' @param newdata Data frame to predict on.
#' @param label Character label prepended to re-emitted warnings.
#'
#' @return Numeric vector of predicted probabilities.
#' @keywords internal
predict_with_warning <- function(model, newdata, label) {
  captured <- character()
  probs <- withCallingHandlers(
    stats::predict(model, newdata = newdata, type = "response"),
    warning = function(w) {
      captured <<- c(captured,
                     paste0(label, " predict: ", conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  for (msg in captured) warning(msg, call. = FALSE)
  probs
}


#' Predict Hazard Under a Counterfactual Treatment Assignment
#'
#' Replaces `treatment_var` in `data` with `treatment_value`, then predicts
#' the hazard with warning capture via [predict_with_warning()]. Returns
#' NULL if `model` is NULL (caller did not request this hazard component).
#' Knows nothing about Y or any specific causal framework.
#'
#' @param model Fitted hazard glm (binomial), or NULL.
#' @param data Person-time data frame.
#' @param treatment_var Character. Name of the column to overwrite.
#' @param treatment_value Numeric. Counterfactual value to assign.
#' @param label Character label prepended to any captured predict warnings.
#' @return Numeric vector of predicted hazards, or NULL.
#' @keywords internal
predict_counterfactual_hazard <- function(model, data, treatment_var,
                                          treatment_value, label) {
  if (is.null(model)) return(NULL)
  newdata <- data
  newdata[[treatment_var]] <- treatment_value
  predict_with_warning(model, newdata, label)
}


#' Fit Hazard Models
#'
#' Fits pooled logistic regression models for the Y-hazard and (optionally)
#' censoring hazard on person-time data. Per-model diagnostics are collected
#' via [check_fitted_positivity()].
#'
#' @param pt_data Data frame in person-time format.
#' @param treatment Character. Treatment column name. Used directly in both
#'   Y- and C-hazard formulas (no working-copy column).
#' @param covariates Character vector. Covariate column names.
#' @param active_methods Character vector. Subset of
#'   `c("g_formula", "ipw")`. Determines which models get fit:
#'   - `"g_formula"`: model_y only
#'   - `"ipw"`:       model_c (only when `ipcw = TRUE`)
#' @param formulas Named list or NULL. User-specified formulas (names `y`,
#'   `c`). Any entry absent falls back to the default.
#' @param ipcw Logical. When FALSE, the censoring model is not fit and
#'   `model_c` stays NULL.
#'
#' @return A list with two entries:
#'   \describe{
#'     \item{models}{Named list: `model_y`, `model_c` (glm objects or NULL).}
#'     \item{checks}{Named list: `y`, `c` (per-model diagnostics or NULL).}
#'   }
#'
#' @details
#' Default formula: `<event> ~ treatment + k + I(k^2) + I(k^3) + covariates`
#' (additive, no interactions). The integer-`k` polynomial in time follows
#' the LOCKED `(0, t_1], ..., (t_{K_max-1}, T_max]` convention from spec
#' §3.0.2 — `k` is the interval index, not its left-edge time.
#'
#' Fit populations follow spec §3.0.6:
#' - Y-hazard (`y_event ~ ...`): rows with `indep_cens == 0 & dep_cens == 0`.
#' - C-hazard (`dep_cens ~ ...`): rows with `indep_cens == 0`.
#'
#' The Y-hazard model fit here is **unweighted**.
#'
#' @keywords internal
fit_hazard_models <- function(pt_data,
                              treatment,
                              covariates,
                              active_methods,
                              formulas,
                              ipcw = TRUE) {

  cov_terms <- if (length(covariates) > 0) {
    paste(covariates, collapse = " + ")
  } else {
    NULL
  }
  time_terms <- "k + I(k^2) + I(k^3)"

  models <- list(model_y = NULL, model_c = NULL)
  checks <- list(y = NULL, c = NULL)

  # --- Y-hazard model (g-formula path; IPW-MSM Y fit is weighted, done later) ---
  # Fit population: rows with indep_cens == 0 AND dep_cens == 0 (spec §3.0.6).
  if ("g_formula" %in% active_methods) {
    y_rows <- pt_data$indep_cens == 0 & pt_data$dep_cens == 0
    fml_y <- formulas$y %||% stats::as.formula(
      paste("y_event ~", paste(c(treatment, time_terms, cov_terms),
                               collapse = " + "))
    )
    fit_result <- fit_logistic(fml_y, pt_data[y_rows, , drop = FALSE],
                               "Y-hazard")
    models$model_y <- fit_result$model
    checks$y <- fit_result$check
  }

  # --- Censoring model (IPW path with ipcw) ---
  # Fit population: rows with indep_cens == 0 (spec §3.0.6). The
  # structural exclusion of indep_cens == 1 rows is what makes the
  # IPCW cliff impossible (no admin-style rows in the fit).
  if ("ipw" %in% active_methods && ipcw) {
    c_rows <- pt_data$indep_cens == 0
    fml_c <- formulas$c %||% stats::as.formula(
      paste("dep_cens ~", paste(c(treatment, time_terms, cov_terms),
                                collapse = " + "))
    )
    fit_result <- fit_logistic(fml_c, pt_data[c_rows, , drop = FALSE],
                               "C-hazard")
    models$model_c <- fit_result$model
    checks$c <- fit_result$check
  }

  list(models = models, checks = checks)
}


#' Per-Subject Cumulative Product of Survival Probability
#'
#' Computes prod_{j <= s}(1 - haz_j) within each subject id. Returned vector
#' is aligned to the input rows. A core survival-analysis primitive.
#'
#' @param haz Numeric vector of hazards.
#' @param id Subject id vector (same length as `haz`). Rows MUST be sorted
#'   by time within each subject id; this function does not re-sort.
#' @return Numeric vector of cumulative survival, aligned to input rows.
#' @keywords internal
cumprod_survival <- function(haz, id) {
  ave(1 - haz, id, FUN = cumprod)
}


#' Weighted Hazard by Interval (Hajek Ratio)
#'
#' Per-interval weighted incidence: numerator is the weighted count of
#' events at each k, denominator is the weighted count at risk at each
#' k. At-risk rows are those with non-NA `event` (the convention is
#' that NA means the subject is no longer at risk for this event type
#' at this k). Framework-agnostic vector primitive used by the
#' `ipw_engine = "km"` survival estimator.
#'
#' Ported from CausalCompetingRisks (`R/ipw_core.R`), un-deferred
#' 2026-05-12 to support `fit_ipw_km()`.
#'
#' @param event Numeric vector. Event indicator at each row: 1 = event,
#'   0 = at risk no event, NA = no longer at risk.
#' @param k Vector of interval indices, same length as `event`.
#' @param weights Numeric vector of per-row weights, same length.
#' @return Named numeric vector: weighted hazard at each unique value
#'   of `k`, names are the k values as character.
#' @keywords internal
weighted_hazard_by_k <- function(event, k, weights) {
  at_risk <- !is.na(event)
  event   <- event[at_risk]
  k       <- k[at_risk]
  weights <- weights[at_risk]

  numer <- tapply(weights * event, k, sum, na.rm = TRUE)
  denom <- tapply(weights,         k, sum, na.rm = TRUE)
  numer / denom
}


#' Cumulative Incidence from Weighted Person-Time (Single-Event)
#'
#' Discrete-time cumulative incidence of `Y` computed from per-row
#' event indicators and Hajek-weighted hazards. Single-event variant
#' (no competing event D), suitable for the CausalSurvival single-
#' outcome setting. Caller is responsible for restricting inputs to a
#' single standing-in arm.
#'
#' Computes weighted hazards by interval via
#' [weighted_hazard_by_k()], aligns to `cut_times`, and applies the
#' standard recursion
#' \deqn{F_Y(t) = \sum_{k \le t} h_Y(k) \, S(k-1)}
#' where \eqn{S(k-1) = \prod_{j < k} (1 - h_Y(j))}.
#'
#' Ported from CausalCompetingRisks (`R/ipw_core.R`), un-deferred
#' 2026-05-12 with the `d_event` arg dropped for the single-event
#' setting.
#'
#' @param y_event Numeric vector of Y event indicators. See
#'   [weighted_hazard_by_k()] for the at-risk encoding.
#' @param k Integer vector of interval indices (1..K_max under spec
#'   §3.0.2). Used as the grouping key for the weighted-hazard fit.
#' @param weights Numeric vector of per-row weights.
#' @param cut_times Numeric vector `c(t_1, ..., T_max)`. Used here for
#'   its length `K_max` (the result is reported at each interval
#'   index `1..K_max`, position-aligned with `cut_times`).
#' @return Numeric vector of cumulative incidence at each interval
#'   index, length `K_max`.
#' @keywords internal
cum_inc_from_weighted <- function(y_event, k, weights, cut_times) {
  K_max      <- length(cut_times)
  haz_y_by_k <- weighted_hazard_by_k(y_event, k, weights)

  # Key by integer interval index 1..K_max (spec §3.0.2). Missing
  # indices (no events / no at-risk rows in that interval for this
  # arm) -> hazard of 0.
  key       <- as.character(seq_len(K_max))
  haz_y_vec <- unname(haz_y_by_k[key])
  haz_y_vec[is.na(haz_y_vec)] <- 0

  # Cumulative survival up to START of each interval (lagged)
  surv <- c(1, cumprod(1 - haz_y_vec)[-length(haz_y_vec)])

  cumsum(haz_y_vec * surv)
}
