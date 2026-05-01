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
#' Used by [fit_hazard_models()] (Y, C hazards on person-time) and by
#' [fit_propensity()] (treatment model on baseline rows).
#'
#' @param formula A fitted model formula.
#' @param data Data frame passed to `glm()`.
#' @param label Human-readable label for warnings and diagnostics.
#'
#' @return A list with `model` (the glm object) and `check` (diagnostics).
#' @keywords internal
fit_logistic <- function(formula, data, label) {
  glm_warnings <- character()

  model <- withCallingHandlers(
    stats::glm(formula, data = data,
               family = stats::binomial(link = "logit")),
    warning = function(w) {
      glm_warnings <<- c(glm_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

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
