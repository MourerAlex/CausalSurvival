#' Propensity Model: P(A | L)
#'
#' Fits and predicts the treatment propensity used by the IPW estimator.
#' The propensity weight `1 / pi(A | L)` extrapolates from the observed
#' treatment mixture to the counterfactual full-population estimand
#' (Hernán & Robins, *Causal Inference: What If*, ch. 12).
#'
#' @keywords internal
NULL


#' Fit Propensity Score Models pi(A | L)
#'
#' Fits a logistic propensity-score model on baseline rows for IPW. With
#' `stabilize = TRUE` (default), also fits a marginal numerator model
#' (`A ~ 1`) for stabilized weights (Robins / Hernán form).
#'
#' Baseline rows only (`k == 0`): the package scope is point treatment
#' (assigned at `k = 0`, frozen thereafter).
#'
#' @param pt_data Person-time data frame.
#' @param treatment Character. Treatment column name.
#' @param covariates Character vector. Baseline covariate column names.
#' @param stabilize Logical. If TRUE (default) also fit the marginal
#'   numerator model A ~ 1.
#' @param formula_full Optional formula override for the conditional
#'   denominator model (default: `A ~ L1 + L2 + ...`).
#' @param formula_num Optional formula override for the numerator model.
#'   Default `A ~ 1` (marginal stabilization). Conditional stabilization
#'   (e.g. `A ~ V_baseline`) is supported by the helper but not yet
#'   exposed by the public API — see `dev/TODO.md`.
#'
#' @return List with elements:
#'   - `model_a` — fitted full propensity glm
#'   - `model_a_num` — fitted numerator glm, or NULL if `stabilize = FALSE`
#'   - `check_a` — diagnostics from [check_fitted_positivity()]
#'   - `check_a_num` — diagnostics for the numerator model, or NULL
#' @keywords internal
fit_propensity <- function(pt_data, treatment, covariates,
                            stabilize = TRUE,
                            formula_full = NULL, formula_num = NULL) {
  baseline <- pt_data[pt_data$k == 0, ]

  # Full conditioning: A ~ L
  if (is.null(formula_full)) {
    rhs <- if (length(covariates) > 0) {
      paste(covariates, collapse = " + ")
    } else {
      "1"
    }
    formula_full <- stats::as.formula(paste(treatment, "~", rhs))
  }
  full <- fit_logistic(formula_full, baseline, "Propensity (full)")

  # Numerator for stabilization: marginal by default, optionally conditional
  num <- NULL
  if (stabilize) {
    if (is.null(formula_num)) {
      formula_num <- stats::as.formula(paste(treatment, "~ 1"))
    }
    num <- fit_logistic(formula_num, baseline, "Propensity (numerator)")
  }

  list(
    model_a     = full$model,
    model_a_num = if (is.null(num)) NULL else num$model,
    check_a     = full$check,
    check_a_num = if (is.null(num)) NULL else num$check
  )
}
