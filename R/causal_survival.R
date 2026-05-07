#' Estimate Counterfactual Survival
#'
#' Discrete-time pooled-logistic estimator of the counterfactual survival
#' (or cumulative incidence) under each arm of a binary treatment.
#'
#' @param pt_data A `person_time` object (from [to_person_time()]).
#' @param method One of `"gformula"` or `"ipw"`. Single method per call.
#' @param formulas Optional named list of model formula overrides. Valid
#'   keys: `y` (Y-hazard), `c` (C-hazard / IPCW denominator), `A`
#'   (propensity denominator), `A_num` (propensity stabilization
#'   numerator). Any absent key falls back to the default linear formula.
#' @param truncate NULL or length-2 numeric `c(lower, upper)` percentile
#'   bounds for IPW weight truncation. NULL = no truncation.
#' @param ipcw NULL or logical. NULL → method-conditional default
#'   (`TRUE` under `"ipw"`, `FALSE` under `"gformula"`).
#' @param stabilize One of `"marginal"` or NULL. Joint switch driving
#'   both treatment and censoring numerators (v1 supports marginal only).
#' @param verbose Logical.
#' @param keep_data Logical. When TRUE, store `pt_data` (and the wide
#'   input it was built from) on the fit for later access.
#'
#' @return S3 object of class `"causal_survival_fit"`.
#' @export
causal_survival <- function(pt_data,
                            method     = "gformula",
                            formulas   = NULL,
                            truncate   = NULL,
                            ipcw       = NULL,
                            stabilize  = "marginal",
                            verbose    = FALSE,
                            keep_data  = TRUE) {

  cl <- match.call()

  # --- pt_data class check (stricter than CCR: classed only) ---
  if (!inherits(pt_data, "person_time")) {
    stop(
      "pt_data must inherit 'person_time'. ",
      "Run to_person_time() on subject-level data first.",
      call. = FALSE
    )
  }

  # --- method (single value, not vector) ---
  valid_methods <- c("gformula", "ipw")
  if (length(method) != 1L || !method %in% valid_methods) {
    stop("method must be one of: ",
         paste(shQuote(valid_methods), collapse = ", "),
         ". Got: ", paste(method, collapse = ", "),
         call. = FALSE)
  }

  # --- ipcw: NULL → method-conditional default ---
  if (is.null(ipcw)) {
    ipcw <- (method == "ipw")
  } else if (!is.logical(ipcw) || length(ipcw) != 1L || is.na(ipcw)) {
    stop("ipcw must be NULL or a single TRUE/FALSE.", call. = FALSE)
  }

  # --- stabilize (v1: NULL or "marginal") ---
  if (!is.null(stabilize) && !identical(stabilize, "marginal")) {
    stop(
      "stabilize must be NULL (no stabilization) or \"marginal\" ",
      "(v1 only allows these two).",
      call. = FALSE
    )
  }

  # --- formulas keys ---
  valid_formula_keys <- c("y", "c", "A", "A_num")
  if (!is.null(formulas)) {
    if (!is.list(formulas) || is.null(names(formulas)) ||
        any(names(formulas) == "")) {
      stop("`formulas` must be a named list (keys: ",
           paste(valid_formula_keys, collapse = ", "), ").",
           call. = FALSE)
    }
    bad <- setdiff(names(formulas), valid_formula_keys)
    if (length(bad) > 0L) {
      stop("Unknown formula key(s): ",
           paste(shQuote(bad), collapse = ", "),
           ". Valid keys: ",
           paste(shQuote(valid_formula_keys), collapse = ", "), ".",
           call. = FALSE)
    }
  }

  # --- truncate (length-2 percentiles in [0,1] with lower < upper) ---
  if (!is.null(truncate)) {
    if (!is.numeric(truncate) || length(truncate) != 2L ||
        any(is.na(truncate)) ||
        truncate[1] < 0 || truncate[2] > 1 ||
        truncate[1] >= truncate[2]) {
      stop("truncate must be NULL or c(lower, upper) percentiles in [0, 1] ",
           "with lower < upper.", call. = FALSE)
    }
  }

  # --- Pull metadata off classed pt_data ---
  cut_times      <- attr(pt_data, "cut_times")
  id_col         <- attr(pt_data, "id_col")
  treatment_col  <- attr(pt_data, "treatment_col")
  covariates_vec <- attr(pt_data, "covariates")

  # [worker dispatch — next chunk]
  stop("worker not yet wired", call. = FALSE)
}


# ----------------------------------------------------------------------------
# Internal helper: clone person-time data for counterfactual prediction
# ----------------------------------------------------------------------------

#' Clone Baseline Across Cut Times
#'
#' Broadcasts each subject's baseline (k = 0) row across all `cut_times`,
#' setting the treatment column to a fixed value `a`. Used by the
#' g-formula worker to predict counterfactual hazards at every (subject, k)
#' regardless of the subject's observed event/censoring time.
#'
#' @param baseline data.frame. One row per subject (typically
#'   `pt_data[pt_data$k == 0, ]`).
#' @param cut_times Numeric vector of interval-end times.
#' @param treatment_col Character. Treatment column name.
#' @param a Numeric (0 or 1). Counterfactual treatment value.
#' @return data.frame with `nrow(baseline) * length(cut_times)` rows.
#' @keywords internal
make_clone <- function(baseline, cut_times, treatment_col, a) {
  n <- nrow(baseline)
  K <- length(cut_times)
  clone <- baseline[rep(seq_len(n), each = K), , drop = FALSE]
  clone$k                <- rep(cut_times, times = n)
  clone[[treatment_col]] <- a
  rownames(clone) <- NULL
  clone
}


# ----------------------------------------------------------------------------
# Internal worker: g-formula path
# ----------------------------------------------------------------------------

#' G-Formula Cumulative Incidence Worker
#'
#' Single-event parametric g-formula. Fits an unweighted Y-hazard pooled
#' logistic model, then for each arm `a in {0, 1}`:
#'
#' 1. Clone baseline covariates across all `cut_times`, setting treatment to `a`.
#' 2. Predict the counterfactual hazard at every (subject, k).
#' 3. Compute per-subject survival `S_i(k) = prod_{j <= k}(1 - h(j | a, L_i))`.
#' 4. Marginalize over L by averaging across subjects: `S(k) = mean_i S_i(k)`.
#' 5. CIF = 1 - S.
#'
#' Equivalent unrolled per-arm form (kept here for explanatory reference;
#' the implementation uses an `lapply` over both arms to avoid duplication):
#'
#' ```
#' clone        <- make_clone(baseline, cut_times, treatment_col, a)
#' haz_a1       <- predict_counterfactual_hazard(model_y, clone,
#'                                               treatment_col, 1, "Y-haz a=1")
#' haz_a0       <- predict_counterfactual_hazard(model_y, clone,
#'                                               treatment_col, 0, "Y-haz a=0")
#' S_a1         <- cumprod_survival(haz_a1, clone[[id_col]])
#' S_a0         <- cumprod_survival(haz_a0, clone[[id_col]])
#' surv_a1_by_k <- tapply(S_a1, clone$k, mean)
#' surv_a0_by_k <- tapply(S_a0, clone$k, mean)
#' ```
#'
#' @keywords internal
fit_gformula <- function(pt_data, id_col, treatment_col, covariates_vec,
                         cut_times, formulas) {

  # 1. Fit unweighted Y-hazard
  fit <- fit_hazard_models(
    pt_data        = pt_data,
    treatment      = treatment_col,
    covariates     = covariates_vec,
    active_methods = "g_formula",
    formulas       = formulas,
    ipcw           = FALSE
  )
  model_y <- fit$models$model_y

  # 2. Baseline per subject (k = 0; left-truncation rejected upstream)
  baseline <- pt_data[pt_data$k == 0, , drop = FALSE]

  # 3. Per-arm CIF: clone -> predict -> cumprod -> mean -> 1 - S
  cif_by_arm <- lapply(c(0, 1), function(a) {
    clone <- make_clone(baseline, cut_times, treatment_col, a)
    haz   <- predict_counterfactual_hazard(
      model_y, clone, treatment_col, a,
      paste0("Y-hazard a=", a)
    )
    if (any(is.na(haz))) {
      warning(
        "G-formula: Y-hazard predictions contain ", sum(is.na(haz)),
        " NA value(s). CIF estimates will be biased.",
        call. = FALSE
      )
    }
    S_i <- cumprod_survival(haz, clone[[id_col]])
    S_k <- as.numeric(tapply(S_i, clone$k, mean))
    1 - S_k
  })

  # 4. Long-format estimates
  estimates <- data.frame(
    treatment = rep(c(0, 1), each = length(cut_times)),
    k         = rep(cut_times, times = 2),
    surv      = c(1 - cif_by_arm[[1]], 1 - cif_by_arm[[2]]),
    inc       = c(    cif_by_arm[[1]],     cif_by_arm[[2]])
  )

  list(
    estimates    = estimates,
    models       = list(y = model_y, c = NULL, A = NULL, A_num = NULL),
    model_checks = fit$checks,
    weights      = NULL
  )
}
