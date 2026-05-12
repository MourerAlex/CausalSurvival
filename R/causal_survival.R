#' Counterfactual risk under a binary point treatment
#'
#' Estimate the counterfactual risk
#' \eqn{E[Y_k^{a, c = 0}] = \Pr(Y_k^{a, c = 0} = 1)} that the event of
#' interest has occurred by the end of discrete interval `k` under
#' arm `a` of a binary treatment assigned at baseline (`t = 0`) and
#' held fixed thereafter, with the treatment-dependent censoring
#' component held at zero. Equivalently, the cumulative incidence
#' \eqn{F^a(t_k) = \Pr(T^a \le t_k, C^{d,a} = 0)}; the survival
#' \eqn{S^a(t_k) = 1 - F^a(t_k)} follows. Time is discretized into a
#' finite grid of intervals `k = 1, …, K_max` under the half-open
#' convention `(t_{k-1}, t_k]`.
#'
#' Identification rests on consistency, exchangeability conditional
#' on the recorded baseline covariates `L_0`, and positivity in both
#' the treatment and the censoring mechanisms. Censoring is split
#' into a treatment-dependent component `dep_cens` (the component
#' `c` set to zero in the estimand) and an administrative or
#' otherwise independent component `indep_cens` (handled by the
#' at-risk set, not by intervention).
#'
#' Two estimators are exposed:
#'
#' - `method = "gformula"`: discrete-time parametric g-formula. Fit
#'   the discrete hazards of `Y` (and, when needed, of `C`) on the
#'   pooled person-time data, then simulate the counterfactual risk
#'   `E[Y_k^{a, c = 0}]` in each arm by recursive substitution.
#' - `method = "ipw"`: inverse-probability weighting. Build the
#'   product weight (treatment and, when applicable, censoring) and
#'   estimate the counterfactual hazard in each arm by a weighted
#'   Hajek pooled-hazard (default, nonparametric in `k`) or by a
#'   weighted pooled-logistic marginal structural model
#'   (`.ipw_engine = "msm"`); the risk is recovered from the
#'   hazards.
#'
#' Standard errors are not produced by this function; obtain them by
#' bootstrap on the subject-level input.
#'
#' @param pt_data A `person_time` object returned by
#'   [to_person_time()].
#' @param method Estimator. One of `"gformula"` or `"ipw"`.
#' @param formulas Optional named list of model formula overrides.
#'   Keys: `y` (Y-hazard), `c` (C-hazard / IPCW denominator), `A`
#'   (propensity denominator), `A_num` (propensity numerator). Any
#'   absent key falls back to the default linear formula.
#' @param truncate `NULL` or length-2 numeric `c(lower, upper)`
#'   percentile bounds for IPW weight truncation. `NULL` leaves the
#'   weights untruncated.
#' @param ipcw `NULL` or logical. `NULL` resolves to the
#'   method-conditional default (`TRUE` under `"ipw"`, `FALSE` under
#'   `"gformula"`).
#' @param stabilize One of `"marginal"` or `NULL`. Drives both the
#'   treatment and the censoring numerators; v1 supports marginal
#'   stabilization only.
#' @param verbose Logical.
#' @param keep_data Logical. When `TRUE`, the `pt_data` (and the
#'   subject-level input it was built from) are retained on the
#'   returned fit.
#' @param .ipw_engine Internal. Pooled-hazard engine under
#'   `method = "ipw"`: `"km"` (default, weighted Hajek hazard,
#'   nonparametric in `k`) or `"msm"` (weighted pooled logistic with
#'   the cubic-in-`k` default). The leading dot flags this as a
#'   developer-facing knob.
#'
#' @return An S3 object of class `"causal_survival_fit"`.
#' @export
causal_survival <- function(pt_data,
                            method      = "gformula",
                            formulas    = NULL,
                            truncate    = NULL,
                            ipcw        = NULL,
                            stabilize   = "marginal",
                            verbose     = FALSE,
                            keep_data   = TRUE,
                            .ipw_engine = "km") {

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

  # --- .ipw_engine (only relevant when method == "ipw") ---
  valid_ipw_engines <- c("km", "msm")
  if (length(.ipw_engine) != 1L || !.ipw_engine %in% valid_ipw_engines) {
    stop(".ipw_engine must be one of: ",
         paste(shQuote(valid_ipw_engines), collapse = ", "),
         ". Got: ", paste(.ipw_engine, collapse = ", "),
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

  if (verbose) message("causal_survival(): fitting method = '", method, "'")

  # --- Worker dispatch with warning collection ---
  # Inner fitters re-emit glm warnings via warning(); the outer handler
  # captures them into `collected_warnings`, muffles propagation, and a
  # single grouped notice is fired at the end so the caller knows to
  # inspect fit$warnings.
  collected_warnings <- character()

  ipw_worker <- switch(.ipw_engine,
    km  = fit_ipw_km,
    msm = fit_ipw_msm
  )
  worker_out <- withCallingHandlers(
    switch(method,
      gformula = fit_gformula(
        pt_data        = pt_data,
        id_col         = id_col,
        treatment_col  = treatment_col,
        covariates_vec = covariates_vec,
        cut_times      = cut_times,
        formulas       = formulas
      ),
      ipw = ipw_worker(
        pt_data        = pt_data,
        id_col         = id_col,
        treatment_col  = treatment_col,
        covariates_vec = covariates_vec,
        cut_times      = cut_times,
        formulas       = formulas,
        ipcw           = ipcw,
        stabilize      = stabilize,
        truncate       = truncate
      )
    ),
    warning = function(w) {
      collected_warnings <<- c(collected_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (length(collected_warnings) > 0L) {
    warning(
      "causal_survival(): ", length(collected_warnings),
      " warning(s) collected during fit. See fit$warnings for details.",
      call. = FALSE
    )
  }

  # --- Assemble S3 fit object ---
  ci_list <- list(gformula = NULL, ipw = NULL)
  ci_list[[method]] <- worker_out$estimates

  fit <- list(
    call                 = cl,
    method               = method,
    ipw_engine           = if (method == "ipw") ipw_engine else NULL,
    cumulative_incidence = ci_list,
    weights              = worker_out$weights,
    models               = worker_out$models,
    model_checks         = worker_out$model_checks,
    model_diagnostics    = NULL,
    warnings             = collected_warnings,
    pt_data              = if (keep_data) pt_data else NULL,
    cut_times            = cut_times,
    treatment_levels     = attr(pt_data, "treatment_levels"),
    id_col               = id_col,
    treatment_col        = treatment_col,
    covariates           = covariates_vec,
    stabilize            = stabilize,
    ipcw                 = ipcw,
    truncate             = truncate
  )
  class(fit) <- "causal_survival_fit"
  fit
}


# ----------------------------------------------------------------------------
# Internal helper: clone person-time data for counterfactual prediction
# ----------------------------------------------------------------------------

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

  # 2. Baseline per subject (k = 1, first analyzable interval under
  # the LOCKED (0, t_1] convention; left-truncation rejected upstream)
  baseline <- pt_data[pt_data$k == 1, , drop = FALSE]

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


# ----------------------------------------------------------------------------
# Internal helper: shared IPW weight construction (steps 1-5)
# ----------------------------------------------------------------------------

#' Build IPW Weights
#'
#' Compute the inverse-probability weight for each person-time row,
#' \deqn{W_i(k) \;=\; \frac{1}{\pi(A_i \mid L_i)} \,
#'                   \prod_{j=1}^{k}\frac{1}{1 - h_C(j;\, A_i, L_i)},}
#' replacing each factor by its stabilized ratio when
#' `stabilize = "marginal"`:
#'
#' | `stabilize`  | propensity numerator | censoring numerator       |
#' |--------------|----------------------|---------------------------|
#' | `"marginal"` | `A ~ 1`              | `dep_cens ~ A` (when ipcw)|
#' | NULL         | none (unstabilized)  | none (unstabilized)       |
#'
#' Truncate at the requested percentile and attach `w_a`, `w_cens`,
#' `w_combined` (and their raw counterparts) to `pt_data`. Shared by
#' [fit_ipw_msm()] and [fit_ipw_km()]; the engines differ only in the
#' downstream survival-estimation step. The censoring numerator is
#' fixed to `dep_cens ~ A` (H&R Technical Point 12.2).
#'
#' @return List with
#'   - `pt_data` — augmented with weight columns.
#'   - `models`  — list of fitted glms (`A`, `A_num`, `c`, `c_num`).
#'   - `checks`  — list of per-model diagnostics.
#'   - `flagged_ids` — ids of subjects with weights flagged by
#'     truncation (or `NULL` when no truncation requested).
#' @keywords internal
fit_ipw_weights <- function(pt_data, id_col, treatment_col,
                            covariates_vec, formulas, ipcw, stabilize,
                            truncate) {

  do_stabilize <- identical(stabilize, "marginal")

  # ---------- 1. Propensity model(s) on the k = 1 row ----------
  prop_fit <- fit_propensity(
    pt_data       = pt_data,
    treatment     = treatment_col,
    covariates    = covariates_vec,
    stabilize     = do_stabilize,
    formula_full  = formulas$A,
    formula_num   = formulas$A_num
  )
  model_a     <- prop_fit$model_a
  model_a_num <- prop_fit$model_a_num

  # ---------- 2. Censoring model(s) for IPCW ----------
  model_c     <- NULL
  model_c_num <- NULL
  check_c     <- NULL
  check_c_num <- NULL
  if (ipcw) {
    haz_fit <- fit_hazard_models(
      pt_data        = pt_data,
      treatment      = treatment_col,
      covariates     = covariates_vec,
      active_methods = "ipw",
      formulas       = formulas,
      ipcw           = TRUE
    )
    model_c <- haz_fit$models$model_c
    check_c <- haz_fit$checks$c

    if (do_stabilize) {
      cnum_fml  <- stats::as.formula(
        paste("dep_cens ~", treatment_col)
      )
      cnum_rows <- pt_data$indep_cens == 0
      cnum_fit  <- fit_logistic(
        cnum_fml, pt_data[cnum_rows, , drop = FALSE],
        "C-hazard (numerator)"
      )
      model_c_num <- cnum_fit$model
      check_c_num <- cnum_fit$check
    }
  }

  # ---------- 3. Build raw weights ----------
  w_a_raw <- ipw_static_trt(
    model_full    = model_a,
    pt_data       = pt_data,
    treatment_col = treatment_col,
    id_col        = id_col,
    model_num     = model_a_num
  )
  w_cens_raw <- if (ipcw) {
    ipw_cens(model_c, pt_data, id_col, model_num = model_c_num)
  } else {
    NULL
  }

  # `*_raw` columns survive truncation so reweight() can re-apply
  # truncation without refitting upstream models.
  pt_data$w_a_raw <- w_a_raw
  pt_data$w_a     <- w_a_raw
  if (ipcw) {
    pt_data$w_cens_raw <- w_cens_raw
    pt_data$w_cens     <- w_cens_raw
  }

  # ---------- 4. Truncation ----------
  trunc_out <- apply_weight_truncation(
    pt_data  = pt_data,
    id_col   = id_col,
    truncate = truncate
  )
  pt_data <- trunc_out$pt_data

  # ---------- 5. Combined per-row weight ----------
  pt_data$w_combined <- if (ipcw) {
    pt_data$w_a * pt_data$w_cens
  } else {
    pt_data$w_a
  }

  list(
    pt_data     = pt_data,
    models      = list(A = model_a, A_num = model_a_num,
                       c = model_c, c_num = model_c_num),
    checks      = list(A = prop_fit$check_a,
                       A_num = prop_fit$check_a_num,
                       c = check_c, c_num = check_c_num),
    flagged_ids = trunc_out$flagged_ids
  )
}


# ----------------------------------------------------------------------------
# Internal worker: IPW path — MSM engine (weighted pooled logistic)
# ----------------------------------------------------------------------------

#' IPW Cumulative Incidence Worker (MSM engine)
#'
#' Build IPW weights via [fit_ipw_weights()], then fit a weighted
#' pooled-logistic Y-MSM and marginalize per arm by
#' clone-predict-marginalize. The MSM is parametric in `k`
#' (`y_event ~ A + k + I(k^2) + I(k^3)` by default).
#'
#' @keywords internal
fit_ipw_msm <- function(pt_data, id_col, treatment_col, covariates_vec,
                        cut_times, formulas, ipcw, stabilize, truncate) {

  w_out   <- fit_ipw_weights(
    pt_data        = pt_data,
    id_col         = id_col,
    treatment_col  = treatment_col,
    covariates_vec = covariates_vec,
    formulas       = formulas,
    ipcw           = ipcw,
    stabilize      = stabilize,
    truncate       = truncate
  )
  pt_data <- w_out$pt_data

  # ---------- 6. Weighted Y-MSM fit ----------
  # Default Y-MSM is marginal in covariates: weights handle confounding,
  # so no covariate adjustment in the outcome model. Users can supply
  # `formulas$y` for a covariate-conditional MSM.
  # Fit population: rows with indep_cens == 0 & dep_cens == 0 (spec §3.0.6).
  time_terms <- "k + I(k^2) + I(k^3)"
  fml_y <- formulas$y %||% stats::as.formula(
    paste("y_event ~",
          paste(c(treatment_col, time_terms), collapse = " + "))
  )
  y_rows  <- pt_data$indep_cens == 0 & pt_data$dep_cens == 0
  msm_fit <- fit_logistic(
    formula = fml_y,
    data    = pt_data[y_rows, , drop = FALSE],
    label   = "Y-MSM (IPW)",
    weights = pt_data$w_combined[y_rows]
  )
  model_y <- msm_fit$model
  check_y <- msm_fit$check

  # ---------- 7. Per-arm CIF: clone -> predict -> cumprod -> mean ----------
  baseline   <- pt_data[pt_data$k == 1, , drop = FALSE]
  cif_by_arm <- lapply(c(0, 1), function(a) {
    clone <- make_clone(baseline, cut_times, treatment_col, a)
    haz   <- predict_counterfactual_hazard(
      model_y, clone, treatment_col, a,
      paste0("Y-MSM a=", a)
    )
    if (any(is.na(haz))) {
      warning(
        "IPW (MSM): Y-MSM predictions contain ", sum(is.na(haz)),
        " NA value(s). CIF estimates will be biased.",
        call. = FALSE
      )
    }
    S_i <- cumprod_survival(haz, clone[[id_col]])
    S_k <- as.numeric(tapply(S_i, clone$k, mean))
    1 - S_k
  })

  estimates <- data.frame(
    treatment = rep(c(0, 1), each = length(cut_times)),
    k         = rep(cut_times, times = 2),
    surv      = c(1 - cif_by_arm[[1]], 1 - cif_by_arm[[2]]),
    inc       = c(    cif_by_arm[[1]],     cif_by_arm[[2]])
  )

  list(
    estimates    = estimates,
    models       = list(
      y     = model_y,
      c     = w_out$models$c,
      A     = w_out$models$A,
      A_num = w_out$models$A_num,
      c_num = w_out$models$c_num
    ),
    model_checks = c(list(y = check_y), w_out$checks),
    weights      = list(
      pt_data_weighted = pt_data,
      weight_summary   = summarize_weights(pt_data),
      truncated_ids    = w_out$flagged_ids,
      truncate         = truncate
    )
  )
}


# ----------------------------------------------------------------------------
# Internal worker: IPW path — KM engine (weighted Hajek pooled hazard)
# ----------------------------------------------------------------------------

#' IPW Cumulative Incidence — Weighted KM Engine
#'
#' Estimate the counterfactual cumulative incidence under each arm by
#' a weighted Hajek pooled-hazard estimator:
#' \deqn{\hat\lambda^a_k \;=\;
#'   \frac{\sum_i W_i\, 1\{Y_{ik}=1,\, A_i = a\}}
#'        {\sum_i W_i\, 1\{\text{at risk at } k,\, A_i = a\}},
#'   \qquad
#'   \hat F^a(k) \;=\; 1 - \prod_{j=1}^{k}(1 - \hat\lambda^a_j).}
#' Weights `W_i` come from [fit_ipw_weights()]. The hazard is
#' nonparametric in `k`; no outcome model is fit.
#'
#' Requires baseline treatment. The arm-specific risk set is undefined
#' under time-varying `A`; use [fit_ipw_msm()] in that case. v1
#' restricts the package to baseline `A`.
#'
#' @keywords internal
fit_ipw_km <- function(pt_data, id_col, treatment_col, covariates_vec,
                       cut_times, formulas, ipcw, stabilize, truncate) {

  w_out   <- fit_ipw_weights(
    pt_data        = pt_data,
    id_col         = id_col,
    treatment_col  = treatment_col,
    covariates_vec = covariates_vec,
    formulas       = formulas,
    ipcw           = ipcw,
    stabilize      = stabilize,
    truncate       = truncate
  )
  pt_data <- w_out$pt_data

  # ---------- 6. Weighted Hajek pooled hazard per arm ----------
  # Fit population (per spec §3.0.6 / §1 KM symmetry with the Y-MSM
  # fit): rows with indep_cens == 0 & dep_cens == 0.
  cif_by_arm <- lapply(c(0, 1), function(a) {
    arm_rows <- pt_data[
      pt_data[[treatment_col]] == a &
        pt_data$dep_cens   == 0 &
        pt_data$indep_cens == 0,
      , drop = FALSE
    ]
    cum_inc_from_weighted(
      y_event   = arm_rows$y_event,
      k         = arm_rows$k,
      weights   = arm_rows$w_combined,
      cut_times = cut_times
    )
  })

  estimates <- data.frame(
    treatment = rep(c(0, 1), each = length(cut_times)),
    k         = rep(cut_times, times = 2),
    surv      = c(1 - cif_by_arm[[1]], 1 - cif_by_arm[[2]]),
    inc       = c(    cif_by_arm[[1]],     cif_by_arm[[2]])
  )

  list(
    estimates    = estimates,
    models       = list(
      y     = NULL,                       # no Y outcome model under KM
      c     = w_out$models$c,
      A     = w_out$models$A,
      A_num = w_out$models$A_num,
      c_num = w_out$models$c_num
    ),
    model_checks = c(list(y = NULL), w_out$checks),
    weights      = list(
      pt_data_weighted = pt_data,
      weight_summary   = summarize_weights(pt_data),
      truncated_ids    = w_out$flagged_ids,
      truncate         = truncate
    )
  )
}
