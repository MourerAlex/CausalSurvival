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

  if (verbose) message("causal_survival(): fitting method = '", method, "'")

  # --- Worker dispatch with warning collection ---
  # Inner fitters re-emit glm warnings via warning(); the outer handler
  # captures them into `collected_warnings`, muffles propagation, and a
  # single grouped notice is fired at the end so the caller knows to
  # inspect fit$warnings.
  collected_warnings <- character()

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
      ipw = fit_ipw(
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


# ----------------------------------------------------------------------------
# Internal worker: IPW path
# ----------------------------------------------------------------------------

#' IPW Cumulative Incidence Worker
#'
#' Inverse-probability-weighted estimator. Builds per-row IPT weights for
#' the static point treatment, optional IPCW for right censoring, applies
#' percentile truncation, fits a weighted Y-MSM, and marginalizes per arm.
#'
#' Joint stabilization (per `stabilize`):
#'
#' | `stabilize`  | propensity numerator | censoring numerator       |
#' |--------------|----------------------|---------------------------|
#' | `"marginal"` | `A ~ 1`              | `c_flag ~ A` (when ipcw)  |
#' | NULL         | none (unstabilized)  | none (unstabilized)       |
#'
#' Censoring numerator is fixed to `c_flag ~ A` per H&R Technical Point
#' 12.2; v1 doesn't expose a `c_num` formula slot.
#'
#' @keywords internal
fit_ipw <- function(pt_data, id_col, treatment_col, covariates_vec,
                    cut_times, formulas, ipcw, stabilize, truncate) {

  do_stabilize <- identical(stabilize, "marginal")

  # ---------- 1. Propensity model(s) on baseline ----------
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
      cnum_fml <- stats::as.formula(
        paste("c_flag ~", treatment_col)
      )
      cnum_fit <- fit_logistic(cnum_fml, pt_data, "C-hazard (numerator)")
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

  # Stash raw + working columns. apply_weight_truncation() clips w_a /
  # w_cens; the *_raw columns are preserved so reweight() can re-apply
  # truncation without refitting the upstream models.
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

  # ---------- 5. Combined per-row weight for the Y-MSM ----------
  combined_w <- if (ipcw) pt_data$w_a * pt_data$w_cens else pt_data$w_a

  # ---------- 6. Weighted Y-MSM fit ----------
  # Default Y-MSM is marginal in covariates: weights handle confounding,
  # so no covariate adjustment in the outcome model. Users can supply
  # `formulas$y` for a covariate-conditional MSM (e.g. for subgroup
  # estimands in v1.1).
  time_terms <- "k + I(k^2) + I(k^3)"
  fml_y <- formulas$y %||% stats::as.formula(
    paste("y_flag ~", paste(c(treatment_col, time_terms), collapse = " + "))
  )
  msm_fit <- fit_logistic(
    formula = fml_y,
    data    = pt_data,
    label   = "Y-MSM (IPW)",
    weights = combined_w
  )
  model_y <- msm_fit$model
  check_y <- msm_fit$check

  # ---------- 7. Per-arm CIF: clone -> predict -> cumprod -> mean ----------
  baseline <- pt_data[pt_data$k == 0, , drop = FALSE]
  cif_by_arm <- lapply(c(0, 1), function(a) {
    clone <- make_clone(baseline, cut_times, treatment_col, a)
    haz   <- predict_counterfactual_hazard(
      model_y, clone, treatment_col, a,
      paste0("Y-MSM a=", a)
    )
    if (any(is.na(haz))) {
      warning(
        "IPW: Y-MSM predictions contain ", sum(is.na(haz)),
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

  stop("fit_ipw: WIP — return list not yet wired", call. = FALSE)
}
