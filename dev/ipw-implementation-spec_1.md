# IPW implementation spec — for coding agents

Status: revised 2026-05-12 after agent grilling. Aligns with the actual
`CausalSurvival` package architecture (`causal_survival.R`,
`propensity.R`, `weights.R`, `hazards.R`) and the LOCKED data convention
in `dev/CAUSAL_SURVIVAL_SPEC.md` §3.0.

`causal_survival(method = "ipw", ...)` is the single user-facing IPW
entry. It selects between two survival estimators via a hidden engine
knob:

- `ipw_engine = "km"` (default) — weighted Hajek-pooled-hazard
  per-interval estimator. Nonparametric in `k`. Newly added in this
  spec; un-defers two CCR helpers (see PORT_DEP_GRAPH).
- `ipw_engine = "msm"` (internal) — weighted pooled-logistic MSM.
  Parametric in `k`. Existing implementation; this is the current
  `fit_ipw()` worker, renamed.

Both share the propensity model fitting and weight construction
(already implemented in the package). They diverge only at the
survival-estimation step.

### Naming clarification

"Weighted KM" here means a discrete-time weighted-hazard product-limit:

$$\hat S^a(k) = \prod_{j=1}^{k}(1 - \hat\lambda^a_j),\quad
\hat\lambda^a_j = \frac{\sum_i W_i \mathbb 1\{Y_{ij}=1,A_i=a\}}{\sum_i W_i \mathbb 1\{\text{at risk at } j, A_i=a\}}.$$

With time-varying $W$ this is not the classical Kaplan-Meier product-
limit estimator — it is a weighted pooled-hazard estimator. The "KM"
label is retained because it is nonparametric in $k$ and reduces to
classical KM under unit weights.

---

## Current package state

`causal_survival.R` dispatches on `method`. Current
`valid_methods <- c("gformula", "ipw")`. The current `"ipw"` worker
(`fit_ipw()`) implements the MSM route. **There is no weighted-KM
implementation yet.**

### Refactor plan

1. Keep `valid_methods <- c("gformula", "ipw")`. The KM/MSM distinction
   is internal — the public method label stays `"ipw"`.
2. Add a new arg `ipw_engine = c("km", "msm")` to `causal_survival()`.
   Default `"km"`. Validated only when `method == "ipw"`.
3. Rename current `fit_ipw()` → `fit_ipw_msm()`. No internal logic
   change beyond the data-schema overhaul (see PORT_DEP_GRAPH Phase 3).
4. Add new `fit_ipw_km()` worker. Re-uses propensity + weight steps
   from `fit_ipw_msm()`, diverges at the survival-estimation step.
5. Extract shared steps 1–5 (propensity + weights) into a helper
   `fit_ipw_weights()` called by both engines.

### Constraints on `ipw_engine = "km"`

`ipw_engine = "km"` is only valid when treatment $A$ is a baseline
point-exposure (single A column, no time-varying assignment). Under
time-varying $A$, the KM hazard ratio is no longer well-defined per
arm; `"msm"` must be used. The dispatch should hard-error if KM is
requested under time-varying $A$.

The bootstrap helper (separate spec) should accept either engine
uniformly.

---

## Shared infrastructure (both engines)

Already implemented:

- `fit_propensity()` — `propensity.R` lines 39–70. Fits denominator
  $A \sim \mathbf L$ glm on baseline rows. With `stabilize = TRUE`,
  also fits numerator $A \sim 1$. Returns both models + diagnostic
  checks.
- `ipw_static_trt()` — `weights.R` lines 103–121. Predicts
  $\hat g(A_i = 1 \mid \mathbf L_i)$ on baseline, picks
  $\hat g(A_i \mid \mathbf L_i)$ per subject's observed treatment,
  computes the stabilized weight ratio, broadcasts to person-time rows
  via `id_col`.
- `apply_weight_truncation()` — `weights.R` lines 178–243. Symmetric
  percentile truncation on `pt_data$w_a` and `pt_data$w_cens`. Logs
  flagged ids.
- `summarize_weights()` — `weights.R` lines 258–279. Distributional
  stats on weight columns.

These are invoked by both `fit_ipw_km()` (new) and `fit_ipw_msm()`
(refactor).

### Stabilization

Package stabilization is **structural** (two-cumprod ratio,
Robins/Hernán form):

```r
sw_i = p_obs_num_i / p_obs_full_i
```

where `p_obs_full_i = P(A = A_i | L_i)` from the denominator model
(`A ~ L`) and `p_obs_num_i = P(A = A_i)` from the numerator model
(`A ~ 1`, marginal). For unstabilized weights, `p_obs_num_i = 1`.

`stabilize` argument is one of `"marginal"` (default) or `NULL`.
Conditional stabilization (e.g., `A ~ V` with $V \subseteq \mathbf L$)
is supported by `fit_propensity()` via `formula_num` but not yet
exposed by the public API.

### Data-schema dependencies

Both engines consume person-time data with the LOCKED schema (spec
§3.0):

- `y_event`, `dep_cens`, `indep_cens` — three-way mutually-exclusive
  event flags.
- `k` is the integer interval index `1..K_end`.
- y-hazard, KM hazard rate, and Y-MSM regression all fit on rows
  with `dep_cens == 0 & indep_cens == 0`. Censoring rows (single
  $k = j$ row per censored subject) are dropped; the C-weight
  up-weights the remaining uncensored subjects.
- The c-hazard model (`dep_cens`) fits on `indep_cens == 0` rows.
- Reporting time grid is `cut_times[1..K_max]`. The `K_end = K_max + 1`
  catch-up interval is excluded from reporting.

---

## Engine 1: `ipw_engine = "km"` — weighted Hajek pooled-hazard (NEW)

Apply the package's own weighted-hazard helper to the weighted person-
time data. **Do not use `survival::survfit()`** — build the helper
on top of the ported CCR primitives to match the rest of the
pipeline.

### Target

At each interval $j$ in arm $a$:

$$\hat\lambda^{a}_j = \frac{\sum_i W_i \, \mathbb 1\{Y_{ij}=1,\, A_i=a,\, \text{dep\_cens}_{ij}=0,\, \text{indep\_cens}_{ij}=0\}}{\sum_i W_i \, \mathbb 1\{\text{at risk at } j,\, A_i=a,\, \text{dep\_cens}_{ij}=0,\, \text{indep\_cens}_{ij}=0\}}.$$

The row filter `dep_cens == 0 & indep_cens == 0` is symmetric with
the MSM Y-glm fit population. Censored-at-$j$ rows are dropped (only
the single $k = j$ row), and the C-weight up-weights remaining
uncensored subjects to absorb dep_cens. Subjects' full pre-censoring
history (rows $k < j$) stays in the regression.

With IPCW (when `ipcw = TRUE`), $W_i = SW^A_i \cdot w^C_i(j)$;
otherwise $W_i = SW^A_i$.

Per-arm survival: $\hat S^a_{\text{IPW-KM}}(k) = \prod_{j=1}^k (1 - \hat\lambda^{a}_j)$.

### Helpers ported from CCR

Per `PORT_DEP_GRAPH.md`, two helpers are un-deferred to support this
engine:

- `weighted_hazard_by_k(event, k, weights)` — Hajek ratio per `k`.
  Vector-input, schema-agnostic. Lands in `R/hazards.R` (CCR origin:
  `R/ipw_core.R`).
- `cum_inc_from_weighted(y_event, k, weights, cut_times)` — discrete-
  time cumulative incidence wrapper. **Single-event variant**: drops
  the `d_event` arg present in CCR (CausalSurvival has no competing
  event D). Lands in `R/hazards.R`.

Sum-of-weights in numerator and denominator means the ratio cancels
the weights' overall scale — works for both stabilized and
unstabilized $w$.

### Worker: `fit_ipw_km()`

Mirrors `fit_ipw_msm()` through step 5 (combined weight). Then:

```r
fit_ipw_km <- function(pt_data, id_col, treatment_col, covariates_vec,
                       cut_times, formulas, ipcw, stabilize, truncate) {

  # Steps 1-5 are shared with fit_ipw_msm() and refactored into
  # fit_ipw_weights():
  #   1. fit_propensity()  -> model_a, model_a_num
  #   2. If ipcw: fit C-hazard on indep_cens == 0 rows (denom + num)
  #   3. ipw_static_trt() -> w_a_raw (broadcast to person-time)
  #   4. apply_weight_truncation() -> clipped w_a (and w_cens if ipcw)
  #   5. pt_data$w_combined <- if (ipcw) pt_data$w_a * pt_data$w_cens
  #                            else        pt_data$w_a

  weights_out <- fit_ipw_weights(...)  # returns list w/ pt_data + models
  pt_data <- weights_out$pt_data

  # 6. Weighted Hajek pooled hazard per arm
  surv_by_arm <- lapply(c(0, 1), function(a) {
    arm_rows <- pt_data[
      pt_data[[treatment_col]] == a &
        pt_data$dep_cens   == 0 &
        pt_data$indep_cens == 0,
      , drop = FALSE
    ]
    inc <- cum_inc_from_weighted(
      y_event   = arm_rows$y_event,
      k         = arm_rows$k,
      weights   = arm_rows$w_combined,
      cut_times = cut_times
    )
    1 - inc  # convert F(t) to S(t)
  })

  estimates <- data.frame(
    treatment = rep(c(0, 1), each = length(cut_times)),
    k         = rep(cut_times, times = 2),
    surv      = c(surv_by_arm[[1]], surv_by_arm[[2]]),
    inc       = c(1 - surv_by_arm[[1]], 1 - surv_by_arm[[2]])
  )

  list(
    estimates    = estimates,
    models       = list(y = NULL, c = weights_out$model_c,
                        A = weights_out$model_a,
                        A_num = weights_out$model_a_num,
                        c_num = weights_out$model_c_num),
    model_checks = weights_out$checks,
    weights      = list(
      pt_data_weighted = pt_data,
      weight_summary   = summarize_weights(pt_data),
      truncated_ids    = weights_out$flagged_ids,
      truncate         = truncate
    )
  )
}
```

**Note:** `models$y = NULL` for `ipw_engine = "km"` because there is
no parametric outcome model. The empty slot is preserved for shape
consistency across engines.

---

## Engine 2: `ipw_engine = "msm"` — weighted pooled-logistic MSM (REFACTOR)

Current `fit_ipw()` in `causal_survival.R` lines 312–459. Behavior is
correct; rename + schema-column updates + shared-helper extraction.

### Refactor

1. Rename `fit_ipw()` → `fit_ipw_msm()`.
2. Extract shared steps 1–5 (propensity + weights) into
   `fit_ipw_weights()` that both `fit_ipw_km()` and `fit_ipw_msm()`
   call. Returns a list: `pt_data` (with `w_a`, `w_cens`,
   `w_combined`), `model_a`, `model_a_num`, `model_c`, `model_c_num`,
   `checks`, `flagged_ids`.
3. `fit_ipw_msm()` retains steps 6–7 (weighted Y-MSM fit +
   clone-predict-marginalize).

### Default Y-MSM formula

Per the LOCKED data schema, the outcome column is `y_event`. Default
formula:

```r
y_event ~ A + k + I(k^2) + I(k^3)
```

Polynomial in the integer interval index `k`. Override via
`formulas$y`. The Y-MSM glm fits on rows with `indep_cens == 0` AND
`dep_cens == 0`.

### Survival estimation

After fitting the weighted MSM, per-arm survival is computed by
**clone-predict-marginalize** — same as the g-formula, but the MSM
lacks $\mathbf L$ as a predictor:

```r
baseline <- pt_data[pt_data$k == 1, , drop = FALSE]
clone <- make_clone(baseline, cut_times, treatment_col, a)
haz   <- predict_counterfactual_hazard(model_y, clone, treatment_col,
                                       a, label)
S_i   <- cumprod_survival(haz, clone[[id_col]])
S_k   <- as.numeric(tapply(S_i, clone$k, mean))
```

The `tapply(S_i, clone$k, mean)` step is the marginalization. For a
marginal MSM (no $\mathbf L$ on RHS), this is a constant over subjects
within $k$ — so the marginalization is trivial. The cloning is kept
for shape consistency with `fit_gformula()`.

Baseline is now `k == 1` (was `k == 0`) per schema overhaul.

---

## Top-level dispatch

`causal_survival.R` lines 44–51 + 113–135 — update to:

```r
valid_methods    <- c("gformula", "ipw")
valid_ipw_engine <- c("km", "msm")

# new arg with engine default
causal_survival <- function(pt_data, method = "gformula",
                            ipw_engine = "km", ...) {
  method     <- match.arg(method, valid_methods)
  ipw_engine <- match.arg(ipw_engine, valid_ipw_engine)
  # ...

  if (method == "ipw" && ipw_engine == "km") {
    # hard error if A is time-varying — see Constraints section above
  }

  worker_out <- switch(method,
    gformula = fit_gformula(...),
    ipw      = switch(ipw_engine,
      km  = fit_ipw_km(pt_data, id_col, treatment_col, covariates_vec,
                       cut_times, formulas, ipcw, stabilize, truncate),
      msm = fit_ipw_msm(pt_data, id_col, treatment_col, covariates_vec,
                        cut_times, formulas, ipcw, stabilize, truncate)
    )
  )
}
```

`ipcw` default conditional: `(method == "ipw")`. Update accordingly.

---

## User-facing API

### Default (weighted KM)

```r
ipw_km_fit <- causal_survival(
  pt_data  = pt_data,
  method   = "ipw",
  formulas = list(A = A ~ age + sex)
)
```

### MSM (engine knob)

```r
ipw_msm_fit <- causal_survival(
  pt_data    = pt_data,
  method     = "ipw",
  ipw_engine = "msm",
  formulas   = list(
    A = A ~ age + sex,
    y = y_event ~ A + k + I(k^2) + I(k^3)
  )
)
```

### MSM internals (for reference)

```r
# Steps 1-5 (shared with ipw_km, via fit_ipw_weights()):
#   propensity fit, stabilized weights, IPCW if requested,
#   truncation, combined per-row weight.

# 6. Weighted Y-MSM fit
fit_rows <- pt_data$indep_cens == 0 & pt_data$dep_cens == 0
msm_fit <- fit_logistic(
  formula = y_event ~ A + k + I(k^2) + I(k^3),
  data    = pt_data[fit_rows, ],
  weights = pt_data$w_combined[fit_rows],
  label   = "Y-MSM (IPW)"
)
model_y <- msm_fit$model

# 7. Per-arm survival via clone-predict-marginalize
baseline <- pt_data[pt_data$k == 1, , drop = FALSE]
surv_by_arm <- lapply(c(0, 1), function(a) {
  clone <- make_clone(baseline, cut_times, treatment_col, a)
  haz   <- predict_counterfactual_hazard(
    model_y, clone, treatment_col, a, paste0("Y-MSM a=", a)
  )
  S_i   <- cumprod_survival(haz, clone[[id_col]])
  as.numeric(tapply(S_i, clone$k, mean))
})
```

---

## Open coding questions

- **Bootstrap.** Full refit on each resample, per blog G6 (i).
  Implement `bootstrap_ci()` that resamples subjects, re-runs
  `causal_survival()`, returns percentile CIs. Separate spec.
- **Diagnostics warnings.** Should the package emit warnings when
  `max(sw) > 10` or `n_eff < 0.5 * n`? Currently
  `check_fitted_positivity()` (in `hazards.R`) emits warnings for
  extreme fitted propensities; no analogous check on the final
  weights.
- **Censoring numerator.** Currently hardcoded to `dep_cens ~ A`
  (Hernán Tech Point 12.2). Consider exposing `formulas$c_num` in a
  future version.
- **Conditional stabilization.** `fit_propensity()` supports
  `formula_num` for $V$-stabilization but the public API doesn't
  expose it. `dev/TODO.md` notes this.
- **Time-varying A guard.** Concrete detection rule for hard-erroring
  on `ipw_engine = "km"` under time-varying treatment.

---

## Cross-references

| Asset | Where used |
|---|---|
| `fit_propensity()` (`propensity.R`) | Shared step 1 — both engines |
| `ipw_static_trt()` (`weights.R`) | Shared step 3 — both engines |
| `apply_weight_truncation()` (`weights.R`) | Shared step 4 — both engines |
| `summarize_weights()` (`weights.R`) | Diagnostics output — both engines |
| `weighted_hazard_by_k()` (`hazards.R`, ported from CCR) | `fit_ipw_km()` only |
| `cum_inc_from_weighted()` (`hazards.R`, ported from CCR, single-event variant) | `fit_ipw_km()` only |
| `fit_logistic()` (`hazards.R`) | `fit_ipw_msm()` only — Y-MSM fit |
| `make_clone()` (`causal_survival.R`) | `fit_ipw_msm()` only (and `fit_gformula()`) |
| `predict_counterfactual_hazard()` (`hazards.R`) | `fit_ipw_msm()` only (and `fit_gformula()`) |
| `cumprod_survival()` (`hazards.R`) | `fit_ipw_msm()` only (and `fit_gformula()`) |
