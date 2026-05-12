# CausalSurvival — Package Specification

Version: 0.1.0 (draft)
Status: pre-implementation
Last updated: 2026-05-09 (interval convention switched from `[a, b)` (left-closed right-open) to `(a, b]` (left-open right-closed). Aligns with survival-package counting-process tradition (Therneau, ABGK), Andersen-Ravn product-integral, gfoRmula, ipw, mstate, pammtools. Events at `time = 0` are now a hard error (no home interval under `(0, t_1]`). Earlier on 2026-05-09: `cens_type` arg renamed to `ipcw` and retyped as logical: scalar logical or per-subject logical vector; TRUE = dependent, FALSE = independent. Consolidates the on/off switch and the per-subject labeling into a single arg; clarified as user's a-priori labeling, never inferred from data. Pre-split mode (`event_cols`) deferred to v1.x.)
Earlier: 2026-05-08 (data convention overhaul: § 3.0 added — splits censoring into dep_cens/indep_cens; structural ordering C_admin → C_dep → Y; interval indexing 1..K_end; T_max arg; cens_type scalar/vector arg; hard error if T_max > max(time); warning if mean(admin_reach) < 0.5; cliff structurally impossible. Supersedes the y_flag/c_flag conventions in §3.1, §4.1, §6.2.)
Earlier: 2026-04-30 (decisions: drop competing-event columns from to_person_time; ipcw orthogonal to method; defer conditional stabilization to v1.1; long-format bootstrap replicates; risk() scale arg; treatment required for v0.1.0; keep_data + person_time/data slots; N3d accessor prefix rule; honest reframing of assumption 5; explicit pre/post-migration in §13; bootstrap re-eval call; coupling trap audit § 12.2 with 6 resolutions)

---

## 3.0 Data convention (LOCKED 2026-05-08)

This section is canonical for the input/output shape of `to_person_time()` and the column schema of the `person_time` class. Where §3.1, §4.1, §6.2 below disagree, §3.0 wins.

### 3.0.1 Structural ordering

Within each interval `k`, the structural event ordering is **C_admin → C_dep → Y** (Young/Stensrud/Tchetgen/Hernán 2020 Stat Med, with admin censoring placed first):

- C_admin removes the subject from observation; nothing else can be observed for that subject in that interval.
- C_dep is informative censoring (LTFU, treatment switch, etc.); fires when C_admin didn't.
- Y is the event of interest; fires when neither censoring did.

Boundary handling at `t = T_max` under the `(a, b]` convention (§3.0.2):
- `time = T_max` is **inside** the last interval `(t_{K_max-1}, T_max]`. Events there fire normally at `k = K_max` (`y_event = 1`); censoring there encodes normally (`dep_cens = 1` or `indep_cens = 1`).
- `time > T_max` is outside all intervals → **admin-censored**: subject contributes at-risk rows up to `k = K_max` with no exit row. No subject is ever dropped (would induce survival bias).

### 3.0.2 Interval indexing

Continuous time over `(0, T_max]` is partitioned into `K_max` analyzable intervals:

- `k = 0`: pre-baseline (not in person-time data).
- `k = 1` covers `(0, t_1]`. `k = K_max` covers `(t_{K_max-1}, T_max]`.
- Left-open right-closed intervals: a subject's event time `t` lands in interval `k` iff `t_{k-1} < t <= t_k`. Boundary value `time = t_k` falls in `k` (right-closed).
- This convention matches survival-package `survSplit` (Therneau), Andersen-Borgan-Gill-Keiding counting processes (`(s, t]` for `dN`), Andersen-Ravn product-integral, gfoRmula (`Yt = 1 if event in (t-1, t]`), ipw, mstate, and pammtools. Justification (Therneau): "covariates apply from start of interval; events occur at end" (càdlàg).
- Events at `time = 0` are not supported (hard error in `to_person_time()`) — `(0, t_1]` excludes 0 by construction.
- Estimates are reported at `k = 1, ..., K_max`. Subjects who reach `T_max` at risk are administratively censored by convention — they contribute at-risk rows through `k = K_max` with no exit row.

### 3.0.3 Column schema

The `person_time` class carries three indicator columns instead of `y_flag` / `c_flag`:

```
columns:    id, k, A, <covariates>, y_event, dep_cens, indep_cens
attributes: cut_times, T_max, K_max, treatment_levels,
            id_col, treatment_col, covariates
```

Per row, at most one of `{y_event, dep_cens, indep_cens}` is `1` (exit row); otherwise all three are `0` (at-risk row, including the final row of admin-censored subjects). Mutual exclusivity is invariant.

### 3.0.4 Row encoding (inferred mode)

| input                                                 | indicator        | row at  |
|-------------------------------------------------------|------------------|---------|
| `status = 1, time <= T_max`                           | `y_event = 1`    | `k`     |
| `status = 1, time > T_max`                            | (no exit row — admin-censored) | —       |
| `status = 0, time <= T_max, ipcw[i] = TRUE` (default) | `dep_cens = 1`   | `k`     |
| `status = 0, time <= T_max, ipcw[i] = FALSE`          | `indep_cens = 1` | `k`     |
| `status = 0, time > T_max`                            | (no exit row — admin-censored) | —       |
| at-risk (incl. admin-censored subjects' final row)    | all `0`          | any `k` |

`ipcw` is honored only when `status = 0` AND `time <= T_max`; otherwise ignored. Per-subject vector form is the user's a-priori labeling of each subject's censoring mechanism, never inferred from data. Events at `time = 0` are a hard error (no home interval under `(0, t_1]`).

### 3.0.5 `to_person_time()` arguments

```r
to_person_time(
  data,
  id           = "id",
  time         = "time",
  status       = "status",     # binary {0, 1}
  ipcw         = TRUE,         # scalar logical OR length-nrow(data)
                               # logical vector. TRUE = dependent
                               # censoring (contributes to c-hazard fit
                               # + gets weighted); FALSE = independent
                               # (weight 1, treated as cause-specific
                               # competitor). User's a-priori labeling,
                               # never inferred from data.
  T_max        = NULL,         # default = max(data[[time]])
  treatment    = "A",
  covariates   = character(),
  cut_points   = NULL,         # NULL → 12 equi-spaced cuts over (0, T_max]
  time_varying = NULL          # reserved
)
```

`T_max` validation:
- `T_max = NULL` → defaults to `max(data[[time]])`.
- `T_max > max(data[[time]])` → **hard error** (would create empty trailing intervals).
- `mean(admin_reach) < 0.5` after encoding → **warning** ("hazard at K_max from thin risk set; CIF at K_max unreliable"). Fit proceeds.
- `event` at `time = 0` → **hard error** listing affected subject id(s) (no home interval under `(0, t_1]`).

### 3.0.6 Hazard fit populations

| model    | rows included                              | response   |
|----------|--------------------------------------------|------------|
| c-hazard | `indep_cens = 0`                           | `dep_cens` |
| y-hazard | `indep_cens = 0` AND `dep_cens = 0`        | `y_event`  |

Implications:
- `y_event = 1` rows have `dep_cens = 0` by mutual exclusivity → included in the c-hazard fit with response 0.
- Admin-truncated subjects (no exit row) contribute at-risk rows up to `k = K_max` — included in both fits with response 0 at every interval.
- Mid-stream `indep_cens = 1` rows (user-supplied via `ipcw[i] = FALSE`) → excluded from both fits.

### 3.0.7 IPCW weight

```
W_C(k) = 1 / cumprod(1 - h_C(j))    for j = 1..k
```

Cumprod through current `k`, using the c-hazard predictions on `dep_cens`. Matches gfoRmula / ipw / ltmle. **No cliff guard needed** — the encoding makes the cliff structurally impossible because `indep_cens = 1` rows never enter the c-hazard fit.

### 3.0.8 What is *not* in the data

- Rows materialized at `k > K_max`.
- `indep_cens = 1` indicators auto-generated mid-stream by the package (only via explicit `ipcw[i] = FALSE` per-subject override).
- A separate `c_flag` column. Replaced by `dep_cens` / `indep_cens` split.
- Synthetic exit rows for admin-truncated subjects. They simply have no exit row.

### 3.0.9 Pre-split mode

Dropped from v1. Pre-split mode (`event_cols` arg accepting subject-level `y_event` / `dep_cens` / `indep_cens` indicator columns directly) was deferred to v1.x to keep the v1 surface small. Users with already-classified censoring must run their classification through `status` + `ipcw` (per-subject vector). Re-add in v1.x if real workflows need it.

Users **cannot** supply already-discretized person-time data; the function always discretizes from subject-level input.

---

## 1. Overview

R package for causal inference on single-event survival outcomes, using discrete-time pooled logistic regression. Two estimators: parametric g-formula (plug-in) and inverse probability weighting. Bootstrap CIs.

Designed as the foundation of a two-package ecosystem. `CausalCompetingRisks` (CCR) imports CausalSurvival's primitives and adds the separable-effects framework for competing events.

### 1.1 Scope (v1)

In:
- Single failure event Y plus right-censoring C
- Binary or multi-level treatment A, baseline-only
- Continuous or categorical baseline covariates L₀
- Discrete time grid (user-supplied or auto-binned)
- G-formula (plug-in, no MC) and IPW (Rep 1 style)
- Cluster bootstrap on subject id
- Identifying assumptions accessor

Out:
- Competing events (→ CCR)
- Time-varying covariates / treatment (v2)
- ICE, MC g-formula (v2)
- Continuous-time / Cox-based estimators
- BCa / studentized bootstrap (v2)
- Targeted learning / doubly-robust estimators

### 1.2 Audience

Epidemiologists and applied statisticians familiar with causal inference. Reading level: Hernán & Robins *Causal Inference: What If*. Not a tidymodels-style ML pkg.

### 1.3 Design philosophy

**Simplicity and auditability above all.** The code should be:

- **Easy to read** — anyone fluent in base R should be able to follow a function end-to-end without consulting external documentation
- **Easy to explain** — each function does one thing; names match what they do; no clever idioms
- **Easy to audit** — a reviewer should be able to trace the math from Hernán & Robins to a specific block of code in one pass

Concrete consequences:

- Pure base R (no dplyr, no data.table) — fewer abstractions between code and computation
- Hardcoded simple defaults (linear formulas, no trim, single method) — user opts into complexity
- Small public surface (~10 exports) — fewer entry points, less to verify
- Explicit error messages — every `stop()` names the offending column / value
- Comments explain the *why* (which equation, which assumption) not the *what*
- No magic — no NSE, no metaprogramming, no global state mutation

Trade-offs accepted:

- More verbose than tidyverse equivalents (e.g., `ave()` calls in g-formula computation)
- User responsible for choices the package could otherwise make implicit (column names, event coding, formula spec)
- Performance not optimized beyond "fast enough for typical epi sample sizes"

Speed, terseness, and feature completeness are subordinate to readability. When in doubt: write more lines.

---

## 2. Package metadata

```
Package:  CausalSurvival
Version:  0.1.0
License:  MIT
Depends:  R (>= 4.1)
Imports:  stats, survival, ggplot2, patchwork, future.apply
Suggests: testthat, knitr, rmarkdown
Author:   Alex Mourer [aut, cre]
```

No dplyr. No data.table. Pure base R for compute. Heavy stat pkgs (`survival`, `gfoRmula`, `WeightIt`) follow the same convention — keep dep tree shallow.

---

## 3. Public API

10 user-facing exports plus auto-exported S3 methods.

### 3.1 Data prep

**Superseded by §3.0.5 as of 2026-05-08.** The signature, column schema, and contract live in §3.0. Highlights:

- `event` / `event_y` / `event_c` args replaced by `status` (binary) + optional `ipcw` (scalar logical or per-subject logical vector).
- Output columns are `y_event` / `dep_cens` / `indep_cens`, not `y_flag` / `c_flag`.
- `T_max` arg added, defaults to `max(time)`.
- `n_intervals` arg dropped; `cut_points = NULL` defaults to 12 equi-spaced cuts.

CCR users: CCR's own fit functions add `A_y <- A; A_d <- A` and `d_event` derivation inline at fit-time. CausalSurvival's `to_person_time` knows nothing about competing events.

### 3.2 Main estimator

```r
causal_survival(
  pt_data,                          # must inherit "person_time"; error otherwise
  method     = "gformula",          # or "ipw"
  formulas   = NULL,                # list(y, c, A, A_num); NULL → linear defaults
  truncate   = NULL,                # NULL or c(lower, upper) percentile bounds
  ipcw       = NULL,                # NULL → method-conditional default: TRUE under "ipw", FALSE under "gformula"
  stabilize  = "marginal",          # "marginal" or NULL
  verbose    = FALSE,
  keep_data  = TRUE                 # store person-time + wide input in fit
)
```

Returns: S3 `c("causal_survival_fit")`. Single method per call (not vector).

### 3.3 Bootstrap

```r
bootstrap(
  fit,
  n_boot   = 500L,
  alpha    = 0.05,
  seed     = NULL,
  parallel = FALSE,
  workers  = NULL,                  # NULL → future::availableCores() - 1 if parallel
  verbose  = FALSE
)
```

Returns: S3 `c("causal_survival_bootstrap")`. Cluster bootstrap on subject id. Internally uses `future.apply::future_lapply` with `future.seed = seed` for reproducible L'Ecuyer streams. Per-replicate seeds discarded. Failed replicates skipped with effective-B reported.

Percentile CIs only in v1.

### 3.4 Accessors

Naming rule (N3d): accessors that **read a quantity off a fit** carry the `causal_` prefix. Helpers/verbs (`bootstrap`, `to_person_time`, `get_*`) stay plain. The prefix resolves real namespace collisions (`emmeans::contrast`, `rms::contrast`) and signals "this projects out of a fit".

```r
causal_risk(fit, scale = "incidence", ci = NULL)   # scale ∈ c("incidence", "survival")
causal_contrast(fit, reference = NULL, contrasts = NULL, ci = NULL, time = NULL)
causal_assumptions(fit)
causal_diagnostic(fit)
causal_risk_table(fit, count = TRUE)
get_person_time(fit)              # returns fit$person_time, errors if keep_data was FALSE
get_data(fit)                     # returns wide validated input, errors if keep_data was FALSE
```

Direct `$` access also supported: `fit$person_time`, `fit$data`. Accessors are stable contract — internal storage may change in future versions, accessors won't.

`causal_risk()`: returns S3 `c("causal_survival_risk")` with `$risk` long-format `data.frame(method, treatment, k, value, lower, upper)`. CIs from `ci` arg (a `causal_survival_bootstrap` object) or NULL.

`causal_contrast()`:
- `reference = NULL, contrasts = NULL` → all pairwise
- `reference = "<level>"` → all vs. reference
- `contrasts = list(name = list(arms = c(...), op = c("-", "/")))` → custom
- `time = NULL` → max cut time
- `time` beyond max → error
- `time` off-grid in range → snap to nearest, message

`causal_assumptions()`: returns S3 `c("causal_survival_assumptions")` with hardcoded baseline list (§4.6).

`causal_diagnostic()`: model checks + weight summary slot.

`causal_risk_table()`: counts at-risk / events / censored per arm at cut_times.

### 3.5 S3 methods (auto-exported)

```
print.causal_survival_fit
print.causal_survival_risk
print.causal_survival_contrast
print.causal_survival_bootstrap
print.causal_survival_assumptions
print.causal_survival_diagnostic
summary.causal_survival_fit
plot.causal_survival_risk
plot.causal_survival_contrast    # placeholder, message-only in v1
plot.causal_survival_diagnostic  # placeholder, message-only in v1
format.causal_survival_assumptions
```

---

## 4. S3 contracts

### 4.1 `person_time`

**Superseded by §3.0.3 as of 2026-05-08.** Class shape:

```
class:      c("person_time", "data.frame")
columns:    id_col, k, A, <covariates>, y_event, dep_cens, indep_cens
attributes: cut_times, T_max, K_max, treatment_levels,
            id_col, treatment_col, covariates
```

Indicator values: `0/1`. Mutual exclusivity invariant — at most one of `{y_event, dep_cens, indep_cens}` is `1` per row (exit row), or all three are `0` (at-risk row). No NAs.

`k` is integer interval index, range `1..K_max`.

CCR users: CCR's fit functions add `A_y`, `A_d`, and `d_event` columns inline before model-fitting. CausalSurvival's `person_time` class itself carries no competing-event columns.

### 4.2 `causal_survival_fit`

```r
list(
  call             = match.call(),
  method           = "gformula" | "ipw",
  estimates        = data.frame(treatment, k, surv, inc),
  models           = list(y = <glm>, c = <glm or NULL>, A = <glm or NULL>, A_num = <glm or NULL>),
  arms             = data.frame(name, label, color),
  weights          = list(W_C = <vector>, W_A = <vector>) | NULL,
  model_checks     = list(),
  warnings         = character(),
  messages         = character(),
  cut_times        = numeric(),
  n                = integer(),
  formulas         = list(),
  args             = list(method, truncate, ipcw, stabilize, keep_data),
  person_time      = <data.frame> | NULL,    # NULL if keep_data = FALSE
  data             = <data.frame> | NULL     # wide validated input; NULL if keep_data = FALSE
)
```

Survival and incidence stored on both scales. Accessor toggles output.

### 4.3 `causal_survival_bootstrap`

```r
list(
  fit_call         = <original fit's call>,
  n_boot_requested = integer(),
  n_boot_effective = integer(),
  alpha            = numeric(),
  replicates       = data.frame(boot_id, treatment, k, value),  # long-format
  ci_lower         = data.frame(treatment, k, lower),
  ci_upper         = data.frame(treatment, k, upper),
  failed_reps      = integer(),
  warnings_count   = integer()
)
```

Raw replicates kept. Per-rep seeds not stored.

### 4.4 `causal_survival_risk`

```r
list(
  risk = data.frame(method, treatment, k, value, lower, upper),
  scale = "survival" | "incidence"
)
```

Long format. CCR's `separable_effects_risk` extends with `arm`, `a_y`, `a_d` columns.

### 4.5 `causal_survival_contrast`

```r
list(
  contrasts = data.frame(name, treatment_a, treatment_b, op, k, estimate, lower, upper),
  time      = numeric(),
  alpha     = numeric()
)
```

### 4.6 `causal_survival_assumptions`

Each assumption: `list(name, statement, status, pointer)`. Four fields, no nesting.

Hardcoded baseline:
1. Consistency — `Y^a = Y` when `A = a`. Untestable.
2. Exchangeability — `Y^a ⊥ A | L`. Untestable.
3. Positivity — `P(A = a | L) > 0` for all l. Checkable → `causal_diagnostic()$weight_summary`.
4. No interference — one subject's `A` doesn't affect another's `Y`. Untestable.
5. Correct model specification. Untestable in general; diagnostic warnings available via `causal_diagnostic()$model_checks` (necessary-not-sufficient checks).
6. Censoring at random (E2) — `T^a ⊥ C | A, L`. Untestable.

CCR appends Δ1, Δ2 isolation conditions and modified-treatment exchangeability.

---

## 5. Internal architecture

```
R/
├── data_prep.R       to_person_time, [helpers for survSplit + flag derivation]
├── validate.R        validate_input_shape, validate_subject_level,
│                     validate_person_time, check_covariate_quality
├── hazards.R         fit_hazard_models, fit_logistic, check_fitted_positivity,
│                     predict_hazard_under, predict_with_warning, cumprod_survival
├── propensity.R      fit_propensity
├── weights.R         ipw, ipw_cens, ipw_static_trt, apply_weight_truncation,
│                     summarize_weights
├── gformula_core.R   gformula_estimate, make_clone, predict_hazards, compute_cum_inc
├── ipw_core.R        weighted_hazard_by_k, cum_inc_from_weighted
├── bootstrap.R       bootstrap (generic on fit_fn callback), bootstrap-internal
├── accessors.R       causal_risk, causal_contrast, causal_diagnostic, snap_time, build_risk_long
├── assumptions.R     causal_assumptions, baseline_assumption_list
├── risk_table.R      causal_risk_table, risk_table_internal
├── print.R           print/summary/format S3 methods
├── plot.R            plot.causal_survival_risk, build_risk_table_plot
├── causal_survival.R causal_survival, fit_causal_survival (worker)
└── zzz.R             package globals, options, .onLoad
```

User-facing exports: ~10 functions. Internals: ~25-30. CCR accesses internals via `@importFrom`.

---

## 6. Defaults

| Argument | Default | Note |
|---|---|---|
| `method` | `"gformula"` | Single string, not vector |
| `formulas` | `NULL` (linear defaults) | `y_flag ~ A + k + cov1 + cov2 + ...` |
| `truncate` | `NULL` (no trim) | Hernán recommends no default trim |
| `ipcw` | `NULL` | Method-conditional default: `TRUE` under `"ipw"`, `FALSE` under `"gformula"`. Orthogonal to `method` — fits W_C regardless of how W_A is handled. |
| `stabilize` | `"marginal"` | Or `NULL`. `"conditional"` deferred to v1.1 (see §15). |
| `verbose` | `FALSE` | |
| `n_intervals` | `12L` | Used if `cut_points = NULL` |
| `n_boot` | `500L` | |
| `alpha` | `0.05` | |
| `seed` | `NULL` | |
| `parallel` | `FALSE` | |
| `workers` | `NULL` | `availableCores() - 1` if parallel |

### 6.1 Estimation modes

`causal_survival()` requires `treatment` in v0.1.0. Two valid modes depending on whether `covariates` are supplied:

| `treatment` | `covariates` | Estimand | Notes |
|---|---|---|---|
| named | empty | Crude survival per arm | No L adjustment |
| named | named | Adjusted survival per arm + contrasts | Standard causal use case |
| `NULL` | any | — | Error in v0.1.0; standardized marginal survival deferred to v1.1 (see §15) |

### 6.2 Default formula construction

When `formulas = NULL`, package builds linear formulas based on what's supplied:

```
# treatment + covariates
y_event  ~ A + k + cov1 + cov2 + ...
dep_cens ~ A + k + cov1 + cov2 + ...      # if ipcw = TRUE
A        ~ cov1 + cov2 + ...              # if method = "ipw"
A_num    ~ 1                              # marginal stabilization

# treatment only (no covariates)
y_event  ~ A + k
dep_cens ~ A + k                          # if ipcw = TRUE
A        ~ 1                              # if method = "ipw" (intercept-only — flags positivity edge case)
```

User opts into polynomials, splines, or interactions via `formulas` arg:

```r
formulas = list(
  y = y_event  ~ A * (k + I(k^2) + I(k^3)) + age + sex,
  c = dep_cens ~ A + k + age,
  A = A ~ age + sex,
  A_num = A ~ 1
)
```

`formulas` is a named list with keys `y`, `c`, `A`, `A_num`. Missing keys fall back to defaults. Extra keys error.

The c-hazard model fits its response (`dep_cens`) on rows with `indep_cens = 0` (per §3.0.6). The y-hazard model fits on rows with `indep_cens = 0` AND `dep_cens = 0`. Mid-stream `indep_cens = 1` rows are excluded from both fits by construction.

---

## 7. Error / warning / silent

### 7.1 Hard errors (`stop()`)

- Missing required column in `data` or `pt_data`
- NA in `id`, `time`, `treatment`, `covariates`
- Non-discrete `time` when discrete required
- `status` column has values outside `{0, 1}`
- `ipcw` not logical, or length not in `{1, nrow(data)}`, or contains `NA`
- Event(s) at `time = 0` (no home interval under `(0, t_1]` convention)
- `T_max > max(data[[time]])` (would create empty trailing intervals)
- Invalid `truncate` range (not length-2, lower ≥ upper, outside [0, 1])
- `pt_data` not of class `person_time`
- `time` arg in `contrast()` beyond `max(cut_times)` or below `min(cut_times)`
- `method` not in `c("gformula", "ipw")`
- `method = "gformula"` AND `ipcw` is a length-`n` vector — per-subject IPCW labels are meaningless when no IPCW step runs.

### 7.2 Warnings (`warning()`)

- NA hazards predicted (model rank-deficient or extrapolation)
- Extreme weights when `truncate = NULL`
- Both `n_intervals` and `cut_points` supplied (cut_points takes precedence)
- Positivity issue (very small fitted hazards or weights > some threshold)
- `mean(admin_reach) < 0.5` after encoding ("hazard at K_max from thin risk set; CIF at K_max unreliable")
- Censoring rate > 50% (combined `dep_cens + indep_cens` rate)
- Zero events in an arm
- Bootstrap replicates failed (printed in summary, captured in `$failed_reps`)
- `method = "gformula"` — always fires: "method='gformula' does not perform IPCW; the `ipcw` arg is inert under this method." Notifies the user that no weighting step occurs, regardless of `ipcw` value.

### 7.3 Silent (collected in `fit$messages`)

- Snap to nearest cut time when `time` is off-grid but in range
- Default formula construction (note which spec was used)

### 7.4 Warning propagation

Top-level `causal_survival()` wraps the fit in `withCallingHandlers`. Captured warnings stored in `fit$warnings`. Re-emitted at console so user notices.

`bootstrap()` suppresses per-replicate warnings, reports count in `$warnings_count`.

`print(fit)` displays "Fit completed with N warnings (see fit$warnings)" when N > 0.

---

## 8. Censoring strategy

| Method | Default | Mechanism |
|---|---|---|
| `gformula` | conditioning under E2 | Hazard models fit on uncensored risk set; covariates must include all predictors of censoring (Hernán & Robins TP 21.10) |
| `gformula` + `ipcw = TRUE` | optional IPCW augmentation | Adds explicit `W_C` for double robustness |
| `ipw` | always with `ipcw = TRUE` | Explicit `W_C` weights, censoring model fitted |

Default `ipcw = TRUE` when `method = "ipw"`. Default `ipcw = FALSE` for `method = "gformula"` (user opts in for double robustness).

---

## 9. Test ladder

10 generic steps (from the user's ladder), all in CausalSurvival's `tests/testthat/`. CCR's tests start at step 11+.

```
1.  Hazards fit correctly (single-arm, no censoring)
2.  Survival from hazards = KM with no censoring
3.  Survival = KM with independent censoring
4.  Cumulative incidence with independent censoring (= 1 - S in single-event)
5.  Two-group: steps 3 + 4 with treatment
6.  Contrasts: ATE, RD, RR
7.  Bootstrap CIs cover truth at nominal rate
8.  Confounding: IPW vs g-formula match (independent censoring)
9.  Confounding: IPW vs g-formula match (dependent censoring)
10. Confounding: IPW + IPCW vs g-formula + IPCW match (dependent censoring)
```

Each step runs as:
- Single-line correctness test (one dataset, point estimate within tolerance)
- `n_rep` simulated test (statistical agreement across replicates)

Test simulator: `tests/testthat/helper-simulate.R` → `simulate_survival_data()`. Internal, not exported.

Reference checks:
- Step 2/3 vs `survival::survfit()` to numerical tolerance
- One published-result reproduction (NHEFS or similar) — TBD

---

## 10. Vignette

Single vignette for v1: `vignettes/getting-started.Rmd`.

Dataset: `survival::lung` (NCCTG lung cancer). Single event (death), ~30% censored, multiple baseline covariates (age, sex, ECOG performance, weight loss).

Workflow demonstrated:
1. Load data
2. `to_person_time()` — convert + bin
3. `causal_survival(method = "gformula")`
4. `causal_risk()` + `plot()`
5. `causal_survival(method = "ipw")` — comparison
6. `bootstrap()` + `causal_contrast()`
7. `causal_assumptions()` — display identifying conditions

README: 1-paragraph description, install, 5-10 line example. Links to vignette and alexmourer.com.

---

## 11. Class hierarchy with CCR

```
CausalSurvival classes (parents):
  person_time
  causal_survival_fit
  causal_survival_risk
  causal_survival_bootstrap
  causal_survival_contrast
  causal_survival_diagnostic
  causal_survival_assumptions

CausalCompetingRisks classes (children):
  separable_effects             extends causal_survival_fit
  separable_effects_risk        extends causal_survival_risk
  separable_effects_bootstrap   extends causal_survival_bootstrap
  separable_effects_contrast    extends causal_survival_contrast
  separable_effects_diagnostic  extends causal_survival_diagnostic
  separable_effects_assumptions extends causal_survival_assumptions
```

Inheritance via `class()` vector, e.g. `class(x) <- c("separable_effects_risk", "causal_survival_risk", "data.frame")`.

S3 dispatch tries CCR's specific method first; falls back to CausalSurvival's parent. CCR's `plot.separable_effects_risk` calls `NextMethod()` for base draw, then layers SE-specific overlays.

---

## 12. Migration plan (CCR → CausalSurvival + CCR)

### 12.1 Sequence

1. Create CausalSurvival skeleton (DESCRIPTION, NAMESPACE, R/, tests/) as sibling repo to CCR
2. **Trap surgery in CCR first** (per §12.2 audit): refactor the 4 functions that need to be ported with a different shape *before* moving anything. CCR's tests must still pass after this step.
3. Move cleanly-generic functions from CCR → CausalSurvival (the ~17 of 23 that have no coupling traps)
4. Build CausalSurvival's forked versions (per §12.2): write fresh, simpler versions of the 5 forked traps in CausalSurvival's R/
5. Update CCR `NAMESPACE`: `Imports: CausalSurvival (>= 0.1.0, < 0.2.0)` (strict pre-1.0)
6. Add `@importFrom CausalSurvival ...` to CCR functions that call moved primitives
7. Update internal calls in CCR's `R/` to reference CausalSurvival functions
8. Audit roxygen `[fn()]` cross-refs — replace with `[CausalSurvival::fn()]` where they cross the boundary
9. Move generic tests to CausalSurvival/tests/; keep SE-specific tests in CCR/tests/
10. `R CMD check` both packages, fix breakage
11. Bump versions: CausalSurvival 0.1.0, CCR 0.2.0 (breaking change due to defaults)

Pinning: strict `(>= 0.1.0, < 0.2.0)` pre-1.0. Floor only post-1.0.

### 12.2 Coupling trap resolutions

Audit (2026-04-30) found 7 coupling traps where generic-looking CCR code silently assumes competing-event structure. Resolutions:

| # | File | Function | Severity | Resolution | Notes |
|---|---|---|---|---|---|
| 1 | `hazards.R` | `fit_hazard_models` | HARD | **Fork** | CausalSurvival fits Y + C; CCR keeps Y + D + C. Both call shared `fit_logistic()` worker. ~10 lines of C-hazard orchestration duplicated. |
| 2 | `validate.R` | `validate_subject_level` | HARD | **Fork** | CausalSurvival validates 2-level event, k-level treatment, `event_y`/`event_c` distinct. CCR keeps current 3-level event / 2-level treatment / `event_y`/`event_d`/`event_c` distinct. Package-specific error messages preserved. |
| 3 | `data_prep.R` | `to_person_time` | HARD | **Wrapper** | CausalSurvival's `to_person_time` does the generic part (Y-flag, C-flag, person-time expansion, 2-slot `event_labels`). CCR ships `to_person_time_ccr()` (or shadow name) that calls CausalSurvival's then derives `d_flag`, adds `A_y`/`A_d`, relabels `event_labels` to 3 slots. |
| 4 | `gformula_core.R` | `gformula_estimate` (+ helpers) | HARD | **Generalize helpers + fork orchestrator** | Helpers (`make_clone`, `predict_hazards`, `compute_cum_inc`) take named-vector treatment assignments (`c(A = 1)` for CausalSurvival, `c(A_y = 1, A_d = 0)` for CCR) and a list of hazard models to predict over. Orchestrators differ only in which assignment vectors they enumerate. CausalSurvival's enumerates `unique(pt_data$A)` (handles k-level treatment per §1.1); CCR's enumerates the 4 fixed `(A_y, A_d)` arms. |
| 5 | `ipw_core.R` | `cum_inc_from_weighted` | SOFT | **Fork** | CausalSurvival ships single-event version (~10 lines) without the `(1 - haz_d)` factor. CCR keeps current. |
| 6 | `risk_table.R` | `risk_table_internal` | SOFT | **Fork** | CausalSurvival has 3 `count` branches (`at_risk`, `events_y`, `censored`). CCR keeps 4 (adds `events_d`). Token vocabulary stays part of API. |
| 7 | `bootstrap.R` | `bootstrap` | HARD | **Fork** | CCR's `bootstrap` builds a 4D array indexed by (k × method × arm × boot_id) — tied to the SE four-arm structure. CausalSurvival's `bootstrap` builds long-format `data.frame(boot_id, treatment, k, value)` per §4.3. Different storage shapes; fresh implementation in CausalSurvival. |

Pattern: when CCR-specific structure is woven through (different methods, different validation, different math), fork. When it's a layered addition (data prep), wrap. When generality is forced anyway by multi-level treatment (helpers in #4), generalize.

Step 2 of §12.1 must complete before step 3: traps #3 and #4 require surgery in CCR (wrapper extraction; helper generalization) so that what's left after the surgery is portable.

### 12.3 Function move list

Inventory of all functions in CCR's `R/` (audit 2026-04-30). Classification:

- **CLEAN-MOVE** — fully generic, ports to CausalSurvival as-is
- **SHARED-HELPER** — small utility used by both; lives in CausalSurvival, CCR imports
- **TRAP-FORK** / **TRAP-WRAPPER** / **TRAP-GENERALIZE** — see §12.2 for trap #
- **CCR-ONLY** — stays in CCR (separable-effects-specific)

| File | Function | Class | Trap# | Note |
|---|---|---|---|---|
| `ipw_core.R` | `weighted_hazard_by_k` | CLEAN-MOVE | — | Hajek hazard from weighted person-time |
| `ipw_core.R` | `cum_inc_from_weighted` | TRAP-FORK | 5 | Single-event version in CausalSurvival |
| `risk_table.R` | `risk_table` | TRAP-FORK | 6 | Becomes `causal_risk_table` (N3d) in CausalSurvival; fresh impl |
| `risk_table.R` | `risk_table_internal` | TRAP-FORK | 6 | Branch on event tokens differs |
| `separable_swap.R` | `swap_d_weights` | CCR-ONLY | — | SE Rep 2 swap weights |
| `separable_swap.R` | `swap_y_weights` | CCR-ONLY | — | SE Rep 2 swap weights |
| `gformula_core.R` | `gformula_estimate` | TRAP-GENERALIZE | 4 | Orchestrator forks; helpers shared |
| `gformula_core.R` | `make_clone` | TRAP-GENERALIZE | 4 | Generalized to named-vector assignments |
| `gformula_core.R` | `predict_hazards` | TRAP-GENERALIZE | 4 | Generalized to model-list |
| `gformula_core.R` | `predict_with_warning` | CLEAN-MOVE | — | Warning-capture wrapper around predict() |
| `gformula_core.R` | `compute_cum_inc` | TRAP-GENERALIZE | 4 | Generalized to hazard list |
| `propensity.R` | `fit_propensity` | SHARED-HELPER | — | Point-treatment propensity |
| `utils.R` | `%\|\|%` | SHARED-HELPER | — | Null-coalescing operator |
| `separable_arms.R` | `separable_arm_hazards` | CCR-ONLY | — | SE four-arm hazard prediction |
| `weights.R` | `ipw` | CLEAN-MOVE | — | Hajek IPW from cumprod ratio |
| `weights.R` | `ipw_cens` | CLEAN-MOVE | — | IPCW wrapper |
| `weights.R` | `ipw_static_trt` | CLEAN-MOVE | — | Baseline-treatment IPTW |
| `weights.R` | `ipw_time_varying_trt` | CCR-ONLY | — | v2 placeholder; stays in CCR until v2 |
| `weights.R` | `apply_weight_truncation` | CLEAN-MOVE | — | Symmetric percentile truncation |
| `weights.R` | `summarize_weights` | CLEAN-MOVE | — | Distribution summary |
| `data_prep.R` | `to_person_time` | TRAP-WRAPPER | 3 | CausalSurvival owns generic; CCR wraps |
| `separable_effects.R` | `separable_effects` | CCR-ONLY | — | SE user-facing orchestrator |
| `separable_effects.R` | `fit_separable_effects` | CCR-ONLY | — | SE worker (Rep 1, Rep 2) |
| `contrasts.R` | `compute_contrasts` | CCR-ONLY | — | SE decomposition (A vs B); CausalSurvival writes fresh `causal_contrast` |
| `separable_ipw.R` | `ipw_estimate` | CCR-ONLY | — | SE IPW orchestrator |
| `separable_ipw.R` | `weighted_arm_cum_inc` | CCR-ONLY | — | SE per-arm CIF |
| `separable_ipw.R` | `estimate_weighted_cum_inc` | CCR-ONLY | — | SE per-rep CIF |
| `arms.R` | `arm_spec` | CCR-ONLY | — | SE four-arm specification |
| `arms.R` | `arm_names` | CCR-ONLY | — | SE four-arm names |
| `print.R` | `print.separable_effects` | CCR-ONLY | — | S3 on CCR class |
| `print.R` | `summary.separable_effects` | CCR-ONLY | — | S3 on CCR class |
| `print.R` | `confint.separable_effects` | CCR-ONLY | — | See §14 re. confint dispatch |
| `print.R` | `print.separable_effects_risk` | CCR-ONLY | — | S3 on CCR class |
| `print.R` | `print.separable_effects_contrast` | CCR-ONLY | — | S3 on CCR class |
| `print.R` | `print.separable_effects_diagnostic` | CCR-ONLY | — | S3 on CCR class |
| `hazards.R` | `fit_hazard_models` | TRAP-FORK | 1 | Y+C in CausalSurvival, Y+D+C in CCR |
| `hazards.R` | `fit_logistic` | SHARED-HELPER | — | Generic GLM fitter w/ warning capture |
| `hazards.R` | `check_fitted_positivity` | CLEAN-MOVE | — | GLM diagnostics |
| `hazards.R` | `predict_hazard_under` | CLEAN-MOVE | — | Counterfactual hazard prediction |
| `hazards.R` | `cumprod_survival` | CLEAN-MOVE | — | Per-subject cumulative survival |
| `assumptions.R` | `assumptions` | CCR-ONLY | — | CausalSurvival writes fresh `causal_assumptions` |
| `assumptions.R` | `isolation_summary` | CCR-ONLY | — | SE-specific |
| `assumptions.R` | `print.separable_effects_assumptions` | CCR-ONLY | — | S3 on CCR class |
| `assumptions.R` | `format.separable_effects_assumptions` | CCR-ONLY | — | S3 on CCR class |
| `accessors.R` | `risk` | CCR-ONLY | — | CausalSurvival writes fresh `causal_risk` (N3d) |
| `accessors.R` | `build_risk_long` | TRAP-FORK | 6 | Pivot logic differs by arm structure |
| `accessors.R` | `contrast` | CCR-ONLY | — | CausalSurvival writes fresh `causal_contrast` (N3d) |
| `accessors.R` | `snap_time` | CLEAN-MOVE | — | Numeric snap-to-grid |
| `accessors.R` | `diagnostic` | CCR-ONLY | — | CausalSurvival writes fresh `causal_diagnostic` (N3d) |
| `bootstrap.R` | `bootstrap` | TRAP-FORK | 7 | Long-format replicates in CausalSurvival |
| `bootstrap.R` | `print.separable_effects_bootstrap` | CCR-ONLY | — | S3 on CCR class |
| `plot.R` | `plot.separable_effects_risk` | CCR-ONLY | — | S3 on CCR class |
| `plot.R` | `build_risk_table_plot` | TRAP-FORK | 6 | CausalSurvival writes simpler version |
| `plot.R` | `build_contrast_annotations` | CCR-ONLY | — | SE decomposition labels |
| `plot.R` | `plot.separable_effects_contrast` | CCR-ONLY | — | S3 on CCR class (placeholder) |
| `plot.R` | `plot.separable_effects_diagnostic` | CCR-ONLY | — | S3 on CCR class (placeholder) |
| `validate.R` | `validate_input_shape` | CLEAN-MOVE | — | data.frame shape checks |
| `validate.R` | `validate_subject_level` | TRAP-FORK | 2 | Different event/treatment-level requirements |
| `validate.R` | `validate_person_time` | CLEAN-MOVE | — | Person-time format checks |
| `validate.R` | `check_covariate_quality` | CLEAN-MOVE | — | NA / cardinality checks |

**Counts (~60 functions total):**
- CLEAN-MOVE: 14
- SHARED-HELPER: 3
- TRAP-FORK: 9 (functions across traps #1, #2, #5, #6, #7)
- TRAP-WRAPPER: 1 (trap #3)
- TRAP-GENERALIZE: 4 (trap #4 helpers + orchestrator)
- CCR-ONLY: ~30

**What actually moves to CausalSurvival's R/ in step 3 of §12.1:** 14 CLEAN-MOVE + 3 SHARED-HELPER = **17 functions ported as-is**. Trap-affected functions (15) are *not* moved — CausalSurvival writes fresh versions per §12.2.

---

## 13. Differences from CCR (defaults that diverge)

"Pre-migration" = CCR's current behavior on `master` as of 2026-04-30. "Post-migration" = the behavior CCR will adopt as part of bumping to 0.2.0 alongside CausalSurvival 0.1.0.

| Default | CausalSurvival 0.1.0 | CCR pre-migration | CCR post-migration |
|---|---|---|---|
| `method` | `"gformula"` (single) | `"all"` (both) | `"all"` (kept; SE workflow always compares g-formula vs IPW) |
| Default formulas | linear (`A + k + cov`) | linear with CCR-specific terms | linear (`A_y + A_d + k + cov`) — matches CausalSurvival shape, SE-specific terms |
| `truncate` | `NULL` (no trim) | `c(0.01, 0.99)` | `NULL` (aligned with CausalSurvival; Hernán recommends no default trim) |
| `time` beyond max in `causal_contrast()` | error | snap to max with message | error (aligned with CausalSurvival; users opt into snapping explicitly) |

Three of four divergences resolve at migration by flipping CCR's defaults to match CausalSurvival. One stays divergent: `method` default. Single-method default is friendlier for the typical CausalSurvival user (one fit, one method); CCR's `"all"` default reflects SE workflow (always compare g-formula vs IPW for sanity).

---

## 14. Open items (TBDs)

- **Coupling trap resolution.** Resolved 2026-04-30; see §12.2.
- **`weighted_hazard_by_k()` orphan check.** Probably moot post-audit: trap #5 (`cum_inc_from_weighted`) is forked, and both forks call `weighted_hazard_by_k`. Confirm during migration; if truly orphan, drop from both packages.
- **Bootstrap re-fit mechanism.** Locked: re-eval `fit$call` with subset `pt_data` (re-runs validation per replicate). Acceptable for v0.1.0; revisit if profiling shows validation cost dominates.
- **Published-result reproduction test.** Pick dataset (NHEFS most likely) and target estimate. Build into test ladder. Not a v0.1.0 release blocker; track as v0.1.0 acceptance criterion.
- **Plot methods for contrast and diagnostic.** Ship as placeholders (message-only) in v0.1.0. Implement in v0.2.
- **`confint.causal_survival_fit` deprecation.** §3.5 lists it as "deprecated, redirect to bootstrap". For a v0.1.0 release with no predecessor, this is odd phrasing — either drop the method entirely or document as "intentionally not provided; use `bootstrap()` for CIs". Decide before v0.1.0 ships.

---

## 15. Roadmap beyond v1

v1.1 candidates:
- BCa bootstrap CIs
- Conditional stabilization weights with V-conditioning
- Subgroup analysis API
- Standardized marginal survival without treatment contrast (closest analogue: `stdReg`)

v2:
- Time-varying covariates and treatment
- ICE g-formula
- MC g-formula (forward simulation)
- Doubly-robust / TMLE estimators

Out of scope indefinitely:
- Continuous-time / Cox-based estimators (use `causalCmprsk`, `concrete`, `mets` instead)
- Mediation analysis (separate package planned)
