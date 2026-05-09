# Port Dependency Graph: CCR → CausalSurvival

Live artifact tracking the dependency-ordered port of 17 functions from
CausalCompetingRisks (CCR) to CausalSurvival, per spec §12.3
(`../../separable_effects/dev/CAUSAL_SURVIVAL_SPEC.md`).

Updated as functions are ported. Status: **Phase 3 in progress**.
Phase 1 + Phase 2 complete — all port-set functions either ported or
deferred (snap_time, weighted_hazard_by_k, cum_inc_from_weighted); Phase 2
trap-fresh-writes landed (`to_person_time`, `validate_subject_level`,
`standardize_treatment`, `fit_hazard_models`).

Phase 3 progress:
- `causal_survival()` skeleton: signature + arg validation written; worker
  dispatch + S3 fit object assembly pending.
- G-formula path: `make_clone` (TRAP-GENERALIZE #4, inlined fresh-write —
  single-event simplification of CCR's two-arm clone) + `fit_gformula`
  worker written. CIF math validated against gfoRmula's `simulate.R`
  (per-subject `poprisk_i = cumsum(h * S(t-1))` ↔ ours `1 - mean_i S_i(t)`,
  algebraically identical for single event).
- IPW path: `fit_ipw()` landed. Smoke test (commit `084fe98`) revealed
  IPCW weight cliff at K_max from the v0.1.0 admin-censoring encoding.
- Accessors / bootstrap / S3 methods: not yet started.

---

## ⚠️ Phase 3 BLOCKED on data-convention overhaul (2026-05-08)

The smoke test exposed a structural problem with the `y_flag` / `c_flag`
encoding: under the v0.1.0 contract, administratively censored subjects
get `c_flag = 1` at the last interval, which collides with the c-hazard
fit (`h_C(K_max) ≈ 1` → IPCW weight `1 / (1 - h_C) ≈ 1e11`).

After methodology review, the data convention has been overhauled. The
canonical spec now lives in `~/Bureau/cowork/separable_effects/dev/CAUSAL_SURVIVAL_SPEC.md` § 3.0 (LOCKED 2026-05-08). Highlights:

- Column schema changes: `y_flag` / `c_flag` → `y_event` / `dep_cens` / `indep_cens`.
- Structural ordering C_admin → C_dep → Y; admin censoring placed first within each interval.
- `k` becomes integer interval index `1..K_max` (was: left-edge interval time).
- Subjects who reach `T_max` at risk are administratively censored by convention — no exit row materialized; their final at-risk row sits at `k = K_max` with all flags `0`.
- `ipcw` arg (scalar logical or per-subject logical vector; TRUE = dependent censoring [contributes to c-hazard fit + gets weighted], FALSE = independent [weight 1, treated as cause-specific competitor]); `T_max` arg with default `max(time)`. Per-subject vector is the user's a-priori labeling of each subject's censoring mechanism, never inferred from data.
- C_admin wins ties at `T_max` (events at `time >= T_max` silently dropped).
- Hard error if `T_max > max(time)`; warning if `mean(admin_reach) < 0.5`.
- Cliff structurally impossible: `indep_cens = 1` rows never enter the c-hazard fit (model fits on `indep_cens = 0` rows only).

### Phase 3 refactor inventory (per-file changes)

**`R/data_prep.R`**
- `to_person_time()`: full rewrite per spec §3.0.5. New args: `ipcw`, `T_max`. Output columns `y_event` / `dep_cens` / `indep_cens`. Admin-truncation handled by dropping events/censoring at `time >= T_max` (no exit row). Intervals follow `(a, b]` convention (matches survSplit native). Hard error on events at `time = 0` (no home interval under `(0, t_1]`).
- `validate_subject_level()`: new checks for `ipcw` shape (scalar logical or length-`n` logical vector); `status ∈ {0, 1}`; `T_max ≤ max(time)` (hard error); `mean(admin_reach) < 0.5` (warning).
- `validate_person_time()`: column set `{y_event, dep_cens, indep_cens}`; mutual-exclusivity invariant (at most one indicator per row); `k ∈ {1, ..., K_max}`.

**`R/hazards.R`**
- `fit_hazard_models()`: c-hazard becomes `glm(dep_cens ~ ...)` on rows with `indep_cens = 0`; y-hazard becomes `glm(y_event ~ ...)` on rows with `indep_cens = 0` AND `dep_cens = 0`. Polynomial-in-`k` interpretation shifts from time-value to integer index; formula form unchanged.

**`R/weights.R`**
- `ipw_cens()`: column references `c_flag` → `dep_cens`. No cliff guard needed — replaced by structural exclusion (admin-truncated subjects have no exit row at `T_max`; `indep_cens = 1` rows excluded from c-hazard fit).
- `ipw()`, `ipw_static_trt()`: column renames only; logic stable.

**`R/causal_survival.R` workers (gformula, ipw)**
- `make_clone()`: clone `k` values become integer `1..K_max` (predictions only reported up to T_max).
- G-formula worker: predicts `y_event` hazard on the clone for `k = 1..K_max`. Cumprod survival over those K_max steps.
- IPW worker: weighted Y-MSM glm fits on rows with `indep_cens = 0` AND `dep_cens = 0` only. Combined weight `w_a * w_cens` as before.
- Baseline-covariate extraction: `pt_data[pt_data$k == 0, ]` → either `pt_data[pt_data$k == 1, ]` or pass the original subject-level data alongside.

**`R/contrast.R` + S3 methods**
- Reporting time grid: `cut_times[1..K_max]`.
- `summary(time=)` and `contrast(time=)` validators reject `time > T_max` (was: silent extrapolation).

**`R/utils.R` / `R/print.person_time.R`**
- `print.person_time()` header: display `K_max`, `T_max` attributes.
- Helpers filtering on `c_flag` → `dep_cens`.

**Tests** (`tests/testthat/`)
- All references to `y_flag` / `c_flag` columns → renamed.
- Tests hard-coding `k = c(0, 2, 4, ...)` → updated to integer `k = 1..K_max`.
- New tests:
  - `ipcw` accepted as scalar logical and as length-n logical vector.
  - Hard error on `T_max > max(time)`.
  - Warning on `mean(admin_reach) < 0.5`.
  - Admin truncation: subjects with `time >= T_max` produce no exit row; final at-risk row sits at `k = K_max`.
  - C_admin wins ties at T_max (`status=1, time=T_max` → no `y_event=1` row).
  - Cliff regression test: smoke run from commit `084fe98` should now complete with finite weights.

**`dev/smoke_run.R`** — DGP updated; verify `summarize_weights()` no longer reports max ~1e11.

**`dev/simulate_separable_data()` helper** (commit e88ab5b) — DGP updated to new column contract.

**Roxygen** — all column references in roxygen blocks across the codebase. `devtools::document()` regenerates `NAMESPACE` and `man/*.Rd`.

### Refactor execution order (next session)

Each chunk is reviewable in isolation; functional state restored only after the whole run.

1. ~~`to_person_time()` signature + roxygen~~ — **done 2026-05-09** (see TODO.md resume).
2. ~~`to_person_time()` body~~ — **done 2026-05-09** (see TODO.md resume).
3. `validate_subject_level()` — extend for `ipcw` shape + T_max range; replace inline placeholder in `to_person_time()`.
4. `validate_person_time()` updates — column set + mutual-exclusivity invariant + integer `k`.
5. `fit_hazard_models()` — c-hazard / y-hazard fit population changes.
6. `fit_ipw()` cliff-guard removal + column references; `dev/smoke_run.R` re-run.
7. Workers (`make_clone`, `fit_gformula`, `fit_ipw`) — `k`-as-integer, baseline extraction.
8. Contrast + S3 methods — reporting grid, time validators.
9. Tests — column renames + new invariant tests.
10. `devtools::document()`; full `R CMD check`.

---

## Red flags found during dep-graph audit (2026-05-01)

Three port-set functions have hidden coupling that contradicts §12.3's
CLEAN-MOVE classification:

1. **`validate_person_time`** hardcodes `d_flag` as a required column
   (`R/validate.R:313` in CCR). CausalSurvival has no D event — this errors
   on every CausalSurvival pt_data. **Resolved during port** as TRAP-FORK:
   `d_flag` dropped from required-cols + flag-validation loop; treatment
   check tightened from `length(unique) == 2` to `setequal(c(0,1))` to
   match downstream IPW contract. Two copies of a small validator are
   honest — the packages have legitimately different domain models.

2. **`summarize_weights` + `apply_weight_truncation`** hardcode SE arm-weight
   names (`w_d_arm_10`, `w_y_arm_10`, …) in their `intersect()` lookup
   (`R/weights.R:239-242, 309-311` in CCR). **Resolved during port:** trimmed
   hardcoded list to `{w_cens, w_a}` (and `_raw` variants for
   `summarize_weights`).

3. **`predict_hazard_under` (CCR)** uses raw `stats::predict()` while
   `predict_with_warning` wraps it. Historical accident — no methodological
   reason. **Resolved during port:** unified into `predict_counterfactual_hazard`
   that routes through `predict_with_warning`. CCR migration step (§12.1)
   inherits a breaking signature change.

---

## Tier 0 — true leaves (no internal port-set deps)

| Function | File (CCR origin) | Status | Note |
|---|---|---|---|
| `%\|\|%` | utils.R | **ported** (R/utils.R) | |
| `snap_time` | accessors.R | **deferred to Phase 3** | only callers (`contrast()`, `summary()`) are Phase-3 surface; re-port when those land |
| `cumprod_survival` | hazards.R | **ported** (R/hazards.R) | doc tightened: time-order precondition |
| `predict_with_warning` | gformula_core.R → **hazards.R** | **ported** (R/hazards.R) | moved file (logical home is hazards.R) |
| `weighted_hazard_by_k` | ipw_core.R | **deferred** | non-parametric Hajek estimator; v1 IPW path is parametric MSM (`fit_logistic` + `predict_counterfactual_hazard`), so no v1 caller. Re-port if v1.x adds a non-parametric (weighted-KM-style) IPW method |
| `cum_inc_from_weighted` | ipw_core.R | **deferred** | thin wrapper over `weighted_hazard_by_k`; v1 parametric-MSM path reuses `predict_counterfactual_hazard` + `cumprod_survival` (no D factor needed for single-event), so no v1 caller. Re-port paired with `weighted_hazard_by_k` if v1.x adds non-parametric Hajek IPW |
| `check_covariate_quality` | validate.R | **ported** (R/validate.R) | NA + unsupported-type promoted to hard errors (was warnings in CCR) |
| `validate_input_shape` | validate.R | **ported** (R/validate.R) | |
| `check_fitted_positivity` | hazards.R | **ported** (R/hazards.R) | |
| `ipw` | weights.R | **ported** (R/weights.R) | `truncate` arg dropped (dead code); `check_prob_vec` helper extracted |
| `ipw_static_trt` | weights.R | **ported** (R/weights.R) | assumes `A ∈ {0,1}` (Option B, standardized in `to_person_time`); `match()` broadcast replaces named-vector lookup |
| `summarize_weights` | weights.R | **ported** (R/weights.R) | red flag #2 resolved; added p001/p01 for symmetric tail diagnostics |
| `apply_weight_truncation` | weights.R | **ported** (R/weights.R) | red flag #2 resolved |

---

## Tier 1 — depends only on Tier 0

| Function | Internal dep | File | Status | Note |
|---|---|---|---|---|
| `fit_logistic` | `check_fitted_positivity` | hazards.R | **ported** (R/hazards.R) | docstring de-CCR'd |
| `predict_counterfactual_hazard` | `predict_with_warning` | hazards.R | **ported** (R/hazards.R) | renamed from `predict_hazard_under`; signature now requires `label` |
| `ipw_cens` | `ipw` | weights.R | **ported** (R/weights.R) | `truncate` arg dropped (consistent with `ipw()`); `model_num` retained for v1.1 stabilized-IPCW path (Phase 3 wires `c ~ A` per H&R §12.6) |
| `validate_person_time` | `check_covariate_quality` | validate.R | **ported** (R/validate.R) | red flag #1 resolved as TRAP-FORK: `d_flag` dropped, treatment tightened to `{0,1}` (matches downstream IPW contract) |

---

## Tier 2

| Function | Internal dep | File | Status |
|---|---|---|---|
| `fit_propensity` | `fit_logistic` | propensity.R | **ported** (R/propensity.R) |

---

## Pipeline view (data flow, runtime order)

Two parallel fit branches (hazard glm vs propensity glm). Both go
through `check_fitted_positivity` as a **post-fit diagnostic**, then
diverge: hazards feed the counterfactual-prediction → survival path,
propensity feeds the IPW path.

```
   ┌─────────────────┐                     ┌─────────────────┐
   │  fit_logistic   │                     │ fit_propensity  │
   │ (Y or C hazard) │                     │   (P(A=1|L))    │
   └────────┬────────┘                     └────────┬────────┘
            │ glm                                   │ glm
            ▼                                       ▼
   ┌─────────────────────────┐         ┌─────────────────────────┐
   │ check_fitted_positivity │         │ check_fitted_positivity │
   │   (post-fit diagnostic) │         │   (post-fit diagnostic) │
   └────────┬────────────────┘         └────────┬────────────────┘
            ▼                                   ▼
   ┌─────────────────────┐            ┌─────────────────────┐
   │ predict_with_       │            │  ipw_static_trt     │
   │   warning           │            │  (consumes p_full   │
   └────────┬────────────┘            │   = predict(prop))  │
            ▼                         └──────────┬──────────┘
   ┌─────────────────────────────┐               │
   │ predict_counterfactual_     │               │ w_a
   │   hazard                    │               ▼
   └────────┬────────────────────┘    [IPW pipeline below]
            │ haz
            ▼
   ┌─────────────────────┐
   │  cumprod_survival   │
   └─────────────────────┘
```

IPW pipeline (parallel branch):

```
   ipw (core, cumprod ratio)
       │
       ├── ipw_cens (IPCW wrapper)
       │
   [ipw_static_trt is independent — inlines its own math, not via ipw()]

   apply_weight_truncation, summarize_weights — applied to assembled w_* columns
```

---

## Convention reminders for the port

- All ported functions stay `@keywords internal` unless §12.3 says exported.
- Roxygen references to `separable_effects()` / `fit$warnings` etc. get
  generalized (CausalSurvival has no SE orchestrator).
- "Y, D, C hazards" → "Y, C hazards" (no competing event D in CausalSurvival).
- Body code unchanged unless we agreed on a refactor. So far:
  - rename + unification of `predict_counterfactual_hazard` (hazards.R)
  - `truncate` arg dropped from `ipw()` and `ipw_cens()` (dead code in v1)
  - `check_prob_vec()` helper extracted at top of `weights.R`
  - `match()` broadcast replaces named-vector lookup in `ipw_static_trt`
  - `summarize_weights` adds p001/p01 for symmetric tail diagnostics
  - `check_covariate_quality`: NA + unsupported-type promoted from
    warnings to hard errors (glm-breaking otherwise)
  - `validate_person_time` treatment check tightened from
    `length(unique) == 2` to `setequal(c(0,1))`
