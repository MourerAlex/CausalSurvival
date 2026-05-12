# Port Dependency Graph: CCR → CausalSurvival

Live artifact tracking the dependency-ordered port of 17 functions from
CausalCompetingRisks (CCR) to CausalSurvival, per spec §12.3
(`./CAUSAL_SURVIVAL_SPEC.md`).

Updated as functions are ported. Status: **Phase 3 in progress**.
Phase 1 + Phase 2 complete — all port-set functions either ported or
deferred (snap_time still deferred until accessors land in Phase 3);
`weighted_hazard_by_k` and `cum_inc_from_weighted` un-deferred
2026-05-12 to support the new `ipw_engine = "km"` (see
`dev/ipw-implementation-spec_1.md`). Phase 2 trap-fresh-writes landed
(`to_person_time`, `validate_subject_level`, `standardize_treatment`,
`fit_hazard_models`).

Phase 3 progress:
- `causal_survival()` skeleton: signature + arg validation written; worker
  dispatch + S3 fit object assembly pending.
- G-formula path: `make_clone` (TRAP-GENERALIZE #4, inlined fresh-write —
  single-event simplification of CCR's two-arm clone) + `fit_gformula`
  worker written. CIF math validated against gfoRmula's `simulate.R`
  (per-subject `poprisk_i = cumsum(h * S(t-1))` ↔ ours `1 - mean_i S_i(t)`,
  algebraically identical for single event).
- IPW path: `fit_ipw()` (MSM) landed. Smoke test (commit `084fe98`)
  revealed IPCW weight cliff at K_max from the v0.1.0 admin-censoring
  encoding → triggered data-convention overhaul (below). New
  `ipw_engine = "km"` engine added per
  `dev/ipw-implementation-spec_1.md`; helpers
  `weighted_hazard_by_k` + `cum_inc_from_weighted` un-deferred + ported
  to `R/hazards.R` (2026-05-12). Worker integration
  (`fit_ipw_km` + shared `fit_ipw_weights`) bundled into step 7 of the
  overhaul below.
- Data-convention overhaul (below): steps 1+2 done 2026-05-12
  (`to_person_time()` signature + roxygen + inferred-mode body).
- Accessors / bootstrap / S3 methods: not yet started.

---

## ⚠️ Phase 3 BLOCKED on data-convention overhaul (2026-05-08)

The smoke test exposed a structural problem with the `y_flag` / `c_flag`
encoding: under the v0.1.0 contract, administratively censored subjects
get `c_flag = 1` at the last interval, which collides with the c-hazard
fit (`h_C(K_max) ≈ 1` → IPCW weight `1 / (1 - h_C) ≈ 1e11`).

After methodology review, the data convention has been overhauled. The
canonical spec now lives in `./CAUSAL_SURVIVAL_SPEC.md` § 3.0 (LOCKED
2026-05-08). Highlights:

- Column schema changes: `y_flag` / `c_flag` → `y_event` / `dep_cens` / `indep_cens`.
- Structural ordering C_admin → C_dep → Y; admin censoring placed first within each interval.
- `k` becomes integer interval index `1..K_max` (was: left-edge interval time).
- Admin-truncated subjects (`time > T_max`) contribute at-risk rows up to `k = K_max` with no exit row. **No K_end interval is materialized** (§3.0.8).
- `ipcw` arg (scalar logical or per-subject logical vector, default `TRUE` = dependent censoring); `T_max` arg with default `max(time)`.
- Boundary at `t = T_max`: events at `time = T_max` fire normally at `k = K_max`; only `time > T_max` is admin-truncated (§3.0.1).
- Hard error if `T_max > max(time)`; hard error on any `time = 0`; warning if `mean(admin_reach) < 0.5`.
- Cliff structurally impossible: `indep_cens = 1` rows never enter the c-hazard fit (model fits on `indep_cens = 0` rows only).

### Phase 3 refactor inventory (per-file changes)

**`R/data_prep.R`**
- `to_person_time()`: full rewrite per spec §3.0.5. New args: `ipcw`, `T_max`. Pre-split mode (`event_cols`) dropped per spec §3.0.9. Output columns `y_event` / `dep_cens` / `indep_cens`. Admin-truncated subjects materialized as at-risk rows up to `k = K_max` with no exit row (no K_end). Hard error on any `time = 0`. **Steps 1+2 complete 2026-05-12.**
- `validate_subject_level()`: new checks for `ipcw` shape (scalar logical or length-`n` logical vector, no NA); `T_max ≤ max(time)` (hard error); `mean(admin_reach) < 0.5` (warning). Currently still has the v0.1 signature; checks live inline in `to_person_time()` pending step 4.
- `validate_person_time()`: column set `{y_event, dep_cens, indep_cens}`; mutual-exclusivity invariant; `k ∈ {1, ..., K_max}`. Left-truncation check (was `k = 0`) updated to `k = 1`.

**`R/hazards.R`**
- `fit_hazard_models()`: c-hazard becomes `glm(dep_cens ~ ...)` on rows with `indep_cens = 0`; y-hazard becomes `glm(y_event ~ ...)` on rows with `indep_cens = 0` AND `dep_cens = 0`. Polynomial-in-`k` interpretation shifts from time-value to integer index; formula form unchanged.

**`R/weights.R`**
- `ipw_cens()`: column references `c_flag` → `dep_cens`. No cliff guard needed — replaced by structural exclusion of `indep_cens = 1` rows from the c-hazard fit.
- `ipw()`, `ipw_static_trt()`: column renames only; logic stable.

**`R/causal_survival.R` workers (gformula, ipw)**
- `make_clone()`: clone `k` values become integer `1..K_max`. Predictions reported over the full `cut_times[1..K_max]` grid.
- G-formula worker: predicts `y_event` hazard on the clone for `k = 1..K_max`. Cumprod survival over those K_max steps.
- IPW worker: weighted Y-MSM glm fits on rows with `indep_cens = 0` AND `dep_cens = 0` only. Combined weight `w_a * w_cens` as before. Also gains `ipw_engine = "km"` engine per `dev/ipw-implementation-spec_1.md` (un-defers `weighted_hazard_by_k` + `cum_inc_from_weighted`, ported 2026-05-12).
- Baseline-covariate extraction: `pt_data[pt_data$k == 0, ]` → `pt_data[pt_data$k == 1, ]`.

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
  - `ipcw` accepted as scalar logical and as length-n logical vector; NA in `ipcw` errors.
  - Hard error on `T_max > max(time)`.
  - Hard error on any `time = 0`.
  - Warning on `mean(admin_reach) < 0.5`.
  - Admin-truncation invariant: subjects with `time > T_max` have at-risk rows up to `k = K_max` with all three flags `= 0` and no exit row.
  - Boundary at `T_max`: `status=1, time=T_max` → `y_event=1` at `k = K_max` (NOT admin-truncated).
  - Cliff regression test: smoke run from commit `084fe98` should now complete with finite weights.

**`dev/smoke_run.R`** — DGP updated; verify `summarize_weights()` no longer reports max ~1e11.

**`dev/simulate_separable_data()` helper** (commit e88ab5b) — DGP updated to new column contract.

**Roxygen** — all column references in roxygen blocks across the codebase. `devtools::document()` regenerates `NAMESPACE` and `man/*.Rd`.

### Refactor execution order

Each chunk is reviewable in isolation; functional state restored only after the whole run.

1. `to_person_time()` signature + roxygen. **Done 2026-05-12.**
2. `to_person_time()` inferred-mode body. **Done 2026-05-12.**
3. ~~`to_person_time()` pre-split mode body + `event_cols` validator.~~ **Dropped** per spec §3.0.9 (pre-split mode deferred to v1.x).
4. `validate_person_time()` updates — column set + mutual-exclusivity + integer `k` starting at 1.
5. `fit_hazard_models()` — c-hazard / y-hazard fit population changes.
6. `fit_ipw()` cliff-guard removal + column references; `dev/smoke_run.R` re-run.
7. Workers (`make_clone`, `fit_gformula`, `fit_ipw`) — `k`-as-integer, baseline extraction. Also wires new `ipw_engine = "km"` worker (`fit_ipw_km`) + shared `fit_ipw_weights` helper per `dev/ipw-implementation-spec_1.md`.
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
| `weighted_hazard_by_k` | ipw_core.R | **scheduled (un-deferred 2026-05-12)** | non-parametric Hajek estimator. Required by new `ipw_engine = "km"` (default IPW engine) per `dev/ipw-implementation-spec_1.md`. Lands in `R/hazards.R` |
| `cum_inc_from_weighted` | ipw_core.R | **scheduled (un-deferred 2026-05-12)** | thin wrapper over `weighted_hazard_by_k`. Single-event variant: `d_event` arg dropped (no competing event D in CausalSurvival). Required by `ipw_engine = "km"`. Lands in `R/hazards.R` |
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
