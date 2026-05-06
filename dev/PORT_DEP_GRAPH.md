# Port Dependency Graph: CCR → CausalSurvival

Live artifact tracking the dependency-ordered port of 17 functions from
CausalCompetingRisks (CCR) to CausalSurvival, per spec §12.3
(`../../separable_effects/dev/CAUSAL_SURVIVAL_SPEC.md`).

Updated as functions are ported. Status: **port in progress (hazards.R + weights.R done).**

---

## Red flags found during dep-graph audit (2026-05-01)

Three port-set functions have hidden coupling that contradicts §12.3's
CLEAN-MOVE classification:

1. **`validate_person_time`** hardcodes `d_flag` as a required column
   (`R/validate.R:313` in CCR). CausalSurvival has no D event — this errors
   on every CausalSurvival pt_data. **Action:** reclassify as TRAP-FORK or
   parameterize the flag list. Spec §12.3 patch needed.

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
| `%\|\|%` | utils.R | pending | |
| `snap_time` | accessors.R | pending | |
| `cumprod_survival` | hazards.R | **ported** (R/hazards.R) | doc tightened: time-order precondition |
| `predict_with_warning` | gformula_core.R → **hazards.R** | **ported** (R/hazards.R) | moved file (logical home is hazards.R) |
| `weighted_hazard_by_k` | ipw_core.R | pending | |
| `check_covariate_quality` | validate.R | pending | |
| `validate_input_shape` | validate.R | pending | |
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
| `validate_person_time` | `check_covariate_quality` | validate.R | pending | ⚠️ red flag #1 (hardcodes `d_flag`) |

---

## Tier 2

| Function | Internal dep | File | Status |
|---|---|---|---|
| `fit_propensity` | `fit_logistic` | propensity.R | pending |

---

## Pipeline view (data flow, not topological)

```
                    ┌─────────────────────┐
                    │ check_fitted_       │
                    │   positivity        │
                    └──────────┬──────────┘
                               │
                               ▼
                    ┌─────────────────────┐
                    │   fit_logistic      │ ──────────────┐
                    └──────────┬──────────┘               │
                               │                          │
                               │ produces glm             │
                               ▼                          ▼
                    ┌─────────────────────┐    ┌──────────────────┐
                    │ predict_with_       │    │  fit_propensity  │
                    │   warning           │    └──────────────────┘
                    └──────────┬──────────┘
                               │
                               ▼
                    ┌─────────────────────────────┐
                    │ predict_counterfactual_     │
                    │   hazard                    │
                    └──────────┬──────────────────┘
                               │ produces haz
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
