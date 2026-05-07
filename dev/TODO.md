# CausalSurvival — TODO / Future scope

Flagged extensions and known limitations. Not in current scope.

## Resume — Phase 3, fit_ipw() mid-implementation

**Branch**: `claude/general-session-1wHyZ` (force-pushed; pull with `--force` if a previous session is checked out).

**Current state**:
- `causal_survival(method = "gformula")` runs end-to-end and returns a complete `causal_survival_fit` S3 object.
- `causal_survival(method = "ipw")` errors at `stop("fit_ipw: WIP — Y-MSM not yet wired")` in `R/causal_survival.R`. Propensity, optional IPCW C-hazard + numerator, raw weights (`w_a`, `w_a_raw`, `w_cens`, `w_cens_raw`), truncation, and `combined_w` are all wired. Y-MSM fit and per-arm marginalization are the remaining pieces.

**Next two chunks** (proposed in chat history; await `ok`):
- **5c**: weighted Y-MSM via `fit_logistic(..., weights = combined_w)`, then per-arm clone → `predict_counterfactual_hazard` → `cumprod_survival` → `tapply(..., mean)`. Default MSM is marginal in L (`y_flag ~ A + k + I(k^2) + I(k^3)`); user override via `formulas$y`. Builds the `estimates` data frame in the same shape as `fit_gformula` returns.
- **5d**: return list — `estimates`, `models = list(y, c, A, A_num, c_num)`, `model_checks` (parallel keys), `weights = list(pt_data_weighted, weight_summary, flagged_ids, flagged_log, truncate)`. Closes the function for real (drops the WIP `stop()`).

**Open questions parked**:
- Symmetrize `models` list shape across `fit_gformula` and `fit_ipw`? (Currently `fit_gformula` has 4 slots, `fit_ipw` will have 5 — user said "I don't know" when asked.)
- Switch the IPW weighted glm to `quasibinomial` to silence "non-integer #successes"? Audit pending — see "Weighted GLM family" section below.

**After fit_ipw closes** (rough order, all per-spec):
1. Accessors — `risk()`, `contrast()`, `summary()`, `print()`, `plot()` (Okabe-Ito palette).
2. `bootstrap(fit)` returning a standalone `causal_survival_bootstrap` object.
3. `reweight()` helper using preserved `*_raw` weight columns.
4. Tests per the "Tests (per file, planned)" checklist below.
5. Build/release checklist (roxygen, `R CMD check`, bundled dataset, `simulate_causal_survival_data()`, vignette).

**Workflow rule in force this session** — Eyeball-10-30 review (chunk-by-chunk approval, no silent edits, commit at file boundaries). Persist this in CLAUDE.md if you want it carried forward across sessions.

## Multi-arm treatment

Currently `A` is binary, standardized to `{0, 1}` in `to_person_time` (hard
error on >2 levels). Both `gfoRmula` and `gfoRmulaICE` support multi-arm
(integer-coded, e.g. `{0, 1, 2, ...}`); the pattern is to enumerate each
regime separately rather than parameterize a single contrast.

Extension path when needed:

1. Relax `to_person_time` standardization from `{0, 1}` to `{0, 1, ..., K-1}`.
2. Treatment model: switch from binomial logistic to multinomial (e.g.
   `nnet::multinom`) — equivalent to `gfoRmula`'s `covtypes = "categorical"`.
3. Regimes: enumerated per-arm rather than a `{0, 1}` contrast (cf.
   `interventionK.A = list(static, rep(level, T))` in gfoRmula).
4. Estimands: pairwise contrasts selectable, with a reference arm.

No architectural lock-in from the current binary-only choice — the
standardization site is the only place that needs to relax.

## Optional pass-through of original `time` / `status` in `to_person_time()`

v1 drops the user's original `time` and `status` columns from the
person-time output (they're redundant with `k` + `y_flag` / `c_flag`
once discretized). Output schema is `id, k, treatment, <covariates>,
y_flag, c_flag`.

**Useful add:** an optional `keep_subject_columns = TRUE` flag (or
similar) that broadcasts the original `time` and `status` to every row.
Lets users cross-check the discretization, diagnose binning artifacts,
and (most importantly) **validate the DGP simulator's discretization**
against known continuous truth when we build
`simulate_causal_survival_data()`.

Defer until the DGP work — that's when we'll actually want it.
Implementation seam already in `to_person_time()` body (the explicit
`pt_data[[time]] <- NULL` / `pt_data[[status]] <- NULL` block).

## Single-arm / natural-course analysis

`validate_subject_level()` and `validate_person_time()` currently hard-error
when treatment has fewer than 2 unique values. This blocks legitimate
single-arm use cases:

- **Natural-course / observational cumulative incidence** — no contrast,
  just the marginal CIF for a single cohort. Useful for descriptive
  reporting or as a baseline against intervened-arm estimates from a
  separate fit.
- **Alternative marginalization / transportability** — standardize a
  single arm to a target covariate distribution that differs from the
  source population.

Extension path:

- Relax the `< 2 unique values` check to a warning, or gate behind an
  explicit `single_arm = TRUE` flag.
- `causal_survival()` returns per-arm CIF only (no contrast object).
- IPW path: `ipw_static_trt` doesn't apply (no propensity contrast);
  only g-formula and censoring weights remain meaningful. Or restrict
  single-arm to `method = "gformula"`.
- Phase 3 design point: how does the orchestrator handle `length(arms)
  == 1`? Probably a separate code path with stripped accessors
  (`contrast()` would error / return NULL).

Defer to v1.x once the binary-arm case is solid.

## Treatment label preservation for display

Treatment is standardized to `{0, 1}` numerically (Option B) — clean
for downstream compute (`ipw_static_trt`, glm positive class,
counterfactual assignments). But the user's original labels
(e.g. `c("ctrl", "trt")`) are lost.

Plan:

- `to_person_time()` stashes `attr(pt_data, "treatment_levels") <-
  c(<level for A=0>, <level for A=1>)` on the returned object.
- Phase 3 `print()` / `summary()` / `plot()` consult the attribute and
  relabel arms for display. Falls back to `"0"` / `"1"` if attribute
  missing.
- Users who bypass `to_person_time()` and supply their own pt_data must
  pre-code `{0, 1}` (validator enforces). v1.1: optional
  `treatment_levels = c("Control", "Treated")` arg on the public entry
  point as a fallback for that case.

Compute path stays plain — labels are presentation-only metadata.

## Phase 3 design points (not blocking the port)

### Stabilization is a single joint switch in v1

`causal_survival()` should expose `stabilize = "marginal"` / `NULL` as a
single switch driving both treatment and censoring numerators jointly:

| `stabilize` | A_num (helper sees) | C_num (helper sees) |
|---|---|---|
| `"marginal"` | `glm(A ~ 1, ...)` (auto) | `glm(c_flag ~ A, ...)` (auto) |
| `NULL` | `model_num = NULL` | `model_num = NULL` |

User customizes Y-hazard, C-hazard (denominator), and treatment
denominator + numerator via `formulas` (`y`, `c`, `A`, `A_num`).
Censoring numerator is internal — there's no `c_num` slot in
`formulas` because v1 doesn't allow conditional stabilization on V,
and `c ~ A` is the only methodologically standard choice (per H&R
§12.6, Technical Point 12.2).

Conditional stabilization (V in numerators + Y-MSM as effect modifier)
is the v1.1 extension, paired with the subgroup / CATE API. The
matching rule (numerator V ⊆ MSM conditioning) means stabilization
choice and Y-MSM shape must move together — this is why a single
joint `stabilize` switch is the right v1 abstraction.

The internal IPW helpers (`ipw_static_trt`, `ipw_cens`) keep the
flexible `model_num` argument so the v1.1 extension lands without
helper rewrites. `causal_survival()` is the place where the joint
construction logic lives.

### Fit object structure (Phase 3)

Inherit CCR's shape (see `../../separable_effects/dev/TODO.md`):

- `fit$cumulative_incidence` as a named list per method run:
  `$gformula` (NULL if not run), `$ipw` (NULL if not run). Plot / print /
  summary iterate over list entries. Slot keys match the user-facing
  `method` arg values (no underscore variant) so the same string
  identifies a method end-to-end.
- `fit$weights` (IPW runs only): `pt_data_weighted` with raw + truncated
  weight columns preserved, `weight_summary`, `truncated_ids`,
  `extreme_weight_adjust`, `extreme_weight_threshold`. Enables
  `reweight()` without refitting.
- `fit$warnings`: collected via `withCallingHandlers()` wrapper around the
  `causal_survival()` body, muffled during run, re-emitted as a single
  group at end. `print(fit)` shows count.
- `fit$model_diagnostics` (planned): post-fit GLM checks consolidated.
- No `fit$contrasts` slot — contrasts computed on demand via
  `contrast(fit, ci = bootstrap(fit))`.

### Bootstrap as standalone object

Same as CCR — `bootstrap(fit)` returns a `causal_survival_bootstrap` (or
similar) object, NOT attached to `fit`. Keeps `fit` pure, makes optional
re-runs explicit. `plot(fit, ci = boot)`, `contrast(fit, ci = boot)`.

### `reweight()` helper

`reweight(fit, extreme_weight_adjust, extreme_weight_threshold)` re-applies
truncation to the preserved raw weights and re-runs the IPW estimator. Fast
sensitivity analysis without refitting. Same plumbing as CCR.

### Post-fit GLM checks (beyond positivity)

`check_fitted_positivity` ported. Still missing — bundle into a single
`check_fitted_quality()` returning a structured diagnostics list:

- Collinearity: glm dropped a term → `NA` coefficient
- Convergence: `!model$converged`
- Separation: extreme standard errors

Surface as `fit$model_diagnostics`.

### Plot color palette

Okabe-Ito subset (vs CCR's 4-arm palette):

- A=1 line: `#D55E00` (vermillion)
- A=0 line: `#0072B2` (blue)
- Per-method overlays use linetype, not color

CI ribbons use the corresponding fill at low alpha.

### Weighted GLM family: `binomial` vs `quasibinomial`

`fit_logistic()` uses `family = binomial("logit")` for both unweighted
and weighted fits. With non-integer IPW weights, `binomial` fires the
"non-integer #successes in a binomial glm!" warning per fit. Point
estimates from `quasibinomial("logit")` are identical; only the
dispersion parameter and naive SE/inference machinery differ — and we
don't use parametric inference (CIs come from the bootstrap).

**Decision deferred.** Before settling, audit how peer packages handle
this:

- `gfoRmula` / `gfoRmulaICE` — what family do they use for weighted
  pooled-logistic MSM fits?
- `ipw` / `ipwErrorY` — same question for the outcome model.
- `survey::svyglm` — uses `quasibinomial` by design for weighted
  regression; is that the convention we want to adopt?
- Robins/Hernán code (e.g. `causaldata`-shipped scripts) — likely
  `binomial` with the warning silenced, but worth confirming.

If the field convention is `quasibinomial` for weighted fits, switch
the IPW path only (keep `binomial` for unweighted Y/C-hazard and
propensity fits to avoid changing existing behavior).

---

## Build / release checklist (v1)

- [ ] `roxygen2::roxygenise()` populates `man/` and `NAMESPACE`
- [ ] `R CMD check` clean (no NOTEs, WARNINGs, ERRORs)
- [ ] Bundled survival example dataset (`data/<name>.rda`) — likely a subset
      of `survival::veteran` or similar, with documentation in `R/data.R`
- [ ] Internal `simulate_causal_survival_data()` helper for unit tests
      (known true effects, edge cases). Not exported.
- [ ] Vignette: walkthrough on the bundled dataset, both methods, with /
      without bootstrap
- [ ] README.md with install + minimal example
- [ ] NEWS.md
- [ ] LICENSE.md
- [ ] `pkgdown` site (optional for v1)
- [ ] CRAN / r-universe submission post-stabilization

## Tests (per file, planned)

- [ ] `test-utils.R` — `%||%`
- [ ] `test-validate.R` — input shape, person-time validator (left-trunc,
      treatment {0,1}, NA hard errors, type hard errors)
- [ ] `test-hazards.R` — `fit_logistic`, `predict_counterfactual_hazard`,
      `cumprod_survival`, `check_fitted_positivity`
- [ ] `test-weights.R` — `ipw`, `ipw_static_trt`, `ipw_cens`, truncation,
      `summarize_weights`
- [ ] `test-propensity.R` — `fit_propensity`
- [ ] `test-to_person_time.R` (Phase 2) — long-format builder, treatment
      standardization to {0,1}, `treatment_levels` attribute
- [ ] `test-causal_survival.R` (Phase 3) — orchestrator, both methods,
      joint stabilization switch
- [ ] `test-accessors.R` (Phase 3) — `risk()`, `contrast()`, `summary()`,
      `print()`, `plot()` smoke tests
- [ ] `test-bootstrap.R` (Phase 3) — resampling, percentile CI

## Lineage note

CCR's `dev/TODO.md` §"Future: shared causal_tools" anticipated factoring
out exactly the functions ported into this package (`to_person_time`,
`validate_*`, pooled-logistic hazard fitting + positivity check, `%||%`,
`weighted_hazard_by_k`-style estimators, cumulative-incidence shape).
CausalSurvival is functionally that shared package for the survival-only
(no D event) case. If a third causal package ever needs the same building
blocks, the CCR-side note becomes the relevant design doc — factor out
properly into a low-level `causal_tools` package at that point.
