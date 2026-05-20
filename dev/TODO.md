# CausalSurvival — TODO / Future scope

Flagged extensions and known limitations. Not in current scope.

## Resume — Phase 3 worker layer complete; accessors next

**Branch**: `claude/general-session-1wHyZ` (force-pushed during this session — pull with `--force` if a previous checkout exists). Last commit: `4511a98 feat: fit_ipw() — return list, closes the IPW worker`.

**Current state**:
- `causal_survival(method = "gformula")` and `causal_survival(method = "ipw")` both run end-to-end and return a `causal_survival_fit` S3 object.
- IPW path: propensity (+ optional `A ~ 1` numerator), optional IPCW C-hazard (+ optional `c_flag ~ A` numerator), `ipw_static_trt` and `ipw_cens` weights, `apply_weight_truncation` (preserves `*_raw` columns for the future `reweight()` seam), weighted Y-MSM (`y_flag ~ A + k + I(k^2) + I(k^3)` by default; user can override via `formulas$y` for a covariate-conditional MSM), per-arm clone → `predict_counterfactual_hazard` → `cumprod_survival` → marginalize.
- Joint stabilization is the single `stabilize` switch (`"marginal"` or `NULL`); both numerators move together.
- Nothing was sanity-checked against R — no R interpreter in the sandbox where this was written. **First action next session: run the package on a tiny synthetic dataset for both methods.**

**Next priorities** (from the Phase 3 design points + checklists below):
1. Accessors: `risk()`, `contrast()`, `summary()`, `print()`, `plot()` (Okabe-Ito palette).
2. `bootstrap(fit)` returning a standalone `causal_survival_bootstrap` object (NOT attached to `fit`).
3. `reweight(fit, truncate = ...)` helper using preserved `*_raw` weight columns.
4. Tests per the "Tests (per file, planned)" checklist below.
5. Build/release: roxygen, `R CMD check`, bundled dataset, `simulate_causal_survival_data()`, vignette.

**Open questions still parked**:
- Symmetrize `models` list shape across workers? `fit_gformula` returns 4 slots (`y, c, A, A_num`), `fit_ipw` returns 5 (adds `c_num`). User said "I don't know" — left asymmetric.
- Switch IPW weighted glm to `quasibinomial` to silence "non-integer #successes"? Audit pending — see "Weighted GLM family" section below.
- `apply_weight_truncation()` internally returns `flagged_ids` and `flagged_log`; the orchestrator only surfaces `truncated_ids` (= `flagged_ids`). The per-row `flagged_log` detail is dropped — re-add if a diagnostics accessor needs it.

**Workflow rule from this session — Eyeball-10-30 review**:
1. Show ~10–30 lines at a time in chat as a fenced code block before writing to disk. One logical unit per chunk.
2. Wait for explicit approval (`ok`, `yes`, `go`, `agree`) before calling Write/Edit. No silent file edits.
3. Flag design choices baked into the chunk in 2–4 bullets after the code.
4. Commit at file boundaries, not at every chunk.
5. No drive-by changes outside the chunk being shown.
6. Split chunks longer than ~30 lines.
7. After approval of chunk N, immediately propose chunk N+1.
8. Preserve user-supplied snippets verbatim if they ask to "keep this somewhere".

Persist this in `CLAUDE.md` to carry it across future sessions automatically.

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
5. `causal_contrast()` accepts general `reference` / `contrasts` per
   spec §3.4 line 266 — reinstate the drafted helper at
   `dev/unused_code/resolve_contrast_pairs_v2.R` (deferred from
   v0.1.0 to keep the binary expansion inline).

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

## v2 — time-varying covariates and treatment

Out of scope for v1 (spec §1.1) but the v1 design should not paint
itself into a corner.

- **Formula RHS / `time_terms`.** The hard-coded `"k + I(k^2) + I(k^3)"`
  in `fit_hazard_models()` is the default time-trend basis. Override
  is already possible via `formulas$y` / `formulas$c`, but the *shape*
  of the default should become a knob — e.g. `time_basis =
  c("poly3", "ns_df4", "factor_k")` — once we have more than one
  realistic option.
- **Time-varying baseline covariates `L_k`.** Today `covariates` is a
  baseline-only channel broadcast to every person-time row. v2 needs
  a separate `time_varying_covariates` argument (per-row, can vary
  with `k`) and a formula RHS that references those columns directly.
  Hazard fits remain pooled-logistic; the row-level RHS just changes.
- **Time-varying treatment `A_k`.** Propensity moves from baseline
  `fit_propensity()` to per-interval `fit_propensity_k()` (or a
  pooled-logistic propensity with `A_k ~ L_k + history`). Weights
  become time-varying per row (cumulative product across `k`).
  `ipw_static_trt()` is replaced by a time-varying weight builder.
- **`ipw_engine = "km"` constraint.** Weighted KM is only well-defined
  for baseline-`A`; under time-varying `A` the engine must hard-error
  and route the user to `ipw_engine = "msm"`. This guard is already
  noted in `dev/ipw-implementation-spec_1.md` open questions.
- **Stabilization for time-varying `A`.** Numerator becomes
  `A_k ~ A_{k-1} + V` (Robins/Hernán §12); current marginal numerator
  collapses to a degenerate case.
- **Fit-population restrictions baked into the hazard fitters.** Today
  `fit_hazard_models()` restricts the Y-hazard fit to
  `indep_cens == 0 & dep_cens == 0` rows and the C-hazard fit to
  `indep_cens == 0` rows; `fit_ipw_msm()` and `fit_ipw_km()` use the
  same masks for the weighted Y-MSM / KM fits. These masks assume the
  three-way censoring split is row-level immutable (the indep/dep
  label travels with each at-risk row through follow-up). Under
  time-varying covariates or treatment, the label may become
  time-varying (a subject's censoring mechanism could switch reasons
  mid-follow-up). The current hard-coded masks would silently
  mis-restrict. v2 needs an explicit feature gate — a `cens_split`
  arg, per-row labels, or a flag that toggles between the v1 baseline-
  immutable convention and a future time-varying one — before the
  time-varying refactor.

Defer all of the above to v2 — flag in roxygen / spec when v1 code
makes a baseline-only assumption that v2 will need to relax.

---

## Bootstrap polish (currently STUB)

`R/bootstrap.R` is a stub: serial loop, `warnings_count = NA_integer_`.
Progress messages are ported from `separable_effects/R/bootstrap.R`
lines 109-136. Polish to match spec §3.3:

- Parallel via `future.apply::future_lapply(future.seed = seed)` for
  reproducible L'Ecuyer streams. (Progress cadence will need to
  switch from in-loop `message()` to a `progressr` handler when this
  happens.)
- Per-replicate `warning()` capture into `warnings_count` (current
  stub suppresses warnings silently inside the replicate fitter).
- Inspect what `inherits(boot_data, "person_time")` actually needs
  preserved across `rbind` — current stub copies all non-built-in
  attributes. May be fragile.

---

## Estimand framing — switch to expectation / risk notation

`causal_survival()` roxygen frames the estimand as the survival
function `S^a(t) = P(T^a > t)`. Correct but not the canonical
causal-inference framing (Hernán & Robins, Stensrud), which leads
with risk / expectation:

- per-interval risk under arm `a` with treatment-dependent censoring
  set to 0: `E[Y_k^{a, c = 0}] = P(Y_k^{a, c = 0} = 1)`
- contrast: `E[Y_k^{a=1, c=0}] - E[Y_k^{a=0, c=0}]` (risk
  difference) and ratios thereof
- under our three-way censoring split, `c` is the
  treatment-dependent component `c^d` (`dep_cens`); the independent
  component `c^i` (`indep_cens`) is handled by the at-risk set, not
  by intervention
- competing-event / separable extensions are expressed on
  `E[Y_k^{a_Y, a_D, c = 0}]`

Rewrite the description + `@details` to lead with
`E[Y_k^{a, c = 0}]` (or equivalently the cumulative incidence
`F^a(t) = P(T^a \le t, C^{d,a} = 0)`) and demote `S^a(t)` to a
secondary view.

---

## Methodology vocabulary audit — survey-sampling vs survival register

The package occasionally uses methodologically-correct but pedagogically
heavy register from survey-sampling literature (e.g. "Hájek pooled-hazard
estimator") where the survival / causal-inference literature uses simpler
canonical names (e.g. "weighted Kaplan-Meier" / "IPW KM" — Hernán & Robins
ch. 17; Cole & Hernán 2004).

**Already done (2026-05-19)**: replaced "Hájek pooled-hazard estimator"
with "weighted pooled-hazard Kaplan-Meier estimator" across `R/*.R` files.
The Hájek qualifier is technically accurate (denominator = sum of weights
vs Horvitz-Thompson's N) but unfamiliar in survival writing. Removal is a
readability win, not a correctness fix. See
`memory/project_methodology_vocabulary.md` for the rationale.

**TODO**: do a fresh scan across `R/*.R` and roxygen for other terms that
may carry the same register mismatch. Candidates to check:

- Mathematical notation choices that lean survey-sampling rather than
  survival (e.g. explicit ratio-estimator notation where a hazard/CIF
  formulation would read more naturally to the target audience).
- Any place where a less-familiar canonical name is used when a standard
  H&R / Andersen-Ravn / Geskus term exists.
- Internal column naming (`y_event`, `dep_cens`, `indep_cens`, etc.) —
  these are package-internal and intentionally explicit; do NOT rename
  without a separate rule-change pass.

Apply replacements `.R`-side only (dev/*.md mentions stay until a
documentation sweep).

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

## Plot UX docs — `cut_times` argument and step-vs-smooth survival

- [ ] Vignette FAQ: "Why is the IPW-MSM survival curve stepped when
  the hazard model is parametric in `k`?" Short answer to embed:
  - The MSM's `y_event ~ A + k + I(k^2) + I(k^3)` is a parametric
    *working model for the hazard*, used internally so the hazard
    has fewer parameters than a fully saturated `factor(k)` MSM
    (matters when `K_max` is large).
  - The *estimand* is `E[Y_k^{a, c=0}]` on the discrete grid
    `k = 1..K_max` per the `(k-1, k]` convention. It is not defined
    between cut times — the at-risk set and event-counting rule are
    interval-based.
  - Reading the hazard at fractional `k` (e.g., 1.5) would be an
    extrapolation/interpolation choice that was *not* part of the
    identification argument. The IPW weights were built to identify
    the discrete-grid value.
  - Visual consistency: the g-formula and IPW-KM estimators are
    genuinely step functions; drawing the MSM smoothly would
    misleadingly imply it estimates a different kind of quantity.
  - Reference: Hernán & Robins ch. 17, discrete-time MSMs are
    conventionally reported on the analysis grid.
  - **TL;DR**: smooth hazard (internal modeling choice), stepped
    survival (output lives on the cut grid).
- [ ] (v0.2 idea) Optional `plot(..., smooth = TRUE)` arg for the MSM
  path only that evaluates the hazard at a fine subgrid and draws a
  continuous curve, with an explicit caveat about extrapolation.
  Off by default.

- [ ] Vignette section explaining `plot(..., cut_times = ...)`:
  - `NULL` (default): show baseline `t = 0` plus every fit cut time.
  - Single positive integer `N`: COUNT mode — show baseline plus `N`
    indices equidistant along `fit$cut_times`.
  - Numeric vector (length >= 2): EXPLICIT SUBSET — values must be in
    `c(0, fit$cut_times)`; `0` is auto-prepended if absent so the
    baseline column is always shown.
  - **Gotcha**: a length-1 numeric is interpreted as a count, not a
    subset. To display a single specific cut time, pair it with `0`,
    e.g. `cut_times = c(0, 5)`. Worth a worked example in the
    vignette since this is the only non-obvious user-facing wrinkle.
- [ ] Optional: add a `risk_table_max_labels` knob that caps label
  density when `cut_points` is large (15+) and uses the same
  equidistant-pick rule as count mode under the hood. Currently the
  user has to pass `cut_times = N` explicitly to thin the labels.

## Lineage note

CCR's `dev/TODO.md` §"Future: shared causal_tools" anticipated factoring
out exactly the functions ported into this package (`to_person_time`,
`validate_*`, pooled-logistic hazard fitting + positivity check, `%||%`,
`weighted_hazard_by_k`-style estimators, cumulative-incidence shape).
CausalSurvival is functionally that shared package for the survival-only
(no D event) case. If a third causal package ever needs the same building
blocks, the CCR-side note becomes the relevant design doc — factor out
properly into a low-level `causal_tools` package at that point.
