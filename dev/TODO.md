# CausalSurvival — TODO / Future scope

Flagged extensions and known limitations. Not in current scope.

## Resume — Phase 3 data-convention overhaul: `to_person_time()` body done; validators next

**Branch**: `claude/gifted-wescoff` (worktree). Last commit at session start: `70bce71 docs: Phase 3 BLOCKED — record data-convention overhaul + cliff smoke run`.

**What landed this session (`R/data_prep.R` only)**:
- `to_person_time()` full rewrite per spec §3.0.5:
  - New args: `ipcw` (logical scalar or per-subject vector; replaces `cens_type`), `T_max` (default `max(time)`).
  - Output columns: `y_event` / `dep_cens` / `indep_cens` (was `y_flag` / `c_flag`).
  - Integer `k = 1..K_max` (was time-value, 0-based). No `K_end` materialization.
  - Inline `ipcw` + `T_max` validation (placeholder until `validate_subject_level()` is extended).
- Interval convention switched from `[a, b)` to **`(a, b]`** (left-open right-closed). Aligns with survSplit native, ABGK / Andersen-Ravn / gfoRmula / ipw / mstate / pammtools. Adversarial check via 3 parallel agents (local resources + papers + packages) refuted the earlier claim that `[a, b)` was the literature standard.
- Boundary handling at `T_max` under `(a, b]`:
  - `time = T_max` is **inside** the last interval → events fire normally at `k = K_max`, censoring encodes normally.
  - `time > T_max` → admin-censored (no exit row). No subject is ever dropped.
- `time = 0` events: hard error (`(0, t_1]` excludes 0). Lists affected subject id(s).

**Next session — refactor execution order**:
1. **`validate_subject_level()`** — extend for `ipcw` shape (scalar logical or length-`n` logical vector); `T_max ≤ max(time)` (hard error); `mean(admin_reach) < 0.5` (warning). Replaces the inline placeholder block in `to_person_time()`.
2. **`validate_person_time()`** — column set `{y_event, dep_cens, indep_cens}`; mutual-exclusivity invariant (at most one indicator per row); `k ∈ {1, ..., K_max}`.
3. **`fit_hazard_models()`** — c-hazard becomes `glm(dep_cens ~ ...)` on rows with `indep_cens = 0`; y-hazard becomes `glm(y_event ~ ...)` on rows with `indep_cens = 0` AND `dep_cens = 0`. Polynomial-in-`k` interpretation shifts from time-value to integer index; formula form unchanged.
4. **`fit_ipw()`** — cliff-guard removal (no longer needed; structural exclusion replaces it). Column references `c_flag` → `dep_cens`. `dev/smoke_run.R` re-run.
5. **Workers** (`make_clone`, `fit_gformula`, `fit_ipw`) — `k` as integer `1..K_max`. Baseline-covariate extraction: `pt_data[pt_data$k == 0, ]` → `pt_data[pt_data$k == 1, ]` or pass original subject-level data alongside.
6. **Contrast + S3 methods** — reporting grid `cut_times[1..K_max]`. `summary(time=)` and `contrast(time=)` reject `time > T_max`. Display layer time-grid notation: `(0, T_max]`.
7. **Tests** — column renames; new invariants (ipcw scalar/vector accepted; admin censoring at `time > T_max`; `time = 0` event hard error; `time = T_max` event fires at `k = K_max`).

**Mostly cascade work**: ~70% from the ipcw consolidation (three-way censoring split), ~20% from integer-`k` change, ~10% from `(a, b]` convention. Same files touched in the same pass.

**Workflow rules from this session**:

*Eyeball-10-30 review*:
1. Show ~10–30 lines at a time in chat as a fenced code block before writing to disk. One logical unit per chunk.
2. Wait for explicit approval (`ok`, `yes`, `go`, `agree`) before calling Write/Edit. No silent file edits.
3. Flag design choices baked into the chunk in 2–4 bullets after the code.
4. Commit at file boundaries, not at every chunk.
5. No drive-by changes outside the chunk being shown.
6. Split chunks longer than ~30 lines.
7. After approval of chunk N, immediately propose chunk N+1.
8. Preserve user-supplied snippets verbatim if they ask to "keep this somewhere".

*Show-only-.R-diffs* (added mid-session): when working on multi-file refactors that touch spec / TODO / PORT_DEP_GRAPH alongside `R/`, show only the `R/` diffs in chat. Spec / TODO / PORT_DEP_GRAPH edits land silently. Do NOT recap them in post-edit summaries either.

*No silent data drops*: in data-prep / validators, anomalous time values (events at `time = 0`, `time > T_max` if structurally disallowed, etc.) raise hard errors with affected subject id(s). Reserve silent behavior only for cases that are structurally well-defined and always-correct (e.g. admin censoring of `time > T_max` subjects via no-exit-row — that's not "silent drop", it's "no exit row by definition").

Persist these in `CLAUDE.md` to carry across future sessions automatically.

**Open questions still parked** (from previous resume):
- Symmetrize `models` list shape across workers? `fit_gformula` returns 4 slots (`y, c, A, A_num`), `fit_ipw` returns 5 (adds `c_num`). User said "I don't know" — left asymmetric.
- Switch IPW weighted glm to `quasibinomial` to silence "non-integer #successes"? Audit pending — see "Weighted GLM family" section below.
- `apply_weight_truncation()` internally returns `flagged_ids` and `flagged_log`; the orchestrator only surfaces `truncated_ids` (= `flagged_ids`). The per-row `flagged_log` detail is dropped — re-add if a diagnostics accessor needs it.

## Display-layer sweep for `(0, T_max]` interval notation

Spec §3.0.2 locks the interval convention as `(a, b]` (left-open
right-closed) — matches survSplit native, ABGK / Andersen-Ravn /
gfoRmula / ipw / mstate / pammtools. Body of `to_person_time()`
already implements this (no `+1` shift, hard error on `time = 0`
events).

Display-layer follow-up still pending:
- `print.causal_survival_fit()` / `summary()` time-grid display:
  show ranges as `(0, T_max]` not `[0, T_max)`.
- `print.causal_survival_risk()` / `print.causal_survival_contrast()`:
  same.
- `causal_risk_table()` row labels: `(t_{k-1}, t_k]` for each interval.
- All roxygen examples and the vignette.

Land when those layers get touched in Phase 3.

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
  - **`ipcw` arg semantics + mixed dep / indep censoring sources**.
    Reader confusion to head off: "if some subjects are labeled `ipcw =
    FALSE` (indep), do their at-risk rows still feed the dep-cens
    hazard fit?" Yes, and there is no bias. The cause-specific hazard
    for dep cens is
    `λ_dep(k | X) = lim P(dep cens at k | at risk at k, X) / dt`.
    An indep-labeled subject IS at risk of dep cens at every interval
    until indep censoring wins the race for them — their `dep_cens = 0`
    rows are correct denominator contributions, not noise. This is the
    same logic that lets admin-truncated subjects donate their
    `dep_cens = 0` rows to the c-hazard fit. The standard "independent
    competing causes given X" assumption (already required for IPCW to
    work) is what makes it valid: indep cens cancels in the weight,
    dep cens is corrected via `1 / (1 - h_dep(k | X))`. Worked example:
    two subjects, one dep-censored at k=3, one indep-censored at k=5 —
    trace each subject's rows through the c-hazard glm fit and through
    the IPCW weight `1 / (1 - h_dep(k | X))`.
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
