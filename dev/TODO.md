# CausalSurvival — TODO / Future scope

Flagged extensions and known limitations. Not in current scope.

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
