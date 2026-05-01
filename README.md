# CausalSurvival

> **Status**: pre-implementation. Skeleton only. See `dev/CAUSAL_SURVIVAL_SPEC.md` in the [CausalCompetingRisks](https://github.com/MourerAlex/CausalCompetingRisks) repository for the full specification.

R package for causal inference on single-event survival outcomes, using discrete-time pooled logistic regression. Provides parametric g-formula and inverse probability weighting estimators for cumulative incidence under static, baseline-only treatment regimes, with bootstrap confidence intervals and identifying-assumption accessors.

Designed as the foundation of a two-package ecosystem: [CausalCompetingRisks](https://github.com/MourerAlex/CausalCompetingRisks) extends to competing events via the separable-effects framework.

## Installation

Not yet on CRAN. Development version:

```r
# install.packages("remotes")
remotes::install_github("MourerAlex/CausalSurvival")
```

## Usage

```r
library(CausalSurvival)

# Convert wide subject-level data to person-time
pt <- to_person_time(
  data        = lung_data,
  id          = "id",
  time        = "time",
  event       = "status",
  treatment   = "trt",
  covariates  = c("age", "sex"),
  event_y     = 1,
  event_c     = 0,
  n_intervals = 12
)

# Fit g-formula
fit <- causal_survival(pt, method = "gformula")

# Risk curves
causal_risk(fit)

# Bootstrap CIs
boot <- bootstrap(fit, n_boot = 500)
causal_risk(fit, ci = boot)

# Contrasts
causal_contrast(fit, ci = boot)

# Identifying assumptions
causal_assumptions(fit)
```

See `vignette("getting-started")` once available.

## License

MIT © Alex Mourer
