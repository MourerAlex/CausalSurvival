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
  data       = lung_data,
  id         = "id",
  time       = "time",
  status     = "status",
  treatment  = "trt",
  covariates = c("age", "sex"),
  cut_points = 12
)

# Fit g-formula
fit <- causal_survival(pt, method = "gformula")

# One-screen fit summary + identifying assumptions
summary(fit)
causal_assumptions(fit)

# Bootstrap CIs (paired with the fit at accessor time)
boot <- bootstrap(fit, n_boot = 500, seed = 1)

# Risk curves, contrasts, and the at-risk table at each cut time
causal_risk(fit, ci = boot)
causal_contrast(fit, ci = boot)
causal_risk_table(fit, count = "at_risk")

# Single plot: cumulative incidence per arm with 95% CI ribbons and a
# stacked "Number at risk" table aligned to the same x-axis
plot(
  causal_risk(fit, scale = "incidence", ci = boot),
  risk_table = "at_risk"
)
```

See `vignette("getting-started")` once available.

## License

MIT © Alex Mourer
