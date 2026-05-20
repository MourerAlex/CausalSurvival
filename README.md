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

# --- Synthetic data ---------------------------------------------------------
set.seed(42); n <- 200
df <- data.frame(id = seq_len(n), L1 = rnorm(n), L2 = rbinom(n, 1, 0.4))
df$A      <- rbinom(n, 1, plogis(-0.2 + 0.3 * df$L1 + 0.5 * df$L2))
df$lambda <- plogis(-2 + 0.3 * df$A + 0.2 * df$L1 + 0.4 * df$L2)
df$time   <- pmin(rgeom(n, df$lambda) + 1L, 10L)
df$status <- as.integer(df$time < 10L)

pt <- to_person_time(df, id = "id", time = "time", status = "status",
                     treatment = "A", covariates = c("L1", "L2"),
                     cut_points = 20)

# --- g-formula --------------------------------------------------------------
fit_g <- causal_survival(pt, method = "gformula")
print(fit_g)
summary(fit_g)
print(causal_assumptions(fit_g))

# --- IPW, default km estimator ----------------------------------------------
fit_i <- causal_survival(pt, method = "ipw", truncate = c(0.01, 0.99))
print(fit_i)

# --- IPW, msm estimator (arg renamed: .ipw_engine -> .ipw_estimator) --------
fit_m <- causal_survival(pt, method = "ipw", truncate = c(0.01, 0.99),
                         .ipw_estimator = "msm")
print(fit_m)
summary(fit_m, ci = boot_g)

# --- Accessors --------------------------------------------------------------
print(causal_risk(fit_g, "incidence"))
print(causal_contrast(fit_g))  # emits the loud ci = NULL warning

# --- Bootstrap + contrast with CI -------------------------------------------
boot_g <- bootstrap(fit_g, n_boot = 500, alpha = 0.05, seed = 1)
print(boot_g)
print(causal_contrast(fit_g, ci = boot_g))
summary(fit_g, ci = boot_g)

# --- Risk-table accessor ----------------------------------------------------
print(causal_risk_table(fit_g, count = "at_risk"))

# --- Plot -------------------------------------------------------------------
# Each fit needs its OWN bootstrap. The replicates are method-specific —
# bootstrap(fit_g) stores g-formula CIFs, so using boot_g with fit_i's plot
# would silently show g-formula bands on an IPW curve.
# plot.1 bootstrap for each fit
boot_i <- bootstrap(fit_i, n_boot = 500, alpha = 0.05, seed = 1)
boot_m <- bootstrap(fit_m, n_boot = 500, alpha = 0.05, seed = 1)

# plot.2 different fits plots with different options
plot(causal_risk(fit_g, "incidence", ci = boot_g), risk_table = "at_risk")
plot(causal_risk(fit_i, "survival", ci = boot_i), risk_table = "events_y")
plot(causal_risk(fit_m, "survival", ci = boot_m))

# plot.3 contrast plot
plot(causal_contrast(fit_g, ci = boot_g))

# plot.4 fit plot with stacked tables and tables only
plot(causal_risk(fit_g),
     risk_table = c("at_risk", "events_y", "censored"))
plot(causal_risk(fit_g),
     risk_table = c("at_risk", "events_y", "censored"),
     curves     = FALSE)
```

See `vignette("getting-started")` once available.

## License

MIT © Alex Mourer
