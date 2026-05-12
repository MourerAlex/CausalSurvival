# dev/smoke_run.R — end-to-end usage on synthetic data, written as a
# normal user would write it. Top-level expressions are wrapped in
# print() / summary() because Rscript does not auto-print at the top
# level (unlike the interactive REPL).

suppressPackageStartupMessages({
  library(stats)
  library(survival)
  library(ggplot2)
})

pkg_root <- "/home/moureralex/Bureau/cowork/CausalSurvival/.claude/worktrees/vibrant-mirzakhani"
for (f in c("utils.R", "validate.R", "hazards.R", "weights.R",
            "propensity.R", "data_prep.R", "causal_survival.R",
            "accessors.R", "assumptions.R", "risk_table.R",
            "bootstrap.R", "print.R", "plot.R")) {
  source(file.path(pkg_root, "R", f))
}

# --- Synthetic data ---------------------------------------------------------
set.seed(42)
n <- 200
df <- data.frame(
  id = seq_len(n),
  L1 = rnorm(n),
  L2 = rbinom(n, 1, 0.4)
)
df$A      <- rbinom(n, 1, plogis(-0.2 + 0.3 * df$L1 + 0.5 * df$L2))
df$lambda <- plogis(-2 + 0.3 * df$A + 0.2 * df$L1 + 0.4 * df$L2)
df$time   <- pmin(rgeom(n, df$lambda) + 1L, 10L)
df$status <- as.integer(df$time < 10L)

pt <- to_person_time(
  df,
  id         = "id",
  time       = "time",
  status     = "status",
  treatment  = "A",
  covariates = c("L1", "L2"),
  cut_points = 5
)

# --- g-formula --------------------------------------------------------------
fit_g <- causal_survival(pt, method = "gformula")
print(fit_g)
summary(fit_g)

# --- IPW, default engine (km) -----------------------------------------------
fit_ipw <- causal_survival(pt, method = "ipw", truncate = c(0.01, 0.99))
print(fit_ipw)
summary(fit_ipw)

# --- IPW, msm engine --------------------------------------------------------
fit_msm <- causal_survival(pt, method = "ipw",
                           truncate = c(0.01, 0.99),
                           .ipw_engine = "msm")
print(fit_msm)

# --- Accessors --------------------------------------------------------------
print(causal_risk(fit_g, scale = "incidence"))
print(causal_risk(fit_g, scale = "survival"))

# Contrast without CI — emits the loud "interpretation discouraged" warning.
print(causal_contrast(fit_g))

# --- Bootstrap + contrast with CI ------------------------------------------
boot <- bootstrap(fit_g, n_boot = 20, alpha = 0.10, seed = 1)
print(boot)
print(causal_contrast(fit_g, ci = boot))
summary(fit_g, ci = boot)

# --- Assumptions accessor ---------------------------------------------------
print(causal_assumptions(fit_g))

# --- Risk-table accessor ----------------------------------------------------
print(causal_risk_table(fit_g, count = "at_risk"))
print(causal_risk_table(fit_g, count = "events_y"))
print(causal_risk_table(fit_g, count = "censored"))

# --- Plots (interactive — top-level renders to active device) ---------------
print(plot(causal_risk(fit_g, scale = "incidence", ci = boot)))
print(plot(causal_risk(fit_g, scale = "survival",  ci = boot)))
print(plot(causal_risk(fit_g, scale = "incidence", ci = boot),
           risk_table = "at_risk"))
