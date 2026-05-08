# dev/smoke_run.R — Phase 3 end-to-end sanity check.
# Sources R/ files in dep order; runs both methods on synthetic data.

suppressPackageStartupMessages({
  library(stats)
  library(survival)
})

pkg_root <- "/home/moureralex/Bureau/cowork/CausalSurvival"
for (f in c("utils.R", "validate.R", "hazards.R", "weights.R",
            "propensity.R", "data_prep.R", "causal_survival.R")) {
  source(file.path(pkg_root, "R", f))
}

set.seed(42)
n <- 200
df <- data.frame(
  id     = seq_len(n),
  L1     = rnorm(n),
  L2     = rbinom(n, 1, 0.4)
)
df$A      <- rbinom(n, 1, plogis(-0.2 + 0.3 * df$L1 + 0.5 * df$L2))
df$lambda <- plogis(-2 + 0.3 * df$A + 0.2 * df$L1 + 0.4 * df$L2)
df$time   <- pmin(rgeom(n, df$lambda), 10L)            # discrete event time
df$status <- as.integer(df$time < 10L)                 # admin censor at t=10

pt <- to_person_time(
  df, id = "id", time = "time", status = "status",
  treatment = "A", covariates = c("L1", "L2"), cut_points = 5
)
cat("\n--- pt_data shape:", nrow(pt), "rows;", ncol(pt), "cols ---\n")
print(head(pt, 3))

cat("\n========== gformula ==========\n")
fit_g <- causal_survival(pt, method = "gformula")
cat("class:", class(fit_g), "\n")
cat("CI head:\n"); print(head(fit_g$cumulative_incidence$gformula))
cat("warnings collected:", length(fit_g$warnings), "\n")

cat("\n========== ipw ==========\n")
fit_i <- causal_survival(pt, method = "ipw", truncate = c(0.01, 0.99))
cat("class:", class(fit_i), "\n")
cat("CI head:\n"); print(head(fit_i$cumulative_incidence$ipw))
cat("weight summary:\n"); print(fit_i$weights$weight_summary)
cat("warnings collected:", length(fit_i$warnings), "\n")
