#' Bootstrap confidence bands (stub)
#'
#' Subject-level percentile bootstrap for a `"causal_survival_fit"`.
#' Stub implementation: serial loop, no progress reporting, no
#' parallelism. Polish per spec §3.3 (future.apply + progress
#' cadence + warnings capture) deferred — tracked in `dev/TODO.md`.
#'
#' Each replicate samples unique subject IDs with replacement, stitches
#' the corresponding person-time rows back together with synthetic IDs,
#' and re-evaluates `fit$call` with the resampled `pt_data` (spec
#' line 767). Failed replicates (errors during refit) are tracked in
#' `$failed_reps`; the effective B is `n_boot - length(failed_reps)`.
#'
#' @param fit A `"causal_survival_fit"` object from
#'   [causal_survival()].
#' @param n_boot Positive integer. Number of bootstrap replicates.
#' @param alpha Two-sided significance level in `(0, 1)`.
#' @param seed Optional integer RNG seed for reproducibility.
#'
#' @return An S3 object of class `"causal_survival_bootstrap"` with
#'   the shape documented in spec §4.3.
#'
#' @export
bootstrap <- function(fit, n_boot = 500, alpha = 0.05, seed = NULL) {
  stopifnot(inherits(fit, "causal_survival_fit"))
  if (!is.numeric(n_boot) || length(n_boot) != 1L ||
      n_boot < 1 || n_boot != round(n_boot)) {
    stop("`n_boot` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be in (0, 1).", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  pt_data <- fit$pt_data
  if (is.null(pt_data)) {
    stop("`fit$pt_data` is NULL; refit `causal_survival(..., keep_data = TRUE)`.",
         call. = FALSE)
  }
  id_col  <- fit$id_col
  method  <- fit$method

  unique_ids <- unique(pt_data[[id_col]])
  pt_by_id   <- split(pt_data, pt_data[[id_col]])
  n          <- length(unique_ids)

  reps_long   <- list()
  failed_reps <- integer()

  tic                <- proc.time()[["elapsed"]]
  estimate_announced <- FALSE

  for (b in seq_len(n_boot)) {
    # Progress messages — ported from separable_effects/R/bootstrap.R
    # lines 109-136. First 50 replicates: every 10. After 50: emit a
    # time estimate once, then every 100.
    if (b == 1 || (b <= 50 && b %% 10 == 0)) {
      message("Bootstrap replicate ", b, "/", n_boot)
    }
    if (!estimate_announced && b == 50 && n_boot > 50) {
      elapsed   <- proc.time()[["elapsed"]] - tic
      per_rep   <- elapsed / 50
      remaining <- (n_boot - 50) * per_rep
      fmt_duration <- function(sec) {
        m <- as.integer(floor(sec / 60))
        s <- as.integer(round(sec - 60 * m))
        if (m > 0) sprintf("%d min %02d sec", m, s) else sprintf("%d sec", s)
      }
      message(sprintf(
        "Bootstrap: 50 replicates done in %s. Estimated remaining: %s (%d more replicates).",
        fmt_duration(elapsed), fmt_duration(remaining), n_boot - 50
      ))
      estimate_announced <- TRUE
    }
    if (b > 50 && b %% 100 == 0) {
      message("Bootstrap replicate ", b, "/", n_boot)
    }

    sampled <- sample(unique_ids, size = n, replace = TRUE)
    boot_data <- do.call(rbind, lapply(seq_along(sampled), function(i) {
      rows <- pt_by_id[[as.character(sampled[i])]]
      rows[[id_col]] <- paste0(sampled[i], "_", i)
      rows
    }))
    # Preserve class + attributes lost by rbind
    class(boot_data) <- class(pt_data)
    for (a in names(attributes(pt_data))) {
      if (a %in% c("names", "row.names", "class")) next
      attr(boot_data, a) <- attr(pt_data, a)
    }

    call_b          <- fit$call
    call_b$pt_data  <- boot_data
    res <- tryCatch(suppressWarnings(eval(call_b)),
                    error = function(e) NULL)
    if (is.null(res)) {
      failed_reps <- c(failed_reps, b)
      next
    }
    est <- res$cumulative_incidence[[method]]
    reps_long[[length(reps_long) + 1L]] <- data.frame(
      boot_id   = b,
      treatment = est$treatment,
      k         = est$k,
      value     = est$inc,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  replicates <- if (length(reps_long) == 0L) {
    data.frame(boot_id = integer(), treatment = numeric(),
               k = integer(), value = numeric(),
               stringsAsFactors = FALSE)
  } else do.call(rbind, reps_long)

  # Percentile bands per (treatment, k)
  q_lower <- alpha / 2
  q_upper <- 1 - alpha / 2
  if (nrow(replicates) > 0L) {
    agg <- stats::aggregate(
      value ~ treatment + k, data = replicates,
      FUN = function(v) {
        c(lower = stats::quantile(v, probs = q_lower, na.rm = TRUE),
          upper = stats::quantile(v, probs = q_upper, na.rm = TRUE))
      }
    )
    ci_lower <- data.frame(treatment = agg$treatment, k = agg$k,
                           lower = agg$value[, "lower"])
    ci_upper <- data.frame(treatment = agg$treatment, k = agg$k,
                           upper = agg$value[, "upper"])
  } else {
    ci_lower <- data.frame(treatment = numeric(), k = integer(),
                           lower = numeric())
    ci_upper <- data.frame(treatment = numeric(), k = integer(),
                           upper = numeric())
  }

  structure(
    list(
      fit_call         = fit$call,
      n_boot_requested = as.integer(n_boot),
      n_boot_effective = as.integer(n_boot - length(failed_reps)),
      alpha            = alpha,
      replicates       = replicates,
      ci_lower         = ci_lower,
      ci_upper         = ci_upper,
      failed_reps      = failed_reps,
      warnings_count   = NA_integer_
    ),
    class = "causal_survival_bootstrap"
  )
}
