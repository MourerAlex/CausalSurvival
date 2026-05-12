# =============================================================================
# resolve_contrast_pairs() — DEFERRED to v2 (multi-arm treatment)
# =============================================================================
#
# Drafted during Phase 3 step 8a (accessors port from separable_effects) on
# 2026-05-12. Removed from the v0.1.0 R/accessors.R because the package
# currently supports binary point treatment only, and `causal_contrast()`
# can expand the (treatment_a, treatment_b) × {"-", "/"} pairs inline for
# the binary case in ~5 lines.
#
# Reinstate this helper when v2 lands multi-arm treatment (see
# dev/TODO.md "Multi-arm treatment"). At that point `causal_contrast()`
# will need to accept the spec §3.4 `reference` / `contrasts` arguments
# in their general form:
#   - `reference = NULL, contrasts = NULL` → all pairwise
#   - `reference = "<level>"`              → all vs. reference
#   - `contrasts = list(name = list(arms = c(...), op = c("-", "/")))` → custom
#
# The helper below covers exactly that translation: it produces the four
# parallel vectors (name, treatment_a, treatment_b, op) the contrast
# engine iterates over.
#
# Do NOT source this file from the package — it lives outside R/ on
# purpose and is not part of the build.
# =============================================================================

#' Resolve the set of arm pairs to contrast
#'
#' Translate the `reference` / `contrasts` arguments of
#' [causal_contrast()] into a flat list of `(name, treatment_a,
#' treatment_b, op)` quadruples. When `contrasts` is `NULL`, all
#' levels other than `reference` are contrasted against it under both
#' operators (`"-"` and `"/"`).
#'
#' @param levels_vec The fit's `treatment_levels` vector.
#' @param reference A single value in `levels_vec`.
#' @param contrasts `NULL` or a named list as documented in
#'   [causal_contrast()].
#' @return List with elements `name`, `treatment_a`, `treatment_b`,
#'   `op`, each a vector of equal length.
#' @family internal
#' @keywords internal
resolve_contrast_pairs <- function(levels_vec, reference, contrasts) {
  if (is.null(contrasts)) {
    others <- setdiff(levels_vec, reference)
    if (length(others) == 0L) {
      stop("No non-reference treatment levels available for contrast.",
           call. = FALSE)
    }
    name <- as.character(rep(others, each = 2L))
    a    <- rep(others, each = 2L)
    b    <- rep(reference, length(others) * 2L)
    op   <- rep(c("-", "/"), times = length(others))
    return(list(name = paste0(name, "_vs_", reference),
                treatment_a = a, treatment_b = b, op = op))
  }
  # Custom contrasts list — validate shape
  if (!is.list(contrasts) || is.null(names(contrasts)) ||
      any(names(contrasts) == "")) {
    stop("`contrasts` must be a named list. See ?causal_contrast.",
         call. = FALSE)
  }
  flat <- list(name = character(), treatment_a = numeric(),
               treatment_b = numeric(), op = character())
  for (nm in names(contrasts)) {
    spec <- contrasts[[nm]]
    if (!is.list(spec) || !all(c("arms", "op") %in% names(spec)) ||
        length(spec$arms) != 2L) {
      stop(sprintf("`contrasts[['%s']]` must be list(arms = c(a, b), op = c(...)).",
                   nm), call. = FALSE)
    }
    if (!all(spec$arms %in% levels_vec)) {
      stop(sprintf("`contrasts[['%s']]$arms` includes unknown level(s).",
                   nm), call. = FALSE)
    }
    if (!all(spec$op %in% c("-", "/"))) {
      stop(sprintf("`contrasts[['%s']]$op` must contain only '-' or '/'.",
                   nm), call. = FALSE)
    }
    flat$name        <- c(flat$name, rep(nm, length(spec$op)))
    flat$treatment_a <- c(flat$treatment_a, rep(spec$arms[1], length(spec$op)))
    flat$treatment_b <- c(flat$treatment_b, rep(spec$arms[2], length(spec$op)))
    flat$op          <- c(flat$op, spec$op)
  }
  flat
}
