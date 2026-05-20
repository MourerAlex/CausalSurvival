#' Identifying assumptions for a causal_survival_fit
#'
#' Returns the hardcoded baseline identification block for the
#' counterfactual risk \eqn{E[Y_k^{a, c = 0}]} targeted by
#' [causal_survival()]. Each entry is a named list with fields
#' `name`, `statement`, `status`, and `pointer` (spec §4.6). The
#' `pointer` field directs the user to the diagnostic accessor that
#' can supply a necessary-but-not-sufficient check, when one exists.
#'
#' Untestable conditions are reported as `status = "untestable"`;
#' testable ones as `status = "testable"`. No data-driven flag is
#' set inside this accessor — the corresponding diagnostic must be
#' invoked separately.
#'
#' @param fit A `"causal_survival_fit"` object from
#'   [causal_survival()].
#'
#' @return An S3 object of class `"causal_survival_assumptions"`
#'   with one list element per assumption.
#'
#' @seealso [causal_diagnostic()], [causal_survival()]
#' @family accessors
#' @export
causal_assumptions <- function(fit) {
  stopifnot(inherits(fit, "causal_survival_fit"))
  structure(
    list(
      list(
        name      = "Consistency",
        statement = paste0(
          "Y^a = Y when A = a — the observed outcome under the ",
          "received treatment equals the counterfactual outcome ",
          "under that same treatment; no hidden version of the ",
          "treatment. Usually understandable from domain knowledge ",
          "(impact of the surgeon, device brand, dosing, ...)."
        ),
        status    = "untestable",
        pointer   = NA_character_
      ),
      list(
        name      = "Exchangeability",
        statement = paste0(
          "Y^a is independent of A given the recorded baseline ",
          "covariates L_0 — no unmeasured baseline confounding."
        ),
        status    = "untestable",
        pointer   = NA_character_
      ),
      list(
        name      = "Positivity",
        statement = paste0(
          "P(A = a | L_0 = l) > 0 for every level a and every l ",
          "with positive density — both arms are reachable for ",
          "every covariate stratum encountered in the data."
        ),
        status    = "testable",
        pointer   = "causal_diagnostic(fit)$weight_summary"
      ),
      list(
        name      = "No interference",
        statement = paste0(
          "One subject's treatment does not affect another ",
          "subject's potential outcomes. Generally can be assumed ",
          "to be true."
        ),
        status    = "untestable",
        pointer   = NA_character_
      ),
      list(
        name      = "Correct model specification",
        statement = paste0(
          "The fitted discrete-time hazard models for Y and (when ",
          "IPCW is used) for C, and the propensity model for A, ",
          "are correctly specified."
        ),
        status    = "untestable",
        pointer   = "causal_diagnostic(fit)$model_checks"
      ),
      list(
        name      = "Censoring at random (E2)",
        statement = paste0(
          "T^a is independent of the treatment-dependent censoring ",
          "C^d conditional on A and L — the censoring mechanism is ",
          "explainable by the recorded variables."
        ),
        status    = "untestable",
        pointer   = NA_character_
      )
    ),
    class = "causal_survival_assumptions"
  )
}


#' Print a causal_survival_assumptions object
#'
#' Renders the baseline identification block as a numbered list with
#' the status tag and (when applicable) the diagnostic pointer.
#'
#' @param x A `"causal_survival_assumptions"` object.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_survival_assumptions <- function(x, ...) {
  cat("Identifying assumptions (causal_survival_assumptions)\n")
  cat("-----------------------------------------------------\n")
  for (i in seq_along(x)) {
    a <- x[[i]]
    tag <- if (a$status == "testable") "[testable]" else "[untestable]"
    cat(sprintf("%d. %s  %s\n", i, a$name, tag))
    cat("   ", a$statement, "\n", sep = "")
    if (!is.na(a$pointer)) {
      cat("   See: ", a$pointer, "\n", sep = "")
    }
    if (i < length(x)) cat("\n")
  }
  invisible(x)
}
