# ============================================================================
# Plot methods for causal_survival_* S3 classes
# ----------------------------------------------------------------------------
# Adapted from `separable_effects/R/plot.R` (627 lines) for the
# CausalSurvival v0.1.0 binary point-treatment scope. The SE-specific
# four-arm display, decomposition annotations
# (`build_contrast_annotations()`), and the contrast / diagnostic plot
# bodies have been dropped or stubbed per spec §3.4 line 769
# ("placeholders for contrast and diagnostic, message-only in v0.1.0").
# The `risk_table` panel arg is preserved in the signature but errors
# with a deferral message until step 8b lands the underlying
# `causal_risk_table()` accessor.
# ============================================================================

#' Default Okabe-Ito palette for binary treatment arms
#'
#' Returns a named character vector of hex colors keyed by the
#' character representation of the two treatment levels in
#' `fit$treatment_levels`. Used as the per-arm color default by
#' [plot.causal_survival_risk()]. Users override individual entries
#' via the `arm_colors` argument.
#'
#' @param levels_vec Length-2 numeric vector (typically `c(0, 1)`).
#' @return Named character vector of hex strings, length 2.
#' @family internal
#' @keywords internal
default_arm_palette <- function(levels_vec) {
  if (length(levels_vec) != 2L) {
    stop("default_arm_palette requires exactly two levels (v0.1.0 ",
         "binary scope).", call. = FALSE)
  }
  # Lower level = blue (#0072B2 control); higher = vermillion
  # (#D55E00 treated). Both are colorblind-safe Okabe-Ito picks.
  out <- c("#0072B2", "#D55E00")
  names(out) <- as.character(c(min(levels_vec), max(levels_vec)))
  out
}


#' Plot cumulative incidence / survival curves
#'
#' Plots the counterfactual curve under each arm of the binary
#' treatment as a step function on the reporting grid `fit$cut_times`,
#' with an optional bootstrap CI ribbon (step-transformed to match the
#' curve). Pairs with [causal_risk()] — pass a `"causal_survival_risk"`
#' object.
#'
#' @param x A `"causal_survival_risk"` object from [causal_risk()].
#' @param arms Numeric vector. Subset of `fit$treatment_levels` to
#'   draw. `NULL` (default) draws all levels.
#' @param arm_colors Named character vector overriding arm hex
#'   colors. Names must be the character representation of treatment
#'   levels. Defaults to Okabe-Ito via [default_arm_palette()]; only
#'   the arms you pass are overridden.
#' @param arm_labels Named character vector overriding the legend
#'   labels for each arm. Defaults to `c("0" = "Control", "1" =
#'   "Treated")` when the levels are `c(0, 1)`; otherwise the level
#'   value is used directly.
#' @param risk_table NULL (default) or one of `"at_risk"`, `"events_y"`,
#'   `"censored"`. Deferred to v0.1.0 polish — currently errors.
#' @param risk_table_height Numeric. Height of the risk table panel
#'   relative to the main plot (which is 1). Default `0.23`.
#' @param title,subtitle Character or NULL. Plot title / subtitle.
#' @param x_label,y_label Axis labels. Defaults reflect the selected
#'   `scale` (incidence vs survival).
#' @param base_size Numeric. Base font size for `theme_minimal()`.
#' @param linewidth Numeric. Width of the step lines. Default 0.8.
#' @param ribbon_alpha Numeric in `[0, 1]`. Transparency of the CI
#'   ribbons. Default 0.15.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object.
#' @family plot
#' @export
plot.causal_survival_risk <- function(x,
                                      arms = NULL,
                                      arm_colors = NULL,
                                      arm_labels = NULL,
                                      risk_table = NULL,
                                      risk_table_height = 0.23,
                                      title = NULL,
                                      subtitle = NULL,
                                      x_label = "Time",
                                      y_label = NULL,
                                      base_size = 11,
                                      linewidth = 0.8,
                                      ribbon_alpha = 0.15,
                                      ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot.causal_survival_risk requires the 'ggplot2' package. ",
         "Install via install.packages('ggplot2').", call. = FALSE)
  }
  if (!is.null(risk_table)) {
    stop("`risk_table` panel is deferred to v0.1.0 polish ",
         "(depends on causal_risk_table(), step 8b not yet ported). ",
         "Track in dev/TODO.md.", call. = FALSE)
  }

  # --- Slice ---
  risk_df    <- x$risk
  have_bands <- any(!is.na(risk_df$lower)) && any(!is.na(risk_df$upper))

  # --- Resolve arms ---
  avail_arms <- sort(unique(risk_df$treatment))
  if (is.null(arms)) arms <- avail_arms
  bad <- setdiff(arms, avail_arms)
  if (length(bad) > 0L) {
    stop("Unknown arm(s): ", paste(bad, collapse = ", "),
         ". Available: ", paste(avail_arms, collapse = ", "),
         call. = FALSE)
  }

  # --- Palette + labels (Okabe-Ito default, user-overridable per-arm) ---
  default_colors <- default_arm_palette(avail_arms)
  default_labels <- stats::setNames(as.character(avail_arms),
                                    as.character(avail_arms))
  if (identical(sort(avail_arms), c(0, 1))) {
    default_labels[["0"]] <- "Control"
    default_labels[["1"]] <- "Treated"
  }
  if (!is.null(arm_colors)) {
    bad_c <- setdiff(names(arm_colors), names(default_colors))
    if (length(bad_c) > 0L) {
      stop("arm_colors has unknown entries: ",
           paste(bad_c, collapse = ", "), call. = FALSE)
    }
    default_colors[names(arm_colors)] <- arm_colors
  }
  if (!is.null(arm_labels)) {
    bad_l <- setdiff(names(arm_labels), names(default_labels))
    if (length(bad_l) > 0L) {
      stop("arm_labels has unknown entries: ",
           paste(bad_l, collapse = ", "), call. = FALSE)
    }
    default_labels[names(arm_labels)] <- arm_labels
  }
  arm_colors_v <- default_colors
  arm_labels_v <- default_labels

  arms_chr <- as.character(arms)

  # --- Shared x-axis ticks (k = 0 origin + cut_times) ---
  k_grid          <- sort(unique(risk_df$time))
  cut_times_full  <- c(0, k_grid)
  shared_x_limits <- c(0, max(cut_times_full))
  pretty_inner    <- pretty(shared_x_limits, n = 5)
  pretty_inner    <- pretty_inner[
    pretty_inner > 0 & pretty_inner < shared_x_limits[2]
  ]
  shared_x_ticks  <- sort(unique(c(shared_x_limits, pretty_inner)))

  # --- Plot data: origin row per arm so curves start at (t = 0, value = 0)
  # for incidence or (t = 0, value = 1) for survival. ---
  origin_value <- if (x$scale == "incidence") 0 else 1
  body_rows <- risk_df[risk_df$treatment %in% arms,
                       c("time", "treatment", "value"), drop = FALSE]
  origin_rows <- data.frame(
    time      = 0,
    treatment = arms,
    value     = origin_value,
    stringsAsFactors = FALSE
  )
  plot_data <- rbind(origin_rows, body_rows)
  plot_data$arm_chr   <- as.character(plot_data$treatment)
  plot_data$arm_label <- arm_labels_v[plot_data$arm_chr]

  # --- Default y-label per scale ---
  if (is.null(y_label)) {
    y_label <- if (x$scale == "incidence")
      "Cumulative incidence" else "Survival"
  }

  # --- Base plot ---
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$time, y = .data$value,
    color = .data$arm_label, group = .data$arm_label
  )) +
    ggplot2::geom_step(linewidth = linewidth) +
    ggplot2::scale_color_manual(
      values = stats::setNames(arm_colors_v[arms_chr],
                               arm_labels_v[arms_chr])
    ) +
    ggplot2::scale_x_continuous(
      breaks = shared_x_ticks,
      limits = shared_x_limits,
      expand = ggplot2::expansion(mult = 0.02)
    ) +
    ggplot2::labs(
      x        = x_label,
      y        = y_label,
      color    = "Arm",
      title    = title %||%
        paste0("Counterfactual ", x$scale,
               " — causal_survival fit"),
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(face = "plain"),
      axis.title.x  = ggplot2::element_text(face = "bold"),
      axis.title.y  = ggplot2::element_text(face = "bold"),
      legend.title  = ggplot2::element_text(face = "bold")
    )

  # --- CI ribbons (step-transformed to match geom_step) ---
  if (have_bands) {
    body_bands <- risk_df[risk_df$treatment %in% arms,
                          c("time", "treatment", "lower", "upper"),
                          drop = FALSE]
    origin_bands <- data.frame(
      time      = 0,
      treatment = arms,
      lower     = origin_value,
      upper     = origin_value,
      stringsAsFactors = FALSE
    )
    ribbon_data <- rbind(origin_bands, body_bands)

    # Step-transform per arm: duplicate each interior point so the ribbon
    # stays flat between cut times and jumps at each cut. Matches the
    # discrete-time hazard step structure.
    ribbon_data <- do.call(rbind, lapply(arms, function(a) {
      d <- ribbon_data[ribbon_data$treatment == a, , drop = FALSE]
      d <- d[order(d$time), , drop = FALSE]
      n <- nrow(d)
      if (n < 2L) return(d)
      time_step  <- c(d$time[1],
                      rep(d$time[2:n], each = 2L))
      lower_step <- c(rep(d$lower[1:(n - 1L)], each = 2L), d$lower[n])
      upper_step <- c(rep(d$upper[1:(n - 1L)], each = 2L), d$upper[n])
      data.frame(
        time      = time_step,
        lower     = lower_step,
        upper     = upper_step,
        treatment = a,
        stringsAsFactors = FALSE
      )
    }))
    ribbon_data$arm_chr   <- as.character(ribbon_data$treatment)
    ribbon_data$arm_label <- arm_labels_v[ribbon_data$arm_chr]

    p <- p + ggplot2::geom_ribbon(
      data = ribbon_data,
      ggplot2::aes(
        x = .data$time, ymin = .data$lower, ymax = .data$upper,
        fill = .data$arm_label
      ),
      alpha = ribbon_alpha, inherit.aes = FALSE
    ) +
      ggplot2::scale_fill_manual(
        values = stats::setNames(arm_colors_v[arms_chr],
                                 arm_labels_v[arms_chr]),
        guide  = "none"
      )
  }

  p
}


#' Plot a causal_survival_contrast object (placeholder)
#'
#' Per spec §3.4 line 769: ships as a message-only placeholder in
#' v0.1.0. Implementation deferred to v0.2.
#'
#' @param x A `"causal_survival_contrast"` object.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`. Emits an informational `message()`.
#' @export
plot.causal_survival_contrast <- function(x, ...) {
  message(
    "plot.causal_survival_contrast() is a placeholder in v0.1.0. ",
    "The contrast table is rendered by print(x); a graphical view ",
    "ships in v0.2 (see dev/TODO.md)."
  )
  invisible(x)
}


#' Plot a causal_survival_diagnostic object (placeholder)
#'
#' Per spec §3.4 line 769: ships as a message-only placeholder in
#' v0.1.0. Implementation deferred to v0.2.
#'
#' @param x A `"causal_survival_diagnostic"` object.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`. Emits an informational `message()`.
#' @export
plot.causal_survival_diagnostic <- function(x, ...) {
  message(
    "plot.causal_survival_diagnostic() is a placeholder in v0.1.0. ",
    "Inspect `x$model_checks` and `x$weight_summary` directly; a ",
    "graphical view ships in v0.2 (see dev/TODO.md)."
  )
  invisible(x)
}
