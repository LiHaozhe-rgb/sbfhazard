# Plot helpers for the main binning experiment outputs.

binplot_family_spec <- function(family) {
  switch(
    family,
    additive_pairwise = list(
      method_labels = c(
        additive_sbf_lc = "SBF LC",
        additive_sbf_ll = "SBF LL",
        additive_pairwise_binning_lc = "Binned LC",
        additive_pairwise_binning_ll = "Binned LL"
      ),
      sbf_methods = c("additive_sbf_lc", "additive_sbf_ll"),
      ylab = c(effect = "Effect", baseline = "Baseline hazard", survival = "Survival at z=0"),
      ylim = list(effect = c(-2, 2), baseline = c(-0.5, 3), survival = c(0, 1))
    ),
    multiplicative_sparse = list(
      method_labels = c(
        multiplicative_sbf = "SBF",
        multiplicative_sparse_binning = "Binned"
      ),
      sbf_methods = "multiplicative_sbf",
      ylab = c(effect = "Log effect", baseline = "Log baseline hazard", survival = "Survival at z=0"),
      ylim = list(effect = c(-2, 2), baseline = c(-4, 4), survival = c(0, 1))
    ),
    stop("family must be 'additive_pairwise' or 'multiplicative_sparse'.", call. = FALSE)
  )
}

binplot_is_blank <- function(x) {
  is.na(x) | !nzchar(as.character(x))
}

binplot_spec_label <- function(df) {
  out <- rep("", nrow(df))
  has_spec <- !binplot_is_blank(df$spec_label)
  method_code <- c(
    equal_width = "EW",
    quantile = "Q",
    equal_width_log = "log-EW"
  )
  time_code <- unname(method_code[as.character(df$time_method)])
  covariate_code <- unname(method_code[as.character(df$covariate_method)])
  time_bins <- ifelse(is.finite(df$time_bins), as.integer(df$time_bins), "")
  covariate_bins <- ifelse(is.finite(df$covariate_bins), as.integer(df$covariate_bins), "")
  representative <- ifelse(
    df$representative == "midpoint",
    "mid",
    as.character(df$representative)
  )
  out[has_spec] <- sprintf(
    "%s%s/%s%s, %s",
    time_code[has_spec],
    time_bins[has_spec],
    covariate_code[has_spec],
    covariate_bins[has_spec],
    representative[has_spec]
  )
  out
}

binplot_kernel_label <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  ifelse(x == "epanechnikov", "Epanechnikov", ifelse(x == "custom", "Custom", x))
}

binplot_show_kernel <- function(df) {
  kernels <- unique(as.character(df$kernel_name[!is.na(df$kernel_name) & nzchar(as.character(df$kernel_name))]))
  length(kernels) > 1L
}

binplot_method_label <- function(method, spec) {
  labels <- spec$method_labels
  out <- unname(labels[as.character(method)])
  out[is.na(out)] <- as.character(method)[is.na(out)]
  out
}

binplot_series_label <- function(df, spec, show_kernel = FALSE) {
  spec_label <- binplot_spec_label(df)
  base <- binplot_method_label(df$method, spec)
  suffix <- ifelse(nzchar(spec_label), spec_label, "")
  if (isTRUE(show_kernel)) {
    kernel_suffix <- binplot_kernel_label(df$kernel_name)
    suffix <- ifelse(nzchar(suffix), paste(suffix, kernel_suffix, sep = "; "), kernel_suffix)
  }
  ifelse(nzchar(suffix), sprintf("%s (%s)", base, suffix), base)
}

binplot_series_key <- function(df) {
  spec_label <- as.character(df$spec_label)
  spec_label[binplot_is_blank(spec_label)] <- "sbf"
  kernel <- as.character(df$kernel_name)
  kernel[is.na(kernel)] <- ""
  fit_type <- as.character(df$fit_type)
  fit_type[is.na(fit_type)] <- ""
  as.character(interaction(
    data.frame(
      method = df$method,
      spec = spec_label,
      kernel = kernel,
      fit_type = fit_type,
      stringsAsFactors = FALSE
    ),
    drop = TRUE,
    lex.order = TRUE
  ))
}

binplot_series_colors <- function(df, family) {
  df$series <- binplot_series_key(df)
  info <- df[!duplicated(df$series), c("series", "method"), drop = FALSE]
  out <- setNames(rep("#4D4D4D", nrow(info)), info$series)
  if (identical(family, "additive_pairwise")) {
    out[info$method == "additive_sbf_lc"] <- "#1B1B1B"
    out[info$method == "additive_sbf_ll"] <- "#7A7A7A"
    for (method in c("additive_pairwise_binning_lc", "additive_pairwise_binning_ll")) {
      idx <- which(info$method == method)
      if (length(idx) > 0L) {
        palette <- if (identical(method, "additive_pairwise_binning_lc")) {
          c("#9ECAE1", "#3182BD", "#08519C")
        } else {
          c("#FC9272", "#DE2D26", "#99000D")
        }
        out[idx] <- if (length(idx) == 1L) palette[2L] else grDevices::colorRampPalette(palette)(length(idx))
      }
    }
  } else {
    out[info$method == "multiplicative_sbf"] <- "#1B1B1B"
    idx <- which(info$method == "multiplicative_sparse_binning")
    if (length(idx) > 0L) {
      palette <- c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#FF7F00", "#B15928")
      out[idx] <- if (length(idx) <= length(palette)) palette[seq_along(idx)] else grDevices::hcl.colors(length(idx), "Dark 3")
    }
  }
  other <- which(!info$method %in% c(
    "additive_sbf_lc", "additive_sbf_ll", "additive_pairwise_binning_lc", "additive_pairwise_binning_ll",
    "multiplicative_sbf", "multiplicative_sparse_binning"
  ))
  if (length(other) > 0L) {
    out[other] <- grDevices::hcl.colors(length(other), "Dark 3")
  }
  out
}

binplot_series_lty <- function(df, spec) {
  method <- as.character(df$method[1L])
  if (method %in% spec$sbf_methods) {
    return(if (grepl("ll$", method)) 2L else 1L)
  }
  spec_label <- as.character(df$spec_label[1L])
  if (grepl("mean", spec_label, fixed = TRUE)) {
    return(3L)
  }
  if (grepl("quantile", spec_label, fixed = TRUE)) {
    return(2L)
  }
  1L
}

binplot_curve_ylim <- function(df, curve_type, spec) {
  ylim <- spec$ylim[[curve_type]]
  truth <- if (identical(curve_type, "baseline")) numeric() else df$truth[is.finite(df$truth)]
  estimate <- df$estimate[is.finite(df$estimate)]
  vals <- c(truth, estimate)
  if (length(vals) > 0L) {
    rng <- range(vals)
    ylim <- c(min(ylim[1L], rng[1L]), max(ylim[2L], rng[2L]))
  }
  ylim
}

binplot_truth_df <- function(df, curve_type) {
  if (identical(curve_type, "baseline")) {
    return(data.frame(grid_value = numeric(), truth = numeric()))
  }
  truth_df <- df[is.finite(df$grid_value) & is.finite(df$truth), c("grid_value", "truth"), drop = FALSE]
  if (nrow(truth_df) == 0L) {
    return(truth_df)
  }
  truth_df <- stats::aggregate(truth ~ grid_value, data = truth_df, FUN = mean)
  truth_df[order(truth_df$grid_value), , drop = FALSE]
}

binplot_converged <- function(runs_tbl) {
  runs_tbl$status == "converged"
}

binplot_open_png <- function(file, width = 2600, height = 1900) {
  grDevices::png(file, width = width, height = height, res = 220, bg = "white")
  par(
    bg = "white",
    mar = c(5.6, 6.2, 4.3, 1.4),
    cex.axis = 2.25,
    cex.lab = 2.4,
    cex.main = 2.35,
    mgp = c(3.7, 1.1, 0)
  )
}

binplot_draw_curve <- function(curves_tbl,
                               family,
                               curve_type,
                               run_id = NULL,
                               include_run_background = TRUE,
                               show_legend = TRUE,
                               main = NULL) {
  spec <- binplot_family_spec(family)
  df <- curves_tbl[curves_tbl$curve_type == curve_type, , drop = FALSE]
  if (!is.null(run_id)) {
    df <- df[df$curve_level == "run" & df$run == run_id, , drop = FALSE]
  }
  df <- df[df$status == "converged", , drop = FALSE]
  df <- df[is.finite(df$grid_value) & is.finite(df$estimate), , drop = FALSE]
  if (nrow(df) == 0L) {
    plot.new()
    title(main = sprintf("No %s data", curve_type))
    return(invisible(NULL))
  }

  summary_df <- if (is.null(run_id)) df[df$curve_level == "summary", , drop = FALSE] else df
  run_df <- if (is.null(run_id) && isTRUE(include_run_background)) df[df$curve_level == "run", , drop = FALSE] else df[FALSE, , drop = FALSE]
  show_kernel <- binplot_show_kernel(df)
  colors <- binplot_series_colors(df, family)
  xlim <- range(df$grid_value, na.rm = TRUE)
  ylim <- binplot_curve_ylim(df, curve_type, spec)
  xlab <- if (identical(curve_type, "effect")) "Covariate value" else "Time"

  plot(xlim, ylim, type = "n", xlab = xlab, ylab = spec$ylab[[curve_type]], main = main)
  if (nrow(run_df) > 0L) {
    run_df$series <- binplot_series_key(run_df)
    for (series in unique(run_df$series)) {
      dd <- run_df[run_df$series == series, , drop = FALSE]
      for (run in unique(dd$run)) {
        rr <- dd[dd$run == run, , drop = FALSE]
        rr <- rr[order(rr$grid_value), , drop = FALSE]
        lines(
          rr$grid_value,
          rr$estimate,
          col = grDevices::adjustcolor(colors[series], alpha.f = 0.13),
          lwd = 1.1,
          lty = binplot_series_lty(rr[1L, , drop = FALSE], spec)
        )
      }
    }
  }

  truth_df <- binplot_truth_df(df, curve_type)
  legend_labels <- character()
  legend_cols <- character()
  legend_lty <- integer()
  legend_lwd <- numeric()
  if (nrow(truth_df) > 0L) {
    lines(truth_df$grid_value, truth_df$truth, col = "black", lwd = 3.4)
    legend_labels <- c(legend_labels, "Truth")
    legend_cols <- c(legend_cols, "black")
    legend_lty <- c(legend_lty, 1L)
    legend_lwd <- c(legend_lwd, 3.4)
  }

  summary_df$series <- binplot_series_key(summary_df)
  for (series in unique(summary_df$series)) {
    ss <- summary_df[summary_df$series == series, , drop = FALSE]
    ss <- ss[order(ss$grid_value), , drop = FALSE]
    method <- as.character(ss$method[1L])
    lwd <- if (method %in% spec$sbf_methods) 3.0 else 2.6
    lty <- binplot_series_lty(ss[1L, , drop = FALSE], spec)
    lines(ss$grid_value, ss$estimate, col = colors[series], lwd = lwd, lty = lty)
    legend_labels <- c(legend_labels, binplot_series_label(ss[1L, , drop = FALSE], spec, show_kernel))
    legend_cols <- c(legend_cols, colors[series])
    legend_lty <- c(legend_lty, lty)
    legend_lwd <- c(legend_lwd, lwd)
  }

  if (isTRUE(show_legend)) {
    legend_ncol <- if (length(legend_labels) >= 9L) 2L else 1L
    legend_cex <- if (legend_ncol == 2L) 1.55 else 1.8
    legend(
      "topright",
      legend = legend_labels,
      col = legend_cols,
      lty = legend_lty,
      lwd = legend_lwd,
      cex = legend_cex,
      ncol = legend_ncol,
      bty = "n"
    )
  }
  invisible(NULL)
}

binplot_plot_curve_file <- function(curves_tbl,
                                    family,
                                    curve_type,
                                    plot_dir,
                                    show_legend = TRUE) {
  out_file <- file.path(plot_dir, sprintf("%s_curve_overlay.png", curve_type))
  binplot_open_png(out_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  binplot_draw_curve(
    curves_tbl,
    family,
    curve_type,
    show_legend = show_legend,
    main = switch(
      curve_type,
      effect = "Effect Curve",
      baseline = "Baseline Curve",
      survival = "Survival Curve"
    )
  )
  out_file
}

binplot_draw_boxplot <- function(groups, xlab, main) {
  groups <- groups[vapply(groups, function(x) any(is.finite(x)), logical(1))]
  if (length(groups) == 0L) {
    plot.new()
    text(0.5, 0.5, sprintf("No finite %s data", xlab), cex = 1.3)
    return(invisible(NULL))
  }
  groups <- groups[order(vapply(groups, stats::median, numeric(1), na.rm = TRUE), decreasing = TRUE)]
  boxplot(
    groups,
    horizontal = TRUE,
    las = 1,
    xlab = xlab,
    main = main,
    col = "#D9EAF7",
    border = "#4D4D4D",
    cex.axis = 1.32,
    cex.lab = 1.48,
    cex.main = 1.42
  )
  abline(v = pretty(unlist(groups, use.names = FALSE)), col = "#ECECEC", lwd = 1)
}

binplot_plot_mse <- function(runs_tbl, family, plot_dir) {
  out_file <- file.path(plot_dir, "mse_boxplot.png")
  binplot_open_png(out_file, width = 3600, height = 2800)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- par(mfrow = c(2, 1), mar = c(5.2, 15.0, 4.1, 1.5), cex.axis = 1.35, cex.lab = 1.5, cex.main = 1.45)
  on.exit(par(op), add = TRUE)
  plot_df <- runs_tbl[binplot_converged(runs_tbl), , drop = FALSE]
  spec <- binplot_family_spec(family)
  show_kernel <- binplot_show_kernel(plot_df)
  effect_ok <- is.finite(plot_df$effect_mse)
  survival_ok <- is.finite(plot_df$survival_mse)
  effect_groups <- split(plot_df$effect_mse[effect_ok], binplot_series_label(plot_df[effect_ok, , drop = FALSE], spec, show_kernel))
  survival_groups <- split(plot_df$survival_mse[survival_ok], binplot_series_label(plot_df[survival_ok, , drop = FALSE], spec, show_kernel))
  binplot_draw_boxplot(effect_groups, "Effect MSE", "Effect MSE Distribution")
  binplot_draw_boxplot(survival_groups, "Survival MSE", "Survival MSE Distribution")
  out_file
}

binplot_plot_runtime <- function(runs_tbl, family, plot_dir) {
  out_file <- file.path(plot_dir, "runtime_boxplot.png")
  binplot_open_png(out_file, width = 3400, height = 1900)
  on.exit(grDevices::dev.off(), add = TRUE)
  runtime <- runs_tbl$fit_runtime_sec
  keep <- binplot_converged(runs_tbl) & is.finite(runtime)
  plot_df <- runs_tbl[keep, , drop = FALSE]
  spec <- binplot_family_spec(family)
  groups <- split(runtime[keep], binplot_series_label(plot_df, spec, binplot_show_kernel(plot_df)))
  binplot_draw_boxplot(groups, "Fit runtime (sec)", "Fit Runtime Distribution")
  out_file
}

sbf_plot_binning_run <- function(run_dir,
                                 family = c("additive_pairwise", "multiplicative_sparse")) {
  family <- match.arg(family)
  runs_file <- file.path(run_dir, "runs.csv")
  curves_file <- file.path(run_dir, "curves.csv")
  if (!file.exists(runs_file) || !file.exists(curves_file)) {
    stop("run_dir must contain runs.csv and curves.csv.", call. = FALSE)
  }
  runs_tbl <- utils::read.csv(runs_file, stringsAsFactors = FALSE)
  curves_tbl <- utils::read.csv(curves_file, stringsAsFactors = FALSE)
  plot_dir <- file.path(run_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  out <- c(
    effect = binplot_plot_curve_file(
      curves_tbl,
      family,
      "effect",
      plot_dir,
      show_legend = FALSE
    ),
    baseline = binplot_plot_curve_file(curves_tbl, family, "baseline", plot_dir),
    survival = binplot_plot_curve_file(curves_tbl, family, "survival", plot_dir),
    mse = binplot_plot_mse(runs_tbl, family, plot_dir),
    runtime = binplot_plot_runtime(runs_tbl, family, plot_dir)
  )
  out
}
