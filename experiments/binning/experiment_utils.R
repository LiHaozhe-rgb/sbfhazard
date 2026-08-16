# Shared helpers for single-scenario binning experiments.

binning_scenario <- function(family, n_train, d, rho, model) {
  prefix <- switch(
    family,
    additive_pairwise = "additive_pairwise_binning",
    multiplicative_sparse = "multiplicative_sparse_binning",
    family
  )
  data.frame(
    scenario_id = sprintf(
      "%s_n%d_rho%s_model%d",
      prefix,
      as.integer(n_train),
      gsub("\\.", "p", as.character(rho)),
      as.integer(model)
    ),
    n_train = as.integer(n_train),
    dimension = as.integer(d),
    rho = as.numeric(rho),
    model = as.integer(model),
    stringsAsFactors = FALSE
  )
}

binning_spec_fields <- function(spec = NULL) {
  if (is.null(spec)) {
    return(data.frame(
      spec_label = NA_character_,
      time_method = NA_character_,
      covariate_method = NA_character_,
      time_bins = NA_integer_,
      covariate_bins = NA_integer_,
      representative = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    spec_label = spec$label,
    time_method = spec$time_method,
    covariate_method = spec$covariate_method,
    time_bins = as.integer(spec$time_bins),
    covariate_bins = as.integer(spec$covariate_bins),
    representative = spec$representative,
    stringsAsFactors = FALSE
  )
}

binning_bind_rows <- function(rows) {
  if (length(rows) == 0L) {
    return(data.frame())
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

binning_rotated_spec_indices <- function(n_specs, run) {
  n_specs <- as.integer(n_specs)
  run <- as.integer(run)
  if (n_specs <= 1L) {
    return(seq_len(max(n_specs, 0L)))
  }
  shift <- (run - 1L) %% n_specs
  c(seq.int(shift + 1L, n_specs), seq_len(shift))
}

binning_sd_or_na <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) <= 1L) {
    return(NA_real_)
  }
  stats::sd(x)
}

binning_mse_or_na <- function(estimate, truth) {
  estimate <- as.numeric(estimate)
  truth <- as.numeric(truth)
  keep <- is.finite(estimate) & is.finite(truth)
  if (!any(keep)) {
    return(NA_real_)
  }
  mean((estimate[keep] - truth[keep])^2)
}

binning_fit_diagnostics <- function(fit, family) {
  out <- fit$fit_diagnostics[1L, , drop = FALSE]
  metadata <- fit$binning$metadata
  time_rows <- metadata$role == "time"
  covariate_rows <- metadata$role == "covariate"
  out$time_effective_bins <- as.numeric(metadata$effective_bins[time_rows][1L])
  out$time_empty_bins <- as.numeric(metadata$empty_bins[time_rows][1L])
  out$mean_covariate_effective_bins <- sbf_mean_or_na(metadata$effective_bins[covariate_rows])
  out$mean_covariate_empty_bins <- sbf_mean_or_na(metadata$empty_bins[covariate_rows])

  if (identical(family, "additive_pairwise")) {
    out$pairwise_nonzero_share <- as.numeric(fit$binning$pairwise_nonzero_share)
    return(out)
  }

  out$sparse_joint_cells <- as.numeric(fit$binning$sparse_joint_cells)
  out$joint_exposure_nonzero_share <- as.numeric(fit$binning$joint_exposure_nonzero_share)
  out$joint_exposure_memory_mb <- as.numeric(fit$binning$joint_exposure_memory_mb)
  out
}

binning_prediction_diagnostics <- function(diagnostics) {
  min_finite <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      return(NA_real_)
    }
    min(x)
  }
  sum_finite <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      return(NA_real_)
    }
    sum(x)
  }

  data.frame(
    prediction_value_denom_p05 = min_finite(diagnostics$value_denom_p05),
    prediction_value_denom_adjusted_count = sum_finite(diagnostics$value_denom_adjusted_count),
    prediction_slope_denom_p05 = min_finite(diagnostics$slope_denom_p05),
    prediction_slope_denom_adjusted_count = sum_finite(diagnostics$slope_denom_adjusted_count),
    prediction_kernel_support_zero_count = sum_finite(diagnostics$kernel_support_zero_count),
    prediction_bin_step_fallback_count = sum_finite(diagnostics$bin_step_fallback_count),
    stringsAsFactors = FALSE
  )
}

binning_truth_curves <- function(family, target_k, d, z_grid, t_grid, model, violate_cox) {
  z0 <- rep(0, d - 1L)
  if (identical(family, "additive")) {
    return(list(
      effect = sbf_true_additive_phi(z_grid, target_k, d - 1L, violate_cox),
      survival = sbf_true_additive_survival(t_grid, z0, model, violate_cox)
    ))
  }
  list(
    effect = sbf_true_multiplicative_phi(z_grid, target_k, violate_cox),
    survival = sbf_true_multiplicative_survival(t_grid, z0, model, violate_cox)
  )
}

binning_predict_curves <- function(fit,
                                   family,
                                   component,
                                   z_grid,
                                   t_grid,
                                   d) {
  component <- as.integer(component)
  z0_values <- rep(0, as.integer(d) - 1L)
  names(z0_values) <- paste0("V", seq_along(z0_values))
  z0 <- as.data.frame(as.list(z0_values), check.names = FALSE)
  effect_raw <- sbf_predict(fit, type = "component", component = component, xout = z_grid, warn_diagnostics = FALSE)
  baseline_raw <- sbf_predict(fit, type = "component", component = 1L, xout = t_grid, warn_diagnostics = FALSE)
  survival_pred <- sbf_predict(fit, type = "survival", data = z0, times = t_grid, warn_diagnostics = FALSE)

  if (identical(family, "multiplicative")) {
    effect <- log(pmax(as.numeric(effect_raw), 1e-8))
    baseline <- log(pmax(as.numeric(baseline_raw), 1e-8))
  } else {
    effect <- as.numeric(effect_raw)
    baseline <- as.numeric(baseline_raw)
  }

  list(
    effect = effect,
    baseline = baseline,
    survival = as.numeric(survival_pred$survival[1L, ]),
    diagnostics = if (!is.null(fit$binning)) {
      binning_prediction_diagnostics(survival_pred$prediction_diagnostics)
    } else {
      data.frame()
    }
  )
}

binning_curve_rows <- function(scenario,
                               family,
                               run,
                               seed,
                               obj,
                               status,
                               z_grid,
                               t_grid,
                               curves,
                               truth) {
  common <- data.frame(
    scenario_id = scenario$scenario_id,
    n_train = scenario$n_train,
    dimension = scenario$dimension,
    rho = scenario$rho,
    model = scenario$model,
    run = as.integer(run),
    seed = as.integer(seed),
    method = obj$method,
    fit_type = obj$fit_type,
    kernel_name = obj$kernel_name,
    status = status,
    stringsAsFactors = FALSE
  )
  common <- cbind(common, binning_spec_fields(obj$spec))
  make_curve <- function(curve_type, grid, estimate, truth_values) {
    data.frame(
      common[rep(1L, length(grid)), , drop = FALSE],
      curve_level = "run",
      curve_type = curve_type,
      grid_value = as.numeric(grid),
      estimate = as.numeric(estimate),
      truth = as.numeric(truth_values),
      stringsAsFactors = FALSE
    )
  }
  rbind(
    make_curve("effect", z_grid, curves$effect, truth$effect),
    make_curve("baseline", t_grid, curves$baseline, rep(NA_real_, length(t_grid))),
    make_curve("survival", t_grid, curves$survival, truth$survival)
  )
}

binning_group_factor <- function(df, cols) {
  factors <- lapply(df[cols], addNA)
  do.call(interaction, c(factors, list(drop = TRUE, lex.order = TRUE)))
}

binning_curve_summary <- function(curve_rows) {
  if (nrow(curve_rows) == 0L) {
    return(curve_rows)
  }
  summary_source <- curve_rows[curve_rows$status == "converged", , drop = FALSE]
  group_cols <- setdiff(names(curve_rows), c("run", "seed", "curve_level", "estimate", "truth", "valid_runs", "estimate_sd"))
  groups <- binning_group_factor(summary_source, group_cols)
  summary_rows <- lapply(split(summary_source, groups), function(df) {
    row <- df[1L, , drop = FALSE]
    row$run <- NA_integer_
    row$seed <- NA_integer_
    row$curve_level <- "summary"
    row$estimate <- sbf_mean_or_na(df$estimate)
    row$truth <- sbf_mean_or_na(df$truth)
    row$valid_runs <- sum(is.finite(df$estimate))
    row$estimate_sd <- binning_sd_or_na(df$estimate)
    row
  })
  curve_rows$curve_level <- "run"
  curve_rows$valid_runs <- NA_integer_
  curve_rows$estimate_sd <- NA_real_
  binning_bind_rows(c(summary_rows, list(curve_rows)))
}

binning_run_row <- function(scenario,
                            family,
                            run,
                            seed,
                            obj,
                            status,
                            curves,
                            truth_effect,
                            truth_survival,
                            paired_speedup_vs_sbf) {
  fit <- obj$fit
  base <- data.frame(
    scenario_id = scenario$scenario_id,
    n_train = scenario$n_train,
    dimension = scenario$dimension,
    rho = scenario$rho,
    model = scenario$model,
    run = as.integer(run),
    seed = as.integer(seed),
    method = obj$method,
    fit_type = obj$fit_type,
    kernel_name = obj$kernel_name,
    status = status,
    fit_runtime_sec = obj$fit_runtime_sec,
    paired_speedup_vs_sbf = paired_speedup_vs_sbf,
    effect_mse = binning_mse_or_na(curves$effect, truth_effect),
    survival_mse = binning_mse_or_na(curves$survival, truth_survival),
    iterations_used = as.numeric(fit$iterations_used),
    final_delta = as.numeric(fit$final_delta),
    stringsAsFactors = FALSE
  )
  if (identical(family, "multiplicative_sparse")) {
    component_floor_count <- as.numeric(
      fit$fit_diagnostics$component_floor_count[[1L]]
    )
    component_update_count <- as.numeric(fit$iterations_used) *
      sum(lengths(fit$x.grid))
    base$component_floor_count <- component_floor_count
    base$component_floor_rate <- if (
      is.finite(component_update_count) && component_update_count > 0
    ) {
      component_floor_count / component_update_count
    } else {
      NA_real_
    }
  }
  cbind(base, binning_spec_fields(obj$spec))
}

binning_diagnostic_row <- function(scenario,
                                   family,
                                   run,
                                   seed,
                                   obj,
                                   status,
                                   fit_diag,
                                   prediction_diag) {
  base <- data.frame(
    scenario_id = scenario$scenario_id,
    n_train = scenario$n_train,
    dimension = scenario$dimension,
    rho = scenario$rho,
    model = scenario$model,
    run = as.integer(run),
    seed = as.integer(seed),
    method = obj$method,
    fit_type = obj$fit_type,
    kernel_name = obj$kernel_name,
    status = status,
    stringsAsFactors = FALSE
  )
  cbind(base, binning_spec_fields(obj$spec), prediction_diag, fit_diag)
}

binning_group_cols <- function(family) {
  base <- c(
    "scenario_id", "n_train", "dimension", "rho", "model",
    "method", "spec_label", "time_method", "covariate_method",
    "time_bins", "covariate_bins", "representative"
  )
  if (identical(family, "additive_pairwise")) {
    return(c(base[1:6], "fit_type", base[7:length(base)], "kernel_name"))
  }
  c(base, "kernel_name")
}

binning_valid_run_rows <- function(df) {
  df$status == "converged" &
    is.finite(as.numeric(df$effect_mse)) &
    is.finite(as.numeric(df$survival_mse))
}

binning_run_summary <- function(df) {
  valid <- binning_valid_run_rows(df)
  out <- data.frame(
    n_total_runs = nrow(df),
    n_converged = sum(df$status == "converged", na.rm = TRUE),
    n_nonconverged = sum(df$status == "nonconverged", na.rm = TRUE),
    valid_runs = sum(valid, na.rm = TRUE),
    mean_fit_runtime_sec = sbf_mean_or_na(as.numeric(df$fit_runtime_sec)[valid]),
    sd_fit_runtime_sec = binning_sd_or_na(as.numeric(df$fit_runtime_sec)[valid]),
    speedup_vs_sbf = NA_real_,
    mean_effect_mse = sbf_mean_or_na(as.numeric(df$effect_mse)[valid]),
    sd_effect_mse = binning_sd_or_na(as.numeric(df$effect_mse)[valid]),
    mean_survival_mse = sbf_mean_or_na(as.numeric(df$survival_mse)[valid]),
    sd_survival_mse = binning_sd_or_na(as.numeric(df$survival_mse)[valid]),
    mean_iterations = sbf_mean_or_na(as.numeric(df$iterations_used)[valid]),
    mean_final_delta = sbf_mean_or_na(as.numeric(df$final_delta)[valid]),
    stringsAsFactors = FALSE
  )
  if ("component_floor_count" %in% names(df)) {
    out$mean_component_floor_count <- sbf_mean_or_na(
      as.numeric(df$component_floor_count)[valid]
    )
  }
  if ("component_floor_rate" %in% names(df)) {
    out$mean_component_floor_rate <- sbf_mean_or_na(
      as.numeric(df$component_floor_rate)[valid]
    )
  }
  out
}

binning_diagnostic_summary <- function(df, family) {
  if (nrow(df) == 0L) {
    out <- data.frame(
      mean_denom_min = NA_real_,
      mean_denom_p05 = NA_real_,
      mean_denom_adjusted_count = NA_real_,
      mean_kernel_norm_min = NA_real_,
      mean_kernel_norm_p05 = NA_real_,
      mean_kernel_norm_zero_count = NA_real_,
      mean_kernel_norm_adjusted_count = NA_real_,
      mean_prediction_kernel_support_zero_count = NA_real_,
      mean_prediction_bin_step_fallback_count = NA_real_,
      mean_time_effective_bins = NA_real_,
      mean_time_empty_bins = NA_real_,
      mean_covariate_effective_bins = NA_real_,
      mean_covariate_empty_bins = NA_real_,
      stringsAsFactors = FALSE
    )
    if (identical(family, "additive_pairwise")) {
      out$mean_pairwise_nonzero_share <- NA_real_
    } else {
      out$mean_sparse_joint_cells <- NA_real_
      out$mean_joint_exposure_nonzero_share <- NA_real_
      out$mean_joint_exposure_memory_mb <- NA_real_
    }
    return(out)
  }

  out <- data.frame(
    mean_denom_min = sbf_mean_or_na(as.numeric(df$denom_min)),
    mean_denom_p05 = sbf_mean_or_na(as.numeric(df$denom_p05)),
    mean_denom_adjusted_count = sbf_mean_or_na(as.numeric(df$denom_adjusted_count)),
    mean_kernel_norm_min = sbf_mean_or_na(as.numeric(df$kernel_norm_min)),
    mean_kernel_norm_p05 = sbf_mean_or_na(as.numeric(df$kernel_norm_p05)),
    mean_kernel_norm_zero_count = sbf_mean_or_na(as.numeric(df$kernel_norm_zero_count)),
    mean_kernel_norm_adjusted_count = sbf_mean_or_na(as.numeric(df$kernel_norm_adjusted_count)),
    mean_prediction_kernel_support_zero_count = sbf_mean_or_na(as.numeric(df$prediction_kernel_support_zero_count)),
    mean_prediction_bin_step_fallback_count = sbf_mean_or_na(as.numeric(df$prediction_bin_step_fallback_count)),
    mean_time_effective_bins = sbf_mean_or_na(as.numeric(df$time_effective_bins)),
    mean_time_empty_bins = sbf_mean_or_na(as.numeric(df$time_empty_bins)),
    mean_covariate_effective_bins = sbf_mean_or_na(as.numeric(df$mean_covariate_effective_bins)),
    mean_covariate_empty_bins = sbf_mean_or_na(as.numeric(df$mean_covariate_empty_bins)),
    stringsAsFactors = FALSE
  )
  if (identical(family, "additive_pairwise")) {
    out$mean_pairwise_nonzero_share <- sbf_mean_or_na(as.numeric(df$pairwise_nonzero_share))
  } else {
    out$mean_sparse_joint_cells <- sbf_mean_or_na(as.numeric(df$sparse_joint_cells))
    out$mean_joint_exposure_nonzero_share <- sbf_mean_or_na(as.numeric(df$joint_exposure_nonzero_share))
    out$mean_joint_exposure_memory_mb <- sbf_mean_or_na(as.numeric(df$joint_exposure_memory_mb))
  }
  out
}

binning_summary_table <- function(runs, diagnostics, family) {
  if (nrow(runs) == 0L) {
    return(data.frame())
  }

  group_cols <- binning_group_cols(family)
  run_groups <- split(runs, binning_group_factor(runs, group_cols))
  diagnostic_groups <- if (nrow(diagnostics) > 0L) {
    split(diagnostics, binning_group_factor(diagnostics, group_cols))
  } else {
    list()
  }

  rows <- lapply(names(run_groups), function(key) {
    df <- run_groups[[key]]
    valid_runs <- df[binning_valid_run_rows(df), c("run", "seed"), drop = FALSE]
    diag <- diagnostic_groups[[key]]
    if (is.null(diag)) {
      diag <- data.frame()
    } else if (nrow(valid_runs) > 0L) {
      diag <- merge(valid_runs, diag, by = c("run", "seed"), sort = FALSE)
    } else {
      diag <- diag[FALSE, , drop = FALSE]
    }
    cbind(
      df[1L, group_cols, drop = FALSE],
      binning_run_summary(df),
      binning_diagnostic_summary(diag, family)
    )
  })
  out <- binning_bind_rows(rows)
  binning_rows <- which(!is.na(out$spec_label))
  for (i in binning_rows) {
    reference_method <- if (identical(family, "additive_pairwise")) {
      paste0("additive_sbf_", out$fit_type[[i]])
    } else {
      "multiplicative_sbf"
    }
    reference <- out[
      out$scenario_id == out$scenario_id[[i]] &
        out$method == reference_method &
        out$kernel_name == out$kernel_name[[i]],
      ,
      drop = FALSE
    ]
    if (nrow(reference) == 1L &&
        is.finite(reference$mean_fit_runtime_sec) &&
        is.finite(out$mean_fit_runtime_sec[[i]]) &&
        out$mean_fit_runtime_sec[[i]] > 0) {
      out$speedup_vs_sbf[[i]] <-
        reference$mean_fit_runtime_sec / out$mean_fit_runtime_sec[[i]]
    }
  }
  out
}

binning_write_outputs <- function(run_dir, summary, runs, diagnostics, curves) {
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list(
    summary_file = file.path(run_dir, "summary.csv"),
    runs_file = file.path(run_dir, "runs.csv"),
    diagnostics_file = file.path(run_dir, "diagnostics.csv"),
    curves_file = file.path(run_dir, "curves.csv")
  )
  sbf_write_csv(summary, files$summary_file)
  sbf_write_csv(runs, files$runs_file)
  sbf_write_csv(diagnostics, files$diagnostics_file)
  sbf_write_csv(curves, files$curves_file)
  files
}

binning_default_run_dir <- function(family, include_sbf, config, n_rep, base_seed) {
  root <- sbf_results_dir(project.root, "runs", "experiments", "binning")
  method_label <- if (isTRUE(include_sbf)) "with_sbf" else "binning_only"
  folder <- sprintf(
    "%s_%s_d%d_n%d_rep%d_m%d_rho%s_seed%d",
    family,
    method_label,
    config$d,
    config$n_train,
    n_rep,
    config$model,
    gsub("\\.", "p", as.character(config$rho)),
    base_seed
  )
  file.path(root, gsub("[^A-Za-z0-9._-]+", "_", folder))
}

binning_finish_outputs <- function(run_dir,
                                   family,
                                   summary,
                                   runs,
                                   diagnostics,
                                   curves) {
  files <- binning_write_outputs(run_dir, summary, runs, diagnostics, curves)
  source(file.path(project.root, "experiments", "binning", "plot_binning.R"))
  plot_files <- sbf_plot_binning_run(run_dir, family = family)
  list(
    run_dir = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
    summary_file = files$summary_file,
    runs_file = files$runs_file,
    diagnostics_file = files$diagnostics_file,
    curves_file = files$curves_file,
    plot_files = plot_files,
    summary = summary,
    runs = runs,
    diagnostics = diagnostics,
    curves = curves
  )
}
