graphics.off()

suppressPackageStartupMessages(library(survival))

source(file.path("experiments", "script_utils.R"))
project.root <- sbf_project_root()
sbf_load_package(project.root)
source(file.path(project.root, "experiments", "binning", "experiment_utils.R"))

default_additive_pairwise_experiment_config <- function(n_train = 640L,
                                                        d = 5L,
                                                        fit_it = 250L,
                                                        target_k = 1L,
                                                        model = 1L,
                                                        violate_cox = TRUE,
                                                        bandwidth = NULL,
                                                        z_grid = seq(-1, 1, length.out = 300),
                                                        t_grid = seq(0, 2, length.out = 300),
                                                        rho = 0.5) {
  d <- as.integer(d)
  if (is.null(bandwidth)) {
    bandwidth <- c(0.7, rep(0.25, d - 1L))
  }
  list(
    n_train = as.integer(n_train),
    d = d,
    fit_it = as.integer(fit_it),
    target_k = as.integer(target_k),
    model = as.integer(model),
    rho = as.numeric(rho),
    violate_cox = isTRUE(violate_cox),
    bandwidth = sbf_expand_bandwidth(bandwidth, d),
    z_grid = as.numeric(z_grid),
    t_grid = as.numeric(t_grid),
    scenario = binning_scenario("additive_pairwise", n_train, d, rho, model)
  )
}

default_additive_pairwise_binning_specs <- function() {
  list(
    list(label = "equal_width_time50_cov30_midpoint",
         time_method = "equal_width", covariate_method = "equal_width",
         time_bins = 50L, covariate_bins = 30L, representative = "midpoint")
    # list(label = "quantile_time50_cov30_midpoint",
    #      time_method = "quantile", covariate_method = "quantile",
    #      time_bins = 50L, covariate_bins = 30L, representative = "midpoint")
  )
}

additive_fit_sbf_obj <- function(data, bandwidth, fit_it, kernel, fit_type) {
  start <- proc.time()[[3L]]
  fit <- sbf_fit(
    data = data,
    bandwidth = bandwidth,
    family = "additive",
    iterations = fit_it,
    kernel = kernel,
    local_constant = identical(fit_type, "lc"),
    warn_nonconvergence = FALSE
  )
  list(
    method = paste0("additive_sbf_", fit_type),
    fit_type = fit_type,
    spec = NULL,
    kernel_name = as.character(fit$fit_settings$kernel_name)[1L],
    fit = fit,
    fit_runtime_sec = as.numeric(proc.time()[[3L]] - start)
  )
}

additive_fit_binning_obj <- function(data, bandwidth, fit_it, kernel, spec, fit_type) {
  start <- proc.time()[[3L]]
  fit <- sbf_fit_binning(
    data = data,
    bandwidth = bandwidth,
    family = "additive",
    time_bins = spec$time_bins,
    covariate_bins = spec$covariate_bins,
    time_binning_method = spec$time_method,
    covariate_binning_method = spec$covariate_method,
    representative = spec$representative,
    local_constant = identical(fit_type, "lc"),
    iterations = fit_it,
    kernel = kernel,
    warn_nonconvergence = FALSE,
    warn_diagnostics = FALSE
  )
  list(
    method = paste0("additive_pairwise_binning_", fit_type),
    fit_type = fit_type,
    spec = spec,
    kernel_name = as.character(fit$fit_settings$kernel_name)[1L],
    fit = fit,
    fit_runtime_sec = as.numeric(proc.time()[[3L]] - start)
  )
}

run_additive_pairwise_binning_experiment <- function(n_rep = 10L,
                                                     base_seed = 2026050111L,
                                                     kernel = "epanechnikov",
                                                     include_sbf = TRUE,
                                                     experiment_config = NULL,
                                                     binning_specs = NULL,
                                                     out_dir = NULL) {
  .sbf_assert_flag(include_sbf, "include_sbf")
  if (is.null(experiment_config)) {
    experiment_config <- default_additive_pairwise_experiment_config()
  }
  if (is.null(binning_specs)) {
    binning_specs <- default_additive_pairwise_binning_specs()
  }

  config <- experiment_config
  scenario <- config$scenario
  truth <- binning_truth_curves(
    family = "additive",
    target_k = config$target_k,
    d = config$d,
    z_grid = config$z_grid,
    t_grid = config$t_grid,
    model = config$model,
    violate_cox = config$violate_cox
  )

  run_rows <- list()
  diagnostic_rows <- list()
  curve_rows <- list()

  for (run in seq_len(as.integer(n_rep))) {
    seed <- as.integer(base_seed) + run
    data <- sbf_simulate_data(
      n = config$n_train,
      d = config$d,
      rho = config$rho,
      family = "additive",
      model = config$model,
      violate_cox = config$violate_cox,
      seed = seed
    )
    fits <- list()
    if (include_sbf) {
      fits[[length(fits) + 1L]] <- additive_fit_sbf_obj(data, config$bandwidth, config$fit_it, kernel, "lc")
      fits[[length(fits) + 1L]] <- additive_fit_sbf_obj(data, config$bandwidth, config$fit_it, kernel, "ll")
    }
    binned_fits <- vector("list", length(binning_specs))
    for (i in binning_rotated_spec_indices(length(binning_specs), run)) {
      spec <- binning_specs[[i]]
      binned_fits[[i]] <- list(
        additive_fit_binning_obj(data, config$bandwidth, config$fit_it, kernel, spec, "lc"),
        additive_fit_binning_obj(data, config$bandwidth, config$fit_it, kernel, spec, "ll")
      )
    }
    fits <- c(fits, unlist(binned_fits, recursive = FALSE))

    sbf_time <- setNames(rep(NA_real_, 2L), c("lc", "ll"))
    for (obj in fits) {
      if (is.null(obj$spec)) {
        sbf_time[[obj$fit_type]] <- obj$fit_runtime_sec
      }
    }

    for (obj in fits) {
      curves <- binning_predict_curves(
        fit = obj$fit,
        family = "additive",
        component = config$target_k + 1L,
        z_grid = config$z_grid,
        t_grid = config$t_grid,
        d = config$d
      )
      status <- if (isTRUE(obj$fit$converged)) "converged" else "nonconverged"
      speedup <- if (!is.null(obj$spec) && is.finite(sbf_time[[obj$fit_type]]) && obj$fit_runtime_sec > 0) {
        sbf_time[[obj$fit_type]] / obj$fit_runtime_sec
      } else {
        NA_real_
      }

      run_rows[[length(run_rows) + 1L]] <- binning_run_row(
        scenario = scenario,
        family = "additive_pairwise",
        run = run,
        seed = seed,
        obj = obj,
        status = status,
        curves = curves,
        truth_effect = truth$effect,
        truth_survival = truth$survival,
        paired_speedup_vs_sbf = speedup
      )

      if (!is.null(obj$spec)) {
        prediction_diag <- curves$diagnostics
        fit_diag <- binning_fit_diagnostics(obj$fit, family = "additive_pairwise")
        diagnostic_rows[[length(diagnostic_rows) + 1L]] <- binning_diagnostic_row(
          scenario = scenario,
          family = "additive_pairwise",
          run = run,
          seed = seed,
          obj = obj,
          status = status,
          fit_diag = fit_diag,
          prediction_diag = prediction_diag
        )
      }
      curve_rows[[length(curve_rows) + 1L]] <- binning_curve_rows(
        scenario = scenario,
        family = "additive_pairwise",
        run = run,
        seed = seed,
        obj = obj,
        status = status,
        z_grid = config$z_grid,
        t_grid = config$t_grid,
        curves = curves,
        truth = truth
      )
    }
  }

  runs <- binning_bind_rows(run_rows)
  diagnostics <- binning_bind_rows(diagnostic_rows)
  curves <- binning_curve_summary(binning_bind_rows(curve_rows))
  summary <- binning_summary_table(runs, diagnostics, family = "additive_pairwise")
  run_dir <- if (is.null(out_dir)) {
    binning_default_run_dir("additive_pairwise", include_sbf, config, n_rep, base_seed)
  } else {
    out_dir
  }
  binning_finish_outputs(
    run_dir = run_dir,
    family = "additive_pairwise",
    summary = summary,
    runs = runs,
    diagnostics = diagnostics,
    curves = curves
  )
}

if (sys.nframe() == 0L) {
  outputs <- run_additive_pairwise_binning_experiment()
  cat("Additive pairwise binning experiment complete.\n")
  cat(sprintf("Summary: %s\n", outputs$summary_file))
  cat(sprintf("Runs: %s\n", outputs$runs_file))
  cat(sprintf("Diagnostics: %s\n", outputs$diagnostics_file))
  cat(sprintf("Curves: %s\n", outputs$curves_file))
}
