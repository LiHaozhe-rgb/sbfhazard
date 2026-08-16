graphics.off()

# Diagnostic only: equal_width_log is kept internal and is not part of the
# public sbf_fit_binning() API. The output schema intentionally matches the
# main multiplicative binning experiment.

suppressPackageStartupMessages(library(survival))

source(file.path("experiments", "binning", "multiplicative_sparse_experiment.R"))

default_multiplicative_log_time_diagnostic_config <- function(n_train = 1000L,
                                                              d = 5L,
                                                              fit_it = 250L,
                                                              target_k = 1L,
                                                              model = 1L,
                                                              violate_cox = TRUE,
                                                              bandwidth = NULL,
                                                              z_grid = seq(-1, 1, length.out = 300),
                                                              t_grid = seq(0, 2, length.out = 300),
                                                              rho = 0.5) {
  default_multiplicative_sparse_experiment_config(
    n_train = n_train,
    d = d,
    fit_it = fit_it,
    target_k = target_k,
    model = model,
    violate_cox = violate_cox,
    bandwidth = bandwidth,
    z_grid = z_grid,
    t_grid = t_grid,
    rho = rho
  )
}

default_multiplicative_log_time_diagnostic_specs <- function() {
  list(
    list(label = "equal_width_time50_cov30_midpoint",
         time_method = "equal_width", covariate_method = "equal_width",
         time_bins = 50L, covariate_bins = 30L, representative = "midpoint"),
    list(label = "log_equal_width_time50_cov30_midpoint",
         time_method = "equal_width_log", covariate_method = "equal_width",
         time_bins = 50L, covariate_bins = 30L, representative = "midpoint")
  )
}

diagnostic_fit_log_time_obj <- function(data, bandwidth, fit_it, kernel, spec) {
  start <- proc.time()[[3L]]
  prepared <- sbfhazard:::.sbf_prepare_fit_inputs(data = data, formula = NULL)
  fit <- sbfhazard:::.sbf_mult_sparse_fit(
    data = prepared$data,
    bandwidth = bandwidth,
    time_bins = spec$time_bins,
    covariate_bins = spec$covariate_bins,
    time_binning_method = "equal_width_log",
    covariate_binning_method = spec$covariate_method,
    representative = spec$representative,
    formula = prepared$formula,
    iterations = fit_it,
    kernel = kernel,
    convergence_tol = 0.001,
    min_component = 1e-8,
    identification = "sample_mean",
    truth_functions = NULL,
    warn_nonconvergence = FALSE,
    warn_diagnostics = FALSE
  )
  list(
    method = "multiplicative_sparse_binning",
    fit_type = NA_character_,
    spec = spec,
    kernel_name = as.character(fit$fit_settings$kernel_name)[1L],
    fit = fit,
    fit_runtime_sec = as.numeric(proc.time()[[3L]] - start)
  )
}

diagnostic_fit_binning_obj <- function(data, bandwidth, fit_it, kernel, spec) {
  if (identical(spec$time_method, "equal_width_log")) {
    return(diagnostic_fit_log_time_obj(data, bandwidth, fit_it, kernel, spec))
  }
  multiplicative_fit_binning_obj(data, bandwidth, fit_it, kernel, spec)
}

run_multiplicative_log_time_diagnostic <- function(n_rep = 30L,
                                                   base_seed = 2026050104L,
                                                   kernel = "epanechnikov",
                                                   config = default_multiplicative_log_time_diagnostic_config(),
                                                   binning_specs = default_multiplicative_log_time_diagnostic_specs(),
                                                   out_dir = NULL) {
  if (is.null(out_dir)) {
    out_dir <- sbf_results_dir(
      project.root,
      "runs",
      "experiments",
      "diagnostics",
      "multiplicative_log_time_diagnostic"
    )
  }

  scenario <- config$scenario
  component <- config$target_k + 1L
  run_rows <- list()
  diagnostic_rows <- list()
  curve_rows <- list()

  for (run in seq_len(as.integer(n_rep))) {
    seed <- as.integer(base_seed) + run
    data <- sbf_simulate_data(
      n = config$n_train,
      d = config$d,
      rho = config$rho,
      family = "multiplicative",
      model = config$model,
      violate_cox = config$violate_cox,
      seed = seed
    )
    truth <- binning_truth_curves(
      family = "multiplicative",
      target_k = config$target_k,
      d = config$d,
      z_grid = config$z_grid,
      t_grid = config$t_grid,
      model = config$model,
      violate_cox = config$violate_cox
    )

    standard_fit <- multiplicative_fit_sbf_obj(data, config$bandwidth, config$fit_it, kernel)
    fit_order <- binning_rotated_spec_indices(length(binning_specs), run)
    binned_fits <- vector("list", length(binning_specs))
    for (i in fit_order) {
      binned_fits[[i]] <- diagnostic_fit_binning_obj(
        data,
        config$bandwidth,
        config$fit_it,
        kernel,
        binning_specs[[i]]
      )
    }
    fits <- c(list(standard_fit), binned_fits)
    sbf_time <- fits[[1L]]$fit_runtime_sec

    for (obj in fits) {
      curves <- binning_predict_curves(
        fit = obj$fit,
        family = "multiplicative",
        component = component,
        z_grid = config$z_grid,
        t_grid = config$t_grid,
        d = config$d
      )
      status <- if (isTRUE(obj$fit$converged)) "converged" else "nonconverged"
      speedup <- if (!is.null(obj$spec) && is.finite(sbf_time) && obj$fit_runtime_sec > 0) {
        sbf_time / obj$fit_runtime_sec
      } else {
        NA_real_
      }

      run_rows[[length(run_rows) + 1L]] <- binning_run_row(
        scenario = scenario,
        family = "multiplicative_sparse",
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
        diagnostic_rows[[length(diagnostic_rows) + 1L]] <- binning_diagnostic_row(
          scenario = scenario,
          family = "multiplicative_sparse",
          run = run,
          seed = seed,
          obj = obj,
          status = status,
          fit_diag = binning_fit_diagnostics(obj$fit, family = "multiplicative_sparse"),
          prediction_diag = curves$diagnostics
        )
      }
      curve_rows[[length(curve_rows) + 1L]] <- binning_curve_rows(
        scenario = scenario,
        family = "multiplicative_sparse",
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
  summary <- binning_summary_table(runs, diagnostics, family = "multiplicative_sparse")

  binning_finish_outputs(
    run_dir = out_dir,
    family = "multiplicative_sparse",
    summary = summary,
    runs = runs,
    diagnostics = diagnostics,
    curves = curves
  )
}

if (sys.nframe() == 0L) {
  outputs <- run_multiplicative_log_time_diagnostic()
  cat("Multiplicative log-time diagnostic complete.\n")
  cat(sprintf("Summary: %s\n", outputs$summary_file))
  cat(sprintf("Runs: %s\n", outputs$runs_file))
  cat(sprintf("Diagnostics: %s\n", outputs$diagnostics_file))
  cat(sprintf("Curves: %s\n", outputs$curves_file))
}
