graphics.off()

source(file.path("experiments", "diagnostics", "multiplicative_log_time_diagnostic.R"))

expect_true <- function(x, message) {
  if (!isTRUE(x)) {
    stop(message, call. = FALSE)
  }
}

config <- default_multiplicative_log_time_diagnostic_config(
  n_train = 50L,
  d = 3L,
  fit_it = 40L,
  bandwidth = rep(0.35, 3),
  z_grid = seq(-0.8, 0.8, length.out = 25),
  t_grid = seq(0, 1.2, length.out = 25),
  rho = 0.2
)
specs <- list(
  list(
    label = "smoke_raw_equal_width",
    time_bins = 6L,
    covariate_bins = 5L,
    time_method = "equal_width",
    covariate_method = "equal_width",
    representative = "midpoint"
  ),
  list(
    label = "smoke_quantile",
    time_bins = 6L,
    covariate_bins = 5L,
    time_method = "quantile",
    covariate_method = "quantile",
    representative = "midpoint"
  ),
  list(
    label = "smoke_log_equal_width",
    time_bins = 6L,
    covariate_bins = 5L,
    time_method = "equal_width_log",
    covariate_method = "equal_width",
    representative = "midpoint"
  )
)
out_dir <- sbf_results_dir(
  project.root,
  "runs",
  "smoke",
  "multiplicative_log_time_diagnostic"
)
outputs <- run_multiplicative_log_time_diagnostic(
  n_rep = 1L,
  base_seed = 20266000L,
  config = config,
  binning_specs = specs,
  out_dir = out_dir
)

files <- c(
  outputs$summary_file,
  outputs$runs_file,
  outputs$diagnostics_file,
  outputs$curves_file,
  outputs$plot_files
)
expect_true(all(file.exists(files)), "Log-time diagnostic did not write every required output.")
expect_true(
  all(c("multiplicative_sbf", "multiplicative_sparse_binning") %in% outputs$runs$method),
  "Log-time diagnostic is missing a fit path."
)
expect_true(
  all(vapply(specs, `[[`, character(1), "label") %in% stats::na.omit(outputs$runs$spec_label)),
  "Log-time diagnostic is missing a specification."
)
expect_true(any(outputs$summary$valid_runs > 0), "Log-time diagnostic has no valid run.")
expect_true(any(is.finite(outputs$runs$fit_runtime_sec)), "Log-time runtime is not finite.")
expect_true(any(is.finite(outputs$runs$effect_mse)), "Log-time effect MSE is not finite.")
expect_true(any(is.finite(outputs$runs$survival_mse)), "Log-time survival MSE is not finite.")
expect_true(any(is.finite(outputs$curves$estimate)), "Log-time curves are not finite.")
expect_true(
  all(is.finite(outputs$runs$component_floor_rate)) &&
    all(outputs$runs$component_floor_rate >= 0) &&
    all(outputs$runs$component_floor_rate <= 1),
  "Log-time component-floor rates are invalid."
)

effect_rows <- outputs$curves[
  outputs$curves$curve_level == "run" & outputs$curves$curve_type == "effect",
  ,
  drop = FALSE
]
expected_truth <- sbf_true_multiplicative_phi(
  effect_rows$grid_value,
  config$target_k,
  config$violate_cox
)
expect_true(
  max(abs(effect_rows$truth - expected_truth), na.rm = TRUE) < 1e-12,
  "Log-time diagnostic truth changed."
)

cat("Multiplicative log-time reproduction check passed.\n")
