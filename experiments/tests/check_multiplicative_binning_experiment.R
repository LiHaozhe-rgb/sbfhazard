graphics.off()

source(file.path("experiments", "binning", "multiplicative_sparse_experiment.R"))

expect_true <- function(x, message) {
  if (!isTRUE(x)) {
    stop(message, call. = FALSE)
  }
}

config <- default_multiplicative_sparse_experiment_config(
  n_train = 60L,
  d = 3L,
  fit_it = 40L,
  bandwidth = rep(0.35, 3),
  z_grid = seq(-0.8, 0.8, length.out = 30),
  t_grid = seq(0, 1.2, length.out = 30),
  rho = 0.2
)
specs <- list(
  list(
    label = "smoke_midpoint",
    time_bins = 6L,
    covariate_bins = 6L,
    time_method = "equal_width",
    covariate_method = "equal_width",
    representative = "midpoint"
  ),
  list(
    label = "smoke_mean",
    time_bins = 6L,
    covariate_bins = 6L,
    time_method = "equal_width",
    covariate_method = "equal_width",
    representative = "mean"
  )
)
out_dir <- sbf_results_dir(project.root, "runs", "smoke", "multiplicative_binning")
outputs <- run_multiplicative_sparse_binning_experiment(
  n_rep = 1L,
  base_seed = 20264000L,
  experiment_config = config,
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
expect_true(all(file.exists(files)), "Multiplicative binning did not write every required output.")

required_run_columns <- c(
  "method", "spec_label", "status", "fit_runtime_sec",
  "effect_mse", "survival_mse", "component_floor_count",
  "component_floor_rate"
)
required_summary_columns <- c(
  "method", "spec_label", "valid_runs", "mean_fit_runtime_sec",
  "speedup_vs_sbf", "mean_effect_mse", "mean_survival_mse",
  "mean_component_floor_rate"
)
expect_true(
  all(required_run_columns %in% names(outputs$runs)),
  "Multiplicative runs are missing required columns."
)
expect_true(
  all(required_summary_columns %in% names(outputs$summary)),
  "Multiplicative summary is missing required columns."
)

expect_true(
  all(c("multiplicative_sbf", "multiplicative_sparse_binning") %in% outputs$runs$method),
  "Multiplicative fit paths are missing."
)
expect_true(
  all(vapply(specs, `[[`, character(1), "label") %in% stats::na.omit(outputs$runs$spec_label)),
  "Multiplicative binning specifications are missing."
)
expect_true(any(outputs$summary$valid_runs > 0), "Multiplicative binning has no valid run.")
expect_true(any(is.finite(outputs$runs$fit_runtime_sec)), "Multiplicative runtime is not finite.")
expect_true(any(is.finite(outputs$runs$effect_mse)), "Multiplicative effect MSE is not finite.")
expect_true(any(is.finite(outputs$runs$survival_mse)), "Multiplicative survival MSE is not finite.")
expect_true(any(is.finite(outputs$curves$estimate)), "Multiplicative curves are not finite.")
expect_true(
  all(is.finite(outputs$runs$component_floor_rate)) &&
    all(outputs$runs$component_floor_rate >= 0) &&
    all(outputs$runs$component_floor_rate <= 1),
  "Component-floor rates must lie in [0, 1]."
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
  "Multiplicative effect truth changed."
)

floor_fit <- list(
  iterations_used = 2L,
  final_delta = 0,
  x.grid = list(1:3, 1:2),
  fit_diagnostics = list(component_floor_count = 2L)
)
floor_row <- binning_run_row(
  scenario = config$scenario,
  family = "multiplicative_sparse",
  run = 1L,
  seed = 1L,
  obj = list(
    method = "multiplicative_sbf",
    fit_type = NA_character_,
    spec = NULL,
    kernel_name = "epanechnikov",
    fit_runtime_sec = 0,
    fit = floor_fit
  ),
  status = "converged",
  curves = list(effect = 0, survival = 1),
  truth_effect = 0,
  truth_survival = 1,
  paired_speedup_vs_sbf = NA_real_
)
expect_true(
  abs(floor_row$component_floor_rate - 0.2) < 1e-12,
  "Component-floor rate does not use iterations and stored grid size."
)

cat("Multiplicative binning experiment smoke passed.\n")
