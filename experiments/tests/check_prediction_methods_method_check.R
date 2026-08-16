graphics.off()

source(file.path("experiments", "method_checks", "prediction_methods.R"))

expect_true <- function(x, message) {
  if (!isTRUE(x)) {
    stop(message, call. = FALSE)
  }
}

config <- default_prediction_methods_method_check_config(
  n_train = 80L,
  d = 3L,
  fit_it = 50L,
  bandwidth = rep(0.35, 3),
  z_grid = seq(-0.8, 0.8, length.out = 30),
  t_grid = seq(0, 1.2, length.out = 30),
  rho = 0.2
)
out_dir <- sbf_results_dir(
  project.root,
  "runs",
  "smoke",
  "prediction_methods_method_check"
)
outputs <- run_prediction_methods_method_check(
  n_rep = 2L,
  base_seed = 2026070301L,
  config = config,
  out_dir = out_dir
)

files <- c(
  outputs$summary_file,
  outputs$runs_file,
  outputs$curves_file,
  outputs$plot_files
)
expect_true(all(file.exists(files)), "Prediction check did not write every required output.")

expected_algorithms <- c(
  "additive_local_linear",
  "additive_local_constant",
  "multiplicative_local_constant"
)
expected_predictions <- c("formula", "interpolation")
expect_true(
  all(expected_algorithms %in% outputs$runs$algorithm),
  "Prediction check is missing an estimator."
)
expect_true(
  all(expected_predictions %in% outputs$runs$prediction),
  "Prediction check is missing a prediction method."
)
expect_true(
  all(c(
    "algorithm", "prediction", "converged", "fit_runtime_sec",
    "effect_mse", "survival_mse"
  ) %in% names(outputs$runs)),
  "Prediction runs are missing required columns."
)
expect_true(any(outputs$runs$converged), "Prediction check has no valid fit.")
expect_true(any(is.finite(outputs$runs$fit_runtime_sec)), "Prediction runtime is not finite.")
expect_true(any(is.finite(outputs$runs$effect_mse)), "Prediction effect MSE is not finite.")
expect_true(any(is.finite(outputs$runs$survival_mse)), "Prediction survival MSE is not finite.")
expect_true(any(is.finite(outputs$curves$estimate)), "Prediction curves are not finite.")
expect_true(
  all(c("effect", "survival") %in% outputs$curves$curve_type),
  "Prediction curve output is incomplete."
)

cat("Prediction methods method-check smoke passed.\n")
