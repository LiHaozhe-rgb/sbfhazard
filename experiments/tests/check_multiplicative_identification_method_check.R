graphics.off()

source(file.path("experiments", "method_checks", "multiplicative_identification.R"))

expect_true <- function(x, message) {
  if (!isTRUE(x)) {
    stop(message, call. = FALSE)
  }
}

config <- default_multiplicative_identification_method_check_config(
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
  "multiplicative_identification_method_check"
)
outputs <- run_multiplicative_identification_method_check(
  n_rep = 2L,
  base_seed = 2026070201L,
  config = config,
  out_dir = out_dir
)

files <- c(
  outputs$summary_file,
  outputs$runs_file,
  outputs$curves_file,
  outputs$plot_file
)
expect_true(all(file.exists(files)), "Identification check did not write every required output.")

expected_rules <- c("sample_mean", "integral", "origin", "jasa")
expect_true(
  all(expected_rules %in% outputs$runs$identification),
  "Identification check is missing a rule."
)
expect_true(
  all(c(
    "identification", "converged", "fit_runtime_sec", "effect_mse",
    "survival_mse", "hazard_max_abs_vs_sample_mean",
    "survival_max_abs_vs_sample_mean"
  ) %in% names(outputs$runs)),
  "Identification runs are missing required columns."
)
expect_true(any(outputs$runs$converged), "Identification check has no valid fit.")
expect_true(any(is.finite(outputs$runs$fit_runtime_sec)), "Identification runtime is not finite.")
expect_true(any(is.finite(outputs$runs$effect_mse)), "Identification effect MSE is not finite.")
expect_true(any(is.finite(outputs$runs$survival_mse)), "Identification survival MSE is not finite.")
expect_true(any(is.finite(outputs$curves$estimate)), "Identification curves are not finite.")

sample_mean <- outputs$runs[
  outputs$runs$identification == "sample_mean",
  ,
  drop = FALSE
]
expect_true(
  max(abs(sample_mean$hazard_max_abs_vs_sample_mean), na.rm = TRUE) < 1e-12,
  "Sample-mean hazard self-difference is not zero."
)
expect_true(
  max(abs(sample_mean$survival_max_abs_vs_sample_mean), na.rm = TRUE) < 1e-12,
  "Sample-mean survival self-difference is not zero."
)

cat("Multiplicative identification method-check smoke passed.\n")
