graphics.off()

source(file.path("experiments", "binning", "additive_pairwise_experiment.R"))

expect_true <- function(x, message) {
  if (!isTRUE(x)) {
    stop(message, call. = FALSE)
  }
}

config <- default_additive_pairwise_experiment_config(
  n_train = 60L,
  d = 3L,
  fit_it = 30L,
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
    representative = "midpoint",
    local_constant = TRUE
  ),
  list(
    label = "smoke_mean",
    time_bins = 6L,
    covariate_bins = 6L,
    time_method = "equal_width",
    covariate_method = "equal_width",
    representative = "mean",
    local_constant = TRUE
  )
)
out_dir <- sbf_results_dir(project.root, "runs", "smoke", "additive_binning")
outputs <- run_additive_pairwise_binning_experiment(
  n_rep = 1L,
  base_seed = 20263000L,
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
expect_true(all(file.exists(files)), "Additive binning did not write every required output.")

required_run_columns <- c(
  "method", "spec_label", "status", "fit_runtime_sec",
  "effect_mse", "survival_mse"
)
required_summary_columns <- c(
  "method", "spec_label", "valid_runs", "mean_fit_runtime_sec",
  "speedup_vs_sbf", "mean_effect_mse", "mean_survival_mse"
)
expect_true(
  all(required_run_columns %in% names(outputs$runs)),
  "Additive runs are missing required columns."
)
expect_true(
  all(required_summary_columns %in% names(outputs$summary)),
  "Additive summary is missing required columns."
)

expected_methods <- c(
  "additive_sbf_lc",
  "additive_sbf_ll",
  "additive_pairwise_binning_lc",
  "additive_pairwise_binning_ll"
)
expected_specs <- vapply(specs, `[[`, character(1), "label")
expect_true(all(expected_methods %in% outputs$runs$method), "Additive fit paths are missing.")
expect_true(
  all(expected_specs %in% stats::na.omit(outputs$runs$spec_label)),
  "Additive binning specifications are missing."
)
expect_true(any(outputs$summary$valid_runs > 0), "Additive binning has no valid run.")
expect_true(any(is.finite(outputs$runs$fit_runtime_sec)), "Additive runtime is not finite.")
expect_true(any(is.finite(outputs$runs$effect_mse)), "Additive effect MSE is not finite.")
expect_true(any(is.finite(outputs$runs$survival_mse)), "Additive survival MSE is not finite.")
expect_true(any(is.finite(outputs$curves$estimate)), "Additive curves are not finite.")
expect_true(
  all(c("effect", "baseline", "survival") %in% outputs$curves$curve_type),
  "Additive curve output is incomplete."
)

sbf_run <- outputs$runs[outputs$runs$method == "additive_sbf_lc", , drop = FALSE][1, ]
binned <- outputs$runs[
  outputs$runs$method == "additive_pairwise_binning_lc",
  ,
  drop = FALSE
][1, ]
summary_input <- rbind(
  sbf_run[rep(1, 2), , drop = FALSE],
  binned[rep(1, 2), , drop = FALSE]
)
summary_input$run <- rep(1:2, 2)
summary_input$status <- rep(c("converged", "nonconverged"), 2)
summary_input$fit_runtime_sec <- c(2, 100, 1, 100)
summary_input$effect_mse <- c(1, 100, 1, 100)
summary_input$survival_mse <- c(1, 100, 1, 100)
summary <- binning_summary_table(summary_input, data.frame(), family = "additive_pairwise")
sbf_summary <- summary[summary$method == "additive_sbf_lc", ]
binned_summary <- summary[summary$method == "additive_pairwise_binning_lc", ]

expect_true(sbf_summary$valid_runs == 1L, "Summary included an invalid SBF run.")
expect_true(binned_summary$valid_runs == 1L, "Summary included an invalid binned run.")
expect_true(binned_summary$speedup_vs_sbf == 2, "Summary speedup is not a ratio of valid mean runtimes.")

cat("Additive binning experiment smoke passed.\n")
