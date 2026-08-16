graphics.off()

source(file.path("experiments", "diagnostics", "binning_layout_diagnostics.R"))

expect_true <- function(x, message) {
  if (!isTRUE(x)) {
    stop(message, call. = FALSE)
  }
}

config <- default_binning_layout_diagnostics_config(
  n_train = 50L,
  d = 3L,
  fit_it = 5L,
  rho = 0.2,
  additive_bandwidth = rep(0.35, 3),
  multiplicative_bandwidth = rep(0.35, 3)
)
specs <- list(
  additive = list(list(
    label = "smoke_additive",
    time_bins = 6L,
    covariate_bins = 5L,
    time_method = "equal_width",
    covariate_method = "equal_width",
    representative = "midpoint",
    local_constant = TRUE
  )),
  multiplicative = list(list(
    label = "smoke_multiplicative",
    time_bins = 6L,
    covariate_bins = 5L,
    time_method = "equal_width",
    covariate_method = "equal_width",
    representative = "midpoint"
  ))
)
out_dir <- sbf_results_dir(
  project.root,
  "runs",
  "smoke",
  "binning_layout_diagnostics"
)
outputs <- run_binning_layout_diagnostics(
  base_seed = 20265000L,
  config = config,
  specs = specs,
  out_dir = out_dir
)

files <- c(
  outputs$paths$data_summary_csv,
  outputs$paths$metadata_csv,
  outputs$paths$fit_summary_csv,
  outputs$paths$plots
)
expect_true(all(file.exists(files)), "Layout diagnostic did not write every required output.")
expect_true(
  setequal(outputs$fit_summary$family, c("additive", "multiplicative")),
  "Layout diagnostic is missing a family."
)
expect_true(all(outputs$fit_summary$status == "ok"), "A layout diagnostic fit failed.")
expect_true(all(is.finite(outputs$fit_summary$runtime_sec)), "Layout runtimes are not finite.")
expect_true(all(outputs$metadata$effective_bins > 0), "Stored bin counts are invalid.")
expect_true(
  all(outputs$metadata$empty_bins >= 0) &&
    all(outputs$metadata$empty_bins <= outputs$metadata$effective_bins),
  "Empty bin counts are invalid."
)

cat("Binning layout reproduction check passed.\n")
