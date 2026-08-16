graphics.off()

source(file.path("experiments", "binning", "additive_pairwise_experiment.R"))
source(file.path("experiments", "binning", "multiplicative_sparse_experiment.R"))
source(file.path("experiments", "diagnostics", "multiplicative_log_time_diagnostic.R"))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 0L) "all" else args[[1L]]
allowed_modes <- c("main", "larger", "diagnostics", "all")
if (!mode %in% allowed_modes) {
  stop(sprintf("mode must be one of: %s", paste(allowed_modes, collapse = ", ")), call. = FALSE)
}

output_root <- sbf_results_dir(
  project.root,
  "runs",
  "experiments",
  "thesis_binning"
)

thesis_output_dir <- function(...) {
  path <- file.path(output_root, ...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

run_thesis_job <- function(label, expr) {
  cat(sprintf("\n== %s ==\n", label))
  start <- proc.time()[[3L]]
  out <- force(expr)
  cat(sprintf("Done in %.1f sec\n", proc.time()[[3L]] - start))
  out
}

additive_main_spec <- list(
  list(
    label = "equal_width_time30_cov30_midpoint",
    time_method = "equal_width",
    covariate_method = "equal_width",
    time_bins = 30L,
    covariate_bins = 30L,
    representative = "midpoint"
  )
)

multiplicative_main_spec <- list(
  list(
    label = "quantile_time50_cov30_midpoint",
    time_method = "quantile",
    covariate_method = "quantile",
    time_bins = 50L,
    covariate_bins = 30L,
    representative = "midpoint"
  )
)

additive_strategy_specs <- list(
  list(label = "equal_width_time50_cov30_midpoint",
       time_method = "equal_width", covariate_method = "equal_width",
       time_bins = 50L, covariate_bins = 30L, representative = "midpoint"),
  list(label = "equal_width_time50_quantile_cov30_midpoint",
       time_method = "equal_width", covariate_method = "quantile",
       time_bins = 50L, covariate_bins = 30L, representative = "midpoint"),
  list(label = "quantile_time50_cov30_midpoint",
       time_method = "quantile", covariate_method = "quantile",
       time_bins = 50L, covariate_bins = 30L, representative = "midpoint"),
  list(label = "quantile_time50_equal_width_cov30_midpoint",
       time_method = "quantile", covariate_method = "equal_width",
       time_bins = 50L, covariate_bins = 30L, representative = "midpoint")
)

additive_representative_specs <- list(
  list(label = "equal_width_time30_cov30_midpoint",
       time_method = "equal_width", covariate_method = "equal_width",
       time_bins = 30L, covariate_bins = 30L, representative = "midpoint"),
  list(label = "equal_width_time30_cov30_mean",
       time_method = "equal_width", covariate_method = "equal_width",
       time_bins = 30L, covariate_bins = 30L, representative = "mean")
)

multiplicative_representative_specs <- list(
  list(label = "quantile_time50_cov30_midpoint",
       time_method = "quantile", covariate_method = "quantile",
       time_bins = 50L, covariate_bins = 30L, representative = "midpoint"),
  list(label = "quantile_time50_cov30_mean",
       time_method = "quantile", covariate_method = "quantile",
       time_bins = 50L, covariate_bins = 30L, representative = "mean")
)

multiplicative_strategy_specs <- list(
  list(label = "equal_width_time50_cov30_midpoint",
       time_method = "equal_width", covariate_method = "equal_width",
       time_bins = 50L, covariate_bins = 30L, representative = "midpoint"),
  list(label = "equal_width_time50_quantile_cov30_midpoint",
       time_method = "equal_width", covariate_method = "quantile",
       time_bins = 50L, covariate_bins = 30L, representative = "midpoint"),
  list(label = "quantile_time50_equal_width_cov30_midpoint",
       time_method = "quantile", covariate_method = "equal_width",
       time_bins = 50L, covariate_bins = 30L, representative = "midpoint"),
  list(label = "quantile_time50_cov30_midpoint",
       time_method = "quantile", covariate_method = "quantile",
       time_bins = 50L, covariate_bins = 30L, representative = "midpoint"),
  list(label = "log_equal_width_time50_cov30_midpoint",
       time_method = "equal_width_log", covariate_method = "equal_width",
       time_bins = 50L, covariate_bins = 30L, representative = "midpoint")
)

outputs <- list()

if (mode %in% c("main", "all")) {
  outputs$main_additive <- run_thesis_job("main/additive_n1000_d5", {
    run_additive_pairwise_binning_experiment(
      n_rep = 30L,
      base_seed = 2026050103L,
      experiment_config = default_additive_pairwise_experiment_config(
        n_train = 1000L,
        d = 5L,
        fit_it = 250L,
        model = 1L,
        rho = 0.5,
        bandwidth = c(0.3, rep(0.25, 4L))
      ),
      binning_specs = additive_main_spec,
      out_dir = thesis_output_dir("main", "additive_n1000_d5")
    )
  })

  outputs$main_multiplicative <- run_thesis_job("main/multiplicative_n1000_d5", {
    run_multiplicative_sparse_binning_experiment(
      n_rep = 30L,
      base_seed = 2026050104L,
      experiment_config = default_multiplicative_sparse_experiment_config(
        n_train = 1000L,
        d = 5L,
        fit_it = 250L,
        model = 1L,
        rho = 0.5,
        bandwidth = c(0.4, rep(0.2, 4L))
      ),
      binning_specs = multiplicative_main_spec,
      out_dir = thesis_output_dir("main", "multiplicative_n1000_d5")
    )
  })
}

if (mode %in% c("larger", "all")) {
  outputs$larger_additive <- run_thesis_job("larger/additive_n2000_d9", {
    run_additive_pairwise_binning_experiment(
      n_rep = 10L,
      base_seed = 2026050103L,
      experiment_config = default_additive_pairwise_experiment_config(
        n_train = 2000L,
        d = 9L,
        fit_it = 250L,
        model = 1L,
        rho = 0.5,
        bandwidth = c(0.3, rep(0.25, 8L))
      ),
      binning_specs = additive_main_spec,
      out_dir = thesis_output_dir("larger", "additive_n2000_d9")
    )
  })

  outputs$larger_multiplicative <- run_thesis_job("larger/multiplicative_n2000_d9", {
    run_multiplicative_sparse_binning_experiment(
      n_rep = 10L,
      base_seed = 2026050104L,
      experiment_config = default_multiplicative_sparse_experiment_config(
        n_train = 2000L,
        d = 9L,
        fit_it = 250L,
        model = 1L,
        rho = 0.5,
        bandwidth = c(0.5, rep(0.2, 8L))
      ),
      binning_specs = multiplicative_main_spec,
      out_dir = thesis_output_dir("larger", "multiplicative_n2000_d9")
    )
  })
}

if (mode %in% c("diagnostics", "all")) {
  outputs$additive_strategy <- run_thesis_job("diagnostics/additive_strategy_n1000_d5", {
    run_additive_pairwise_binning_experiment(
      n_rep = 30L,
      base_seed = 2026050103L,
      include_sbf = TRUE,
      experiment_config = default_additive_pairwise_experiment_config(
        n_train = 1000L,
        d = 5L,
        fit_it = 250L,
        model = 1L,
        rho = 0.5,
        bandwidth = c(0.7, rep(0.25, 4L))
      ),
      binning_specs = additive_strategy_specs,
      out_dir = thesis_output_dir("diagnostics", "additive_strategy_n1000_d5")
    )
  })

  outputs$multiplicative_strategy <- run_thesis_job("diagnostics/multiplicative_strategy_n1000_d5", {
    run_multiplicative_log_time_diagnostic(
      n_rep = 30L,
      base_seed = 2026050104L,
      config = default_multiplicative_log_time_diagnostic_config(
        n_train = 1000L,
        d = 5L,
        fit_it = 250L,
        model = 1L,
        rho = 0.5,
        bandwidth = c(0.4, rep(0.2, 4L))
      ),
      binning_specs = multiplicative_strategy_specs,
      out_dir = thesis_output_dir("diagnostics", "multiplicative_strategy_n1000_d5")
    )
  })

  outputs$additive_representatives <- run_thesis_job("diagnostics/additive_representatives_n1000_d5", {
    run_additive_pairwise_binning_experiment(
      n_rep = 30L,
      base_seed = 2026050103L,
      include_sbf = FALSE,
      experiment_config = default_additive_pairwise_experiment_config(
        n_train = 1000L,
        d = 5L,
        fit_it = 250L,
        model = 1L,
        rho = 0.5,
        bandwidth = c(0.3, rep(0.25, 4L))
      ),
      binning_specs = additive_representative_specs,
      out_dir = thesis_output_dir("diagnostics", "additive_representatives_n1000_d5")
    )
  })

  outputs$multiplicative_representatives <- run_thesis_job("diagnostics/multiplicative_representatives_n1000_d5", {
    run_multiplicative_sparse_binning_experiment(
      n_rep = 30L,
      base_seed = 2026050104L,
      include_sbf = FALSE,
      experiment_config = default_multiplicative_sparse_experiment_config(
        n_train = 1000L,
        d = 5L,
        fit_it = 250L,
        model = 1L,
        rho = 0.5,
        bandwidth = c(0.4, rep(0.2, 4L))
      ),
      binning_specs = multiplicative_representative_specs,
      out_dir = thesis_output_dir("diagnostics", "multiplicative_representatives_n1000_d5")
    )
  })
}

cat(sprintf("\nThesis binning mode '%s' complete.\n", mode))
cat(sprintf("Outputs: %s\n", output_root))
