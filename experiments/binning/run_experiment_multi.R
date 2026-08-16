graphics.off()

source(file.path("experiments", "binning", "multiplicative_sparse_experiment.R"))

experiment_config <- default_multiplicative_sparse_experiment_config(
  n_train = 1000L,
  d = 5L,
  fit_it = 250L,
  target_k = 1L,
  model = 2L,
  violate_cox = TRUE,
  bandwidth = NULL,
  z_grid = seq(-1, 1, length.out = 300),
  t_grid = seq(0, 2, length.out = 300),
  rho = 0.5
)

outputs <- run_multiplicative_sparse_binning_experiment(
  n_rep = 20L,
  base_seed = 2026050104L,
  experiment_config = experiment_config,
  binning_specs = default_multiplicative_sparse_binning_specs()
)

cat(sprintf("Summary: %s\n", outputs$summary_file))
