graphics.off()

source(file.path("experiments", "binning", "additive_pairwise_experiment.R"))

experiment_config <- default_additive_pairwise_experiment_config(
  n_train = 1000L,
  d = 5L,
  fit_it = 250L,
  target_k = 1L,
  model = 1L,
  violate_cox = TRUE,
  bandwidth = c(0.2, rep(0.25, 4L)),
  z_grid = seq(-1, 1, length.out = 300),
  t_grid = seq(0, 2, length.out = 300),
  rho = 0.5
)

outputs <- run_additive_pairwise_binning_experiment(
  n_rep = 20L,
  base_seed = 2026050103L,
  experiment_config = experiment_config,
  binning_specs = default_additive_pairwise_binning_specs(),
  include_sbf = FALSE
)

cat(sprintf("Summary: %s\n", outputs$summary_file))
