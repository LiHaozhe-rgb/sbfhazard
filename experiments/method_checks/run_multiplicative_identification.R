graphics.off()

source(file.path("experiments", "method_checks", "multiplicative_identification.R"))

config <- default_multiplicative_identification_method_check_config(
  n_train = 600L,
  d = 5L,
  fit_it = 250L,
  target_k = 1L,
  model = 1L,
  violate_cox = TRUE,
  bandwidth = c(0.2, rep(0.3, 4L)),
  z_grid = seq(-1, 1, length.out = 150),
  t_grid = seq(0, 2, length.out = 150),
  rho = 0.5
)

outputs <- run_multiplicative_identification_method_check(
  n_rep = 20L,
  base_seed = 1020L,
  config = config
)

cat(sprintf("Summary: %s\n", outputs$summary_file))
cat(sprintf("Runs: %s\n", outputs$runs_file))
cat(sprintf("Curves: %s\n", outputs$curves_file))
cat(sprintf("Plot: %s\n", outputs$plot_file))
