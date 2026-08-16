graphics.off()

suppressPackageStartupMessages({
  library(survival)
  library(prodlim)
})

source(file.path("experiments", "script_utils.R"))
project.root <- sbf_project_root()
sbf_load_package(project.root)
source(file.path(project.root, "experiments", "benchmarks", "benchmark_origin_api.R"))
source(file.path(project.root, "experiments", "benchmarks", "benchmark_algorithms.R"))
source(file.path(project.root, "experiments", "benchmarks", "benchmark_outputs.R"))
source(file.path(project.root, "baselines", "SBF_MH_LL_Additive_Origin.R"))
source(file.path(project.root, "baselines", "SBF_sim_0.1_origin.R"))

expect_true <- function(x, message) {
  if (!isTRUE(x)) {
    stop(message, call. = FALSE)
  }
}

expect_equal <- function(x, y, message, tolerance = 1e-8) {
  if (!isTRUE(all.equal(x, y, tolerance = tolerance))) {
    stop(message, call. = FALSE)
  }
}

add_origin_fun <- SBF.MH.CLL
mult_origin_fun <- SBF.MH.LC
z_grid <- seq(-0.8, 0.8, length.out = 30)
t_grid <- seq(0, 1.2, length.out = 30)
results <- list(
  run_additive_local_linear(
    add_origin_fun = add_origin_fun,
    n = 50L,
    d = 3L,
    rho = 0.2,
    model = 1L,
    violate_cox = TRUE,
    iterations = 6L,
    n_rep = 1L,
    seed_base = 20261000L,
    target_k = 1L,
    bandwidth = rep(0.35, 3),
    z_grid = z_grid,
    t_grid = t_grid
  ),
  run_additive_local_constant(
    add_origin_fun = add_origin_fun,
    n = 50L,
    d = 3L,
    rho = 0.2,
    model = 1L,
    violate_cox = TRUE,
    iterations = 6L,
    n_rep = 1L,
    seed_base = 20261010L,
    target_k = 1L,
    bandwidth = rep(0.35, 3),
    z_grid = z_grid,
    t_grid = t_grid
  ),
  run_multiplicative_local_constant(
    mul_origin_fun = mult_origin_fun,
    n = 40L,
    d = 3L,
    rho = 0.2,
    model = 1L,
    iterations = 6L,
    n_rep = 1L,
    seed_base = 20262000L,
    target_k = 1L,
    bandwidth = rep(0.35, 3),
    z_grid = z_grid,
    t_grid = t_grid
  )
)
out_dir <- sbf_results_dir(
  project.root,
  "runs",
  "smoke",
  "benchmark_sbf_vs_origin"
)
outputs <- benchmark_write_outputs(results, out_dir)

for (path in c(
  outputs$paths$summary_csv,
  outputs$paths$runs_csv,
  outputs$paths$pointwise_csv,
  outputs$paths$plots
)) {
  expect_true(file.exists(path), sprintf("Missing benchmark output: %s", path))
}

expected_algorithms <- c(
  "additive_local_linear",
  "additive_local_constant",
  "multiplicative_local_constant"
)
expect_true(
  all(expected_algorithms %in% outputs$runs$algorithm),
  "Benchmark runs are missing an estimator."
)
expect_true(
  all(expected_algorithms %in% outputs$pointwise$algorithm),
  "Benchmark curves are missing an estimator."
)
expect_true(
  all(c(
    "algorithm", "time_origin", "time_sbf", "paired_converged",
    "effect_mse_sbf_vs_origin", "survival_mse_sbf_vs_origin"
  ) %in% names(outputs$runs)),
  "Benchmark runs are missing required columns."
)
expect_true(
  all(c("curve_type", "truth", "origin", "sbf") %in% names(outputs$pointwise)),
  "Benchmark curves are missing required columns."
)
expect_true(any(outputs$runs$paired_converged), "No valid benchmark comparison was produced.")
expect_true(all(is.finite(outputs$runs$time_origin)), "Origin runtime is not finite.")
expect_true(all(is.finite(outputs$runs$time_sbf)), "Optimized runtime is not finite.")
expect_true(any(is.finite(outputs$pointwise$origin)), "Origin curves are not finite.")
expect_true(any(is.finite(outputs$pointwise$sbf)), "Optimized curves are not finite.")

mult_effect <- outputs$pointwise[
  outputs$pointwise$algorithm == "multiplicative_local_constant" &
    outputs$pointwise$curve_type == "effect",
  ,
  drop = FALSE
]
expect_equal(
  mult_effect$truth,
  sbf_true_multiplicative_phi(mult_effect$grid_value, 1L, violate_cox = TRUE),
  "Multiplicative benchmark truth changed.",
  tolerance = 1e-12
)

additive_data <- sbf_simulate_data(
  n = 30L,
  d = 3L,
  family = "additive",
  seed = 20263000L
)
additive_fit <- sbf_fit_origin(
  data = additive_data,
  bandwidth = rep(0.4, 3),
  family = "additive",
  origin_engine = add_origin_fun,
  iterations = 1L,
  kernel = "epanechnikov",
  local_constant = TRUE
)
additive_reference <- add_origin_fun(
  formula = survival::Surv(time, status) ~ V1 + V2,
  data = additive_fit$data,
  bandwidth = rep(0.4, 3),
  weight = "sw",
  x.grid = NULL,
  n.grid.additional = 0L,
  integral.approx = "midd",
  it = 2L,
  kern = sbfhazard:::.sbf_kernel_spec("epanechnikov")$kernel,
  initial = NULL,
  convergence_tol = 0.001,
  kcorr = TRUE,
  LC = TRUE,
  wrong = FALSE,
  classic.backfit = FALSE,
  print = FALSE
)
expect_true(additive_fit$iterations_used == 1L, "Additive adapter sweep count changed.")
expect_equal(
  additive_fit$alpha_backfit,
  additive_reference$alpha_backfit,
  "Additive adapter no longer matches one original sweep."
)

multiplicative_data <- sbf_simulate_data(
  n = 30L,
  d = 3L,
  family = "multiplicative",
  seed = 20263001L
)
truth_functions <- lapply(seq_len(2L), function(k) {
  force(k)
  function(z) sbf_true_multiplicative_phi(z, k, violate_cox = TRUE)
})
multiplicative_fit <- sbf_fit_origin(
  data = multiplicative_data,
  bandwidth = rep(0.4, 3),
  family = "multiplicative",
  origin_engine = mult_origin_fun,
  iterations = 1L,
  kernel = "epanechnikov",
  truth_functions = truth_functions
)
multiplicative_reference <- mult_origin_fun(
  data = multiplicative_fit$data,
  b.grid = rep(0.4, 3),
  x.grid = NULL,
  n.grid.additional = 0L,
  integral.approx = "midd",
  it = 2L,
  kern = sbfhazard:::.sbf_kernel_spec("epanechnikov")$kernel,
  initial = NULL,
  convergence_tol = 0.001,
  verbose = FALSE,
  phi = truth_functions
)
expect_true(multiplicative_fit$iterations_used == 1L, "Multiplicative adapter sweep count changed.")
expect_equal(
  multiplicative_fit$alpha_backfit,
  multiplicative_reference$alpha_backfit,
  "Multiplicative adapter no longer matches one original sweep."
)

cat("Benchmark sbf-vs-origin smoke passed.\n")
