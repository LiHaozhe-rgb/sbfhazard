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

out_dir <- sbf_results_dir(
  project.root,
  "runs",
  "experiments",
  "benchmarks",
  "sbf_vs_origin"
)

source(file.path(project.root, "baselines", "SBF_MH_LL_Additive_Origin.R"))
source(file.path(project.root, "baselines", "SBF_sim_0.1_origin.R"))

cat("Running SBF-vs-origin benchmark...\n")

benchmark_results <- list(
  run_additive_local_linear(
    add_origin_fun = SBF.MH.CLL,
    n = 600L,
    d = 5L,
    rho = 0.5,
    model = 1L,
    violate_cox = TRUE,
    iterations = 300L,
    n_rep = 30L,
    seed_base = 1020L,
    target_k = 1L,
    bandwidth = c(0.2, 0.3, 0.3, 0.3, 0.3),
    z_grid = seq(-1, 1, length.out = 150),
    t_grid = seq(0, 2.0, length.out = 150)
  ),
  run_additive_local_constant(
    add_origin_fun = SBF.MH.CLL,
    n = 600L,
    d = 5L,
    rho = 0.5,
    model = 1L,
    violate_cox = TRUE,
    iterations = 300L,
    n_rep = 30L,
    seed_base = 1020L,
    target_k = 1L,
    bandwidth = c(0.2, 0.3, 0.3, 0.3, 0.3),
    z_grid = seq(-1, 1, length.out = 150),
    t_grid = seq(0, 2.0, length.out = 150)
  ),
  run_multiplicative_local_constant(
    mul_origin_fun = SBF.MH.LC,
    n = 600L,
    d = 5L,
    rho = 0.5,
    model = 1L,
    iterations = 300L,
    n_rep = 30L,
    seed_base = 1020L,
    target_k = 1L,
    bandwidth = c(0.2, 0.3, 0.3, 0.3, 0.3),
    z_grid = seq(-1, 1, length.out = 150),
    t_grid = seq(0, 2.0, length.out = 150)
  )
)

benchmark_outputs <- benchmark_write_outputs(benchmark_results, out_dir)

cat("Benchmark complete.\n")
cat(sprintf("Output directory: %s\n", out_dir))
cat(sprintf("Summary: %s\n", benchmark_outputs$paths$summary_csv))
cat(sprintf("Runs: %s\n", benchmark_outputs$paths$runs_csv))
cat(sprintf("Pointwise: %s\n", benchmark_outputs$paths$pointwise_csv))
