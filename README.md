# SBF Hazard

This repository contains the alpha `sbfhazard` R package and the experiments
used to check its SBF fits and binned approximations. It also keeps the
historical source and baseline code used for implementation comparisons.

## Install

Install the package from GitHub with:

```r
remotes::install_github("LiHaozhe-rgb/sbfhazard", subdir = "sbfhazard")
```

## Main folders

- `sbfhazard/` — the R package for fitting, prediction, simulation and binning.
- `experiments/` — benchmark, binning, diagnostic and smoke-test scripts.
- `baselines/` — historical implementations used in the original-code
  comparison.
- `source/` — older source files kept for reference.
- `test_results/` — selected reference results and generated run outputs.

## Main experiments

Run these commands from the repository root:

```sh
Rscript experiments/benchmarks/benchmark_sbf_vs_origin.R
Rscript experiments/binning/additive_pairwise_experiment.R
Rscript experiments/binning/multiplicative_sparse_experiment.R
Rscript experiments/real_data/trace_additive_hazard.R
```

Larger binning presets are available in:

```sh
Rscript experiments/binning/run_experiment_additive.R
Rscript experiments/binning/run_experiment_multi.R
```

See [`experiments/README.md`](experiments/README.md) for the other checks and
script details. See [`sbfhazard/README.md`](sbfhazard/README.md) for package
installation and usage.

## Checks

Run the package tests with:

```sh
Rscript -e "library(testthat); pkgload::load_all('sbfhazard', export_all = FALSE, quiet = TRUE); testthat::test_dir('sbfhazard/tests/testthat', reporter = 'summary', env = globalenv())"
```

The experiment smoke suite checks the main benchmark, binning experiments and
method checks:

```sh
Rscript experiments/tests/run_smoke_suite.R
```

The experiment scripts use `sbf_simulate_data()`, `sbf_fit()`,
`sbf_fit_binning()` and `sbf_predict()`. Historical code under `source/` and
`baselines/` is kept for comparison and is not part of the package runtime.
