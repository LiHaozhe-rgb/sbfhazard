# Experiments

The scripts in this directory use the `sbfhazard` package for benchmarks,
binning comparisons, diagnostics and the TRACE example. Run them from the
repository root.

## Main scripts

- `benchmarks/benchmark_sbf_vs_origin.R` compares package SBF fits with the
  historical implementations in `baselines/`.
- `binning/additive_pairwise_experiment.R` runs the additive binning study.
- `binning/multiplicative_sparse_experiment.R` runs the multiplicative binning
  study.
- `binning/run_experiment_additive.R` runs a larger additive setting.
- `binning/run_experiment_multi.R` runs a larger multiplicative setting.
- `real_data/trace_additive_hazard.R` compares SBF and binned additive LL fits
  on the TRACE cohort.

## Method checks

- `method_checks/multiplicative_identification.R` checks the multiplicative
  identification rules.
- `method_checks/prediction_methods.R` compares formula and interpolation
  prediction.
- `method_checks/run_multiplicative_identification.R` runs the larger
  identification check.
- `method_checks/run_prediction_methods.R` runs the larger prediction check.

## Diagnostics

- `diagnostics/binning_layout_diagnostics.R` summarizes simulated data and bin
  layouts.
- `diagnostics/multiplicative_log_time_diagnostic.R` checks sensitivity to
  log-spaced time bins. It is diagnostic-only, not a main benchmark.

## Checks

The smoke suite checks the main benchmark, both binning experiments and the
method checks:

```sh
Rscript experiments/tests/run_smoke_suite.R
```

Run the TRACE and diagnostic checks separately when their inputs or outputs
change:

```sh
Rscript experiments/tests/check_trace_additive_hazard.R
Rscript experiments/tests/check_binning_layout_diagnostics.R
Rscript experiments/tests/check_multiplicative_log_time_diagnostic.R
```

The experiments use `sbf_simulate_data()`, `sbf_fit()`,
`sbf_fit_binning()` and `sbf_predict()`. Results are written below
`test_results/runs/`, with smoke outputs kept separate from full experiment
outputs.
