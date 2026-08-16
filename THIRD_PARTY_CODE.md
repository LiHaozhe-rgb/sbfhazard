# Historical and Third-Party Code

The `sbfhazard/` package is the maintained implementation in this repository.
Files under `baselines/` and `source/` are retained for benchmark
reproducibility and historical reference. They are not loaded by the package at
runtime.

## Benchmark baselines

- `baselines/SBF_MH_LL_Additive_Origin.R` is adapted from the public replication
  code for *Smooth Backfitting for Additive Hazard Rates*. The upstream
  repository is
  <https://github.com/MHiabu/Replicate-Smooth-Backfitting-for-Additive-Hazard-Rates>.
- `baselines/SBF_sim_0.1_origin.R` is adapted from the historical simulation
  implementation used for the multiplicative smooth-backfitting comparison.

The benchmark copies retain the original component-update calculations. Local
changes are limited to convergence results, diagnostics, and return values used
by the common benchmark interface.

## Historical source

The files under `source/` are working and replication files collected during
the project. They include earlier estimator, simulation, and prediction code.
They are provided to document provenance and are not maintained as a public API.

The repository's MIT license applies to the maintained project code. Historical
or third-party files remain subject to any terms attached to their original
sources.
