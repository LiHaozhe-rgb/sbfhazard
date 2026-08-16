# sbfhazard

Alpha R package for smooth backfitting hazard estimators.

`sbfhazard` provides additive and multiplicative SBF estimators, binned
approximations, prediction, and simulation. The package is intended for the
examples and experiments in this repository.

## Install

From GitHub:

```r
remotes::install_github("LiHaozhe-rgb/sbfhazard", subdir = "sbfhazard")
```

From a local clone, run `R CMD INSTALL sbfhazard` from the repository root.

## Data and formulas

The fitting functions accept a numeric data frame with `time` and `status`
columns, or a survival formula. With a formula, the right-hand side is used
to build the fitted covariates. For example:

```r
survival::Surv(time, status) ~ V1 + log(V2 + 2)
```

Without a formula, columns other than `time` and `status` are treated as
numeric fitted features. At least one observed event is needed after missing
rows are removed.

For formula prediction, `data` can contain the original right-hand-side
variables. Matrix and vector inputs are already on the fitted feature scale.
Component prediction also uses this fitted scale: component 1 is the
baseline/time component and later components are the fitted covariates.

Factor levels and contrasts from the training data are stored with the fit.
New prediction data must use levels that were present during fitting.

The prediction help pages describe interpolation, support boundaries, numerical
integration, floors and diagnostic messages. These details are kept with the
function documentation rather than repeated here.

## Additive SBF fit

```r
library(sbfhazard)

dat <- sbf_simulate_data(n = 100, d = 3, family = "additive", seed = 1)
fit_add <- sbf_fit(
  data = dat,
  bandwidth = c(0.3, 0.3, 0.3),
  family = "additive",
  formula = survival::Surv(time, status) ~ V1 + V2,
  local_constant = TRUE,
  iterations = 20,
  convergence_tol = 1e-3,
  warn_nonconvergence = FALSE
)

sbf_predict(
  fit_add,
  type = "component",
  component = 2,
  xout = seq(-0.8, 0.8, length.out = 25)
)
```

## Multiplicative SBF fit

```r
dat <- sbf_simulate_data(n = 100, d = 3, family = "multiplicative", seed = 2)
fit_mult <- sbf_fit(
  data = dat,
  bandwidth = c(0.3, 0.3, 0.3),
  family = "multiplicative",
  formula = survival::Surv(time, status) ~ V1 + V2,
  iterations = 20,
  warn_nonconvergence = FALSE
)

sbf_predict(
  fit_mult,
  type = "survival",
  data = data.frame(V1 = 0, V2 = 0),
  times = c(0.25, 0.5, 0.75)
)
```

## Binned fits

`sbf_fit_binning()` evaluates a binned approximation to the SBF update
equations. The SBF fit is the reference when comparing speed and accuracy.

```r
dat <- sbf_simulate_data(n = 100, d = 3, family = "additive", seed = 3)
fit_bin_add <- sbf_fit_binning(
  dat,
  bandwidth = c(0.3, 0.3, 0.3),
  family = "additive",
  time_bins = 8,
  covariate_bins = 8,
  iterations = 20,
  warn_nonconvergence = FALSE,
  warn_diagnostics = FALSE
)

dat <- sbf_simulate_data(n = 100, d = 3, family = "multiplicative", seed = 4)
fit_bin_mult <- sbf_fit_binning(
  dat,
  bandwidth = c(0.3, 0.3, 0.3),
  family = "multiplicative",
  time_bins = 8,
  covariate_bins = 8,
  iterations = 20,
  warn_nonconvergence = FALSE,
  warn_diagnostics = FALSE
)
```

## Examples and returned objects

The example scripts simulate data, fit SBF and binned models, and make simple
component and survival predictions:

```sh
Rscript sbfhazard/examples/additive_demo.R
Rscript sbfhazard/examples/multiplicative_demo.R
```

Fit objects contain fitted components, support grids, settings and
diagnostics. Binned fits also store binning metadata.

## Public API

- `sbf_fit()` — fit an additive or multiplicative SBF estimator.
- `sbf_fit_binning()` — fit an additive or multiplicative binned approximation.
- `sbf_predict()` — return component, hazard or survival predictions.
- `sbf_simulate_data()` — generate data for examples and experiments.
