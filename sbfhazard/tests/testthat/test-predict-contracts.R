linear_interpolation_reference <- function(grid, values, xout) {
  keep <- is.finite(grid) & is.finite(values)
  grid <- as.numeric(grid[keep])
  values <- as.numeric(values[keep])
  order_index <- order(grid)
  grid <- grid[order_index]
  values <- values[order_index]
  unique_index <- !duplicated(grid)

  stats::approx(
    x = grid[unique_index],
    y = values[unique_index],
    xout = xout,
    method = "linear",
    rule = 2
  )$y
}

test_that("formula prediction is the default", {
  dat <- sbf_simulate_data(
    n = 60,
    d = 3,
    family = "additive",
    seed = 901
  )
  fit <- sbf_fit(
    dat,
    bandwidth = rep(0.4, 3),
    family = "additive",
    local_constant = TRUE,
    iterations = 4,
    warn_nonconvergence = FALSE
  )

  default <- sbf_predict(
    fit,
    type = "component",
    component = 2,
    xout = c(-0.1, 0.1),
    warn_diagnostics = FALSE
  )
  formula <- sbf_predict(
    fit,
    type = "component",
    component = 2,
    xout = c(-0.1, 0.1),
    prediction = "formula",
    warn_diagnostics = FALSE
  )

  expect_equal(default, formula)
})

test_that("formula survival preserves unsorted and duplicate times", {
  dat <- sbf_simulate_data(
    n = 60,
    d = 3,
    family = "additive",
    seed = 902
  )
  fit <- sbf_fit(
    dat,
    bandwidth = rep(0.4, 3),
    family = "additive",
    local_constant = TRUE,
    iterations = 4,
    warn_nonconvergence = FALSE
  )
  ordered_times <- as.numeric(stats::quantile(fit$x.grid[[1]], c(0.25, 0.5, 0.75)))
  query_times <- ordered_times[c(3, 1, 3, 2)]
  newdata <- matrix(c(0, 0), nrow = 1)

  prediction <- sbf_predict(
    fit,
    type = "survival",
    data = newdata,
    times = query_times,
    warn_diagnostics = FALSE
  )
  reference <- sbf_predict(
    fit,
    type = "survival",
    data = newdata,
    times = ordered_times,
    warn_diagnostics = FALSE
  )

  expect_equal(prediction$times, query_times)
  expect_equal(
    prediction$survival,
    reference$survival[, match(query_times, ordered_times), drop = FALSE]
  )
})

test_that("SBF interpolation uses fitted-grid values and endpoints", {
  dat <- sbf_simulate_data(
    n = 60,
    d = 3,
    family = "additive",
    seed = 903
  )
  fit <- sbf_fit(
    dat,
    bandwidth = rep(0.4, 3),
    family = "additive",
    local_constant = TRUE,
    iterations = 4,
    warn_nonconvergence = FALSE
  )
  grid <- fit$x.grid[[2]]
  xout <- c(min(grid) - 1, 0, max(grid) + 1)

  prediction <- sbf_predict(
    fit,
    type = "component",
    component = 2,
    xout = xout,
    prediction = "interpolation",
    warn_diagnostics = FALSE
  )

  expect_equal(
    as.numeric(prediction),
    linear_interpolation_reference(grid, fit$alpha_backfit[[2]], xout)
  )
  expect_equal(prediction[1], fit$alpha_backfit[[2]][1])
  expect_equal(prediction[3], tail(fit$alpha_backfit[[2]], 1))
})

test_that("additive interpolation gives the assembled hazard and survival", {
  dat <- sbf_simulate_data(
    n = 70,
    d = 3,
    family = "additive",
    seed = 904
  )
  fit <- sbf_fit(
    dat,
    bandwidth = rep(0.4, 3),
    family = "additive",
    local_constant = TRUE,
    iterations = 4,
    warn_nonconvergence = FALSE
  )
  newdata <- matrix(c(0.1, -0.1, -0.2, 0.2), nrow = 2, byrow = TRUE)
  times <- c(0.2, 0.4)

  hazard <- sbf_predict(
    fit,
    type = "hazard",
    data = newdata,
    times = times,
    prediction = "interpolation",
    warn_diagnostics = FALSE
  )
  baseline <- linear_interpolation_reference(
    fit$x.grid[[1]],
    fit$alpha_backfit[[1]],
    times
  )
  covariates <- cbind(
    linear_interpolation_reference(fit$x.grid[[2]], fit$alpha_backfit[[2]], newdata[, 1]),
    linear_interpolation_reference(fit$x.grid[[3]], fit$alpha_backfit[[3]], newdata[, 2])
  )
  survival <- sbf_predict(
    fit,
    type = "survival",
    data = newdata,
    times = times,
    prediction = "interpolation",
    warn_diagnostics = FALSE
  )
  support_start <- fit$component_info$support_min[1]
  positive <- which(times > support_start)
  integration_grid <- sort(unique(c(support_start, fit$x.grid[[1]], times[positive])))
  integration_hazard <- sbf_predict(
    fit,
    type = "hazard",
    data = newdata,
    times = integration_grid,
    prediction = "interpolation",
    warn_diagnostics = FALSE
  )$hazard
  cumulative_hazard <- t(apply(integration_hazard, 1, function(values) {
    c(0, cumsum(0.5 * (head(values, -1) + tail(values, -1)) * diff(integration_grid)))
  }))
  expected_survival <- matrix(1, nrow = nrow(newdata), ncol = length(times))
  expected_survival[, positive] <- exp(
    -cumulative_hazard[, match(times[positive], integration_grid), drop = FALSE]
  )

  expect_equal(
    unname(hazard$hazard),
    unname(outer(rowSums(covariates), baseline, "+"))
  )
  expect_equal(
    survival$survival,
    expected_survival
  )
})

test_that("multiplicative interpolation multiplies fitted components", {
  dat <- sbf_simulate_data(
    n = 70,
    d = 3,
    family = "multiplicative",
    seed = 905
  )
  fit <- sbf_fit(
    dat,
    bandwidth = rep(0.4, 3),
    family = "multiplicative",
    iterations = 4,
    warn_nonconvergence = FALSE
  )
  newdata <- matrix(c(0.1, -0.1), nrow = 1)
  times <- c(0.2, 0.4)

  prediction <- sbf_predict(
    fit,
    type = "hazard",
    data = newdata,
    times = times,
    prediction = "interpolation",
    warn_diagnostics = FALSE
  )
  baseline <- pmax(1e-8, linear_interpolation_reference(
    fit$x.grid[[1]], fit$alpha_backfit[[1]], times
  ))
  covariate_term <- prod(c(
    pmax(1e-8, linear_interpolation_reference(
      fit$x.grid[[2]], fit$alpha_backfit[[2]], newdata[1, 1]
    )),
    pmax(1e-8, linear_interpolation_reference(
      fit$x.grid[[3]], fit$alpha_backfit[[3]], newdata[1, 2]
    ))
  ))

  expect_equal(as.numeric(prediction$hazard), covariate_term * baseline)
})

test_that("binned fits support formula prediction but not interpolation", {
  dat <- sbf_simulate_data(
    n = 60,
    d = 3,
    family = "additive",
    seed = 906
  )
  fit <- sbf_fit_binning(
    dat,
    bandwidth = rep(0.4, 3),
    family = "additive",
    time_bins = 5,
    covariate_bins = 5,
    local_constant = TRUE,
    iterations = 4,
    warn_nonconvergence = FALSE,
    warn_diagnostics = FALSE
  )

  prediction <- sbf_predict(
    fit,
    type = "hazard",
    data = matrix(c(0, 0), nrow = 1),
    times = c(0.2, 0.4),
    warn_diagnostics = FALSE
  )
  expect_true(all(is.finite(prediction$hazard)))
  expect_error(
    sbf_predict(
      fit,
      type = "component",
      component = 2,
      xout = 0,
      prediction = "interpolation",
      warn_diagnostics = FALSE
    ),
    "formula"
  )
})

test_that("binned component prediction uses bin-step fallback at zero support", {
  dat <- sbf_simulate_data(
    n = 60,
    d = 3,
    family = "additive",
    seed = 907
  )
  fit <- sbf_fit_binning(
    dat,
    bandwidth = rep(0.05, 3),
    family = "additive",
    time_bins = 6,
    covariate_bins = 6,
    local_constant = TRUE,
    iterations = 3,
    warn_nonconvergence = FALSE,
    warn_diagnostics = FALSE
  )
  xout <- max(fit$component_info$support_max[2]) + 10
  prediction <- sbf_predict(
    fit,
    type = "component",
    component = 2,
    xout = xout,
    warn_diagnostics = FALSE
  )
  diagnostics <- attr(prediction, "prediction_diagnostics")

  expect_equal(diagnostics$kernel_support_zero_count, 1L)
  expect_equal(diagnostics$bin_step_fallback_count, 1L)
  expect_equal(as.numeric(prediction), tail(fit$alpha_backfit[[2]], 1))
})

test_that("multiplicative component protection replaces invalid values", {
  protected <- sbfhazard:::.sbf_clamp_component_values(
    c(NA_real_, -1, 0, 2),
    lower_bound = 1e-8
  )

  expect_equal(protected, c(1e-8, 1e-8, 1e-8, 2))
})

test_that("hazard prediction applies the requested lower bound", {
  dat <- sbf_simulate_data(
    n = 50,
    d = 3,
    family = "additive",
    seed = 908
  )
  fit <- sbf_fit(
    dat,
    bandwidth = rep(0.4, 3),
    family = "additive",
    local_constant = TRUE,
    iterations = 3,
    warn_nonconvergence = FALSE
  )
  prediction <- sbf_predict(
    fit,
    type = "hazard",
    data = matrix(c(0, 0), nrow = 1),
    times = c(0.2, 0.4),
    clamp_nonnegative = TRUE,
    min_hazard = 1e6,
    warn_diagnostics = FALSE
  )

  expect_equal(as.numeric(prediction$hazard), rep(1e6, 2))
  expect_equal(prediction$clamp_count, 2L)
})

test_that("prediction requires arguments for the selected output", {
  dat <- sbf_simulate_data(
    n = 40,
    d = 3,
    family = "additive",
    seed = 909
  )
  fit <- sbf_fit(
    dat,
    bandwidth = rep(0.4, 3),
    family = "additive",
    local_constant = TRUE,
    iterations = 3,
    warn_nonconvergence = FALSE
  )

  expect_error(
    sbf_predict(fit, type = "component", component = 2),
    "xout"
  )
  expect_error(sbf_predict(fit, type = "hazard"), "data")
})
