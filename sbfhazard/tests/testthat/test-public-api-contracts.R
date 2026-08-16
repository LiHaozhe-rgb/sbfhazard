test_that("the package exports its four public functions", {
  expect_setequal(
    getNamespaceExports("sbfhazard"),
    c("sbf_fit", "sbf_fit_binning", "sbf_predict", "sbf_simulate_data")
  )
})

test_that("public functions reject invalid method choices", {
  dat <- sbf_simulate_data(
    n = 40,
    d = 3,
    family = "additive",
    seed = 401
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
    sbf_simulate_data(n = 20, d = 3, family = "unknown"),
    "arg"
  )
  expect_error(
    sbf_fit(dat, bandwidth = rep(0.4, 3), family = "unknown"),
    "arg"
  )
  expect_error(sbf_predict(fit, type = "unknown"), "arg")
})

test_that("multiplicative binning rejects the additive local_constant argument", {
  dat <- sbf_simulate_data(
    n = 40,
    d = 3,
    family = "multiplicative",
    seed = 402
  )

  expect_error(
    sbf_fit_binning(
      dat,
      bandwidth = rep(0.4, 3),
      family = "multiplicative",
      time_bins = 5,
      covariate_bins = 5,
      local_constant = TRUE,
      iterations = 3,
      warn_nonconvergence = FALSE,
      warn_diagnostics = FALSE
    ),
    "local_constant"
  )
})
