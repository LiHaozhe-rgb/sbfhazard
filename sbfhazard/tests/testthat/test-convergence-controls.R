test_that("fit change includes every fitted component", {
  old <- list(c(0, 0), c(0, 0))
  new <- list(c(0.1, 0), c(0, 0.5))

  expect_equal(sbfhazard:::.sbf_fit_delta(old, new), 0.5)
  expect_true(is.na(sbfhazard:::.sbf_fit_delta(
    list(c(1, Inf)),
    list(c(1, Inf))
  )))
})

test_that("local-linear convergence includes values and slopes", {
  dat <- sbf_simulate_data(
    n = 50,
    d = 3,
    family = "additive",
    seed = 601
  )
  sbf_fit_obj <- sbf_fit(
    dat,
    bandwidth = rep(0.4, 3),
    family = "additive",
    local_constant = FALSE,
    iterations = 1,
    convergence_tol = .Machine$double.xmin,
    warn_nonconvergence = FALSE
  )
  binned <- sbf_fit_binning(
    dat,
    bandwidth = rep(0.4, 3),
    family = "additive",
    time_bins = 5,
    covariate_bins = 5,
    local_constant = FALSE,
    iterations = 1,
    convergence_tol = .Machine$double.xmin,
    warn_nonconvergence = FALSE,
    warn_diagnostics = FALSE
  )

  for (fit in list(sbf_fit_obj, binned)) {
    expect_equal(fit$iterations_used, 1L)
    expect_equal(fit$iterations_max, 1L)
    expect_equal(fit$convergence_tol, .Machine$double.xmin)
    expect_true(is.finite(fit$final_delta))
    expect_true(all(is.finite(unlist(fit$alpha_backfit))))
    expect_true(all(is.finite(unlist(fit$alpha.1_backfit))))
  }
})

test_that("integral approximation choices remain public", {
  dat <- sbf_simulate_data(
    n = 40,
    d = 3,
    family = "additive",
    seed = 602
  )

  for (rule in c("midd", "left", "right")) {
    fit <- sbf_fit(
      dat,
      bandwidth = rep(0.4, 3),
      family = "additive",
      local_constant = TRUE,
      integral_approx = rule,
      iterations = 2,
      warn_nonconvergence = FALSE
    )
    expect_equal(fit$fit_settings$integral_approx, rule)
  }
})

test_that("iteration limit and tolerance are validated", {
  dat <- sbf_simulate_data(
    n = 40,
    d = 3,
    family = "additive",
    seed = 603
  )

  expect_error(
    sbf_fit(
      dat,
      bandwidth = rep(0.4, 3),
      family = "additive",
      iterations = 0,
      warn_nonconvergence = FALSE
    ),
    "iterations"
  )
  expect_error(
    sbf_fit(
      dat,
      bandwidth = rep(0.4, 3),
      family = "additive",
      convergence_tol = 0,
      warn_nonconvergence = FALSE
    ),
    "convergence_tol"
  )
})
