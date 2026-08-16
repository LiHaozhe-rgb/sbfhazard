test_that("multiplicative identification scales follow the four rules", {
  alpha <- c(1, 2, 4, 8)
  grid <- c(-1, -0.5, 0.5, 1)
  widths <- c(1, 2, 2, 1)
  truth <- list(function(x) 0.2 * x)

  sample_mean <- sbfhazard:::.sbf_multiplicative_identification_scale(
    alpha_component = alpha,
    grid = grid,
    widths = widths,
    identification = "sample_mean",
    min_component = 1e-8,
    sample_weights = c(3, 1, 1, 1)
  )
  integral <- sbfhazard:::.sbf_multiplicative_identification_scale(
    alpha_component = alpha,
    grid = grid,
    widths = widths,
    identification = "integral",
    min_component = 1e-8
  )
  origin <- sbfhazard:::.sbf_multiplicative_identification_scale(
    alpha_component = alpha,
    grid = grid,
    widths = widths,
    identification = "origin",
    min_component = 1e-8,
    truth_functions = truth,
    covariate_index = 1
  )
  jasa <- sbfhazard:::.sbf_multiplicative_identification_scale(
    alpha_component = alpha,
    grid = grid,
    widths = widths,
    identification = "jasa",
    min_component = 1e-8,
    truth_functions = truth,
    covariate_index = 1
  )

  expect_equal(sample_mean, exp(stats::weighted.mean(log(alpha), c(3, 1, 1, 1))))
  expect_equal(integral, stats::weighted.mean(alpha, widths))
  expect_equal(origin, exp(log(alpha[2]) - truth[[1]](grid[2])))

  anchor <- c(2L, 3L)
  truth_x <- c(min(abs(grid + 0.5)), min(abs(grid - 0.5)))
  expect_equal(
    jasa,
    exp(mean(log(alpha[anchor]) - truth[[1]](truth_x)))
  )
})

test_that("multiplicative identification preserves the assembled hazard", {
  baseline <- c(0.5, 0.8)
  components <- list(
    c(1, 2, 4, 8),
    c(1.5, 1, 2, 3)
  )
  grids <- list(
    c(-1, -0.5, 0.5, 1),
    c(-1, -0.5, 0.5, 1)
  )
  widths <- list(rep(1, 4), rep(1, 4))
  truth <- list(function(x) 0 * x, function(x) 0.1 * x)

  for (rule in c("sample_mean", "integral", "origin", "jasa")) {
    scales <- vapply(seq_along(components), function(k) {
      sbfhazard:::.sbf_multiplicative_identification_scale(
        alpha_component = components[[k]],
        grid = grids[[k]],
        widths = widths[[k]],
        identification = rule,
        min_component = 1e-8,
        sample_weights = rep(1, 4),
        truth_functions = if (rule %in% c("origin", "jasa")) truth else NULL,
        covariate_index = k
      )
    }, numeric(1))

    before <- outer(baseline, components[[1]] * components[[2]])
    after <- outer(
      baseline * prod(scales),
      (components[[1]] / scales[1]) * (components[[2]] / scales[2])
    )
    expect_equal(after, before)
  }
})

test_that("SBF sample-mean identification normalizes covariates", {
  dat <- sbf_simulate_data(
    n = 60,
    d = 3,
    family = "multiplicative",
    seed = 801
  )
  fit <- sbf_fit(
    dat,
    bandwidth = rep(0.4, 3),
    family = "multiplicative",
    identification = "sample_mean",
    iterations = 4,
    warn_nonconvergence = FALSE
  )

  for (component in 2:3) {
    expect_equal(
      mean(log(fit$alpha_backfit[[component]])),
      0,
      tolerance = 1e-8
    )
  }
})

test_that("binned integral identification normalizes covariates", {
  dat <- sbf_simulate_data(
    n = 60,
    d = 3,
    family = "multiplicative",
    seed = 802
  )
  fit <- sbf_fit_binning(
    dat,
    bandwidth = rep(0.4, 3),
    family = "multiplicative",
    time_bins = 6,
    covariate_bins = 5,
    identification = "integral",
    iterations = 4,
    warn_nonconvergence = FALSE,
    warn_diagnostics = FALSE
  )
  state <- fit$multiplicative_sparse_predict_state

  for (component in 2:3) {
    average <- stats::weighted.mean(
      fit$alpha_backfit[[component]],
      state$widths[[component]]
    )
    expect_equal(average, 1, tolerance = 1e-8)
  }
})

test_that("weak-information multiplicative fits keep components positive", {
  dat <- data.frame(
    time = seq_len(8),
    x1 = seq(-1, 1, length.out = 8),
    x2 = seq(1, -1, length.out = 8),
    status = c(1, rep(0, 7))
  )
  sbf_fit_obj <- sbf_fit(
    dat,
    bandwidth = rep(0.15, 3),
    family = "multiplicative",
    iterations = 1,
    warn_nonconvergence = FALSE
  )
  binned <- sbf_fit_binning(
    dat,
    bandwidth = rep(0.15, 3),
    family = "multiplicative",
    time_bins = 4,
    covariate_bins = 4,
    iterations = 1,
    warn_nonconvergence = FALSE,
    warn_diagnostics = FALSE
  )

  for (fit in list(sbf_fit_obj, binned)) {
    values <- unlist(fit$alpha_backfit, use.names = FALSE)
    expect_true(all(is.finite(values)))
    expect_true(all(values > 0))
    expect_gt(fit$fit_diagnostics$component_floor_count, 0L)
  }
})

test_that("multiplicative identification arguments are validated", {
  multiplicative_data <- sbf_simulate_data(
    n = 40,
    d = 3,
    family = "multiplicative",
    seed = 803
  )
  additive_data <- sbf_simulate_data(
    n = 40,
    d = 3,
    family = "additive",
    seed = 804
  )

  expect_error(
    sbf_fit(
      multiplicative_data,
      bandwidth = rep(0.4, 3),
      family = "multiplicative",
      identification = "unknown"
    ),
    "should be one of"
  )
  expect_error(
    sbf_fit(
      multiplicative_data,
      bandwidth = rep(0.4, 3),
      family = "multiplicative",
      identification = "origin"
    ),
    "truth_functions"
  )
  expect_error(
    sbf_fit(
      multiplicative_data,
      bandwidth = rep(0.4, 3),
      family = "multiplicative",
      identification = "origin",
      truth_functions = list(
        function(x) rep(NA_real_, length(x)),
        function(x) x
      )
    ),
    "truth_functions"
  )
  expect_error(
    sbf_fit(
      additive_data,
      bandwidth = rep(0.4, 3),
      family = "additive",
      identification = "integral"
    ),
    "identification"
  )
})
