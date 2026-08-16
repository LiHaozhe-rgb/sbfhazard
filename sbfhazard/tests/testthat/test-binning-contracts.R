test_that("equal-width and quantile grids follow their definitions", {
  x <- 0:11
  equal_width <- sbfhazard:::.sbf_binning_bin_numeric_vector(
    x,
    bins = 3,
    rule = "equal_width",
    representative = "midpoint"
  )
  quantile <- sbfhazard:::.sbf_binning_bin_numeric_vector(
    x,
    bins = 3,
    rule = "quantile",
    representative = "midpoint"
  )

  expect_equal(equal_width$breaks, seq(0, 11, length.out = 4))
  expect_equal(
    equal_width$representatives,
    0.5 * (equal_width$breaks[-1] + head(equal_width$breaks, -1))
  )
  expect_equal(
    quantile$breaks,
    as.numeric(stats::quantile(x, seq(0, 1, length.out = 4), type = 8))
  )
  expect_lte(diff(range(tabulate(quantile$bin_index, nbins = 3))), 1)
})

test_that("mean representatives use observed values and midpoint empty bins", {
  x <- c(0, 0.1, 0.2, 9.8, 9.9, 10)
  midpoint <- sbfhazard:::.sbf_binning_bin_numeric_vector(
    x,
    bins = 4,
    rule = "equal_width",
    representative = "midpoint"
  )
  mean_rep <- sbfhazard:::.sbf_binning_bin_numeric_vector(
    x,
    bins = 4,
    rule = "equal_width",
    representative = "mean"
  )
  counts <- tabulate(mean_rep$bin_index, nbins = 4)

  expect_equal(mean_rep$representatives[1], mean(x[x < midpoint$breaks[2]]))
  expect_equal(
    mean_rep$representatives[counts == 0],
    midpoint$representatives[counts == 0]
  )
})

test_that("repeated quantile cut points reduce the stored grid", {
  result <- sbfhazard:::.sbf_binning_bin_numeric_vector(
    c(rep(0, 30), 1:7),
    bins = 6,
    rule = "quantile",
    representative = "midpoint"
  )

  expect_equal(result$bins_requested, 6L)
  expect_equal(result$effective_bins, 2L)
  expect_equal(result$empty_bins, 0L)
})

test_that("pairwise exposure is grouped by both covariate bins", {
  overlap <- rbind(c(1, 0), c(0, 2), c(1, 1))
  bin_index <- list(
    c(1L, 2L, 2L),
    c(1L, 1L, 2L),
    c(1L, 2L, 3L)
  )

  result <- sbfhazard:::.sbf_add_pairwise_pairwise_exposure(
    bin_index = bin_index,
    overlap = overlap,
    n_grid = c(2L, 2L, 3L)
  )
  expected <- matrix(
    c(1, 0, 2, 0, 0, 2),
    nrow = 2,
    ncol = 3
  )

  expect_equal(result$pairwise_sum[[2]][[3]], expected)
  expect_equal(result$pairwise_sum[[3]][[2]], t(expected))
})

test_that("multiplicative exposure is grouped by observed covariate-bin pattern", {
  dat <- data.frame(
    time = c(1, 2, 3, 4),
    z1 = c(0, 0, 1, 1),
    z2 = c(0, 1, 0, 1),
    status = c(1, 0, 1, 0)
  )
  precompute <- sbfhazard:::.sbf_mult_sparse_build_fit_precompute(
    data = dat,
    bandwidth = c(1, 1, 1),
    time_bins = 2,
    covariate_bins = 2,
    time_binning_method = "equal_width",
    covariate_binning_method = "equal_width",
    representative = "midpoint"
  )

  covariate_bins <- lapply(2:3, function(component) {
    cut(
      dat[[component]],
      breaks = precompute$breaks[[component]],
      include.lowest = TRUE,
      labels = FALSE
    )
  })
  pattern <- sbfhazard:::.sbf_mult_sparse_joint_key(covariate_bins)
  overlap <- sbfhazard:::.sbf_add_pairwise_time_overlap(
    dat$time,
    head(precompute$breaks[[1]], -1),
    tail(precompute$breaks[[1]], -1)
  )
  expected <- rowsum(overlap, group = as.integer(pattern), reorder = TRUE)

  expect_equal(precompute$joint_exposure, expected)
  expect_equal(precompute$sparse_joint_cells, 4L)
})

test_that("binned fitting protects zero denominators", {
  dat <- data.frame(
    time = seq_len(8),
    z = seq(-1, 1, length.out = 8),
    status = rep(c(1, 0), 4)
  )
  fit <- sbf_fit_binning(
    dat,
    bandwidth = c(1, 1),
    family = "multiplicative",
    time_bins = 4,
    covariate_bins = 4,
    kernel = function(u) rep(0, length(u)),
    iterations = 2,
    convergence_tol = .Machine$double.xmin,
    warn_nonconvergence = FALSE,
    warn_diagnostics = FALSE
  )
  update_count <- fit$iterations_used * sum(lengths(fit$x.grid))

  expect_equal(fit$fit_diagnostics$denom_zero_count, update_count)
  expect_equal(fit$fit_diagnostics$denom_adjusted_count, update_count)
  expect_true(all(is.finite(unlist(fit$alpha_backfit))))
})

test_that("bin metadata records stored and empty bins", {
  dat <- data.frame(
    time = seq_len(12),
    z = c(0:5, 100:105),
    status = rep(c(1, 0), 6)
  )
  fit <- sbf_fit_binning(
    dat,
    bandwidth = c(4, 20),
    family = "additive",
    time_bins = 4,
    covariate_bins = 6,
    time_binning_method = "equal_width",
    covariate_binning_method = "equal_width",
    local_constant = TRUE,
    iterations = 2,
    warn_nonconvergence = FALSE,
    warn_diagnostics = FALSE
  )
  metadata <- fit$binning$metadata
  z_metadata <- metadata[metadata$dimension == "z", ]

  expect_equal(z_metadata$bins_requested, 6L)
  expect_equal(z_metadata$effective_bins, 6L)
  expect_equal(z_metadata$empty_bins, 4L)
  expect_equal(length(fit$x.grid[[2]]), z_metadata$effective_bins)
})

test_that("distinct-value grids retain all distinct values", {
  x <- rep(c(0, 1, 2, 3), each = 2)
  result <- sbfhazard:::.sbf_binning_bin_numeric_vector(
    x,
    bins = 5,
    rule = "equal_width",
    representative = "midpoint"
  )

  expect_equal(result$grid_mode, "unique")
  expect_equal(result$representatives, 0:3)
  expect_equal(result$empty_bins, 0L)
})

test_that("binning requires two distinct values", {
  expect_error(
    sbfhazard:::.sbf_binning_bin_numeric_vector(
      rep(1, 8),
      bins = 4,
      rule = "equal_width",
      representative = "midpoint"
    ),
    "at least two distinct values"
  )
})

test_that("log-time equal-width bins use the log1p scale", {
  x <- c(0.01, 0.02, 0.03, 0.1, 1, 3, 8, 12)
  midpoint <- sbfhazard:::.sbf_binning_bin_numeric_vector(
    x,
    bins = 4,
    rule = "equal_width_log",
    representative = "midpoint"
  )

  expect_equal(
    diff(log1p(midpoint$breaks)),
    rep(diff(log1p(midpoint$breaks))[1], 4)
  )
  expect_equal(
    midpoint$representatives,
    expm1(0.5 * (head(log1p(midpoint$breaks), -1) + tail(log1p(midpoint$breaks), -1)))
  )
  expect_error(
    sbfhazard:::.sbf_binning_bin_numeric_vector(
      c(-1, x),
      bins = 4,
      rule = "equal_width_log",
      representative = "midpoint"
    ),
    "nonnegative"
  )
})

test_that("log-time binning remains an internal diagnostic option", {
  dat <- sbf_simulate_data(
    n = 40,
    d = 3,
    family = "additive",
    seed = 701
  )

  expect_error(
    sbf_fit_binning(
      dat,
      bandwidth = rep(0.4, 3),
      family = "additive",
      time_bins = 5,
      covariate_bins = 5,
      time_binning_method = "equal_width_log",
      iterations = 2,
      warn_nonconvergence = FALSE,
      warn_diagnostics = FALSE
    ),
    "arg"
  )
})

test_that("public bin counts must be at least two", {
  dat <- sbf_simulate_data(
    n = 40,
    d = 3,
    family = "additive",
    seed = 702
  )

  expect_error(
    sbf_fit_binning(
      dat,
      bandwidth = rep(0.4, 3),
      family = "additive",
      time_bins = 1,
      covariate_bins = 5,
      iterations = 2,
      warn_nonconvergence = FALSE,
      warn_diagnostics = FALSE
    ),
    "time_bins"
  )
})
