test_that("simulation returns fit-ready data", {
  for (family in c("additive", "multiplicative")) {
    dat <- sbf_simulate_data(
      n = 30,
      d = 3,
      family = family,
      model = 1,
      seed = 100
    )

    expect_s3_class(dat, "data.frame")
    expect_equal(nrow(dat), 30L)
    expect_true(all(c("time", "status", "V1", "V2") %in% names(dat)))
    expect_true(all(is.finite(dat$time)))
    expect_true(all(dat$status %in% c(0L, 1L)))
  }
})

test_that("Makeham simulation returns additive and multiplicative data", {
  skip_if_not_installed("VGAM")

  for (family in c("additive", "multiplicative")) {
    dat <- sbf_simulate_data(
      n = 20,
      d = 3,
      family = family,
      model = 2,
      seed = 200
    )

    expect_s3_class(dat, "data.frame")
    expect_equal(nrow(dat), 20L)
    expect_true(all(is.finite(dat$time)))
  }
})

test_that("the four fitting paths return finite component predictions", {
  additive_data <- sbf_simulate_data(
    n = 60,
    d = 3,
    family = "additive",
    seed = 301
  )
  multiplicative_data <- sbf_simulate_data(
    n = 60,
    d = 3,
    family = "multiplicative",
    seed = 302
  )
  bandwidth <- rep(0.4, 3)

  fits <- list(
    sbf_fit(
      additive_data,
      bandwidth = bandwidth,
      family = "additive",
      local_constant = TRUE,
      iterations = 4,
      warn_nonconvergence = FALSE
    ),
    sbf_fit_binning(
      additive_data,
      bandwidth = bandwidth,
      family = "additive",
      time_bins = 5,
      covariate_bins = 5,
      local_constant = TRUE,
      iterations = 4,
      warn_nonconvergence = FALSE,
      warn_diagnostics = FALSE
    ),
    sbf_fit(
      multiplicative_data,
      bandwidth = bandwidth,
      family = "multiplicative",
      iterations = 4,
      warn_nonconvergence = FALSE
    ),
    sbf_fit_binning(
      multiplicative_data,
      bandwidth = bandwidth,
      family = "multiplicative",
      time_bins = 5,
      covariate_bins = 5,
      iterations = 4,
      warn_nonconvergence = FALSE,
      warn_diagnostics = FALSE
    )
  )
  expected_classes <- c(
    "sbf_additive_fit",
    "sbf_additive_pairwise_binning_fit",
    "sbf_multiplicative_fit",
    "sbf_multiplicative_sparse_binning_fit"
  )

  for (i in seq_along(fits)) {
    expect_s3_class(fits[[i]], expected_classes[i])
    prediction <- sbf_predict(
      fits[[i]],
      type = "component",
      component = 2,
      xout = 0,
      warn_diagnostics = FALSE
    )
    expect_true(is.finite(prediction))
  }
})
