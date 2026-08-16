test_that("named time and status columns may appear in any order", {
  dat <- sbf_simulate_data(
    n = 50,
    d = 3,
    family = "additive",
    seed = 501
  )
  dat <- dat[, c("V2", "status", "time", "V1")]

  fit <- sbf_fit(
    dat,
    bandwidth = rep(0.4, 3),
    family = "additive",
    local_constant = TRUE,
    iterations = 3,
    warn_nonconvergence = FALSE
  )

  expect_equal(fit$feature_names, c("V2", "V1"))
  expect_true(is.finite(sbf_predict(
    fit,
    type = "component",
    component = 2,
    xout = 0,
    warn_diagnostics = FALSE
  )))
})

test_that("formula fits reuse their feature transformation for prediction", {
  raw <- sbf_simulate_data(
    n = 60,
    d = 3,
    family = "additive",
    seed = 502
  )
  dat <- data.frame(
    observed_time = raw$time,
    event = raw$status,
    z1 = raw$V1,
    z2 = raw$V2
  )
  fit <- sbf_fit(
    dat,
    formula = survival::Surv(observed_time, event) ~ log(z1 + 2) + I(z2^2),
    bandwidth = rep(0.45, 3),
    family = "additive",
    local_constant = TRUE,
    iterations = 3,
    warn_nonconvergence = FALSE
  )

  raw_prediction <- sbf_predict(
    fit,
    type = "hazard",
    data = data.frame(z1 = 0, z2 = 0),
    times = c(0.25, 0.5),
    warn_diagnostics = FALSE
  )
  transformed_prediction <- sbf_predict(
    fit,
    type = "hazard",
    data = matrix(c(log(2), 0), nrow = 1),
    times = c(0.25, 0.5),
    warn_diagnostics = FALSE
  )

  expect_equal(raw_prediction$hazard, transformed_prediction$hazard)
  expect_error(
    sbf_predict(
      fit,
      type = "hazard",
      data = data.frame(z1 = 0, other = 0),
      times = 0.5,
      warn_diagnostics = FALSE
    ),
    "columns"
  )
})

test_that("formula prediction reuses training factor levels", {
  raw <- sbf_simulate_data(
    n = 60,
    d = 2,
    family = "additive",
    seed = 503
  )
  dat <- data.frame(
    time = raw$time,
    status = raw$status,
    group = factor(rep(c("A", "B", "C"), length.out = nrow(raw)))
  )
  fit <- sbf_fit(
    dat,
    formula = survival::Surv(time, status) ~ group,
    bandwidth = rep(0.5, 3),
    family = "additive",
    local_constant = TRUE,
    iterations = 3,
    warn_nonconvergence = FALSE
  )

  raw_prediction <- sbf_predict(
    fit,
    type = "hazard",
    data = data.frame(group = "B"),
    times = 0.5,
    warn_diagnostics = FALSE
  )
  design_prediction <- sbf_predict(
    fit,
    type = "hazard",
    data = data.frame(groupB = 1, groupC = 0),
    times = 0.5,
    warn_diagnostics = FALSE
  )

  expect_equal(raw_prediction$hazard, design_prediction$hazard)
  expect_error(
    sbf_predict(
      fit,
      type = "hazard",
      data = data.frame(group = "D"),
      times = 0.5,
      warn_diagnostics = FALSE
    ),
    "new level"
  )
})

test_that("formula prediction reuses training contrasts", {
  raw <- sbf_simulate_data(
    n = 60,
    d = 2,
    family = "additive",
    seed = 504
  )
  dat <- data.frame(
    time = raw$time,
    status = raw$status,
    group = factor(rep(c("A", "B", "C"), length.out = nrow(raw)))
  )

  old_contrasts <- getOption("contrasts")
  on.exit(options(contrasts = old_contrasts), add = TRUE)
  options(contrasts = c("contr.sum", "contr.poly"))
  fit <- sbf_fit(
    dat,
    formula = survival::Surv(time, status) ~ group,
    bandwidth = rep(0.5, 3),
    family = "additive",
    local_constant = TRUE,
    iterations = 3,
    warn_nonconvergence = FALSE
  )

  options(contrasts = c("contr.helmert", "contr.poly"))
  raw_prediction <- sbf_predict(
    fit,
    type = "hazard",
    data = data.frame(group = "B"),
    times = 0.5,
    warn_diagnostics = FALSE
  )
  design_prediction <- sbf_predict(
    fit,
    type = "hazard",
    data = fit$data[2, fit$feature_names, drop = FALSE],
    times = 0.5,
    warn_diagnostics = FALSE
  )

  expect_equal(raw_prediction$hazard, design_prediction$hazard)
})

test_that("fitting requires at least one observed event", {
  dat <- sbf_simulate_data(
    n = 40,
    d = 3,
    family = "additive",
    seed = 505
  )
  dat$status <- 0L

  expect_error(
    sbf_fit(
      dat,
      bandwidth = rep(0.4, 3),
      family = "additive",
      local_constant = TRUE,
      warn_nonconvergence = FALSE
    ),
    "at least one observed event"
  )
})

test_that("invalid fit inputs fail at the public boundary", {
  dat <- sbf_simulate_data(
    n = 40,
    d = 3,
    family = "additive",
    seed = 506
  )
  nonnumeric <- dat
  nonnumeric$V1 <- as.character(nonnumeric$V1)

  expect_error(
    sbf_fit(as.matrix(dat), bandwidth = rep(0.4, 3), family = "additive"),
    "data.frame"
  )
  expect_error(
    sbf_fit(nonnumeric, bandwidth = rep(0.4, 3), family = "additive"),
    "numeric"
  )
  expect_error(
    sbf_fit(dat, bandwidth = rep(0.4, 2), family = "additive"),
    "bandwidth"
  )
})
