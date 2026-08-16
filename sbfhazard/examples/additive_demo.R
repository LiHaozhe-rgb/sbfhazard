library(sbfhazard)

dat <- sbf_simulate_data(
  n = 160,
  d = 3,
  rho = 0.2,
  family = "additive",
  seed = 101
)
model <- survival::Surv(time, status) ~ V1 + V2
bandwidth <- rep(0.4, 3)

sbf_fit_obj <- sbf_fit(
  data = dat,
  bandwidth = bandwidth,
  family = "additive",
  formula = model,
  local_constant = TRUE,
  iterations = 100
)
binned <- sbf_fit_binning(
  data = dat,
  bandwidth = bandwidth,
  family = "additive",
  formula = model,
  time_bins = 12,
  covariate_bins = 12,
  local_constant = TRUE,
  iterations = 100
)

xout <- seq(min(dat$V1), max(dat$V1), length.out = 100)
sbf_component <- sbf_predict(
  sbf_fit_obj,
  type = "component",
  component = 2,
  xout = xout
)
binned_component <- sbf_predict(
  binned,
  type = "component",
  component = 2,
  xout = xout
)

times <- seq(min(dat$time), stats::quantile(dat$time, 0.9), length.out = 100)
profile <- data.frame(V1 = 0, V2 = 0)
sbf_survival <- sbf_predict(
  sbf_fit_obj,
  type = "survival",
  data = profile,
  times = times
)$survival[1, ]
binned_survival <- sbf_predict(
  binned,
  type = "survival",
  data = profile,
  times = times
)$survival[1, ]

print(head(data.frame(xout, sbf_component, binned_component)))
print(head(data.frame(times, sbf_survival, binned_survival)))

if (interactive()) {
  old_par <- par(mfrow = c(1, 2))

  matplot(
    xout,
    cbind(sbf_component, binned_component),
    type = "l",
    lty = c(1, 2),
    col = c("black", "#0072B2"),
    xlab = "V1",
    ylab = "Component"
  )
  legend(
    "topright",
    c("SBF", "Binned"),
    lty = c(1, 2),
    col = c("black", "#0072B2"),
    bty = "n"
  )
  matplot(
    times,
    cbind(sbf_survival, binned_survival),
    type = "l",
    lty = c(1, 2),
    col = c("black", "#0072B2"),
    xlab = "Time",
    ylab = "Survival"
  )
  par(old_par)
}
