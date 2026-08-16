library(sbfhazard)

dat <- sbf_simulate_data(
  n = 160,
  d = 3,
  rho = 0.2,
  family = "multiplicative",
  seed = 201
)
model <- survival::Surv(time, status) ~ V1 + V2
bandwidth <- rep(0.4, 3)

sbf_fit_obj <- sbf_fit(
  data = dat,
  bandwidth = bandwidth,
  family = "multiplicative",
  formula = model,
  iterations = 100
)
binned <- sbf_fit_binning(
  data = dat,
  bandwidth = bandwidth,
  family = "multiplicative",
  formula = model,
  time_bins = 12,
  covariate_bins = 12,
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

print(head(data.frame(xout, sbf_component, binned_component)))

if (interactive()) {
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
}
