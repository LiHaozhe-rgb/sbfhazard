source(file.path("experiments", "script_utils.R"))
project.root <- sbf_project_root()
sbf_load_package(project.root)

if (!requireNamespace("timereg", quietly = TRUE)) {
  stop("Package 'timereg' is required for the TRACE analysis.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0L) {
  out_dir <- args[[1L]]
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
} else {
  out_dir <- sbf_results_dir(
    project.root, "runs", "experiments", "real_data", "trace_additive_hazard"
  )
}

bandwidth <- c(time = 0.1, age = 15, wmi = 0.8)
iterations <- 100L
convergence_tol <- 0.001
kernel <- "epanechnikov"
time_bins <- 100L
covariate_bins <- 30L

utils::data("TRACE", package = "timereg", envir = environment())
raw <- get("TRACE", envir = environment())

# Follow-up starts at AMI. Patients younger than 40 are excluded rather than
# entering the risk set later. Follow-up is censored at five years or age 85.
observed_time <- as.numeric(raw$time)
age <- as.numeric(raw$age)
wmi <- as.numeric(raw$wmi)
status <- as.numeric(raw$status)

keep <- complete.cases(observed_time, age, wmi, status) &
  is.finite(observed_time) & is.finite(age) &
  is.finite(wmi) & is.finite(status) &
  age >= 40 & age < 85
censoring_time <- pmin(5, 85 - age)
retained_time <- pmin(observed_time, censoring_time)
event <- status != 0 & observed_time <= censoring_time
keep <- keep & is.finite(retained_time) & retained_time > 0

dat <- data.frame(
  time = retained_time[keep],
  status = as.integer(event[keep]),
  age = age[keep],
  wmi = wmi[keep]
)
stopifnot(nrow(dat) == 1786L, sum(dat$status) == 725L)
cat(sprintf("TRACE analysis data: n=%d, events=%d.\n", nrow(dat), sum(dat$status)))

start <- proc.time()[[3L]]
sbf_fit_obj <- sbf_fit(
  data = dat,
  formula = survival::Surv(time, status) ~ age + wmi,
  family = "additive",
  bandwidth = bandwidth,
  integral_approx = "midd",
  local_constant = FALSE,
  iterations = iterations,
  convergence_tol = convergence_tol,
  kernel = kernel,
  warn_nonconvergence = FALSE
)
sbf_runtime <- proc.time()[[3L]] - start

start <- proc.time()[[3L]]
binned_fit <- sbf_fit_binning(
  data = dat,
  formula = survival::Surv(time, status) ~ age + wmi,
  family = "additive",
  bandwidth = bandwidth,
  time_bins = time_bins,
  covariate_bins = covariate_bins,
  time_binning_method = "equal_width",
  covariate_binning_method = "equal_width",
  representative = "midpoint",
  local_constant = FALSE,
  iterations = iterations,
  convergence_tol = convergence_tol,
  kernel = kernel,
  warn_nonconvergence = FALSE,
  warn_diagnostics = FALSE
)
binned_runtime <- proc.time()[[3L]] - start

fits <- list(sbf = sbf_fit_obj, pairwise_binned = binned_fit)
stopifnot(all(vapply(fits, function(fit) isTRUE(fit$converged), logical(1))))
runtimes <- c(
  sbf = as.numeric(sbf_runtime),
  pairwise_binned = as.numeric(binned_runtime)
)

component_grids <- list(
  baseline = seq(0, 5, length.out = 200),
  age = seq(40, 85, length.out = 200),
  wmi = seq(0.3, 3, length.out = 200)
)
component_ids <- c(baseline = 1L, age = 2L, wmi = 3L)
curve_rows <- list()
diagnostic_rows <- list()
row <- 1L

for (method in names(fits)) {
  for (component_name in names(component_grids)) {
    prediction <- sbf_predict(
      fits[[method]],
      type = "component",
      component = component_ids[[component_name]],
      xout = component_grids[[component_name]],
      prediction = "formula",
      warn_diagnostics = FALSE
    )
    prediction_diagnostics <- attr(prediction, "prediction_diagnostics")

    curve_rows[[row]] <- data.frame(
      method = method,
      component = component_name,
      grid_value = component_grids[[component_name]],
      estimate = as.numeric(prediction)
    )
    diagnostic_rows[[row]] <- data.frame(
      method = method,
      component = component_name,
      zero_kernel_support = sum(
        as.numeric(prediction_diagnostics$kernel_support_zero_count),
        na.rm = TRUE
      ),
      fallback_count = sum(
        as.numeric(prediction_diagnostics$bin_step_fallback_count),
        na.rm = TRUE
      )
    )
    row <- row + 1L
  }
}

component_curves <- do.call(rbind, curve_rows)
prediction_diagnostics <- do.call(rbind, diagnostic_rows)

component_mse <- data.frame(
  component = names(component_grids),
  mse = NA_real_
)
for (i in seq_len(nrow(component_mse))) {
  component_name <- component_mse$component[[i]]
  sbf_values <- component_curves$estimate[
    component_curves$method == "sbf" &
      component_curves$component == component_name
  ]
  binned_values <- component_curves$estimate[
    component_curves$method == "pairwise_binned" &
      component_curves$component == component_name
  ]
  component_mse$mse[[i]] <- mean((binned_values - sbf_values)^2)
}

followup_quantiles <- stats::quantile(
  dat$time, c(0.25, 0.75), names = FALSE, type = 8
)
age_quantiles <- stats::quantile(
  dat$age, c(0.25, 0.75), names = FALSE, type = 8
)
wmi_quantiles <- stats::quantile(
  dat$wmi, c(0.25, 0.75), names = FALSE, type = 8
)
cohort_summary <- data.frame(
  n = nrow(dat),
  n_events = sum(dat$status),
  censoring_percent = 100 * mean(dat$status == 0),
  followup_min = min(dat$time),
  followup_q25 = followup_quantiles[[1L]],
  followup_median = stats::median(dat$time),
  followup_q75 = followup_quantiles[[2L]],
  followup_max = max(dat$time),
  followup_distinct = length(unique(dat$time)),
  followup_at_five = sum(dat$time == 5),
  age_min = min(dat$age),
  age_q25 = age_quantiles[[1L]],
  age_median = stats::median(dat$age),
  age_q75 = age_quantiles[[2L]],
  age_max = max(dat$age),
  age_distinct = length(unique(dat$age)),
  wmi_min = min(dat$wmi),
  wmi_q25 = wmi_quantiles[[1L]],
  wmi_median = stats::median(dat$wmi),
  wmi_q75 = wmi_quantiles[[2L]],
  wmi_max = max(dat$wmi),
  wmi_distinct = length(unique(dat$wmi))
)

fit_summary_rows <- list()
for (method in names(fits)) {
  fit <- fits[[method]]
  method_diagnostics <- prediction_diagnostics[
    prediction_diagnostics$method == method, , drop = FALSE
  ]
  is_binned <- method == "pairwise_binned"

  fit_summary_rows[[method]] <- data.frame(
    method = method,
    runtime_sec = runtimes[[method]],
    speedup = runtimes[["sbf"]] / runtimes[[method]],
    iterations_used = fit$iterations_used,
    converged = fit$converged,
    final_delta = fit$final_delta,
    denominator_adjustments = sum(
      as.numeric(fit$fit_diagnostics$denom_adjusted_count),
      na.rm = TRUE
    ),
    kernel_normalizer_adjustments = sum(
      as.numeric(fit$fit_diagnostics$kernel_norm_adjusted_count),
      na.rm = TRUE
    ),
    component_zero_kernel_support = sum(
      method_diagnostics$zero_kernel_support
    ),
    component_fallback_count = sum(method_diagnostics$fallback_count),
    time_grid_size = length(fit$x.grid[[1L]]),
    age_grid_size = length(fit$x.grid[[2L]]),
    wmi_grid_size = length(fit$x.grid[[3L]]),
    time_bins_requested = if (is_binned) time_bins else NA_integer_,
    covariate_bins_requested = if (is_binned) {
      covariate_bins
    } else {
      NA_integer_
    }
  )
}
fit_summary <- do.call(rbind, fit_summary_rows)
rownames(fit_summary) <- NULL

cohort_file <- sbf_write_csv(
  cohort_summary, file.path(out_dir, "cohort_summary.csv")
)
fit_summary_file <- sbf_write_csv(
  fit_summary, file.path(out_dir, "fit_summary.csv")
)
component_curves_file <- sbf_write_csv(
  component_curves, file.path(out_dir, "component_curves.csv")
)
component_mse_file <- sbf_write_csv(
  component_mse, file.path(out_dir, "component_mse.csv")
)

plot_file <- file.path(out_dir, "trace_additive_curves.png")
method_color <- c(sbf = "#111111", pairwise_binned = "#0072B2")
method_line <- c(sbf = 1, pairwise_binned = 2)
method_label <- c(sbf = "SBF", pairwise_binned = "Binned")
component_title <- c(
  baseline = "Baseline",
  age = "Age component",
  wmi = "WMI component"
)
component_xlab <- c(
  baseline = "Time (years)",
  age = "Age at AMI (years)",
  wmi = "WMI"
)

grDevices::png(plot_file, width = 3000, height = 1000, res = 220)
par(
  mfrow = c(1, 3),
  mar = c(5.6, 5.8, 4.2, 1.2),
  oma = c(0, 0, 3.8, 0),
  cex.axis = 1.9,
  cex.lab = 2.2,
  cex.main = 2.1,
  mgp = c(3.7, 1.1, 0)
)

for (component_name in names(component_title)) {
  values <- component_curves[
    component_curves$component == component_name, , drop = FALSE
  ]
  plot(
    range(values$grid_value), range(values$estimate), type = "n",
    xlab = component_xlab[[component_name]],
    ylab = "Component estimate",
    main = component_title[[component_name]]
  )
  abline(h = 0, col = "#BBBBBB", lty = 3)

  for (method in names(method_label)) {
    curve <- values[values$method == method, , drop = FALSE]
    lines(
      curve$grid_value,
      curve$estimate,
      col = method_color[[method]],
      lty = method_line[[method]],
      lwd = 3
    )
  }

  if (component_name == "wmi") {
    legend(
      "topright",
      legend = method_label,
      col = method_color,
      lty = method_line,
      lwd = 3,
      cex = 1.75,
      bty = "n"
    )
  }
}
mtext(
  "TRACE additive local linear fits",
  outer = TRUE,
  font = 2,
  cex = 2.05
)
invisible(grDevices::dev.off())
plot_file <- normalizePath(plot_file, winslash = "/", mustWork = TRUE)

cat(sprintf("Cohort summary: %s\n", cohort_file))
cat(sprintf("Fit summary: %s\n", fit_summary_file))
cat(sprintf("Component curves: %s\n", component_curves_file))
cat(sprintf("Component MSE: %s\n", component_mse_file))
cat(sprintf("Plot: %s\n", plot_file))
