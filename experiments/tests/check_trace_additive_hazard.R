if (!requireNamespace("timereg", quietly = TRUE)) {
  cat("TRACE reproduction check skipped: package 'timereg' is not installed.\n")
} else {
  source(file.path("experiments", "script_utils.R"))
  project.root <- sbf_project_root()

  expect_true <- function(x, message) {
    if (!isTRUE(x)) {
      stop(message, call. = FALSE)
    }
  }

  out_dir <- file.path(
    project.root, "test_results", "runs", "smoke", "trace_additive_hazard"
  )
  unlink(out_dir, recursive = TRUE)

  script <- file.path(
    project.root, "experiments", "real_data", "trace_additive_hazard.R"
  )
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = shQuote(c(script, out_dir))
  )
  expect_true(status == 0L, "TRACE analysis script failed.")

  files <- c(
    cohort = file.path(out_dir, "cohort_summary.csv"),
    fit_summary = file.path(out_dir, "fit_summary.csv"),
    components = file.path(out_dir, "component_curves.csv"),
    component_mse = file.path(out_dir, "component_mse.csv"),
    plot = file.path(out_dir, "trace_additive_curves.png")
  )
  expect_true(all(file.exists(files)), "TRACE did not write every required output.")

  cohort <- utils::read.csv(files[["cohort"]])
  fit_summary <- utils::read.csv(files[["fit_summary"]])
  component_curves <- utils::read.csv(files[["components"]])
  component_mse <- utils::read.csv(files[["component_mse"]])

  utils::data("TRACE", package = "timereg", envir = environment())
  raw <- get("TRACE", envir = environment())
  observed_time <- as.numeric(raw$time)
  age <- as.numeric(raw$age)
  wmi <- as.numeric(raw$wmi)
  raw_status <- as.numeric(raw$status)
  keep <- complete.cases(observed_time, age, wmi, raw_status) &
    is.finite(observed_time) & is.finite(age) &
    is.finite(wmi) & is.finite(raw_status) &
    age >= 40 & age < 85
  censoring_time <- pmin(5, 85 - age)
  retained_time <- pmin(observed_time, censoring_time)
  retained_status <- as.integer(
    raw_status != 0 & observed_time <= censoring_time
  )
  keep <- keep & is.finite(retained_time) & retained_time > 0

  expect_true(sum(keep) == 1786L, "TRACE row count changed.")
  expect_true(
    sum(retained_status[keep]) == 725L,
    "TRACE event count changed."
  )
  expect_true(cohort$n == 1786L, "Reported TRACE row count changed.")
  expect_true(cohort$n_events == 725L, "Reported TRACE event count changed.")
  expect_true(
    all(fit_summary$converged),
    "A TRACE fit did not converge."
  )
  expect_true(
    all(is.finite(component_curves$estimate)),
    "TRACE component curves are not finite."
  )
  expect_true(
    all(is.finite(component_mse$mse)),
    "TRACE component MSE is not finite."
  )
  expect_true(
    setequal(fit_summary$method, c("sbf", "pairwise_binned")),
    "TRACE fit summary has unexpected method identifiers."
  )
  expect_true(
    setequal(component_curves$method, c("sbf", "pairwise_binned")),
    "TRACE component curves have unexpected method identifiers."
  )

  for (component_name in component_mse$component) {
    curves <- component_curves[
      component_curves$component == component_name, , drop = FALSE
    ]
    sbf_values <- curves$estimate[curves$method == "sbf"]
    binned_values <- curves$estimate[
      curves$method == "pairwise_binned"
    ]
    reported_mse <- component_mse$mse[
      component_mse$component == component_name
    ]
    expect_true(
      isTRUE(all.equal(
        reported_mse,
        mean((binned_values - sbf_values)^2)
      )),
      sprintf("TRACE %s MSE changed.", component_name)
    )
  }

  cat("TRACE reproduction check passed.\n")
}
