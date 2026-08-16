#' Fit smooth backfitting hazard models
#'
#' Fit additive or multiplicative smooth backfitting hazard models, or their
#' binning-based approximations.
#'
#' @param data A numeric data frame with `time` and `status` columns, or a data
#'   frame used with `formula`. Non-formula fits use all columns other than
#'   `time` and `status` as fitted-scale features. The complete fitting data
#'   must contain at least one observed event with `status = 1`.
#' @param bandwidth Numeric bandwidth vector with one value per model component,
#'   including time.
#' @param family Model family: `"additive"` or `"multiplicative"`.
#' @param formula Optional formula of the form
#'   `survival::Surv(time, status) ~ features`.
#'   Right-hand side terms are converted to fitted-scale numeric features with
#'   `prodlim::model.design(maxOrder = 1, dropIntercept = TRUE)`.
#' @param integral_approx SBF-grid integration rule.
#' @param iterations Maximum number of completed backfitting iterations. One
#'   iteration updates every model component once.
#' @param kernel Kernel name or custom kernel function.
#' @param convergence_tol Positive convergence tolerance for backfitting.
#' @param local_constant For additive fits, whether to use local constant
#'   updates. Multiplicative SBF fits require `FALSE`.
#' @param identification Multiplicative component identification rule.
#'   The default `"sample_mean"` centers the log component over the training
#'   sample.
#' @param truth_functions Truth log-effect functions for truth-based
#'   multiplicative identification; one function per covariate component.
#' @param warn_nonconvergence Whether to warn when the fit does not reach the
#'   convergence tolerance within `iterations`.
#'
#' @return A fitted additive or multiplicative smooth backfitting object.
#'
#' @details
#' `sbf_fit()` evaluates the SBF estimators on observation-value grids.
#' `sbf_fit_binning()` returns speed-oriented binned approximations with bin
#' layout metadata in the fitted object.
#' Common fit fields include `alpha_backfit`, `x.grid`, `component_info`,
#' `fit_settings`, `fit_diagnostics`, and `feature_names`; binned fits also
#' include `binning`.
#' Multiplicative diagnostics include `component_floor_count`, the total number
#' of baseline and covariate raw-update values replaced by the fixed positive
#' floor across completed backfitting iterations.
#' Multiplicative `fit_settings` record this fixed value as `min_component`.
#' `iterations_used` records completed iterations and `iterations_max` records
#' the requested iteration limit for both SBF and binned fits.
#'
#' The multiplicative identification rules are: Sample-mean identification
#' (`"sample_mean"`) for training-sample log-mean centering, Integral
#' identification (`"integral"`) for fitted-grid integral mean one,
#' Single-point truth-based identification (`"origin"`) for midpoint
#' truth-anchor alignment, and Two-point truth-based identification (`"jasa"`)
#' for the source-script convention at -0.5 and 0.5. The `"jasa"` rule evaluates
#' the fitted component at the nearest grid points and evaluates truth at the
#' corresponding absolute distances, following the source script. The
#' two truth-based rules require `truth_functions` and are mainly for simulation
#' benchmarks.
#'
#' @examples
#' dat <- sbf_simulate_data(n = 50, d = 3, family = "additive", seed = 1)
#' fit <- sbf_fit(
#'   dat,
#'   bandwidth = c(0.35, 0.35, 0.35),
#'   family = "additive",
#'   local_constant = TRUE,
#'   iterations = 5,
#'   warn_nonconvergence = FALSE
#' )
#'
#' bin_fit <- sbf_fit_binning(
#'   dat,
#'   bandwidth = c(0.35, 0.35, 0.35),
#'   family = "additive",
#'   time_bins = 5,
#'   covariate_bins = 5,
#'   iterations = 5,
#'   warn_nonconvergence = FALSE,
#'   warn_diagnostics = FALSE
#' )
#'
#' mdat <- sbf_simulate_data(n = 50, d = 3, family = "multiplicative", seed = 2)
#' mfit <- sbf_fit(
#'   mdat,
#'   bandwidth = c(0.35, 0.35, 0.35),
#'   family = "multiplicative",
#'   iterations = 5,
#'   warn_nonconvergence = FALSE
#' )
#'
#' mbin_fit <- sbf_fit_binning(
#'   mdat,
#'   bandwidth = c(0.35, 0.35, 0.35),
#'   family = "multiplicative",
#'   time_bins = 5,
#'   covariate_bins = 5,
#'   iterations = 5,
#'   warn_nonconvergence = FALSE,
#'   warn_diagnostics = FALSE
#' )
#'
#' @rdname sbf_fit
#' @export

sbf_fit <- function(data,
                    bandwidth,
                    family = c("additive", "multiplicative"),
                    formula = NULL,
                    integral_approx = c("midd", "left", "right"),
                    iterations = 100L,
                    kernel = "epanechnikov",
                    convergence_tol = 0.001,
                    local_constant = FALSE,
                    identification = c("sample_mean", "integral", "origin", "jasa"),
                    truth_functions = NULL,
                    warn_nonconvergence = TRUE) {
  family <- match.arg(family)
  identification <- match.arg(identification)
  prepared <- .sbf_prepare_fit_inputs(data = data, formula = formula)
  data <- prepared$data
  d <- ncol(data) - 1L
  bandwidth <- .sbf_assert_numeric_vector(bandwidth, "bandwidth", expected_length = d, positive = TRUE)
  integral_approx <- match.arg(integral_approx)
  iterations <- .sbf_assert_scalar_whole_number(iterations, "iterations", min_value = 1L)
  convergence_tol <- .sbf_assert_scalar_numeric(convergence_tol, "convergence_tol", positive = TRUE)
  .sbf_kernel_spec(kernel, arg_name = "kernel")
  .sbf_assert_flag(local_constant, "local_constant")
  .sbf_assert_flag(warn_nonconvergence, "warn_nonconvergence")

  if (identical(family, "additive")) {
    if (!identical(identification, "sample_mean")) {
      .sbf_stop_bad_arg("identification", "\"sample_mean\" when family = \"additive\"")
    }
    if (!is.null(truth_functions)) {
      .sbf_stop_bad_arg("truth_functions", "NULL when family = \"additive\"")
    }
    fit <- .sbf_add_fit(
      data = data,
      formula = prepared$formula,
      bandwidth = bandwidth,
      integral_approx = integral_approx,
      iterations = iterations,
      kernel = kernel,
      convergence_tol = convergence_tol,
      kernel_correction = TRUE,
      local_constant = local_constant,
      wrong_taylor = FALSE,
      classic_backfit = FALSE,
      print_trace = FALSE,
      warn_nonconvergence = warn_nonconvergence
    )
    fit$formula_design <- prepared$formula_design
    return(fit)
  }

  if (!identical(local_constant, FALSE)) {
    .sbf_stop_bad_arg("local_constant", "FALSE when family = \"multiplicative\"")
  }
  if (identification %in% c("origin", "jasa")) {
    .sbf_assert_multiplicative_truth_functions(truth_functions, p = d - 1L)
  } else if (!is.null(truth_functions)) {
    .sbf_stop_bad_arg("truth_functions", "NULL unless identification is \"origin\" or \"jasa\"")
  }
  fit <- .sbf_mult_fit(
    data = data,
    bandwidth = bandwidth,
    formula = prepared$formula,
    integral_approx = integral_approx,
    iterations = iterations,
    kernel = kernel,
    convergence_tol = convergence_tol,
    identification = identification,
    truth_functions = truth_functions,
    warn_nonconvergence = warn_nonconvergence
  )
  fit$formula_design <- prepared$formula_design
  fit
}
