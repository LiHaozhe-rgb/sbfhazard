#' @param time_bins,covariate_bins Binning resolution.
#' @param time_binning_method,covariate_binning_method Binning method for time
#'   and covariate components.
#' @param representative Bin representative: `"midpoint"` or `"mean"`.
#' @param warn_diagnostics Whether to emit diagnostic warnings for numerical or
#'   binning issues.
#'
#' @rdname sbf_fit
#' @export

sbf_fit_binning <- function(data,
                            bandwidth,
                            time_bins,
                            covariate_bins,
                            family = c("additive", "multiplicative"),
                            time_binning_method = c("quantile", "equal_width"),
                            covariate_binning_method = c("quantile", "equal_width"),
                            representative = c("midpoint", "mean"),
                            formula = NULL,
                            local_constant = NULL,
                            iterations = 100L,
                            kernel = "epanechnikov",
                            convergence_tol = 0.001,
                            identification = c("sample_mean", "integral", "origin", "jasa"),
                            truth_functions = NULL,
                            warn_nonconvergence = TRUE,
                            warn_diagnostics = TRUE) {
  family <- match.arg(family)
  identification <- match.arg(identification)
  prepared <- .sbf_prepare_fit_inputs(data = data, formula = formula)
  data <- prepared$data
  d <- ncol(data) - 1L
  bandwidth <- .sbf_assert_numeric_vector(bandwidth, "bandwidth", expected_length = d, positive = TRUE)
  iterations <- .sbf_assert_scalar_whole_number(iterations, "iterations", min_value = 1L)
  convergence_tol <- .sbf_assert_scalar_numeric(convergence_tol, "convergence_tol", positive = TRUE)
  .sbf_assert_flag(warn_nonconvergence, "warn_nonconvergence")
  .sbf_assert_flag(warn_diagnostics, "warn_diagnostics")
  time_binning_method <- match.arg(time_binning_method)
  covariate_binning_method <- match.arg(covariate_binning_method)
  representative <- match.arg(representative)
  time_bins <- .sbf_assert_scalar_whole_number(time_bins, "time_bins", min_value = 2L)
  covariate_bins <- .sbf_assert_scalar_whole_number(covariate_bins, "covariate_bins", min_value = 2L)
  .sbf_kernel_spec(kernel, arg_name = "kernel")

  if (identical(family, "additive")) {
    if (!identical(identification, "sample_mean")) {
      .sbf_stop_bad_arg("identification", "\"sample_mean\" when family = \"additive\"")
    }
    if (!is.null(truth_functions)) {
      .sbf_stop_bad_arg("truth_functions", "NULL when family = \"additive\"")
    }
    local_constant <- if (is.null(local_constant)) TRUE else .sbf_assert_flag(local_constant, "local_constant")
    fit <- .sbf_add_pairwise_fit(
      data = data,
      formula = prepared$formula,
      bandwidth = bandwidth,
      time_bins = time_bins,
      covariate_bins = covariate_bins,
      time_binning_method = time_binning_method,
      covariate_binning_method = covariate_binning_method,
      representative = representative,
      local_constant = local_constant,
      kernel = kernel,
      kernel_correction = TRUE,
      iterations = iterations,
      convergence_tol = convergence_tol,
      warn_nonconvergence = warn_nonconvergence,
      warn_diagnostics = warn_diagnostics
    )
    fit$formula_design <- prepared$formula_design
    return(fit)
  }

  if (!is.null(local_constant)) {
    .sbf_stop_bad_arg("local_constant", "NULL when family = \"multiplicative\"")
  }
  if (identification %in% c("origin", "jasa")) {
    .sbf_assert_multiplicative_truth_functions(truth_functions, p = d - 1L)
  } else if (!is.null(truth_functions)) {
    .sbf_stop_bad_arg("truth_functions", "NULL unless identification is \"origin\" or \"jasa\"")
  }
  fit <- .sbf_mult_sparse_fit(
    data = data,
    bandwidth = bandwidth,
    time_bins = time_bins,
    covariate_bins = covariate_bins,
    time_binning_method = time_binning_method,
    covariate_binning_method = covariate_binning_method,
    representative = representative,
    formula = prepared$formula,
    iterations = iterations,
    kernel = kernel,
    convergence_tol = convergence_tol,
    min_component = 1e-8,
    identification = identification,
    truth_functions = truth_functions,
    warn_nonconvergence = warn_nonconvergence,
    warn_diagnostics = warn_diagnostics
  )
  fit$formula_design <- prepared$formula_design
  fit
}
