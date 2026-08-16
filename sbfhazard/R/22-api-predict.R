#' Predict from smooth backfitting hazard models
#'
#' Predict components, hazards, or survival curves from fitted smooth
#' backfitting hazard models.
#'
#' @param result A fit object returned by `sbf_fit()` or
#'   `sbf_fit_binning()`.
#' @param type Prediction type: `"component"`, `"hazard"`, or `"survival"`.
#' @param component Component index to evaluate when `type = "component"`.
#' @param xout Evaluation points when `type = "component"`. These are on the
#'   fitted feature scale for formula fits.
#' @param data New covariate data when predicting hazards or survival. For
#'   formula fits, a data frame containing the raw right-hand-side variables is
#'   transformed automatically using the factor levels and contrasts actually
#'   used by the training design; unseen factor levels are rejected.
#'   Matrix/vector inputs are treated as fitted-scale covariates.
#' @param times Optional evaluation times for hazard or survival prediction.
#' @param clamp_nonnegative Whether to replace non-finite hazards and hazards
#'   below `min_hazard` by `min_hazard`. If `FALSE`, finite hazards are left
#'   unchanged and non-finite hazards raise an error.
#' @param min_hazard Positive lower bound used when clamping hazards.
#' @param prediction SBF-fit prediction method. `"formula"` evaluates the
#'   fitted estimating equations at new points. `"interpolation"` linearly
#'   interpolates fitted SBF-grid component values and uses the nearest fitted
#'   value outside the grid. Binning fits only support `"formula"`.
#' @param warn_diagnostics Whether to emit diagnostic warnings.
#'
#' @return A component vector, hazard prediction object, or survival prediction
#'   object.
#'
#' @details
#' `component = 1` is the baseline/time component. Components `2` and above
#' refer to the fitted feature components in `component_info`. Component
#' prediction `xout` is always on the fitted feature scale.
#'
#' Formula and interpolation survival prediction use an integration grid formed
#' from the fitted support start, the baseline grid, and the requested times,
#' then apply trapezoidal integration. Interpolation continues with the
#' nearest-boundary hazard after the last fitted baseline point. At or before
#' the support start, survival is one.
#'
#' @examples
#' dat <- sbf_simulate_data(n = 50, d = 3, family = "additive", seed = 3)
#' fit <- sbf_fit(
#'   dat,
#'   bandwidth = c(0.35, 0.35, 0.35),
#'   family = "additive",
#'   formula = survival::Surv(time, status) ~ V1 + V2,
#'   local_constant = TRUE,
#'   iterations = 5,
#'   warn_nonconvergence = FALSE
#' )
#'
#' sbf_predict(fit, type = "component", component = 2, xout = c(-0.2, 0, 0.2))
#' sbf_predict(
#'   fit,
#'   type = "hazard",
#'   data = data.frame(V1 = 0, V2 = 0),
#'   times = c(0.25, 0.5),
#'   warn_diagnostics = FALSE
#' )
#' sbf_predict(
#'   fit,
#'   type = "survival",
#'   data = data.frame(V1 = 0, V2 = 0),
#'   times = c(0.25, 0.5),
#'   warn_diagnostics = FALSE
#' )
#'
#' @export

sbf_predict <- function(result,
                        type = c("component", "hazard", "survival"),
                        component = NULL,
                        xout = NULL,
                        data = NULL,
                        times = NULL,
                        clamp_nonnegative = TRUE,
                        min_hazard = 1e-8,
                        prediction = c("formula", "interpolation"),
                        warn_diagnostics = TRUE) {
  type <- match.arg(type)
  prediction <- match.arg(prediction)
  if (!inherits(result, "sbf_additive_fit") && !inherits(result, "sbf_multiplicative_fit")) {
    .sbf_stop_bad_arg("result", "an sbfhazard fit object")
  }
  .sbf_assert_required_fields(result, fields = c("component_info", "x.grid"))
  if (!is.data.frame(result$component_info) || nrow(result$component_info) < 1L) {
    .sbf_stop_bad_arg("result", "an sbfhazard fit object with component_info")
  }

  n_components <- nrow(result$component_info)
  p <- n_components - 1L
  is_binning <- inherits(result, "sbf_additive_pairwise_binning_fit") ||
    inherits(result, "sbf_multiplicative_sparse_binning_fit")

  if (is_binning && identical(prediction, "interpolation")) {
    .sbf_stop_bad_arg("prediction", "\"formula\" for binned fits")
  }

  .sbf_assert_flag(warn_diagnostics, "warn_diagnostics")

  if (identical(type, "component")) {
    if (is.null(component)) {
      .sbf_stop_bad_arg("component", "a component index when type = \"component\"")
    }
    if (is.null(xout)) {
      .sbf_stop_bad_arg("xout", "evaluation points when type = \"component\"")
    }
    component <- .sbf_assert_scalar_whole_number(component, "component", min_value = 1L)
    if (component > n_components) {
      .sbf_stop_bad_arg(
        "component",
        sprintf("an index between 1 and %d for this result", n_components),
        actual = as.character(component)
      )
    }
    
    xout <- .sbf_assert_numeric_vector(xout, "xout", min_length = 1L)
    if (inherits(result, "sbf_additive_fit")) {
      if (inherits(result, "sbf_additive_pairwise_binning_fit")) {
        return(.sbf_add_pairwise_predict_component(
          result = result,
          component = component,
          xout = xout,
          warn_diagnostics = warn_diagnostics
        ))
      }
      return(.sbf_add_predict_component(
        result = result,
        component = component,
        xout = xout,
        prediction = prediction,
        warn_diagnostics = warn_diagnostics
      ))
    }
    if (inherits(result, "sbf_multiplicative_sparse_binning_fit")) {
      return(.sbf_mult_sparse_predict_component(
        result = result,
        component = component,
        xout = xout,
        warn_diagnostics = warn_diagnostics
      ))
    }
    return(.sbf_mult_predict_component(
      result = result,
      component = component,
      xout = xout,
      min_component = 1e-8,
      prediction = prediction,
      warn_diagnostics = warn_diagnostics
    ))
  }

  if (is.null(data)) {
    .sbf_stop_bad_arg("data", sprintf("new covariate data when type = \"%s\"", type))
  }
  .sbf_assert_flag(clamp_nonnegative, "clamp_nonnegative")
  min_hazard <- .sbf_assert_scalar_numeric(min_hazard, "min_hazard", positive = TRUE)
  if (!is.null(times)) {
    times <- .sbf_assert_numeric_vector(times, "times", min_length = 1L)
  }
  cov_data <- .sbf_prepare_prediction_feature_matrix(
    result = result,
    data = data,
    p = p
  )

  if (identical(type, "hazard")) {
    eval_times <- if (is.null(times)) as.numeric(result$x.grid[[1L]]) else times
    if (inherits(result, "sbf_additive_fit")) {
      if (inherits(result, "sbf_additive_pairwise_binning_fit")) {
        return(.sbf_add_pairwise_predict_hazard(
          result = result,
          cov_data = cov_data,
          eval_times = eval_times,
          clamp_nonnegative = clamp_nonnegative,
          min_hazard = min_hazard,
          warn_diagnostics = warn_diagnostics
        ))
      }
      return(.sbf_add_predict_hazard(
        result = result,
        cov_data = cov_data,
        eval_times = eval_times,
        clamp_nonnegative = clamp_nonnegative,
        min_hazard = min_hazard,
        prediction = prediction,
        warn_diagnostics = warn_diagnostics
      ))
    }
    if (inherits(result, "sbf_multiplicative_sparse_binning_fit")) {
      return(.sbf_mult_sparse_predict_hazard(
        result = result,
        cov_data = cov_data,
        eval_times = eval_times,
        clamp_nonnegative = clamp_nonnegative,
        min_hazard = min_hazard,
        warn_diagnostics = warn_diagnostics
      ))
    }
    return(.sbf_mult_predict_hazard(
      result = result,
      cov_data = cov_data,
      eval_times = eval_times,
      clamp_nonnegative = clamp_nonnegative,
      min_hazard = min_hazard,
      min_component = 1e-8,
      prediction = prediction,
      warn_diagnostics = warn_diagnostics
    ))
  }

  query_times <- if (is.null(times)) {
    baseline_grid <- .sbf_finite_baseline_grid(result$x.grid)
    if (identical(prediction, "interpolation") && !is_binning) {
      baseline_grid[-1L]
    } else {
      baseline_grid
    }
  } else {
    times
  }
  if (inherits(result, "sbf_additive_fit")) {
    if (inherits(result, "sbf_additive_pairwise_binning_fit")) {
      return(.sbf_add_pairwise_predict_survival(
        result = result,
        cov_data = cov_data,
        query_times = query_times,
        clamp_nonnegative = clamp_nonnegative,
        min_hazard = min_hazard,
        warn_diagnostics = warn_diagnostics
      ))
    }
    return(.sbf_add_predict_survival(
      result = result,
      cov_data = cov_data,
      query_times = query_times,
      clamp_nonnegative = clamp_nonnegative,
      min_hazard = min_hazard,
      prediction = prediction,
      warn_diagnostics = warn_diagnostics
    ))
  }
  if (inherits(result, "sbf_multiplicative_sparse_binning_fit")) {
    return(.sbf_mult_sparse_predict_survival(
      result = result,
      cov_data = cov_data,
      query_times = query_times,
      clamp_nonnegative = clamp_nonnegative,
      min_hazard = min_hazard,
      warn_diagnostics = warn_diagnostics
    ))
  }
  .sbf_mult_predict_survival(
    result = result,
    cov_data = cov_data,
    query_times = query_times,
    clamp_nonnegative = clamp_nonnegative,
    min_hazard = min_hazard,
    min_component = 1e-8,
    prediction = prediction,
    warn_diagnostics = warn_diagnostics
  )
}
