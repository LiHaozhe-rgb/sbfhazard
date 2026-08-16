# Multiplicative direct-prediction numerical engine.

.sbf_multiplicative_get_predict_state <- function(result) {
  .sbf_assert_required_fields(
    result,
    fields = c("x.grid", "alpha_backfit", "multiplicative_predict_state")
  )
  state <- result$multiplicative_predict_state
  state$alpha <- lapply(result$alpha_backfit, as.numeric)
  state$x_grid <- result$x.grid
  state
}

.sbf_clamp_component_values <- function(x, lower_bound) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- lower_bound
  pmax(lower_bound, x)
}

.sbf_predict_formula_component_multiplicative_engine <- function(state,
                                                                 component_index,
                                                                 xout) {
  xout <- as.numeric(xout)

  K.new <- .sbf_eval_kernel(
    outer(state$train_x[, component_index], xout, FUN = "-"),
    state$bandwidth[component_index],
    state$kernel
  )
  K.new.norm <- sweep(K.new, 1, state$k.X.b[[component_index]], "/")
  O <- colSums(K.new.norm * state$status)

  if (component_index == 1L) {
    K.time.new <- .sbf_eval_kernel(
      outer(as.numeric(state$x_grid[[1]]), xout, FUN = "-"),
      state$bandwidth[1],
      state$kernel
    )
    K.time.new.norm <- sweep(K.time.new, 1, state$k.b, "/")
    D <- as.numeric(crossprod(as.numeric(state$baseline.exposure), K.time.new.norm))
  } else {
    scalar <- state$covariate.scalar[[component_index]]
    D <- as.numeric(crossprod(scalar, K.new.norm))
  }

  O / D
}

.sbf_predict_hazard_multiplicative_engine <- function(state,
                                                      cov_data,
                                                      eval_times,
                                                      clamp_nonnegative = TRUE,
                                                      min_hazard = 1e-8,
                                                      min_component = 1e-8) {
  eval_times <- as.numeric(eval_times)
  n <- nrow(cov_data)
  p <- length(state$alpha) - 1L

  baseline <- .sbf_predict_formula_component_multiplicative_engine(state, component_index = 1L, xout = eval_times)
  baseline <- .sbf_clamp_component_values(baseline, min_component)

  covariate_components <- matrix(1, nrow = n, ncol = p)
  if (p > 0L) {
    for (j in seq_len(p)) {
      vals <- .sbf_predict_formula_component_multiplicative_engine(state, component_index = j + 1L, xout = cov_data[, j])
      covariate_components[, j] <- .sbf_clamp_component_values(vals, min_component)
    }
  }

  covariate_term <- if (p > 0L) apply(covariate_components, 1, prod) else rep(1, n)
  hazard <- tcrossprod(covariate_term, baseline)
  raw_hazard <- hazard
  if (!clamp_nonnegative && any(!is.finite(raw_hazard))) {
    stop(
      "Prediction produced non-finite hazards with clamp_nonnegative = FALSE.",
      call. = FALSE
    )
  }
  if (clamp_nonnegative) {
    hazard[!is.finite(hazard)] <- min_hazard
    hazard <- matrix(pmax(min_hazard, as.numeric(hazard)), nrow = n, ncol = length(eval_times), byrow = FALSE)
  }
  clamp_count <- if (clamp_nonnegative) {
    sum(!is.finite(raw_hazard) | raw_hazard < min_hazard)
  } else {
    0L
  }

  colnames(hazard) <- paste0("t", seq_along(eval_times))

  list(
    hazard = hazard,
    times = eval_times,
    baseline = baseline,
    covariate_term = covariate_term,
    covariate_components = covariate_components,
    clamp_count = clamp_count
  )
}

.sbf_predict_interpolation_hazard_multiplicative_engine <- function(state,
                                                                    cov_data,
                                                                    eval_times,
                                                                    clamp_nonnegative = TRUE,
                                                                    min_hazard = 1e-8,
                                                                    min_component = 1e-8) {
  eval_times <- as.numeric(eval_times)
  n <- nrow(cov_data)
  p <- length(state$alpha) - 1L

  baseline <- .sbf_interpolate_grid_values(
    grid = state$x_grid[[1L]],
    values = state$alpha[[1L]],
    xout = eval_times
  )
  baseline <- .sbf_clamp_component_values(baseline, min_component)

  covariate_components <- matrix(1, nrow = n, ncol = p)
  if (p > 0L) {
    for (j in seq_len(p)) {
      vals <- .sbf_interpolate_grid_values(
        grid = state$x_grid[[j + 1L]],
        values = state$alpha[[j + 1L]],
        xout = cov_data[, j]
      )
      covariate_components[, j] <- .sbf_clamp_component_values(vals, min_component)
    }
  }

  covariate_term <- if (p > 0L) apply(covariate_components, 1L, prod) else rep(1, n)
  hazard <- tcrossprod(covariate_term, baseline)
  raw_hazard <- hazard
  if (!clamp_nonnegative && any(!is.finite(raw_hazard))) {
    stop(
      "Prediction produced non-finite hazards with clamp_nonnegative = FALSE.",
      call. = FALSE
    )
  }
  if (clamp_nonnegative) {
    hazard[!is.finite(hazard)] <- min_hazard
    hazard <- matrix(pmax(min_hazard, as.numeric(hazard)), nrow = n, ncol = length(eval_times), byrow = FALSE)
  }
  clamp_count <- if (clamp_nonnegative) {
    sum(!is.finite(raw_hazard) | raw_hazard < min_hazard)
  } else {
    0L
  }
  colnames(hazard) <- paste0("t", seq_along(eval_times))

  list(
    hazard = hazard,
    times = eval_times,
    baseline = baseline,
    covariate_term = covariate_term,
    covariate_components = covariate_components,
    clamp_count = clamp_count
  )
}

.sbf_predict_survival_multiplicative_engine <- function(state,
                                                        cov_data,
                                                        query_times,
                                                        baseline_grid,
                                                        integration_start,
                                                        clamp_nonnegative = TRUE,
                                                        min_hazard = 1e-8,
                                                        min_component = 1e-8) {
  query_times <- as.numeric(query_times)
  n <- nrow(cov_data)
  survival <- matrix(1, nrow = n, ncol = length(query_times))

  positive_idx <- which(query_times > integration_start)
  if (length(positive_idx) == 0L) {
    return(list(
      survival = survival,
      times = query_times,
      covariate_term = rep(1, n),
      covariate_components = NULL,
      clamp_count = 0L
    ))
  }

  integration_grid <- sort(unique(c(integration_start, baseline_grid, query_times[positive_idx])))
  hazard_pred <- .sbf_predict_hazard_multiplicative_engine(
    state = state,
    cov_data = cov_data,
    eval_times = integration_grid,
    clamp_nonnegative = clamp_nonnegative,
    min_hazard = min_hazard,
    min_component = min_component
  )

  delta <- diff(integration_grid)
  query_index <- match(query_times[positive_idx], integration_grid)
  for (i in seq_len(n)) {
    h <- as.numeric(hazard_pred$hazard[i, ])
    # Union-grid right-endpoint alternative:
    # cumhaz <- c(0, cumsum(h[-1L] * delta))
    # Active implementation: union-grid trapezoidal integration.
    cumhaz <- c(0, cumsum(0.5 * (h[-length(h)] + h[-1L]) * delta))
    survival[i, positive_idx] <- exp(-cumhaz[query_index])
  }

  list(
    survival = survival,
    times = query_times,
    covariate_term = hazard_pred$covariate_term,
    covariate_components = hazard_pred$covariate_components,
    clamp_count = hazard_pred$clamp_count
  )
}

.sbf_predict_interpolation_survival_multiplicative_engine <- function(state,
                                                                      cov_data,
                                                                      query_times,
                                                                      baseline_grid,
                                                                      integration_start,
                                                                      clamp_nonnegative = TRUE,
                                                                      min_hazard = 1e-8,
                                                                      min_component = 1e-8) {
  query_times <- as.numeric(query_times)
  n <- nrow(cov_data)
  survival <- matrix(1, nrow = n, ncol = length(query_times))
  positive_idx <- which(query_times > integration_start)
  integration_grid <- sort(unique(c(
    integration_start,
    baseline_grid,
    query_times[positive_idx]
  )))
  hazard_pred <- .sbf_predict_interpolation_hazard_multiplicative_engine(
    state = state,
    cov_data = cov_data,
    eval_times = integration_grid,
    clamp_nonnegative = clamp_nonnegative,
    min_hazard = min_hazard,
    min_component = min_component
  )

  delta <- diff(integration_grid)
  query_index <- match(query_times[positive_idx], integration_grid)
  for (i in seq_len(n)) {
    h <- as.numeric(hazard_pred$hazard[i, ])
    # Union-grid right-endpoint alternative:
    # cumhaz <- c(0, cumsum(h[-1L] * delta))
    # Active implementation: union-grid trapezoidal integration.
    cumhaz <- c(0, cumsum(0.5 * (h[-length(h)] + h[-1L]) * delta))
    survival[i, positive_idx] <- exp(-cumhaz[query_index])
  }

  list(
    survival = survival,
    times = query_times,
    covariate_term = hazard_pred$covariate_term,
    covariate_components = hazard_pred$covariate_components,
    clamp_count = hazard_pred$clamp_count
  )
}


# Multiplicative SBF prediction adapter layer -------------------------------

.sbf_mult_predict_component <- function(result,
                                              component,
                                              xout,
                                              min_component = 1e-8,
                                              prediction = "formula",
                                              warn_diagnostics = TRUE) {
  state <- .sbf_multiplicative_get_predict_state(result)

  support_info <- .sbf_component_support(result, component_index = component)
  values <- if (identical(prediction, "interpolation")) {
    .sbf_interpolate_grid_values(
      grid = state$x_grid[[component]],
      values = state$alpha[[component]],
      xout = xout
    )
  } else {
    .sbf_predict_formula_component_multiplicative_engine(state, component_index = component, xout = xout)
  }
  clamp_count <- sum(!is.finite(values) | values < min_component)
  values <- .sbf_clamp_component_values(values, min_component)

  out <- .sbf_attach_component_support(values, xout = xout, support_info = support_info)
  attr(out, "prediction_diagnostics") <- .sbf_prediction_diagnostics(
    out,
    extrapolated_count = sum(as.logical(attr(out, "extrapolated")), na.rm = TRUE),
    clamp_count = clamp_count
  )
  .sbf_warn_fit_prediction_diagnostics(
    out,
    label = "Multiplicative SBF component prediction",
    warn_diagnostics = warn_diagnostics
  )
  out
}

.sbf_mult_predict_hazard <- function(result,
                                           cov_data,
                                           eval_times,
                                           clamp_nonnegative = TRUE,
                                           min_hazard = 1e-8,
                                           min_component = 1e-8,
                                           prediction = "formula",
                                           warn_diagnostics = TRUE) {
  state <- .sbf_multiplicative_get_predict_state(result)

  pred <- if (identical(prediction, "interpolation")) {
    .sbf_predict_interpolation_hazard_multiplicative_engine(
      state = state,
      cov_data = cov_data,
      eval_times = eval_times,
      clamp_nonnegative = clamp_nonnegative,
      min_hazard = min_hazard,
      min_component = min_component
    )
  } else {
    .sbf_predict_hazard_multiplicative_engine(
      state = state,
      cov_data = cov_data,
      eval_times = eval_times,
      clamp_nonnegative = clamp_nonnegative,
      min_hazard = min_hazard,
      min_component = min_component
    )
  }

  out <- .sbf_build_hazard_output(
    hazard = pred$hazard,
    times = pred$times,
    baseline = pred$baseline,
    covariate_term = pred$covariate_term,
    covariate_components = pred$covariate_components,
    clamp_count = pred$clamp_count
  )
  support_count <- .sbf_count_prediction_support_violations(result, eval_times, cov_data)
  out$prediction_diagnostics <- .sbf_prediction_diagnostics(
    out$hazard,
    extrapolated_count = support_count,
    clamp_count = pred$clamp_count
  )
  .sbf_warn_fit_prediction_diagnostics(
    out,
    label = "Multiplicative SBF hazard prediction",
    warn_diagnostics = warn_diagnostics
  )
  out
}

.sbf_mult_predict_survival <- function(result,
                                             cov_data,
                                             query_times,
                                             clamp_nonnegative = TRUE,
                                             min_hazard = 1e-8,
                                             min_component = 1e-8,
                                             prediction = "formula",
                                             warn_diagnostics = TRUE) {
  state <- .sbf_multiplicative_get_predict_state(result)
  baseline_grid <- .sbf_finite_baseline_grid(state$x_grid)
  baseline_support <- .sbf_component_support(result, component_index = 1L)

  pred <- if (identical(prediction, "interpolation")) {
    .sbf_predict_interpolation_survival_multiplicative_engine(
      state = state,
      cov_data = cov_data,
      query_times = query_times,
      baseline_grid = baseline_grid,
      integration_start = baseline_support$support_min,
      clamp_nonnegative = clamp_nonnegative,
      min_hazard = min_hazard,
      min_component = min_component
    )
  } else {
    .sbf_predict_survival_multiplicative_engine(
      state = state,
      cov_data = cov_data,
      query_times = query_times,
      baseline_grid = baseline_grid,
      integration_start = baseline_support$support_min,
      clamp_nonnegative = clamp_nonnegative,
      min_hazard = min_hazard,
      min_component = min_component
    )
  }

  out <- .sbf_build_survival_output(
    survival = pred$survival,
    times = pred$times,
    covariate_term = pred$covariate_term,
    covariate_components = pred$covariate_components,
    clamp_count = pred$clamp_count,
    integration_start = baseline_support$support_min,
    baseline_support = baseline_support
  )
  support_count <- .sbf_count_prediction_support_violations(result, query_times, cov_data)
  out$prediction_diagnostics <- .sbf_prediction_diagnostics(
    out$survival,
    extrapolated_count = support_count,
    clamp_count = pred$clamp_count
  )
  .sbf_warn_fit_prediction_diagnostics(
    out,
    label = "Multiplicative SBF survival prediction",
    warn_diagnostics = warn_diagnostics
  )
  out
}
