# Additive direct-prediction numerical engine.

.sbf_additive_get_predict_state <- function(result) {
  .sbf_assert_required_fields(
    result,
    fields = c("x.grid", "alpha_backfit", "alpha.1_backfit", "additive_predict_state")
  )
  state <- result$additive_predict_state
  state$alpha <- lapply(result$alpha_backfit, as.numeric)
  state$alpha1 <- lapply(result$alpha.1_backfit, as.numeric)
  state$x_grid <- result$x.grid
  state
}

.sbf_safe_divide_additive <- function(num, den) {
  out <- as.numeric(num) / as.numeric(den)
  out[!is.finite(out)] <- NA_real_
  out
}

.sbf_predict_smoothed_state_additive <- function(state, component_index, xout) {
  xout <- as.numeric(xout)
  grid <- as.numeric(state$x_grid[[component_index]])
  alpha <- as.numeric(state$alpha[[component_index]])
  alpha1 <- as.numeric(state$alpha1[[component_index]])
  if (length(alpha1) != length(alpha)) {
    alpha1 <- rep(0, length(alpha))
  }

  basis.minus.eval <- outer(grid, xout, FUN = "-") # length(grid) x length(xout) matrix
  K.raw <- .sbf_eval_kernel(basis.minus.eval, state$bandwidth[component_index], state$kernel)
  K.raw[!is.finite(K.raw)] <- 0
  dx <- as.numeric(state$dx[[component_index]])
  if (isTRUE(state$kernel_correction)) {
    denom <- colSums(sweep(K.raw, 1, dx, "*"))
    denom <- pmax(as.numeric(denom), 1e-12)
  } else {
    denom <- rep(1, length(xout))
  }
  weights <- sweep(sweep(K.raw, 1, dx, "*"), 2, denom, "/")
  weights[!is.finite(weights)] <- 0

  value <- as.numeric(crossprod(alpha, weights) + crossprod(alpha1, basis.minus.eval * weights))
  slope <- as.numeric(crossprod(alpha1, weights))
  grid_idx <- match(signif(xout, digits = 12), signif(grid, digits = 12))
  grid_match <- !is.na(grid_idx)
  if (any(grid_match)) {
    value[grid_match] <- alpha[grid_idx[grid_match]]
    slope[grid_match] <- alpha1[grid_idx[grid_match]]
  }
  list(value = value, slope = slope)
}

.sbf_predict_formula_component_additive_engine <- function(state, component_index, xout) {
  xout <- as.numeric(xout)
  old_state <- .sbf_predict_smoothed_state_additive(state, component_index = component_index, xout = xout)

  if (component_index == 1L) {
    delta.obs <- outer(xout, state$X[, 1], FUN = "-")
    K.obs.raw <- .sbf_eval_kernel(delta.obs, state$bandwidth[1], state$kernel)
    K.obs <- sweep(K.obs.raw, 2, state$k.X.b[[1]], "/")
    K.obs[!is.finite(K.obs)] <- 0
    O0 <- as.numeric(K.obs %*% state$status)
    O1 <- as.numeric((delta.obs * K.obs) %*% state$status)

    grid.minus.xout <- outer(as.numeric(state$x_grid[[1]]), xout, FUN = "-")
    xout.minus.grid <- -grid.minus.xout
    K.time.raw <- .sbf_eval_kernel(grid.minus.xout, state$bandwidth[1], state$kernel)
    time.scale <- as.numeric(state$dx[[1]]) / state$k.b
    time.weights <- sweep(K.time.raw, 1, time.scale, "*")
    time.weights[!is.finite(time.weights)] <- 0
    slope.time.weights <- xout.minus.grid * time.weights
    slope.denom.weights <- xout.minus.grid^2 * time.weights

    D0 <- as.numeric(crossprod(state$Y.colsum, time.weights))
    D1 <- as.numeric(crossprod(state$Y.colsum, slope.denom.weights))
    O3 <- as.numeric(crossprod(state$Y.colsum, slope.time.weights))
    N0 <- as.numeric(crossprod(state$nuisance.baseline, time.weights))
    N1 <- as.numeric(crossprod(state$nuisance.baseline, slope.time.weights))
  } else {
    delta <- outer(xout, state$X[, component_index], FUN = "-")
    K.raw <- .sbf_eval_kernel(delta, state$bandwidth[component_index], state$kernel)
    K.new <- sweep(K.raw, 2, state$k.X.b[[component_index]], "/")
    K.new[!is.finite(K.new)] <- 0

    O0 <- as.numeric(K.new %*% state$status)
    O1 <- as.numeric((delta * K.new) %*% state$status)
    D0 <- as.numeric(K.new %*% state$Y.dx.1)
    D1 <- as.numeric((delta^2 * K.new) %*% state$Y.dx.1)
    O3 <- as.numeric((delta * K.new) %*% state$Y.dx.1)
    N0 <- as.numeric(K.new %*% state$nuisance.covariate[[component_index]])
    N1 <- as.numeric((delta * K.new) %*% state$nuisance.covariate[[component_index]])
  }

  value <- .sbf_safe_divide_additive(O0 - old_state$slope * O3 - N0, D0)
  if (!isTRUE(state$local_constant)) {
    slope <- .sbf_safe_divide_additive(O1 - value * O3 - N1, D1)
    attr(value, "slope") <- slope
  }
  value
}

.sbf_predict_hazard_additive_engine <- function(state,
                                                cov_data,
                                                eval_times,
                                                clamp_nonnegative = FALSE,
                                                min_hazard = 1e-8) {
  eval_times <- as.numeric(eval_times)
  n <- nrow(cov_data)
  p <- ncol(cov_data)

  baseline <- .sbf_predict_formula_component_additive_engine(state, component_index = 1L, xout = eval_times)
  covariate_components <- matrix(0, nrow = n, ncol = p)
  if (p > 0L) {
    for (j in seq_len(p)) {
      covariate_components[, j] <- .sbf_predict_formula_component_additive_engine(
        state,
        component_index = j + 1L,
        xout = cov_data[, j]
      )
    }
  }
  covariate_term <- if (p > 0L) rowSums(covariate_components) else rep(0, n)

  hazard <- outer(covariate_term, baseline, "+")
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

.sbf_predict_interpolation_hazard_additive_engine <- function(state,
                                                              cov_data,
                                                              eval_times,
                                                              clamp_nonnegative = FALSE,
                                                              min_hazard = 1e-8) {
  eval_times <- as.numeric(eval_times)
  n <- nrow(cov_data)
  p <- ncol(cov_data)

  baseline <- .sbf_interpolate_grid_values(
    grid = state$x_grid[[1L]],
    values = state$alpha[[1L]],
    xout = eval_times
  )
  covariate_components <- matrix(0, nrow = n, ncol = p)
  if (p > 0L) {
    for (j in seq_len(p)) {
      covariate_components[, j] <- .sbf_interpolate_grid_values(
        grid = state$x_grid[[j + 1L]],
        values = state$alpha[[j + 1L]],
        xout = cov_data[, j]
      )
    }
  }
  covariate_term <- if (p > 0L) rowSums(covariate_components) else rep(0, n)

  hazard <- outer(covariate_term, baseline, "+")
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

.sbf_predict_survival_additive_engine <- function(state,
                                                  cov_data,
                                                  query_times,
                                                  baseline_grid,
                                                  integration_start,
                                                  clamp_nonnegative = TRUE,
                                                  min_hazard = 1e-8) {
  query_times <- as.numeric(query_times)
  n <- nrow(cov_data)
  survival <- matrix(1, nrow = n, ncol = length(query_times))

  positive_idx <- which(query_times > integration_start)
  if (length(positive_idx) == 0L) {
    return(list(
      survival = survival,
      times = query_times,
      covariate_term = rep(0, n),
      covariate_components = NULL,
      clamp_count = 0L
    ))
  }

  integration_grid <- sort(unique(c(integration_start, baseline_grid, query_times[positive_idx])))
  hazard_pred <- .sbf_predict_hazard_additive_engine(
    state = state,
    cov_data = cov_data,
    eval_times = integration_grid,
    clamp_nonnegative = clamp_nonnegative,
    min_hazard = min_hazard
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

.sbf_predict_interpolation_survival_additive_engine <- function(state,
                                                                cov_data,
                                                                query_times,
                                                                baseline_grid,
                                                                integration_start,
                                                                clamp_nonnegative = TRUE,
                                                                min_hazard = 1e-8) {
  query_times <- as.numeric(query_times)
  n <- nrow(cov_data)
  survival <- matrix(1, nrow = n, ncol = length(query_times))
  positive_idx <- which(query_times > integration_start)
  integration_grid <- sort(unique(c(
    integration_start,
    baseline_grid,
    query_times[positive_idx]
  )))
  hazard_pred <- .sbf_predict_interpolation_hazard_additive_engine(
    state = state,
    cov_data = cov_data,
    eval_times = integration_grid,
    clamp_nonnegative = clamp_nonnegative,
    min_hazard = min_hazard
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


# Additive SBF prediction adapter layer ------------------------------------

.sbf_add_predict_component <- function(result,
                                       component,
                                       xout,
                                       prediction = "formula",
                                       warn_diagnostics = TRUE) {
  state <- .sbf_additive_get_predict_state(result)

  support_info <- .sbf_component_support(result, component_index = component)
  values <- if (identical(prediction, "interpolation")) {
    .sbf_interpolate_grid_values(
      grid = state$x_grid[[component]],
      values = state$alpha[[component]],
      xout = xout
    )
  } else {
    .sbf_predict_formula_component_additive_engine(state, component_index = component, xout = xout)
  }
  out <- .sbf_attach_component_support(values, xout = xout, support_info = support_info)
  attr(out, "prediction_diagnostics") <- .sbf_prediction_diagnostics(
    out,
    extrapolated_count = sum(as.logical(attr(out, "extrapolated")), na.rm = TRUE),
    clamp_count = 0L
  )
  .sbf_warn_fit_prediction_diagnostics(
    out,
    label = "Additive SBF component prediction",
    warn_diagnostics = warn_diagnostics
  )
  out
}

.sbf_add_predict_hazard <- function(result,
                                    cov_data,
                                    eval_times,
                                    clamp_nonnegative = FALSE,
                                    min_hazard = 1e-8,
                                    prediction = "formula",
                                    warn_diagnostics = TRUE) {
  state <- .sbf_additive_get_predict_state(result)
  pred <- if (identical(prediction, "interpolation")) {
    .sbf_predict_interpolation_hazard_additive_engine(
      state = state,
      cov_data = cov_data,
      eval_times = eval_times,
      clamp_nonnegative = clamp_nonnegative,
      min_hazard = min_hazard
    )
  } else {
    .sbf_predict_hazard_additive_engine(
      state = state,
      cov_data = cov_data,
      eval_times = eval_times,
      clamp_nonnegative = clamp_nonnegative,
      min_hazard = min_hazard
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
    label = "Additive SBF hazard prediction",
    warn_diagnostics = warn_diagnostics
  )
  out
}

.sbf_add_predict_survival <- function(result,
                                      cov_data,
                                      query_times,
                                      clamp_nonnegative = TRUE,
                                      min_hazard = 1e-8,
                                      prediction = "formula",
                                      warn_diagnostics = TRUE) {
  state <- .sbf_additive_get_predict_state(result)
  baseline_grid <- .sbf_finite_baseline_grid(state$x_grid)
  baseline_support <- .sbf_component_support(result, component_index = 1L)

  pred <- if (identical(prediction, "interpolation")) {
    .sbf_predict_interpolation_survival_additive_engine(
      state = state,
      cov_data = cov_data,
      query_times = query_times,
      baseline_grid = baseline_grid,
      integration_start = baseline_support$support_min,
      clamp_nonnegative = clamp_nonnegative,
      min_hazard = min_hazard
    )
  } else {
    .sbf_predict_survival_additive_engine(
      state = state,
      cov_data = cov_data,
      query_times = query_times,
      baseline_grid = baseline_grid,
      integration_start = baseline_support$support_min,
      clamp_nonnegative = clamp_nonnegative,
      min_hazard = min_hazard
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
    label = "Additive SBF survival prediction",
    warn_diagnostics = warn_diagnostics
  )
  out
}
