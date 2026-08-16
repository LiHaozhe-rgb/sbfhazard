# Multiplicative binning prediction engine.

# State validation ------------------------------------------------------------

.sbf_mult_sparse_get_predict_state <- function(result) {
  .sbf_assert_required_fields(
    result,
    fields = c("x.grid", "alpha_backfit", "multiplicative_sparse_predict_state")
  )
  result$multiplicative_sparse_predict_state
}

# Fallback helper -------------------------------------------------------------

# Return stored component values at the bin containing each xout point.
.sbf_mult_sparse_bin_step_values <- function(result,
                                     state,
                                     component,
                                     xout) {
  min_component <- as.numeric(state$min_component)
  breaks <- as.numeric(state$breaks[[component]])
  values <- pmax(min_component, as.numeric(result$alpha_backfit[[component]]))
  idx <- cut(as.numeric(xout), breaks = breaks, include.lowest = TRUE, labels = FALSE)
  idx[is.na(idx) & as.numeric(xout) < breaks[1L]] <- 1L
  idx[is.na(idx) & as.numeric(xout) > breaks[length(breaks)]] <- length(values)
  idx <- pmax(1L, pmin(length(values), as.integer(idx)))
  values[idx]
}

.sbf_mult_sparse_predict_component <- function(result,
                                               component,
                                               xout,
                                               warn_diagnostics = TRUE) {
  state <- .sbf_mult_sparse_get_predict_state(result)
  min_component <- as.numeric(state$min_component)
  xout <- as.numeric(xout)
  kernel_info <- .sbf_mult_sparse_kernel_source_target(
    source = state$x_grid[[component]],
    target = xout,
    widths = state$widths[[component]],
    bandwidth = state$bandwidth[component],
    kernel = state$kernel
  )
  kernel_norm <- kernel_info$weights
  O <- .sbf_mult_sparse_component_occurrence(state, component, kernel_norm = kernel_norm)
  D_raw <- .sbf_mult_sparse_component_denominator(
    precomp = state,
    alpha = NULL,
    component = component,
    kernel_norm = kernel_norm,
    min_component = min_component,
    smoothed_components = state$smoothed_components
  )
  D <- .sbf_mult_sparse_guard_positive(
    D_raw,
    floor_value = as.numeric(state$denominator_floor)
  )
  values <- O / D$value
  values[!is.finite(values)] <- min_component
  values <- pmax(min_component, values)

  support_empty <- colSums(abs(kernel_norm) > 0) == 0
  if (any(support_empty)) {
    values[support_empty] <- .sbf_mult_sparse_bin_step_values(
      result = result,
      state = state,
      component = component,
      xout = xout[support_empty]
    )
  }
  attr(values, "prediction_diagnostics") <- c(
    as.list(.sbf_binning_denominator_diagnostics(
      value_raw = D_raw,
      floor_value = as.numeric(state$denominator_floor),
      value_adjusted_count = D$adjusted_count
    )),
    list(
      kernel_support_zero_count = sum(support_empty),
      bin_step_fallback_count = sum(support_empty)
    )
  )
  support_info <- .sbf_component_support(result, component_index = component)
  out <- .sbf_attach_component_support(values, xout = xout, support_info = support_info)
  .sbf_warn_binning_prediction_diagnostics(
    out,
    label = "Multiplicative binned component prediction",
    warn_diagnostics = warn_diagnostics
  )
  out
}

.sbf_mult_sparse_predict_hazard <- function(result,
                                            cov_data,
                                            eval_times,
                                            clamp_nonnegative = TRUE,
                                            min_hazard = 1e-8,
                                            warn_diagnostics = TRUE) {
  state <- .sbf_mult_sparse_get_predict_state(result)
  min_component <- as.numeric(state$min_component)
  n <- nrow(cov_data)
  p <- ncol(cov_data)

  baseline <- .sbf_mult_sparse_predict_component(
    result = result,
    component = 1L,
    xout = eval_times,
    warn_diagnostics = FALSE
  )
  diagnostics <- vector("list", p + 1L)
  diagnostics[[1L]] <- attr(baseline, "prediction_diagnostics")
  baseline <- pmax(min_component, as.numeric(baseline))

  covariate_components <- matrix(1, nrow = n, ncol = p)
  if (p > 0L) {
    for (j in seq_len(p)) {
      comp <- .sbf_mult_sparse_predict_component(
        result = result,
        component = j + 1L,
        xout = cov_data[, j],
        warn_diagnostics = FALSE
      )
      covariate_components[, j] <- pmax(min_component, as.numeric(comp))
      diagnostics[[j + 1L]] <- attr(comp, "prediction_diagnostics")
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

  out <- .sbf_build_hazard_output(
    hazard = hazard,
    times = eval_times,
    baseline = baseline,
    covariate_term = covariate_term,
    covariate_components = covariate_components,
    clamp_count = clamp_count
  )
  out$prediction_diagnostics <- .sbf_binning_prediction_diagnostics_table(
    diagnostics,
    result = result
  )
  .sbf_warn_binning_prediction_diagnostics(
    out,
    label = "Multiplicative binned hazard prediction",
    warn_diagnostics = warn_diagnostics
  )
  out
}

.sbf_mult_sparse_predict_survival <- function(result,
                                              cov_data,
                                              query_times,
                                              clamp_nonnegative = TRUE,
                                              min_hazard = 1e-8,
                                              warn_diagnostics = TRUE) {
  state <- .sbf_mult_sparse_get_predict_state(result)
  baseline_grid <- .sbf_finite_baseline_grid(state$x_grid)
  baseline_support <- .sbf_component_support(result, component_index = 1L)

  n <- nrow(cov_data)
  survival <- matrix(1, nrow = n, ncol = length(query_times))
  t0 <- baseline_support$support_min
  positive_idx <- which(query_times > t0)
  if (length(positive_idx) == 0L) {
    out <- .sbf_build_survival_output(
      survival = survival,
      times = query_times,
      covariate_term = rep(1, n),
      clamp_count = 0L,
      integration_start = t0,
      baseline_support = baseline_support
    )
    out$prediction_diagnostics <- .sbf_binning_prediction_diagnostics_table(
      list(),
      result = result
    )
    return(out)
  }

  integration_grid <- sort(unique(c(t0, baseline_grid, query_times[positive_idx])))
  hazard_pred <- .sbf_mult_sparse_predict_hazard(
    result = result,
    cov_data = cov_data,
    eval_times = integration_grid,
    clamp_nonnegative = clamp_nonnegative,
    min_hazard = min_hazard,
    warn_diagnostics = FALSE
  )
  delta <- diff(integration_grid)
  query_index <- match(query_times[positive_idx], integration_grid)
  for (i in seq_len(n)) {
    h <- as.numeric(hazard_pred$hazard[i, ])
    cumhaz <- c(0, cumsum(0.5 * (h[-length(h)] + h[-1L]) * delta))
    survival[i, positive_idx] <- exp(-cumhaz[query_index])
  }

  out <- .sbf_build_survival_output(
    survival = survival,
    times = query_times,
    covariate_term = hazard_pred$covariate_term,
    covariate_components = hazard_pred$covariate_components,
    clamp_count = hazard_pred$clamp_count,
    integration_start = t0,
    baseline_support = baseline_support
  )
  out$prediction_diagnostics <- hazard_pred$prediction_diagnostics
  .sbf_warn_binning_prediction_diagnostics(
    out,
    label = "Multiplicative binned survival prediction",
    warn_diagnostics = warn_diagnostics
  )
  out
}
