# Additive pairwise prediction engine
# State validation ------------------------------------------------------------

# Validate the fit object and return the stored prediction state.
.sbf_add_pairwise_get_predict_state <- function(result) {
  .sbf_assert_required_fields(result, fields = c("x.grid", "alpha_backfit", "additive_pairwise_predict_state"))
  state <- result$additive_pairwise_predict_state
  state
}

# Kernel and fallback helpers -------------------------------------------------

# Build K_{j,0}(x,r), K_{j,1}(x,r), and K_{j,2}(x,r) from
# eq:add-predict-kernel for arbitrary query points x.
.sbf_add_pairwise_predict_kernel_rows <- function(xout,
                                        grid,
                                        bandwidth,
                                        kernel,
                                        correction_constants) {
  diff_mat <- outer(as.numeric(xout), as.numeric(grid), FUN = "-")
  K0 <- .sbf_eval_kernel(diff_mat, bandwidth, kernel)
  K0[!is.finite(K0)] <- 0
  K0 <- sweep(K0, 2L, as.numeric(correction_constants), "/")
  K0[!is.finite(K0)] <- 0
  list(
    K0 = K0,
    K1 = diff_mat * K0,
    K2 = diff_mat^2 * K0
  )
}

# Map evaluation points to the nearest/stored bin value for zero-support
# fallback outside the effective kernel support.
.sbf_add_pairwise_bin_step_state <- function(result, state, component, xout) {
  xout <- as.numeric(xout)
  alpha <- as.numeric(result$alpha_backfit[[component]])
  alpha1 <- as.numeric(result$alpha.1_backfit[[component]])
  breaks <- as.numeric(state$breaks[[component]])

  idx <- cut(xout, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  idx[is.na(idx) & xout < breaks[1L]] <- 1L
  idx[is.na(idx) & xout > breaks[length(breaks)]] <- length(alpha)
  idx <- as.integer(idx)
  idx[is.na(idx)] <- 1L
  idx <- pmax(1L, pmin(length(alpha), idx))

  list(value = alpha[idx], slope = alpha1[idx])
}

.sbf_add_pairwise_predict_guard_denominator <- function(x, floor_value) {
  x <- as.numeric(x)
  adjusted <- !is.finite(x) | x < floor_value
  x[adjusted] <- floor_value
  list(value = x, adjusted_count = sum(adjusted))
}

# Smoothing and component evaluation -----------------------------------------

# Evaluate the fitted local-linear state at x:
#   fitted value = S0 alpha + S1 alpha1,
#   fitted slope = S0 alpha1.
.sbf_add_pairwise_predict_smoothed_state <- function(result, state, component, xout) {
  xout <- as.numeric(xout)
  grid <- as.numeric(state$x_grid[[component]])
  widths <- as.numeric(state$widths[[component]])
  alpha <- as.numeric(result$alpha_backfit[[component]])
  alpha1 <- as.numeric(result$alpha.1_backfit[[component]])

  diff_mat <- outer(xout, grid, FUN = "-")
  K0 <- .sbf_eval_kernel(diff_mat, state$bandwidth[component], state$kernel)
  K0[!is.finite(K0)] <- 0
  weighted <- sweep(K0, 2L, widths, "*")
  denom <- rowSums(weighted)
  denom[!is.finite(denom) | abs(denom) < 1e-12] <- 1e-12
  S0 <- sweep(weighted, 1L, denom, "/")
  basis <- outer(xout, grid, function(eval_x, basis_x) basis_x - eval_x)

  value <- as.numeric(S0 %*% alpha + (S0 * basis) %*% alpha1)
  slope <- as.numeric(S0 %*% alpha1)
  grid_idx <- match(signif(xout, digits = 12), signif(grid, digits = 12))
  grid_match <- !is.na(grid_idx)
  if (any(grid_match)) {
    value[grid_match] <- alpha[grid_idx[grid_match]]
    slope[grid_match] <- alpha1[grid_idx[grid_match]]
  }
  list(value = value, slope = slope)
}

# Predict one LC component by evaluating the LC report equation at x:
#   alpha_j^pr(x) = (O_j(x) - B_j(x; alpha_-j)) / D_j(x).
.sbf_add_pairwise_predict_component_lc <- function(result, state, component, xout) {
  kernels <- .sbf_add_pairwise_predict_kernel_rows(
    xout = xout,
    grid = state$x_grid[[component]],
    bandwidth = state$bandwidth[component],
    kernel = state$kernel,
    correction_constants = state$kernel_correction_constants[[component]]
  )
  exposure <- if (component == 1L) {
    as.numeric(state$time_exposure)
  } else {
    rowSums(state$exposure_moments[[component]]$R0)
  }
  floor_value <- as.numeric(state$denominator_floor)

  D_raw <- as.numeric(kernels$K0 %*% exposure)
  D <- .sbf_add_pairwise_predict_guard_denominator(D_raw, floor_value)
  O <- as.numeric(kernels$K0 %*% as.numeric(state$occurrence[[component]]))
  B <- as.numeric(kernels$K0 %*% as.numeric(state$component_nuisance[[component]]))

  values <- (O - B) / D$value
  values[!is.finite(values)] <- NA_real_
  support_empty <- rowSums(kernels$K0 != 0) == 0
  if (any(support_empty)) {
    fallback <- .sbf_add_pairwise_bin_step_state(
      result = result,
      state = state,
      component = component,
      xout = as.numeric(xout)[support_empty]
    )
    values[support_empty] <- fallback$value
  }
  attr(values, "prediction_diagnostics") <- c(
    as.list(.sbf_binning_denominator_diagnostics(
      value_raw = D_raw,
      floor_value = floor_value,
      value_adjusted_count = D$adjusted_count
    )),
    list(
      kernel_support_zero_count = sum(support_empty),
      bin_step_fallback_count = sum(support_empty)
    )
  )
  values
}

# Predict one LL component by evaluating D0/C/D2/O0/O1/B0/B1 at x, then taking
# one report-style value update followed by one slope update.
.sbf_add_pairwise_predict_component_ll <- function(result, state, component, xout) {
  kernels <- .sbf_add_pairwise_predict_kernel_rows(
    xout = xout,
    grid = state$x_grid[[component]],
    bandwidth = state$bandwidth[component],
    kernel = state$kernel,
    correction_constants = state$kernel_correction_constants[[component]]
  )
  exposure <- if (component == 1L) {
    as.numeric(state$time_exposure)
  } else {
    rowSums(state$exposure_moments[[component]]$R0)
  }
  floor_value <- as.numeric(state$denominator_floor)
  old_state <- .sbf_add_pairwise_predict_smoothed_state(result, state, component, xout)

  D0_raw <- as.numeric(kernels$K0 %*% exposure)
  C <- as.numeric(kernels$K1 %*% exposure)
  D2_raw <- as.numeric(kernels$K2 %*% exposure)
  D0 <- .sbf_add_pairwise_predict_guard_denominator(D0_raw, floor_value)
  D2 <- .sbf_add_pairwise_predict_guard_denominator(D2_raw, floor_value)

  occurrence <- as.numeric(state$occurrence[[component]])
  O0 <- as.numeric(kernels$K0 %*% occurrence)
  O1 <- as.numeric(kernels$K1 %*% occurrence)

  nuisance <- as.numeric(state$component_nuisance[[component]])
  B0 <- as.numeric(kernels$K0 %*% nuisance)
  B1 <- as.numeric(kernels$K1 %*% nuisance)

  rhs0 <- O0 - B0
  rhs1 <- O1 - B1
  value_raw <- rhs0 - old_state$slope * C
  value <- value_raw / D0$value
  slope <- (rhs1 - value * C) / D2$value
  value[!is.finite(value)] <- NA_real_
  slope[!is.finite(slope)] <- NA_real_
  support_empty <- rowSums(kernels$K0 != 0) == 0
  if (any(support_empty)) {
    fallback <- .sbf_add_pairwise_bin_step_state(
      result = result,
      state = state,
      component = component,
      xout = as.numeric(xout)[support_empty]
    )
    value[support_empty] <- fallback$value
    slope[support_empty] <- fallback$slope
  }
  attr(value, "prediction_diagnostics") <- c(
    as.list(.sbf_binning_denominator_diagnostics(
      value_raw = D0_raw,
      slope_raw = D2_raw,
      floor_value = floor_value,
      value_adjusted_count = D0$adjusted_count,
      slope_adjusted_count = D2$adjusted_count
    )),
    list(
      kernel_support_zero_count = sum(support_empty),
      bin_step_fallback_count = sum(support_empty)
    )
  )
  attr(value, "slope") <- slope
  value
}


# Additive binned prediction adapter layer -----------------------------------

.sbf_add_pairwise_predict_component <- function(result,
                                                component,
                                                xout,
                                                warn_diagnostics = TRUE) {
  state <- .sbf_add_pairwise_get_predict_state(result)

  values <- if (isTRUE(state$local_constant)) {
    .sbf_add_pairwise_predict_component_lc(result, state, component, as.numeric(xout))
  } else {
    .sbf_add_pairwise_predict_component_ll(result, state, component, as.numeric(xout))
  }
  support_info <- .sbf_component_support(result, component_index = component)
  out <- .sbf_attach_component_support(values, xout = xout, support_info = support_info)
  .sbf_warn_binning_prediction_diagnostics(
    out,
    label = "Additive binned component prediction",
    warn_diagnostics = warn_diagnostics
  )
  out
}

.sbf_add_pairwise_predict_hazard <- function(result,
                                             cov_data,
                                             eval_times,
                                             clamp_nonnegative = FALSE,
                                             min_hazard = 1e-8,
                                             warn_diagnostics = TRUE) {
  state <- .sbf_add_pairwise_get_predict_state(result)
  n <- nrow(cov_data)
  p <- ncol(cov_data)

  covariate_components <- matrix(0, nrow = n, ncol = p)
  diagnostics <- vector("list", p + 1L)
  if (p > 0L) {
    for (j in seq_len(p)) {
      pred_j <- .sbf_add_pairwise_predict_component(
        result = result,
        component = j + 1L,
        xout = cov_data[, j],
        warn_diagnostics = FALSE
      )
      covariate_components[, j] <- as.numeric(pred_j)
      diagnostics[[j + 1L]] <- attr(pred_j, "prediction_diagnostics")
    }
  }
  covariate_term <- if (p > 0L) rowSums(covariate_components) else rep(0, n)

  baseline <- .sbf_add_pairwise_predict_component(
    result = result,
    component = 1L,
    xout = eval_times,
    warn_diagnostics = FALSE
  )
  diagnostics[[1L]] <- attr(baseline, "prediction_diagnostics")
  baseline <- as.numeric(baseline)

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
    label = "Additive binned hazard prediction",
    warn_diagnostics = warn_diagnostics
  )
  out
}

.sbf_add_pairwise_predict_survival <- function(result,
                                               cov_data,
                                               query_times,
                                               clamp_nonnegative = TRUE,
                                               min_hazard = 1e-8,
                                               warn_diagnostics = TRUE) {
  state <- .sbf_add_pairwise_get_predict_state(result)
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
      covariate_term = rep(0, n),
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
  hazard_pred <- .sbf_add_pairwise_predict_hazard(
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
    label = "Additive binned survival prediction",
    warn_diagnostics = warn_diagnostics
  )
  out
}
