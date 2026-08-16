# Multiplicative binning fit solver.

# Numerical guards ------------------------------------------------------------

# Floor non-positive or non-finite denominators/components before division.
.sbf_mult_sparse_guard_positive <- function(x,
                                    floor_value) {
  x <- as.numeric(x)
  bad <- !is.finite(x) | x < floor_value
  x[bad] <- floor_value
  list(value = x, adjusted_count = sum(bad))
}

# Smoothed occurrence numerator for one component.
.sbf_mult_sparse_component_occurrence <- function(precomp, component, kernel_norm = NULL) {
  if (is.null(kernel_norm)) {
    kernel_norm <- precomp$kernel_norm[[component]]
  }
  as.numeric(crossprod(as.numeric(precomp$occurrence[[component]]), kernel_norm))
}

# Smooth component values over the bin grid.
.sbf_mult_sparse_smooth_component <- function(precomp, alpha, component) {
  as.numeric(
    precomp$kernel_norm[[component]] %*%
      (as.numeric(alpha[[component]]) * as.numeric(precomp$widths[[component]]))
  )
}

# Denominator construction ----------------------------------------------------

# For each observed joint cell g, compute prod_{j in covariates} alpha_j(r_j(g)).
.sbf_mult_sparse_cov_product_by_joint <- function(precomp,
                                          alpha,
                                          covariates,
                                          min_component = 1e-8,
                                          smoothed_components = NULL) {
  product <- rep(1, nrow(precomp$joint_bins))
  for (j in covariates) {
    component <- j + 1L
    smooth_j <- if (is.null(smoothed_components)) {
      .sbf_mult_sparse_smooth_component(precomp, alpha, component)
    } else {
      as.numeric(smoothed_components[[component]])
    }
    smooth_j <- pmax(min_component, smooth_j)
    product <- product * smooth_j[precomp$joint_bins[, j]]
  }
  product
}

# Compute the denominator in the update.
# - component 1: sum_s K_0(l,s) sum_g E_g(s) prod_j \bar alpha_j(r_j(g)).
# - component k>1: group joint cells by r_k(g)=r, multiply baseline-smoothed
#   exposure by the product over all other covariates, then smooth over r.
.sbf_mult_sparse_component_denominator <- function(precomp,
                                           alpha,
                                           component,
                                           kernel_norm = NULL,
                                           min_component = 1e-8,
                                           smoothed_components = NULL) {
  if (is.null(kernel_norm)) {
    kernel_norm <- precomp$kernel_norm[[component]]
  }
  d <- length(precomp$n_grid)
  p <- d - 1L

  if (component == 1L) {
    cov_product <- .sbf_mult_sparse_cov_product_by_joint(
      precomp = precomp,
      alpha = alpha,
      covariates = seq_len(p),
      min_component = min_component,
      smoothed_components = smoothed_components
    )
    weighted_time <- colSums(sweep(precomp$joint_exposure, 1L, cov_product, "*"))
    return(as.numeric(crossprod(weighted_time, kernel_norm)))
  }

  k_cov <- component - 1L
  baseline_smooth <- if (is.null(smoothed_components)) {
    .sbf_mult_sparse_smooth_component(precomp, alpha, component = 1L)
  } else {
    as.numeric(smoothed_components[[1L]])
  }
  baseline_smooth <- pmax(min_component, baseline_smooth)
  other_product <- .sbf_mult_sparse_cov_product_by_joint(
    precomp = precomp,
    alpha = alpha,
    covariates = setdiff(seq_len(p), k_cov),
    min_component = min_component,
    smoothed_components = smoothed_components
  )
  exposure_with_baseline <- sweep(precomp$joint_exposure, 2L, baseline_smooth, "*")
  scalar_by_joint <- rowSums(exposure_with_baseline) * other_product
  source_scalar <- numeric(precomp$n_grid[component])
  grouped <- rowsum(
    scalar_by_joint,
    group = precomp$joint_bins[, k_cov],
    reorder = FALSE
  )
  source_scalar[as.integer(rownames(grouped))] <- as.numeric(grouped[, 1L])
  as.numeric(crossprod(source_scalar, kernel_norm))
}

# Normalization and convergence ----------------------------------------------

# Move each covariate component's scale into alpha_0.
.sbf_mult_sparse_normalize_alpha <- function(precomp,
                                             alpha,
                                             min_component = 1e-8,
                                             identification = "sample_mean",
                                             truth_functions = NULL) {
  d <- length(alpha)
  if (d <= 1L) {
    return(alpha)
  }

  for (component in 2:d) {
    scale <- .sbf_multiplicative_identification_scale(
      alpha_component = alpha[[component]],
      grid = precomp$x_grid[[component]],
      widths = precomp$widths[[component]],
      identification = identification,
      min_component = min_component,
      sample_weights = precomp$observation_counts[[component]],
      truth_functions = truth_functions,
      covariate_index = component - 1L
    )
    alpha[[component]] <- as.numeric(alpha[[component]]) / scale
    alpha[[1L]] <- as.numeric(alpha[[1L]]) * scale
  }
  alpha
}

# Pick a stable denominator floor from denominators.
.sbf_mult_sparse_denominator_floor <- function(denominators) {
  positive <- denominators[is.finite(denominators) & denominators > 0]
  if (length(positive) == 0L) {
    return(1e-10)
  }
  max(1e-12, 1e-10 * stats::median(positive))
}

# Diagnostics -----------------------------------------------------------------

.sbf_mult_sparse_fit_diagnostics <- function(precomp,
                                             raw_denominators,
                                             denominator_floor,
                                             denominator_floor_count,
                                             component_floor_count) {
  denominator_diag <- .sbf_binning_denominator_diagnostics(
    value_raw = raw_denominators,
    floor_value = denominator_floor,
    value_adjusted_count = denominator_floor_count
  )
  diagnostics <- cbind(denominator_diag, precomp$kernel_norm_diagnostics)
  diagnostics$component_floor_count <- as.integer(component_floor_count)
  diagnostics
}

# Prediction state ------------------------------------------------------------

.sbf_mult_sparse_build_predict_state <- function(precomp,
                                                 alpha,
                                                 denominator_floor,
                                                 min_component) {
  list(
    occurrence = precomp$occurrence,
    joint_bins = precomp$joint_bins,
    joint_exposure = precomp$joint_exposure,
    breaks = precomp$breaks,
    widths = precomp$widths,
    x_grid = precomp$x_grid,
    n_grid = precomp$n_grid,
    bandwidth = precomp$bandwidth,
    kernel = precomp$kernel,
    kernel_norm = precomp$kernel_norm,
    smoothed_components = lapply(seq_along(alpha), function(k) {
      .sbf_mult_sparse_smooth_component(precomp, alpha, k)
    }),
    denominator_floor = denominator_floor,
    min_component = min_component
  )
}

# Fixed-point iteration -------------------------------------------------------

# Run the multiplicative binning iteration:
# initialize positive alpha, update components 0..p sequentially by the ratio
# formulas, normalize all covariate scales, then check absolute change.
.sbf_mult_sparse_fit_engine <- function(precomp,
                                        iterations,
                                        convergence_tol,
                                        min_component,
                                        identification = "sample_mean",
                                        truth_functions = NULL,
                                        initial = NULL) {
  d <- length(precomp$n_grid)
  alpha <- if (is.null(initial)) {
    lapply(precomp$n_grid, function(n) rep(1, n))
  } else {
    initial
  }

  warm_D <- unlist(lapply(seq_len(d), function(k) {
    .sbf_mult_sparse_component_denominator(
      precomp = precomp,
      alpha = alpha,
      component = k,
      min_component = min_component
    )
  }), use.names = FALSE)
  denominator_floor <- .sbf_mult_sparse_denominator_floor(warm_D)

  converged <- FALSE
  final_delta <- NA_real_
  iterations_used <- 0L
  component_floor_count <- 0L
  raw_denominators <- list()
  denominator_floor_count <- 0L

  for (iter in seq_len(iterations)) {
    old <- alpha
    for (component in seq_len(d)) {
      O <- .sbf_mult_sparse_component_occurrence(precomp, component)
      D_raw <- .sbf_mult_sparse_component_denominator(
        precomp = precomp,
        alpha = alpha,
        component = component,
        min_component = min_component
      )
      D <- .sbf_mult_sparse_guard_positive(D_raw, denominator_floor)
      raw_denominators[[length(raw_denominators) + 1L]] <- D_raw
      denominator_floor_count <- denominator_floor_count + D$adjusted_count
      values <- O / D$value
      floor_idx <- !is.finite(values) | values < min_component
      component_floor_count <- component_floor_count + sum(floor_idx)
      values[floor_idx] <- min_component
      alpha[[component]] <- values
    }

    alpha <- .sbf_mult_sparse_normalize_alpha(
      precomp = precomp,
      alpha = alpha,
      min_component = min_component,
      identification = identification,
      truth_functions = truth_functions
    )
    final_delta <- .sbf_fit_delta(old, alpha)
    iterations_used <- iter
    if (is.finite(final_delta) && final_delta <= convergence_tol) {
      converged <- TRUE
      break
    }
  }

  diagnostics <- .sbf_mult_sparse_fit_diagnostics(
    precomp = precomp,
    raw_denominators = raw_denominators,
    denominator_floor = denominator_floor,
    denominator_floor_count = denominator_floor_count,
    component_floor_count = component_floor_count
  )

  list(
    alpha = alpha,
    converged = converged,
    final_delta = final_delta,
    iterations_used = iterations_used,
    denominator_floor = denominator_floor,
    diagnostics = diagnostics
  )
}


# Internal multiplicative binned-fit adapter used by the public API wrapper.
.sbf_mult_sparse_fit <- function(data,
                                 bandwidth,
                                 time_bins,
                                 covariate_bins,
                                 time_binning_method = "quantile",
                                 covariate_binning_method = "quantile",
                                 representative = "midpoint",
                                 formula = NULL,
                                 iterations = 200L,
                                 kernel = "epanechnikov",
                                 initial = NULL,
                                 convergence_tol = 0.001,
                                 min_component = 1e-8,
                                 identification = "sample_mean",
                                 truth_functions = NULL,
                                 warn_nonconvergence = TRUE,
                                 warn_diagnostics = TRUE) {
  precomp <- .sbf_mult_sparse_build_fit_precompute(
    data = data,
    bandwidth = bandwidth,
    time_bins = time_bins,
    covariate_bins = covariate_bins,
    time_binning_method = time_binning_method,
    covariate_binning_method = covariate_binning_method,
    representative = representative,
    kernel = kernel
  )
  engine <- .sbf_mult_sparse_fit_engine(
    precomp = precomp,
    iterations = iterations,
    initial = initial,
    convergence_tol = convergence_tol,
    min_component = min_component,
    identification = identification,
    truth_functions = truth_functions
  )
  alpha <- engine$alpha

  fit <- list(
    alpha_backfit = alpha,
    x.grid = precomp$x_grid,
    converged = isTRUE(engine$converged),
    iterations_used = engine$iterations_used,
    iterations_max = iterations,
    final_delta = engine$final_delta,
    convergence_tol = convergence_tol,
    model_type = "multiplicative",
    formula = formula,
    feature_names = precomp$feature_names,
    data = data,
    component_info = precomp$component_info,
    fit_settings = list(
      bandwidth = precomp$bandwidth,
      iterations = iterations,
      convergence_tol = convergence_tol,
      kernel = precomp$kernel,
      kernel_name = precomp$kernel_name,
      min_component = min_component,
      identification = identification
    ),
    binning = list(
      spec = precomp$binning_spec,
      n_fit_rows = nrow(data),
      sparse_joint_cells = precomp$sparse_joint_cells,
      joint_exposure_nonzero_share = precomp$joint_exposure_nonzero_share,
      joint_exposure_memory_mb = precomp$joint_exposure_memory_mb,
      metadata = precomp$bin_metadata
    ),
    fit_diagnostics = engine$diagnostics,
    multiplicative_sparse_predict_state = .sbf_mult_sparse_build_predict_state(
      precomp = precomp,
      alpha = alpha,
      denominator_floor = engine$denominator_floor,
      min_component = min_component
    )
  )
  class(fit) <- unique(c(
    "sbf_multiplicative_sparse_binning_fit",
    "sbf_multiplicative_fit",
    class(fit)
  ))
  if (warn_nonconvergence && !isTRUE(fit$converged)) {
    .sbf_emit_nonconvergence_warning(
      model_type = "multiplicative binned",
      iterations_used = fit$iterations_used,
      iterations_max = fit$iterations_max,
      final_delta = fit$final_delta,
      convergence_tol = fit$convergence_tol
    )
  }
  .sbf_warn_binning_fit_diagnostics(
    fit,
    label = "Multiplicative binned fit",
    warn_diagnostics = warn_diagnostics
  )
  fit
}
