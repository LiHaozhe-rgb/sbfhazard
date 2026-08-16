# Additive pairwise fit solver for binning.
#
# Numerical guards and identifiability ---------------------------------------

.sbf_add_pairwise_guard_denominator <- function(D) {
  raw <- as.numeric(unlist(D, use.names = FALSE))
  positive <- raw[is.finite(raw) & raw > 0]
  D_floor <- if (length(positive) == 0L) {
    1e-10
  } else {
    max(1e-12, 1e-10 * stats::median(positive))
  }

  guard_denominator <- function(x) {
    x <- as.numeric(x)
    x[!is.finite(x) | x < D_floor] <- D_floor
    x
  }

  list(
    D = lapply(D, guard_denominator),
    D_floor = D_floor
  )
}

# Apply the additive identifiability step from the LC algorithm
.sbf_add_pairwise_center_alpha <- function(alpha, k, widths) {
  w <- as.numeric(widths[[k]])
  shift <- as.numeric(crossprod(alpha[[k]], w) / sum(w))
  alpha[[1L]] <- as.numeric(alpha[[1L]]) + shift
  alpha[[k]] <- as.numeric(alpha[[k]]) - shift
  alpha
}

# Prediction state ------------------------------------------------------------

.sbf_add_pairwise_build_component_nuisance <- function(precomp, alpha, alpha1) {
  smoothed_components <- lapply(seq_along(alpha), function(k) {
    S <- precomp$smoother_moments[[k]]
    as.numeric(S$S0 %*% as.numeric(alpha[[k]]) +
      S$S1 %*% as.numeric(alpha1[[k]]))
  })

  d <- length(precomp$n_grid)
  nuisance <- lapply(precomp$n_grid, numeric)

  for (k in seq_len(d)) {
    for (j in seq_len(d)) {
      if (j == k) {
        next
      }
      contribution <- if (k == 1L && j > 1L) {
        t(precomp$exposure_moments[[j]]$R0) %*% smoothed_components[[j]]
      } else if (k > 1L && j == 1L) {
        precomp$exposure_moments[[k]]$R0 %*% smoothed_components[[1L]]
      } else {
        precomp$pair_exposure0[[k]][[j]] %*% smoothed_components[[j]]
      }
      nuisance[[k]] <- nuisance[[k]] + as.numeric(contribution)
    }
  }
  nuisance
}

# Store the O/R/R^+/K/S summaries needed by eq:add-predict-kernel and the
# one-step prediction formulas, without keeping the full subject-level data.
.sbf_add_pairwise_build_predict_state <- function(precomp,
                                                  alpha,
                                                  alpha1,
                                                  local_constant,
                                                  terms,
                                                  bandwidth,
                                                  kernel) {
  list(
    occurrence = precomp$occurrence,
    time_exposure = precomp$time_exposure,
    exposure_moments = precomp$exposure_moments,
    widths = precomp$widths,
    x_grid = precomp$x_grid,
    n_grid = precomp$n_grid,
    breaks = precomp$breaks,
    bandwidth = as.numeric(bandwidth),
    kernel = kernel,
    kernel_correction_constants = precomp$kernel_correction_constants,
    local_constant = isTRUE(local_constant),
    denominator_floor = terms$denominator$D_floor,
    component_nuisance = .sbf_add_pairwise_build_component_nuisance(
      precomp = precomp,
      alpha = alpha,
      alpha1 = alpha1
    )
  )
}

# Diagnostics -----------------------------------------------------------------

.sbf_add_pairwise_fit_diagnostics <- function(precomp,
                                              fit_result) {
  terms <- fit_result$terms
  alpha <- fit_result$alpha
  alpha1 <- fit_result$alpha1
  denominator_floor <- terms$denominator$D_floor
  value_raw <- if (isTRUE(fit_result$local_constant)) terms$D_raw else terms$D0_raw
  slope_raw <- if (isTRUE(fit_result$local_constant)) NULL else terms$D2_raw
  denominator_diag <- .sbf_binning_denominator_diagnostics(
    value_raw = value_raw,
    slope_raw = slope_raw,
    floor_value = as.numeric(denominator_floor)
  )
  cbind(denominator_diag, precomp$kernel_norm_diagnostics)
}

# Formula-term construction ---------------------------------------------------

# Cross-component exposure block used in both LC and LL B terms:
# - k=1, j>1: baseline update uses R_j^q(I_0^s).
# - k>1, j=1: covariate update uses R_k^r(I_0^s) against baseline alpha_0.
# - k>1, j>1: covariate-covariate update uses R_{k,j}^{r,q,+}.
.sbf_add_pairwise_cross_exposure_block <- function(precomp,
                                                   k,
                                                   j,
                                                   moment) {
  K <- precomp$update_moments[[k]][[moment]]
  if (k == 1L && j > 1L) {
    return(K %*% t(precomp$exposure_moments[[j]]$R0))
  }
  if (k > 1L && j == 1L) {
    return(K %*% precomp$exposure_moments[[k]]$R0)
  }
  if (k > 1L && j > 1L) {
    return(K %*% precomp$pair_exposure0[[k]][[j]])
  }
  stop("Invalid cross-exposure block request.", call. = FALSE)
}

# Build the fixed LC short notation:
#   O_k = K_{k,0} O_k^bin,
#   D_k = K_{k,0} exposure_k,
#   B_{k,0} = cross_exposure_kj S_{j,0} alpha_j.
# Per-iteration work only evaluates alpha_k = (O_k - B_{k,0}) / D_k.
.sbf_add_pairwise_build_lc_update_terms <- function(precomp) {
  d <- length(precomp$n_grid)
  K0 <- lapply(precomp$update_moments, function(moment) moment$K0)
  S0 <- lapply(precomp$smoother_moments, function(moment) moment$S0)

  exposure <- vector("list", d)
  exposure[[1L]] <- precomp$time_exposure
  if (d >= 2L) {
    for (k in 2:d) {
      exposure[[k]] <- rowSums(precomp$exposure_moments[[k]]$R0)
    }
  }
  D_raw <- lapply(seq_len(d), function(k) {
    as.numeric(K0[[k]] %*% exposure[[k]])
  })
  den <- .sbf_add_pairwise_guard_denominator(D_raw)
  D <- den$D

  O <- lapply(seq_len(d), function(k) {
    as.numeric(K0[[k]] %*% precomp$occurrence[[k]])
  })

  B_blocks <- vector("list", d)
  for (k in seq_len(d)) {
    B_blocks[[k]] <- vector("list", d)
  }
  if (d >= 2L) {
    for (k in seq_len(d)) {
      for (j in seq_len(d)) {
        if (j == k) {
          next
        }
        block <- .sbf_add_pairwise_cross_exposure_block(
          precomp = precomp,
          k = k,
          j = j,
          moment = "K0"
        ) %*% S0[[j]]
        B_blocks[[k]][[j]] <- block
      }
    }
  }

  list(
    O = O,
    B_blocks = B_blocks,
    D = D,
    D_raw = D_raw,
    denominator = den
  )
}

# Build the fixed LL short notation:
#   value: alpha_k = (O_{k,0} - C_k alpha_{k,1} - B_{k,0}) / D_{k,0}
#   slope: alpha_{k,1} = (O_{k,1} - C_k alpha_k - B_{k,1}) / D_{k,2}.
# The B blocks multiply smoothed values S0 alpha_j + S1 alpha_{j,1}.
.sbf_add_pairwise_build_ll_update_terms <- function(precomp) {
  d <- length(precomp$n_grid)
  K0 <- lapply(precomp$update_moments, function(moment) moment$K0)
  K1 <- lapply(precomp$update_moments, function(moment) moment$K1)
  K2 <- lapply(precomp$update_moments, function(moment) moment$K2)
  S0 <- lapply(precomp$smoother_moments, function(moment) moment$S0)
  S1 <- lapply(precomp$smoother_moments, function(moment) moment$S1)

  exposure <- vector("list", d)
  exposure[[1L]] <- precomp$time_exposure
  if (d >= 2L) {
    for (k in 2:d) {
      exposure[[k]] <- rowSums(precomp$exposure_moments[[k]]$R0)
    }
  }

  D0_raw <- D2_raw <- C <- vector("list", d)
  O0 <- O1 <- vector("list", d)
  for (k in seq_len(d)) {
    D0_raw[[k]] <- as.numeric(K0[[k]] %*% exposure[[k]])
    D2_raw[[k]] <- as.numeric(K2[[k]] %*% exposure[[k]])
    C[[k]] <- as.numeric(K1[[k]] %*% exposure[[k]])
    O0[[k]] <- as.numeric(K0[[k]] %*% precomp$occurrence[[k]])
    O1[[k]] <- as.numeric(K1[[k]] %*% precomp$occurrence[[k]])
  }

  # Boundary bins can have zero second moments in LL; floor them for a stable sweep.
  den <- .sbf_add_pairwise_guard_denominator(c(D0_raw, D2_raw))
  D0 <- den$D[seq_len(d)]
  D2 <- den$D[d + seq_len(d)]

  B_blocks <- list(
    value_from_value = vector("list", d),
    value_from_slope = vector("list", d),
    slope_from_value = vector("list", d),
    slope_from_slope = vector("list", d)
  )
  for (k in seq_len(d)) {
    B_blocks$value_from_value[[k]] <- vector("list", d)
    B_blocks$value_from_slope[[k]] <- vector("list", d)
    B_blocks$slope_from_value[[k]] <- vector("list", d)
    B_blocks$slope_from_slope[[k]] <- vector("list", d)
  }
  if (d >= 2L) {
    for (k in seq_len(d)) {
      for (j in seq_len(d)) {
        if (j == k) {
          next
        }
        block0 <- .sbf_add_pairwise_cross_exposure_block(
          precomp = precomp,
          k = k,
          j = j,
          moment = "K0"
        )
        block1 <- .sbf_add_pairwise_cross_exposure_block(
          precomp = precomp,
          k = k,
          j = j,
          moment = "K1"
        )
        B_blocks$value_from_value[[k]][[j]] <- block0 %*% S0[[j]]
        B_blocks$value_from_slope[[k]][[j]] <- block0 %*% S1[[j]]
        B_blocks$slope_from_value[[k]][[j]] <- block1 %*% S0[[j]]
        B_blocks$slope_from_slope[[k]][[j]] <- block1 %*% S1[[j]]
      }
    }
  }

  list(
    O0 = O0,
    O1 = O1,
    C = C,
    B_blocks = B_blocks,
    D0 = D0,
    D2 = D2,
    D0_raw = D0_raw,
    D2_raw = D2_raw,
    denominator = den
  )
}

# Iteration entry points ------------------------------------------------------

# Run LC backfitting until convergence or the iteration cap.
.sbf_add_pairwise_fit_lc_engine <- function(precomp,
                                            iterations,
                                            convergence_tol,
                                            initial = NULL) {
  terms <- .sbf_add_pairwise_build_lc_update_terms(precomp)

  alpha <- if (is.null(initial)) {
    lapply(precomp$n_grid, function(n) rep(0, n))
  } else {
    initial
  }
  final_delta <- NA_real_
  iterations_used <- 0L
  for (iter in seq_len(iterations)) {
    old_alpha <- alpha
    for (k in seq_along(alpha)) {
      B <- numeric(length(terms$D[[k]]))
      for (j in seq_along(alpha)) {
        if (j == k) {
          next
        }
        B <- B + as.numeric(terms$B_blocks[[k]][[j]] %*% as.numeric(alpha[[j]]))
      }

      alpha[[k]] <- (terms$O[[k]] - B) / terms$D[[k]]
      alpha[[k]][!is.finite(alpha[[k]])] <- NA_real_

      if (k > 1L) {
        alpha <- .sbf_add_pairwise_center_alpha(alpha, k, precomp$widths)
      }
    }
    iterations_used <- iter
    final_delta <- .sbf_fit_delta(old_alpha, alpha)
    if (is.finite(final_delta) && final_delta <= convergence_tol) {
      break
    }
  }

  alpha1 <- lapply(precomp$n_grid, function(n) rep(0, n))
  finite_alpha <- all(vapply(alpha, function(x) all(is.finite(x)), logical(1)))
  converged <- finite_alpha && is.finite(final_delta) && final_delta <= convergence_tol
  list(
    local_constant = TRUE,
    terms = terms,
    alpha = alpha,
    alpha1 = alpha1,
    iterations_used = iterations_used,
    final_delta = final_delta,
    converged = converged
  )
}

# Run LL backfitting until convergence or the iteration cap.
.sbf_add_pairwise_fit_ll_engine <- function(precomp,
                                            iterations,
                                            convergence_tol,
                                            initial = NULL) {
  terms <- .sbf_add_pairwise_build_ll_update_terms(precomp)

  parts <- list(
    alpha = if (is.null(initial)) lapply(precomp$n_grid, function(n) rep(0, n)) else initial,
    alpha1 = lapply(precomp$n_grid, function(n) rep(0, n))
  )
  final_delta <- NA_real_
  iterations_used <- 0L
  for (iter in seq_len(iterations)) {
    old_parts <- parts
    for (k in seq_along(parts$alpha)) {
      B0 <- numeric(length(terms$D0[[k]]))
      B1 <- numeric(length(terms$D2[[k]]))
      for (j in seq_along(parts$alpha)) {
        if (j == k) {
          next
        }
        alpha_j <- as.numeric(parts$alpha[[j]])
        alpha1_j <- as.numeric(parts$alpha1[[j]])
        B0 <- B0 +
          as.numeric(terms$B_blocks$value_from_value[[k]][[j]] %*% alpha_j) +
          as.numeric(terms$B_blocks$value_from_slope[[k]][[j]] %*% alpha1_j)
        B1 <- B1 +
          as.numeric(terms$B_blocks$slope_from_value[[k]][[j]] %*% alpha_j) +
          as.numeric(terms$B_blocks$slope_from_slope[[k]][[j]] %*% alpha1_j)
      }

      value_update <- (terms$O0[[k]] -
        terms$C[[k]] * parts$alpha1[[k]] -
        B0) / terms$D0[[k]]
      value_update[!is.finite(value_update)] <- NA_real_
      parts$alpha[[k]] <- value_update

      slope_update <- (terms$O1[[k]] -
        terms$C[[k]] * parts$alpha[[k]] -
        B1) / terms$D2[[k]]
      slope_update[!is.finite(slope_update)] <- NA_real_
      parts$alpha1[[k]] <- slope_update

      if (k > 1L) {
        parts$alpha <- .sbf_add_pairwise_center_alpha(parts$alpha, k, precomp$widths)
      }
    }
    iterations_used <- iter
    final_delta <- .sbf_fit_delta(
      c(old_parts$alpha, old_parts$alpha1),
      c(parts$alpha, parts$alpha1)
    )
    if (is.finite(final_delta) && final_delta <= convergence_tol) {
      break
    }
  }

  finite_parts <- all(vapply(parts$alpha, function(x) all(is.finite(x)), logical(1))) &&
    all(vapply(parts$alpha1, function(x) all(is.finite(x)), logical(1)))
  converged <- finite_parts && is.finite(final_delta) && final_delta <= convergence_tol
  list(
    local_constant = FALSE,
    terms = terms,
    alpha = parts$alpha,
    alpha1 = parts$alpha1,
    iterations_used = iterations_used,
    final_delta = final_delta,
    converged = converged
  )
}


# Internal additive binned-fit adapter used by the public API wrapper.
.sbf_add_pairwise_fit <- function(data,
                                  bandwidth,
                                  time_bins,
                                  covariate_bins,
                                  time_binning_method = "quantile",
                                  covariate_binning_method = "quantile",
                                  representative = "midpoint",
                                  formula = NULL,
                                  local_constant = TRUE,
                                  kernel = "epanechnikov",
                                  kernel_correction = TRUE,
                                  initial = NULL,
                                  iterations = 200L,
                                  convergence_tol = 0.001,
                                  warn_nonconvergence = TRUE,
                                  warn_diagnostics = TRUE) {
  precomp <- .sbf_add_pairwise_build_fit_precompute(
    data = data,
    bandwidth = bandwidth,
    time_bins = time_bins,
    covariate_bins = covariate_bins,
    time_binning_method = time_binning_method,
    covariate_binning_method = covariate_binning_method,
    representative = representative,
    kernel = kernel,
    kernel_correction = kernel_correction
  )
  fit_result <- if (isTRUE(local_constant)) {
    .sbf_add_pairwise_fit_lc_engine(
      precomp = precomp,
      initial = initial,
      iterations = iterations,
      convergence_tol = convergence_tol
    )
  } else {
    .sbf_add_pairwise_fit_ll_engine(
      precomp = precomp,
      initial = initial,
      iterations = iterations,
      convergence_tol = convergence_tol
    )
  }

  alpha <- fit_result$alpha
  alpha1 <- fit_result$alpha1
  terms <- fit_result$terms
  diagnostics <- .sbf_add_pairwise_fit_diagnostics(precomp, fit_result)
  kernel <- precomp$kernel
  fit_suffix <- if (isTRUE(fit_result$local_constant)) "lc" else "ll"

  fit <- list(
    alpha_backfit = alpha,
    alpha.1_backfit = alpha1,
    x.grid = precomp$x_grid,
    data = data,
    converged = isTRUE(fit_result$converged),
    iterations_used = fit_result$iterations_used,
    iterations_max = iterations,
    final_delta = fit_result$final_delta,
    convergence_tol = convergence_tol,
    model_type = "additive",
    formula = formula,
    feature_names = precomp$feature_names,
    component_info = precomp$component_info,
    binning = list(
      spec = precomp$binning_spec,
      n_fit_rows = nrow(data),
      pairwise_total_cells = precomp$pairwise_diagnostics$total_cells,
      pairwise_nonzero_cells = precomp$pairwise_diagnostics$nonzero_cells,
      pairwise_nonzero_share = precomp$pairwise_diagnostics$nonzero_share,
      pairwise_memory_mb = precomp$pairwise_diagnostics$memory_mb,
      metadata = precomp$bin_metadata
    ),
    fit_settings = list(
      bandwidth = as.numeric(bandwidth),
      iterations = iterations,
      convergence_tol = convergence_tol,
      kernel = kernel,
      kernel_name = precomp$kernel_name,
      local_constant = isTRUE(fit_result$local_constant)
    ),
    fit_diagnostics = diagnostics,
    additive_pairwise_predict_state = .sbf_add_pairwise_build_predict_state(
      precomp = precomp,
      alpha = alpha,
      alpha1 = alpha1,
      local_constant = fit_result$local_constant,
      terms = terms,
      bandwidth = bandwidth,
      kernel = kernel
    )
  )
  class(fit) <- unique(c(
    sprintf("sbf_additive_pairwise_binning_%s_fit", fit_suffix),
    "sbf_additive_pairwise_binning_fit",
    "sbf_additive_fit",
    class(fit)
  ))
  if (warn_nonconvergence && !isTRUE(fit$converged)) {
    .sbf_emit_nonconvergence_warning(
      model_type = "additive binned",
      iterations_used = fit$iterations_used,
      iterations_max = fit$iterations_max,
      final_delta = fit$final_delta,
      convergence_tol = fit$convergence_tol
    )
  }
  .sbf_warn_binning_fit_diagnostics(
    fit,
    label = "Additive binned fit",
    warn_diagnostics = warn_diagnostics
  )
  fit
}
