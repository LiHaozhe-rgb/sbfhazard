# Additive pairwise fit precompute builders for binning.

# Build W_i(I_0^s): each subject's at-risk overlap with time bin s.
.sbf_add_pairwise_time_overlap <- function(time, lower, upper) {
  time <- as.numeric(time)
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  time_mat <- matrix(time, nrow = length(time), ncol = length(lower))
  upper_mat <- matrix(upper, nrow = length(time), ncol = length(upper), byrow = TRUE)
  clipped_upper <- pmin(time_mat, upper_mat)
  lower_mat <- matrix(lower, nrow = length(time), ncol = length(lower), byrow = TRUE)
  R0 <- clipped_upper - lower_mat
  R0[!is.finite(R0) | R0 < 0] <- 0
  R0
}

# Build R_k^r(I_0^s) by summing W_i(I_0^s) within component-k bin r.
.sbf_add_pairwise_marginal_exposure_matrix <- function(bin_index,
                                                       overlap,
                                                       n_bins) {
  overlap <- as.matrix(overlap)
  out <- matrix(0, nrow = n_bins, ncol = ncol(overlap))
  grouped <- rowsum(overlap, group = as.integer(bin_index), reorder = FALSE)
  idx <- as.integer(rownames(grouped))
  out[idx, ] <- grouped
  out
}

# Sum total at-risk time for each pair of covariate bins.
.sbf_add_pairwise_pairwise_exposure <- function(bin_index, overlap, n_grid) {
  d <- length(n_grid)
  pairwise_sum <- vector("list", d)
  total_cells <- 0
  nonzero_cells <- 0

  for (k in seq_len(d)) {
    pairwise_sum[[k]] <- vector("list", d)
  }

  subject_exposure <- rowSums(overlap, na.rm = TRUE)
  active_subjects <- is.finite(subject_exposure) & subject_exposure > 0

  if (d >= 3L) {
    for (k in 2:d) {
      for (j in 2:d) {
        if (j == k) {
          next
        }
        # for (r in seq_len(n_grid[k])) {
        #   for (q in seq_len(n_grid[j])) {
        #     pair_exposure[r, q] <- sum(subject_exposure[k_bin == r & j_bin == q])
        #   }
        # }
        
        k_bin <- bin_index[[k]]
        j_bin <- bin_index[[j]]
        n_pair_cells <- n_grid[k] * n_grid[j]
        pair_id <- k_bin + (j_bin - 1L) * n_grid[k]

        pair_exposure <- numeric(n_pair_cells)
        grouped <- rowsum(
          matrix(subject_exposure[active_subjects], ncol = 1L),
          group = pair_id[active_subjects],
          reorder = FALSE
        )
        pair_exposure[as.integer(rownames(grouped))] <- as.numeric(grouped[, 1L])
        pairwise_sum[[k]][[j]] <- matrix(pair_exposure, nrow = n_grid[k], ncol = n_grid[j])

        total_cells <- total_cells + n_pair_cells
        nonzero_cells <- nonzero_cells + sum(pair_exposure > 0)
      }
    }
  }

  list(
    pairwise_sum = pairwise_sum,
    total_cells = total_cells,
    nonzero_cells = nonzero_cells,
    nonzero_share = if (total_cells > 0L) nonzero_cells / total_cells else NA_real_,
    memory_mb = total_cells * 8 / 1024^2
  )
}

# Kernel moment construction --------------------------------------------------

.sbf_add_pairwise_kernel_correction_constants <- function(grid,
                                                widths,
                                                bandwidth,
                                                kernel,
                                                kernel_correction = TRUE) {
  eps <- 1e-12
  if (!isTRUE(kernel_correction)) {
    constants <- rep(1, length(grid))
    return(list(
      constants = constants,
      raw_constants = constants,
      zero_count = 0L
    ))
  }

  diff_mat <- outer(grid, grid, FUN = "-")
  K0 <- .sbf_eval_kernel(diff_mat, bandwidth, kernel) # source-bin-by-source-bin kernel matrix for the component grid
  K0[!is.finite(K0)] <- 0
  raw_constants <- colSums(sweep(K0, 1L, widths, "*")) # source-bin vector of raw correction constants
  zero <- !is.finite(raw_constants) | raw_constants < eps
  constants <- raw_constants
  constants[zero] <- eps
  list(
    constants = constants,
    raw_constants = raw_constants,
    zero_count = sum(zero)
  )
}

# Collapse per-component boundary normalizers into one fit diagnostics row.
.sbf_add_pairwise_kernel_norm_diagnostics <- function(correction_info) {
  .sbf_kernel_norm_diagnostics(
    kernel_norm = lapply(correction_info, function(info) info$raw_constants),
    floor_value = 1e-12
  )
}

# Build the two fixed kernel moment blocks used by additive binning:
# update moments K0/K1/K2 for backfitting, and source smoothers S0/S1 for
# other-component terms.  Both blocks share the same raw grid kernel matrix.
.sbf_add_pairwise_kernel_moments <- function(grid,
                                             widths,
                                             bandwidth,
                                             kernel,
                                             correction_constants) {
  diff_mat <- outer(grid, grid, FUN = "-")
  K_raw <- .sbf_eval_kernel(diff_mat, bandwidth, kernel)
  K_raw[!is.finite(K_raw)] <- 0

  K0 <- sweep(K_raw, 2L, correction_constants, "/")
  K0[!is.finite(K0)] <- 0
  update <- list(
    K0 = K0,
    K1 = diff_mat * K0,
    K2 = diff_mat^2 * K0
  )

  weighted <- sweep(K_raw, 2L, widths, "*")
  denom <- rowSums(weighted)
  denom[!is.finite(denom) | abs(denom) < 1e-12] <- 1e-12
  S0 <- sweep(weighted, 1L, denom, "/")
  smoother <- list(
    S0 = S0,
    S1 = S0 * (-diff_mat)
  )

  list(update = update, smoother = smoother)
}

# Precompute -------------------------------------------------------------

# Main precompute builder.  It constructs all fixed objects
# fitting: O, R, R_{k,j}^{+}, K, S, bin metadata, and support diagnostics.  It

.sbf_add_pairwise_build_fit_precompute <- function(data,
                                                 bandwidth,
                                                 time_bins = NULL,
                                                 covariate_bins,
                                                 time_binning_method = "quantile",
                                                 covariate_binning_method = "quantile",
                                                 representative = "midpoint",
                                                 kernel = "epanechnikov",
                                                 kernel_correction = TRUE) {
  
  # Normalize the binning spec. The public API has already prepared numeric fit data.
  feature_names <- setdiff(names(data), c("time", "status"))
  spec <- .sbf_binning_prepare_spec(
    time_bins = time_bins,
    covariate_bins = covariate_bins,
    time_binning_method = time_binning_method,
    covariate_binning_method = covariate_binning_method,
    representative = representative
  )

  # Prepare variables for summary construction.
  mm <- data
  X <- as.matrix(mm[, seq_len(ncol(mm) - 1L), drop = FALSE])
  status <- as.numeric(mm$status)
  d <- ncol(X)
  component_names <- c("time", feature_names)
  roles <- c("time", rep("covariate", length(feature_names)))

  kernel_spec <- .sbf_kernel_spec(kernel, arg_name = "kernel")
  kernel <- kernel_spec$kernel

  bin_results <- vector("list", d)
  bin_results[[1L]] <- .sbf_binning_bin_numeric_vector(
    X[, 1L],
    bins = spec$time_bins,
    rule = spec$time_method,
    representative = spec$representative
  )
  if (d >= 2L) {
    for (k in 2:d) {
      bin_results[[k]] <- .sbf_binning_bin_numeric_vector(
        X[, k],
        bins = spec$covariate_bins,
        rule = spec$covariate_method,
        representative = spec$representative
      )
    }
  }

  x_grid <- lapply(bin_results, function(x) as.numeric(x$representatives))
  bin_index <- lapply(bin_results, function(x) as.integer(x$bin_index))
  breaks <- lapply(bin_results, function(result) result$breaks)
  widths <- lapply(bin_results, function(result) result$widths)
  n_grid <- vapply(x_grid, length, integer(1))

  time_breaks <- breaks[[1L]]
  lower <- head(time_breaks, -1L)
  upper <- tail(time_breaks, -1L)
  overlap <- .sbf_add_pairwise_time_overlap(as.numeric(X[, 1L]), lower, upper)

  # Event occurrence summaries O0s/Ojr 
  occurrence <- vector("list", d)
  occurrence[[1L]] <- as.numeric(tapply(
    status,
    factor(bin_index[[1L]], levels = seq_len(n_grid[1L])),
    sum,
    default = 0
  ))

  # Marginal exposure summaries R_j used by denominator terms D_j,a.
  exposure_moments <- vector("list", d)
  exposure_moments[[1L]] <- list(
    R0 = as.numeric(colSums(overlap))
  )

  # Compute marginal exposure and occurrence summaries for each covariate.
  if (d >= 2L) {
    for (k in 2:d) {
      occurrence[[k]] <- as.numeric(tapply(
        status,
        factor(bin_index[[k]], levels = seq_len(n_grid[k])),
        sum,
        default = 0
      ))
      R0k <- .sbf_add_pairwise_marginal_exposure_matrix(bin_index[[k]], overlap, n_grid[k])
      exposure_moments[[k]] <- list(R0 = R0k)
    }
  }

  pairwise <- .sbf_add_pairwise_pairwise_exposure(bin_index = bin_index, overlap = overlap, n_grid = n_grid)
  
  # Compute K_j,a update moments and S_j
  correction_info <- lapply(seq_len(d), function(k) {
    .sbf_add_pairwise_kernel_correction_constants(
      grid = x_grid[[k]],
      widths = widths[[k]],
      bandwidth = bandwidth[k],
      kernel = kernel,
      kernel_correction = isTRUE(kernel_correction)
    )
  })
  kernel_correction_constants <- lapply(correction_info, function(info) info$constants)
  kernel_norm_diagnostics <- .sbf_add_pairwise_kernel_norm_diagnostics(correction_info)
  kernel_moments <- lapply(seq_len(d), function(k) {
    .sbf_add_pairwise_kernel_moments(
      x_grid[[k]],
      widths[[k]],
      bandwidth[k],
      kernel,
      correction_constants = kernel_correction_constants[[k]]
    )
  })
  update_moments <- lapply(kernel_moments, function(moment) moment$update)
  smoother_moments <- lapply(kernel_moments, function(moment) moment$smoother)

  pairwise_diagnostics <- list(
    total_cells = pairwise$total_cells,
    nonzero_cells = pairwise$nonzero_cells,
    nonzero_share = pairwise$nonzero_share,
    memory_mb = pairwise$memory_mb
  )

  support_min <- vapply(bin_results, function(result) result$support_min, numeric(1))
  support_max <- vapply(bin_results, function(result) result$support_max, numeric(1))
  effective_support_min <- vapply(bin_results, function(result) result$effective_support_min, numeric(1))
  effective_support_max <- vapply(bin_results, function(result) result$effective_support_max, numeric(1))
  component_info <- .sbf_component_info_from_grid(
    x_grid = x_grid,
    component_names = component_names,
    roles = roles,
    support_min = support_min,
    support_max = support_max,
    effective_support_min = effective_support_min,
    effective_support_max = effective_support_max
  )
  bin_metadata <- do.call(rbind, lapply(seq_len(d), function(k) {
    data.frame(
      dimension = component_names[k],
      role = roles[k],
      binned = TRUE,
      grid_mode = bin_results[[k]]$grid_mode,
      rule = if (k == 1L) bin_results[[k]]$rule else spec$covariate_method,
      representative = spec$representative,
      bins_requested = if (k == 1L) spec$time_bins else spec$covariate_bins,
      effective_bins = bin_results[[k]]$effective_bins,
      support_min = bin_results[[k]]$support_min,
      support_max = bin_results[[k]]$support_max,
      effective_support_min = bin_results[[k]]$effective_support_min,
      effective_support_max = bin_results[[k]]$effective_support_max,
      empty_bins = bin_results[[k]]$empty_bins,
      stringsAsFactors = FALSE
    )
  }))
  rownames(bin_metadata) <- NULL
  time_exposure <- as.numeric(exposure_moments[[1L]]$R0)

  list(
    feature_names = feature_names,
    binning_spec = spec,
    kernel = kernel,
    kernel_name = kernel_spec$name,
    kernel_correction = isTRUE(kernel_correction),
    kernel_correction_constants = kernel_correction_constants,
    kernel_norm_diagnostics = kernel_norm_diagnostics,
    x_grid = x_grid,
    n_grid = n_grid,
    breaks = breaks,
    widths = widths,
    occurrence = occurrence,
    exposure_moments = exposure_moments,
    pair_exposure0 = pairwise$pairwise_sum,
    pairwise_diagnostics = pairwise_diagnostics,
    time_exposure = time_exposure,
    update_moments = update_moments,
    smoother_moments = smoother_moments,
    component_info = component_info,
    bin_metadata = bin_metadata
  )
}
