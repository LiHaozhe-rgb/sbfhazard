# Multiplicative binning fit precompute builders.

# Kernel and joint-grid helpers ----------------------------------------------

# Build normalized kernel weights from source bins to target points.
.sbf_mult_sparse_kernel_source_target <- function(source, target, widths, bandwidth, kernel) {
  source <- as.numeric(source)
  target <- as.numeric(target)
  widths <- as.numeric(widths)
  K_grid <- .sbf_eval_kernel(outer(source, source, FUN = "-"), bandwidth, kernel)
  normalizer_raw <- rowSums(sweep(K_grid, 2L, widths, "*"))
  normalizer <- normalizer_raw
  normalizer[!is.finite(normalizer) | normalizer <= 1e-12] <- 1e-12
  K_new <- .sbf_eval_kernel(outer(source, target, FUN = "-"), bandwidth, kernel)
  list(
    weights = sweep(K_new, 1L, normalizer, "/"),
    normalizer_raw = normalizer_raw
  )
}

# Encode each subject's observed covariate-bin pattern as a sparse joint-cell id.
.sbf_mult_sparse_joint_key <- function(cov_bin_index) {
  if (length(cov_bin_index) == 1L) {
    idx <- as.integer(cov_bin_index[[1L]])
    return(factor(idx, levels = sort(unique(idx))))
  }
  cov_codes <- lapply(cov_bin_index, as.integer)
  do.call(interaction, c(cov_codes, list(drop = TRUE, lex.order = TRUE)))
}

# Build joint-pattern summaries used by the multiplicative fit and predict code.
.sbf_mult_sparse_build_fit_precompute <- function(data,
                                                      bandwidth,
                                                      time_bins,
                                                      covariate_bins,
                                                      time_binning_method = "quantile",
                                                      covariate_binning_method = "quantile",
                                                      representative = "midpoint",
                                                      kernel = "epanechnikov") {

  feature_names <- setdiff(names(data), c("time", "status"))
  spec <- .sbf_binning_prepare_spec(
    time_bins = time_bins,
    covariate_bins = covariate_bins,
    time_binning_method = time_binning_method,
    covariate_binning_method = covariate_binning_method,
    representative = representative
  )

  d <- ncol(data) - 1L
  p <- d - 1L
  component_names <- c("time", feature_names)
  roles <- c("time", rep("covariate", p))
  kernel_spec <- .sbf_kernel_spec(kernel, arg_name = "kernel")
  kernel <- kernel_spec$kernel

  X <- as.matrix(data[, seq_len(d), drop = FALSE])
  status <- as.numeric(data$status)
  bin_results <- vector("list", d)
  for (k in seq_len(d)) {
    bin_results[[k]] <- .sbf_binning_bin_numeric_vector(
      X[, k],
      bins = if (k == 1L) spec$time_bins else spec$covariate_bins,
      rule = if (k == 1L) spec$time_method else spec$covariate_method,
      representative = spec$representative
    )
  }

  x_grid <- lapply(bin_results, function(x) as.numeric(x$representatives))
  bin_index <- lapply(bin_results, function(x) as.integer(x$bin_index))
  breaks <- lapply(bin_results, function(result) result$breaks)
  widths <- lapply(bin_results, function(result) result$widths)
  n_grid <- vapply(x_grid, length, integer(1))

  time_breaks <- breaks[[1L]]
  lower <- head(time_breaks, -1L)
  upper <- tail(time_breaks, -1L)
  overlap <- outer(as.numeric(X[, 1L]), seq_along(lower), function(ti, s) {
    pmax(0, pmin(ti, upper[s]) - lower[s])
  })
  overlap[!is.finite(overlap)] <- 0

  occurrence <- lapply(seq_len(d), function(k) {
    as.numeric(tapply(
      status,
      factor(bin_index[[k]], levels = seq_len(n_grid[k])),
      sum,
      default = 0
    ))
  })

  cov_bin_index <- bin_index[-1L]

  # Observed sparse joint cells g.  joint_id maps subjects to g, while
  # joint_bins[g, j] stores the marginal bin r_j(g) used in the update formulas.
  joint_factor <- .sbf_mult_sparse_joint_key(cov_bin_index)
  joint_levels <- levels(joint_factor)
  joint_id <- as.integer(joint_factor)
  if (p == 1L) {
    joint_bins <- matrix(as.integer(joint_levels), ncol = 1L)
  } else {
    joint_bins <- do.call(rbind, strsplit(joint_levels, ".", fixed = TRUE))
    joint_bins <- matrix(as.integer(joint_bins), nrow = length(joint_levels))
  }
  colnames(joint_bins) <- feature_names
  n_joint <- nrow(joint_bins)

  # E_g(s): total at-risk time in time bin s for subjects in joint cell g.
  joint_exposure <- rowsum(overlap, group = joint_id, reorder = TRUE)

  kernel_info <- lapply(seq_len(d), function(k) {
    .sbf_mult_sparse_kernel_source_target(
      source = x_grid[[k]],
      target = x_grid[[k]],
      widths = widths[[k]],
      bandwidth = bandwidth[k],
      kernel = kernel
    )
  })
  kernel_norm <- lapply(kernel_info, function(info) info$weights)
  kernel_norm_diagnostics <- .sbf_kernel_norm_diagnostics(
    kernel_norm = lapply(kernel_info, function(info) info$normalizer_raw),
    floor_value = 1e-12
  )
  # Counts weight the empirical log-mean normalization in the fit engine.
  observation_counts <- lapply(seq_len(d), function(k) {
    tabulate(bin_index[[k]], nbins = n_grid[k])
  })

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
      rule = if (k == 1L) spec$time_method else spec$covariate_method,
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

  list(
    feature_names = feature_names,
    binning_spec = spec,
    x_grid = x_grid,
    n_grid = n_grid,
    breaks = breaks,
    widths = widths,
    occurrence = occurrence,
    observation_counts = observation_counts,
    joint_bins = joint_bins,
    joint_exposure = joint_exposure,
    sparse_joint_cells = n_joint,
    joint_exposure_nonzero_share = sum(joint_exposure > 0) / length(joint_exposure),
    joint_exposure_memory_mb = length(joint_exposure) * 8 / 1024^2,
    kernel_norm = kernel_norm,
    kernel_norm_diagnostics = kernel_norm_diagnostics,
    bandwidth = as.numeric(bandwidth),
    kernel = kernel,
    kernel_name = kernel_spec$name,
    component_info = component_info,
    bin_metadata = bin_metadata
  )
}
