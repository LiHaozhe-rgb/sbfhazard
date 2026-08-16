# Shared kernel helpers for SBF and active binning methods.
#
# Methods share one kernel interface: validate a named/custom kernel once, then
# evaluate bandwidth-scaled kernel matrices consistently.

# Normalize a kernel argument into display metadata plus the callable kernel.
.sbf_kernel_spec <- function(kernel = "epanechnikov",
                            arg_name = "kernel") {
  if (is.function(kernel)) {
    return(list(
      name = "custom",
      kernel = kernel
    ))
  }

  if (!is.character(kernel) || length(kernel) != 1L || is.na(kernel)) {
    .sbf_stop_bad_arg(arg_name, "the string 'epanechnikov' or a kernel function")
  }
  kernel_name <- match.arg(tolower(trimws(kernel)), choices = "epanechnikov")

  list(
    name = kernel_name,
    kernel = function(u) 0.75 * (1 - u^2) * (abs(u) < 1)
  )
}

# Evaluate K(u / h) / h for matrix-like offsets, including custom kernels that
# may not be vectorized over matrices.
.sbf_eval_kernel <- function(u, bandwidth, kernel) {
  ku <- kernel(u / bandwidth)
  if (length(ku) != length(u)) {
    ku <- apply(u / bandwidth, 1:2, kernel)
  }
  matrix(ku, nrow = nrow(u), ncol = ncol(u)) / bandwidth
}

# Record fit kernel normalization summaries for fit_diagnostics.
.sbf_kernel_norm_diagnostics <- function(kernel_norm,
                                         floor_value = 1e-12) {
  vals <- as.numeric(unlist(kernel_norm, use.names = FALSE))
  finite <- vals[is.finite(vals)]
  adjusted <- !is.finite(vals) | vals < floor_value
  data.frame(
    kernel_norm_min = if (length(finite) > 0L) min(finite) else NA_real_,
    kernel_norm_p05 = if (length(finite) > 0L) as.numeric(stats::quantile(finite, probs = 0.05, names = FALSE, type = 8)) else NA_real_,
    kernel_norm_floor = as.numeric(floor_value),
    kernel_norm_zero_count = as.integer(sum(is.finite(vals) & vals <= 0)),
    kernel_norm_adjusted_count = as.integer(sum(adjusted)),
    stringsAsFactors = FALSE
  )
}
