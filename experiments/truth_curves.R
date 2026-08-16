# Truth-curve definitions for the simulation models (component phi and survival).
# Shared by benchmark and binning experiment scripts.

sbf_true_additive_phi <- function(z, k, d_cov, violate_cox = TRUE) {
  if (isTRUE(violate_cox)) {
    sign_k <- if ((as.integer(k) %% 2L) == 1L) 1 else -1
    return(sign_k * (2 / sqrt(d_cov)) * sin(pi * z))
  }
  if ((as.integer(k) %% 2L) == 1L) {
    return(-2 * z)
  }
  2 * z
}

sbf_true_multiplicative_phi <- function(z, k, violate_cox = TRUE) {
  if (isTRUE(violate_cox)) {
    if ((as.integer(k) %% 2L) == 1L) {
      return(2 * sin(pi * z))
    }
    return(2 * z)
  }
  if ((as.integer(k) %% 2L) == 1L) {
    return(-z)
  }
  z
}

sbf_true_additive_survival <- function(times, z0, model = 1L, violate_cox = TRUE) {
  top <- 0
  for (k in seq_along(z0)) {
    top <- top + sbf_true_additive_phi(z0[[k]], k, length(z0), violate_cox)
  }
  event_rate <- top + 1
  if (as.integer(model) == 1L) {
    return(stats::pexp(times, rate = event_rate, lower.tail = FALSE))
  }
  if (!requireNamespace("VGAM", quietly = TRUE)) {
    stop("Package 'VGAM' is required for additive Makeham truth.", call. = FALSE)
  }
  VGAM::pmakeham(times, scale = 5, shape = 1, epsilon = event_rate, lower.tail = FALSE)
}

sbf_true_multiplicative_survival <- function(times, z0, model = 1L, violate_cox = TRUE) {
  top <- 0
  for (k in seq_along(z0)) {
    top <- top + sbf_true_multiplicative_phi(z0[[k]], k, violate_cox)
  }
  if (as.integer(model) == 1L) {
    return(stats::pexp(times, rate = exp(top), lower.tail = FALSE))
  }
  if (!requireNamespace("VGAM", quietly = TRUE)) {
    stop("Package 'VGAM' is required for multiplicative Makeham truth.", call. = FALSE)
  }
  VGAM::pmakeham(times, scale = 1, shape = exp(top), lower.tail = FALSE)
}
