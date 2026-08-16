# Multiplicative simulation DGP internals.

.sbf_simulate_multiplicative_covariates <- function(n, d = 5L, rho = 0, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for multiplicative simulation.", call. = FALSE)
  }
  stddev <- rep(1, d - 1)
  cor_mat <- matrix(rho, nrow = d - 1, ncol = d - 1)
  cor_mat[col(cor_mat) == row(cor_mat)] <- 1
  cov_mat <- stddev %*% t(stddev) * cor_mat
  z <- MASS::mvrnorm(n = n, mu = rep(0, d - 1), Sigma = cov_mat, empirical = FALSE)
  2.5 * atan(z) / pi
}

.sbf_multiplicative_rate <- function(z, model = 1L, violate_cox = TRUE) {
  z <- .sbf_as_numeric_matrix(z, arg_name = "z")
  d_cov <- ncol(z)
  phi <- vector("list", d_cov)

  if (violate_cox) {
    for (k in seq_len(d_cov)) {
      if ((k %% 2L) == 1L) {
        phi[[k]] <- function(x) 2 * sin(pi * x)
      } else {
        phi[[k]] <- function(x) 2 * x
      }
    }
  } else {
    for (k in seq_len(d_cov)) {
      if ((k %% 2L) == 1L) {
        phi[[k]] <- function(x) -x
      } else {
        phi[[k]] <- function(x) x
      }
    }
  }

  top <- rep(0, nrow(z))
  for (k in seq_len(d_cov)) {
    top <- top + phi[[k]](z[, k])
  }

  surv_function <- switch(
    as.character(model),
    "1" = function(t) stats::pexp(t, rate = exp(top), lower.tail = FALSE),
    "2" = function(t) VGAM::pmakeham(t, scale = 1, shape = exp(top), lower.tail = FALSE),
    stop("Unsupported model. Use 1 or 2.")
  )

  list(top = top, phi = phi, surv_function = surv_function)
}

.sbf_normalize_multiplicative_effect_truth <- function(phi_fun,
                                                       xout,
                                                       observed_x,
                                                       observation_weights = NULL) {
  .sbf_assert_function(phi_fun, "phi_fun")
  xout <- as.numeric(xout)
  observed_x <- as.numeric(observed_x)
  if (any(!is.finite(xout))) {
    .sbf_stop_bad_arg("xout", "a numeric vector with only finite values")
  }
  if (any(!is.finite(observed_x))) {
    .sbf_stop_bad_arg("observed_x", "a numeric vector with only finite values")
  }
  if (is.null(observation_weights)) {
    observation_weights <- rep(1, length(observed_x))
  } else {
    .sbf_assert_numeric_vector(
      observation_weights,
      "observation_weights",
      expected_length = length(observed_x),
      positive = TRUE
    )
  }

  effect <- as.numeric(phi_fun(xout))
  shift <- stats::weighted.mean(as.numeric(phi_fun(observed_x)), w = observation_weights)
  effect <- effect - shift

  effect
}

.sbf_simulate_multiplicative_data <- function(n = 200L,
                                              d = 5L,
                                              rho = 0,
                                              model = 1L,
                                              violate_cox = TRUE,
                                              seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  z <- .sbf_simulate_multiplicative_covariates(n = n, d = d, rho = rho)
  top <- .sbf_multiplicative_rate(z, model = model, violate_cox = violate_cox)$top

  if (model == 1L) {
    time <- stats::rexp(n, exp(top))
    censor <- stats::rexp(n, exp(top) / 1.75)
  } else if (model == 2L) {
    if (!requireNamespace("VGAM", quietly = TRUE)) {
      stop("Package 'VGAM' is required for multiplicative Makeham simulation.", call. = FALSE)
    }
    beta <- 1
    alpha <- exp(top)
    time <- VGAM::rmakeham(n, beta, alpha)
    censor <- VGAM::rmakeham(n, beta, alpha / 1.75)
  } else {
    stop("Unsupported model. Use 1 or 2.")
  }

  observed_time <- time * (time <= censor) + censor * (time > censor)
  status <- as.integer(time <= censor)

  data.frame(time = observed_time, status = status, as.data.frame(z))
}
