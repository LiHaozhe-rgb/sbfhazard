# Additive simulation DGP internals.

.sbf_simulate_additive_covariates <- function(n, d = 5L, rho = 0, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for additive simulation.", call. = FALSE)
  }
  stddev <- rep(1, d - 1)
  cor_mat <- matrix(rho, nrow = d - 1, ncol = d - 1)
  cor_mat[col(cor_mat) == row(cor_mat)] <- 1
  cov_mat <- stddev %*% t(stddev) * cor_mat
  z <- MASS::mvrnorm(n = n, mu = rep(0, d - 1), Sigma = cov_mat, empirical = FALSE)
  2.5 * atan(z) / pi
}

.sbf_additive_rate <- function(z, model = 1L, violate_cox = TRUE) {
  z <- .sbf_as_numeric_matrix(z, arg_name = "z")
  d_cov <- ncol(z)
  phi <- vector("list", d_cov)

  if (violate_cox) {
    for (k in seq_len(d_cov)) {
      if ((k %% 2L) == 1L) {
        phi[[k]] <- function(x) (2 / sqrt(d_cov)) * sin(pi * x)
      } else {
        phi[[k]] <- function(x) -(2 / sqrt(d_cov)) * sin(pi * x)
      }
    }
  } else {
    for (k in seq_len(d_cov)) {
      if ((k %% 2L) == 1L) {
        phi[[k]] <- function(x) -2 * x
      } else {
        phi[[k]] <- function(x) 2 * x
      }
    }
  }

  top <- rep(0, nrow(z))
  for (k in seq_len(d_cov)) {
    top <- top + phi[[k]](z[, k])
  }

  event_rate <- top + 1
  surv_function <- switch(
    as.character(model),
    "1" = function(t) stats::pexp(t, rate = event_rate, lower.tail = FALSE),
    "2" = function(t) VGAM::pmakeham(t, scale = 5, shape = 1, epsilon = event_rate, lower.tail = FALSE),
    stop("Unsupported model. Use 1 or 2.")
  )

  list(top = top, phi = phi, surv_function = surv_function)
}

.sbf_simulate_additive_data <- function(n = 200L,
                                        d = 5L,
                                        rho = 0,
                                        model = 1L,
                                        violate_cox = TRUE,
                                        seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  z <- .sbf_simulate_additive_covariates(n = 10L * n, d = d, rho = rho)
  top <- .sbf_additive_rate(z, model = model, violate_cox = violate_cox)$top

  top <- top + 1
  z <- z[top > 1, , drop = FALSE]
  top <- top[top > 1]

  top <- top[seq_len(n)]
  z <- z[seq_len(n), , drop = FALSE]

  if (model == 1L) {
    time <- stats::rexp(length(top), top)
    censor <- stats::rexp(length(top), top / 1.75)
  } else if (model == 2L) {
    if (!requireNamespace("VGAM", quietly = TRUE)) {
      stop("Package 'VGAM' is required for additive Makeham simulation.", call. = FALSE)
    }
    beta <- 5
    alpha <- 1
    time <- VGAM::rmakeham(n, beta, alpha, epsilon = top)
    censor <- VGAM::rmakeham(n, beta, alpha / 1.75, epsilon = top)
  } else {
    stop("Unsupported model. Use 1 or 2.")
  }

  observed_time <- time * (time <= censor) + censor * (time > censor)
  status <- as.integer(time <= censor)

  data.frame(time = observed_time, status = status, as.data.frame(z))
}
