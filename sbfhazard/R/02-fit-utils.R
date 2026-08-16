# Shared fitting helpers.

# Grid and integration 

.sbf_prepare_grid_spec <- function(observed_grid,
                                   x_grid = NULL,
                                   n_grid_additional = 0L,
                                   x_min = NULL,
                                   x_max = NULL,
                                   integral_approx = "midd") {
  if (!is.list(observed_grid) || length(observed_grid) == 0L) {
    .sbf_stop_bad_arg("observed_grid", "a non-empty list of numeric vectors")
  }
  d <- length(observed_grid)

  if (is.null(x_grid)) {
    x_grid <- observed_grid
  }
  if (!is.list(x_grid) || length(x_grid) != d) {
    .sbf_stop_bad_arg("x_grid", sprintf("NULL or a list with length %d", d))
  }

  x_grid <- lapply(seq_len(d), function(k) {
    x <- as.numeric(x_grid[[k]])
    if (length(x) == 0L || any(!is.finite(x))) {
      .sbf_stop_bad_arg(sprintf("x_grid[[%d]]", k), "finite numeric grid values")
    }
    sort(x)
  })

  if (is.null(x_min)) {
    x_min <- vapply(x_grid, function(x) x[1L], numeric(1))
  } else {
    x_min <- as.numeric(x_min)
  }
  if (is.null(x_max)) {
    x_max <- vapply(x_grid, function(x) x[length(x)], numeric(1))
  } else {
    x_max <- as.numeric(x_max)
  }

  if (length(x_min) != d || length(x_max) != d ||
      any(!is.finite(c(x_min, x_max))) ||
      any(x_min >= x_max)) {
    .sbf_stop_bad_arg("x_min/x_max", sprintf("finite vectors with length %d and x_min < x_max", d))
  }

  n_grid_additional <- .sbf_assert_scalar_whole_number(
    n_grid_additional,
    "n_grid_additional",
    min_value = 0L
  )

  x_grid_additional <- lapply(seq_len(d), function(k) {
    seq(x_min[k], x_max[k], length.out = n_grid_additional)
  })
  x_grid <- lapply(seq_len(d), function(k) {
    sort(c(x_grid[[k]], x_grid_additional[[k]]))
  })
  n_grid <- lengths(x_grid)
  if (any(n_grid < 2L)) {
    .sbf_stop_bad_arg("x_grid", "at least two grid values per component")
  }

  dx <- lapply(seq_len(d), function(k) {
    grid <- x_grid[[k]]
    gaps <- diff(grid)
    n_k <- length(grid)

    if (identical(integral_approx, "midd")) {
      first <- (gaps[1L] / 2) + (grid[1L] - x_min[k])
      last <- (gaps[n_k - 1L] / 2) + (x_max[k] - grid[n_k])
      if (n_k == 2L) {
        weights <- c(first, last)
      } else {
        middle <- (gaps[-(n_k - 1L)] + gaps[-1L]) / 2
        weights <- c(first, middle, last)
      }
    } else if (identical(integral_approx, "left")) {
      weights <- c(gaps, x_max[k] - grid[n_k])
    } else if (identical(integral_approx, "right")) {
      weights <- c(grid[1L] - x_min[k], gaps)
    } else {
      .sbf_stop_bad_arg("integral_approx", "\"midd\", \"left\", or \"right\"")
    }

    if (length(weights) != n_k || any(!is.finite(weights)) || any(weights < 0)) {
      .sbf_stop_bad_arg(
        sprintf("SBF grid integration weights for component %d", k),
        "finite nonnegative weights"
      )
    }
    weights
  })

  list(
    x_grid = x_grid,
    n_grid = n_grid,
    dx = dx
  )
}

# Build component metadata from stored grids and optional support bounds.
.sbf_component_info_from_grid <- function(x_grid,
                                          component_names,
                                          roles,
                                          support_min = NULL,
                                          support_max = NULL,
                                          effective_support_min = NULL,
                                          effective_support_max = NULL) {
  if (!is.list(x_grid) || length(x_grid) == 0L) {
    .sbf_stop_bad_arg("x_grid", "a non-empty list of numeric grids")
  }
  d <- length(x_grid)
  if (length(component_names) != d || length(roles) != d) {
    .sbf_stop_bad_arg("component_names/roles", sprintf("vectors with length %d", d))
  }

  grid_range <- function(k) {
    x <- as.numeric(x_grid[[k]])
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      c(NA_real_, NA_real_)
    } else {
      range(x)
    }
  }
  defaults <- t(vapply(seq_len(d), grid_range, numeric(2)))

  resolve_bound <- function(values, fallback) {
    if (is.null(values)) {
      return(fallback)
    }
    values <- as.numeric(values)
    if (length(values) != d) {
      .sbf_stop_bad_arg("component support bounds", sprintf("vectors with length %d", d))
    }
    finite <- is.finite(values)
    fallback[finite] <- values[finite]
    fallback
  }

  support_min <- resolve_bound(support_min, defaults[, 1])
  support_max <- resolve_bound(support_max, defaults[, 2])
  effective_support_min <- resolve_bound(effective_support_min, support_min)
  effective_support_max <- resolve_bound(effective_support_max, support_max)

  data.frame(
    index = seq_len(d),
    name = as.character(component_names),
    role = as.character(roles),
    support_min = as.numeric(support_min),
    support_max = as.numeric(support_max),
    effective_support_min = as.numeric(effective_support_min),
    effective_support_max = as.numeric(effective_support_max),
    stringsAsFactors = FALSE
  )
}

# Convergence -----------------------------------------------------------------

.sbf_fit_delta <- function(old, new) {
  delta <- abs(
    as.numeric(unlist(new, use.names = FALSE)) -
      as.numeric(unlist(old, use.names = FALSE))
  )
  if (length(delta) == 0L || any(!is.finite(delta))) {
    return(NA_real_)
  }
  max(delta)
}

# Multiplicative identification ----------------------------------------------

.sbf_assert_multiplicative_truth_functions <- function(truth_functions, p) {
  if (!is.list(truth_functions) || length(truth_functions) != p ||
      !all(vapply(truth_functions, is.function, logical(1)))) {
    .sbf_stop_bad_arg("truth_functions", sprintf("a list of %d truth functions", p))
  }
  invisible(truth_functions)
}

.sbf_multiplicative_identification_scale <- function(alpha_component,
                                                      grid,
                                                      widths,
                                                      identification,
                                                      min_component,
                                                      sample_weights = NULL,
                                                      truth_functions = NULL,
                                                      covariate_index = NULL) {
  alpha_component <- pmax(min_component, as.numeric(alpha_component))

  if (identical(identification, "sample_mean")) {
    eta <- log(alpha_component)
    if (is.null(sample_weights)) {
      return(exp(mean(eta)))
    }
    return(exp(stats::weighted.mean(eta, w = as.numeric(sample_weights))))
  }

  if (identical(identification, "integral")) {
    widths <- as.numeric(widths)
    return(sum(alpha_component * widths) / sum(widths))
  }

  grid <- as.numeric(grid)
  if (identical(identification, "origin")) {
    anchor_idx <- as.integer(length(grid) / 2)
    truth_x <- grid[anchor_idx]
  } else if (identical(identification, "jasa")) {
    anchor_idx <- c(
      which.min(abs(grid + 0.5)),
      which.min(abs(grid - 0.5))
    )
    truth_x <- c(
      min(abs(grid + 0.5)),
      min(abs(grid - 0.5))
    )
  } else {
    .sbf_stop_bad_arg("identification", "\"sample_mean\", \"integral\", \"origin\", or \"jasa\"")
  }
  truth_values <- as.numeric(truth_functions[[covariate_index]](truth_x))
  if (length(truth_values) != length(truth_x) || any(!is.finite(truth_values))) {
    .sbf_stop_bad_arg("truth_functions", "functions returning finite numeric values at identification anchors")
  }
  scale <- exp(mean(log(alpha_component[anchor_idx]) - truth_values))
  if (!is.finite(scale) || scale <= 0) {
    .sbf_stop_bad_arg("identification", "a finite positive scale")
  }
  scale
}
