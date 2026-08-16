# Shared internal helpers for prediction modules.

.sbf_assert_required_fields <- function(x,
                                        fields,
                                        label = "result") {
  if (!is.list(x)) {
    .sbf_stop_bad_arg(label, "a list", actual = sprintf("class=%s", paste(class(x), collapse = "/")))
  }
  object_names <- names(x)
  missing_fields <- fields[vapply(fields, function(field) {
    !(field %in% object_names) || is.null(x[[field]])
  }, logical(1))]
  if (length(missing_fields) > 0L) {
    .sbf_stop_bad_arg(label, sprintf("an object containing %s", paste(fields, collapse = ", ")), actual = sprintf("missing %s", paste(missing_fields, collapse = ", ")))
  }
  invisible(x)
}

# Support and extrapolation ---------------------------------------------------

# Return the baseline time grid used when prediction times are omitted.
.sbf_finite_baseline_grid <- function(x_grid) {
  if (!is.list(x_grid) || length(x_grid) < 1L) {
    .sbf_stop_bad_arg("x_grid", "a non-empty list of prediction grids")
  }
  baseline_grid <- as.numeric(x_grid[[1L]])
  baseline_grid <- sort(unique(baseline_grid[is.finite(baseline_grid)]))
  if (length(baseline_grid) == 0L) {
    stop("No finite baseline time grid.", call. = FALSE)
  }
  baseline_grid
}

# Return support bounds for one component from fit metadata.
.sbf_component_support <- function(result,
                                  component_index) {
  row_idx <- which(result$component_info$index == component_index)
  if (length(row_idx) == 0L) {
    .sbf_stop_bad_arg("result", "component_info covering the requested component")
  }
  row <- result$component_info[row_idx[1L], , drop = FALSE]
  list(
    support_min = as.numeric(row$support_min[1L]),
    support_max = as.numeric(row$support_max[1L]),
    effective_support_min = as.numeric(row$effective_support_min[1L]),
    effective_support_max = as.numeric(row$effective_support_max[1L])
  )
}

# Count queries outside effective fitted support.
.sbf_count_prediction_support_violations <- function(result,
                                                    times,
                                                    cov_data) {
  count_component <- function(component_index, xout) {
    info <- .sbf_component_support(result, component_index = component_index)
    xout <- as.numeric(xout)
    tol <- max(.Machine$double.eps^0.5, 1e-8)
    sum(xout < (info$effective_support_min - tol) |
          xout > (info$effective_support_max + tol), na.rm = TRUE)
  }

  support_count <- count_component(1L, times)
  if (ncol(cov_data) > 0L) {
    for (j in seq_len(ncol(cov_data))) {
      support_count <- support_count + count_component(j + 1L, cov_data[, j])
    }
  }
  support_count
}

# Attach support metadata and extrapolation flags.
.sbf_attach_component_support <- function(values,
                                         xout,
                                         support_info) {
  xout <- as.numeric(xout)
  tol <- max(.Machine$double.eps^0.5, 1e-8)
  extrapolated <- xout < (support_info$effective_support_min - tol) |
    xout > (support_info$effective_support_max + tol)
  attr(values, "support_min") <- support_info$support_min
  attr(values, "support_max") <- support_info$support_max
  attr(values, "effective_support_min") <- support_info$effective_support_min
  attr(values, "effective_support_max") <- support_info$effective_support_max
  attr(values, "extrapolated") <- extrapolated
  values
}

# Interpolate fitted grid values at new evaluation points.
.sbf_interpolate_grid_values <- function(grid,
                                         values,
                                         xout) {
  grid <- as.numeric(grid)
  values <- as.numeric(values)
  xout <- as.numeric(xout)

  keep <- is.finite(grid) & is.finite(values)
  grid <- grid[keep]
  values <- values[keep]
  ord <- order(grid)
  grid <- grid[ord]
  values <- values[ord]
  unique_idx <- !duplicated(grid)
  grid <- grid[unique_idx]
  values <- values[unique_idx]

  if (length(grid) < 2L) {
    stop("Interpolation prediction requires at least two finite grid values.", call. = FALSE)
  }

  out <- numeric(length(xout))
  for (i in seq_along(xout)) {
    x <- xout[i]
    if (x <= grid[1L]) {
      out[i] <- values[1L]
      next
    }
    if (x >= grid[length(grid)]) {
      out[i] <- values[length(values)]
      next
    }

    index <- which.min(abs(x - grid))
    error <- grid[index] - x
    if (error == 0) {
      out[i] <- values[index]
      next
    }
    index2 <- if (error > 0) index - 1L else index + 1L
    error2 <- grid[index2] - x
    # out[i] <- (abs(error) * values[index] + abs(error2) * values[index2]) /
    #   (abs(error) + abs(error2)) # original source weighting
    # With the endpoint rules above, this equals stats::approx(..., rule = 2).
    out[i] <- (abs(error2) * values[index] + abs(error) * values[index2]) /
      (abs(error) + abs(error2))
  }

  out
}

# Prediction diagnostics ------------------------------------------------------

# Summarize nonfinite, extrapolation, and clamp counts for SBF predictions.
.sbf_prediction_diagnostics <- function(values,
                                       extrapolated_count = 0L,
                                       clamp_count = 0L) {
  vals <- as.numeric(unlist(values, use.names = FALSE))
  list(
    extrapolated_count = as.integer(extrapolated_count),
    nonfinite_count = as.integer(sum(!is.finite(vals))),
    clamp_count = as.integer(clamp_count)
  )
}

# Prediction output builders --------------------------------------------------

# Build the shared list shape returned by hazard prediction.
.sbf_build_hazard_output <- function(hazard,
                                    times,
                                    baseline,
                                    covariate_term,
                                    covariate_components,
                                    clamp_count) {
  list(
    hazard = hazard,
    times = times,
    baseline = baseline,
    covariate_term = covariate_term,
    covariate_components = covariate_components,
    clamp_count = as.integer(clamp_count)
  )
}

# Build the shared list shape returned by survival prediction.
.sbf_build_survival_output <- function(survival,
                                      times,
                                      covariate_term,
                                      clamp_count,
                                      integration_start,
                                      baseline_support,
                                      covariate_components = NULL) {
  list(
    survival = survival,
    times = times,
    covariate_term = covariate_term,
    covariate_components = covariate_components,
    clamp_count = as.integer(clamp_count),
    integration_start = integration_start,
    baseline_support_min = baseline_support$support_min,
    baseline_support_max = baseline_support$support_max,
    baseline_effective_support_min = baseline_support$effective_support_min,
    baseline_effective_support_max = baseline_support$effective_support_max
  )
}
