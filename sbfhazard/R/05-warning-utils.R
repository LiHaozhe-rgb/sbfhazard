# Shared diagnostic warning helpers.
#
# These helpers turn existing fit/prediction diagnostics into concise warnings.
# They do not recompute denominator, kernel, exposure, or model quantities.

# Sum one numeric diagnostics field, ignoring missing or non-finite entries.
.sbf_diag_sum <- function(x, name) {
  vals <- as.numeric(x[[name]])
  sum(vals[is.finite(vals)])
}

# Return the first finite value from a scalar diagnostics field.
.sbf_first_finite <- function(x) {
  vals <- as.numeric(x)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) {
    return(NA_real_)
  }
  vals[1L]
}

# Emit one combined warning message from the collected diagnostics messages.
.sbf_emit_diagnostic_warnings <- function(messages) {
  if (length(messages) > 0L) {
    warning(paste(messages, collapse = " "), call. = FALSE)
  }
  invisible(messages)
}

# Warn when a fit stops without satisfying the convergence tolerance.
.sbf_emit_nonconvergence_warning <- function(model_type,
                                            iterations_used,
                                            iterations_max,
                                            final_delta,
                                            convergence_tol) {
  warning(
    sprintf(
      paste(
        "%s fit did not converge within %d iterations.",
        "iterations_used=%d, final_delta=%s, tolerance=%s.",
        "Consider increasing `iterations`, checking bandwidth settings"
      ),
      model_type,
      iterations_max,
      iterations_used,
      format(signif(final_delta, 6), trim = TRUE, scientific = FALSE),
      format(signif(convergence_tol, 6), trim = TRUE, scientific = FALSE)
    ),
    call. = FALSE
  )
}

# Warn about fit-time binning issues that affect support or denominators.
.sbf_warn_binning_fit_diagnostics <- function(fit,
                                             label = "Binned fit",
                                             warn_diagnostics = TRUE) {
  if (!isTRUE(warn_diagnostics)) {
    return(invisible(character()))
  }

  diagnostics <- fit$fit_diagnostics

  messages <- character()
  denom_adjusted <- .sbf_first_finite(diagnostics$denom_adjusted_count)
  if (is.finite(denom_adjusted) && denom_adjusted > 0) {
    messages <- c(messages, sprintf(
      "%s used denominator flooring for %d entries; see fit diagnostics.",
      label,
      as.integer(denom_adjusted)
    ))
  }

  metadata <- fit$binning$metadata
  grid_mode <- if ("grid_mode" %in% names(metadata)) {
    as.character(metadata$grid_mode)
  } else {
    rep("binned", nrow(metadata))
  }

  unique_grid <- which(grid_mode == "unique")
  if (length(unique_grid) > 0L) {
    messages <- c(messages, sprintf(
      paste(
        "%s used distinct-value grids because requested bins reached",
        "the number of distinct observed values for: %s."
      ),
      label,
      paste(as.character(metadata$dimension[unique_grid]), collapse = ", ")
    ))
  }

  binned_grid <- grid_mode == "binned" | is.na(grid_mode)
  reduced_grid <- which(
    binned_grid &
      metadata$effective_bins < metadata$bins_requested
  )
  if (length(reduced_grid) > 0L) {
    messages <- c(messages, sprintf(
      "%s stored fewer bins than requested after repeated boundaries were removed for: %s.",
      label,
      paste(as.character(metadata$dimension[reduced_grid]), collapse = ", ")
    ))
  }

  empty_grid <- which(binned_grid & metadata$empty_bins > 0L)
  empty_bins <- sum(as.numeric(metadata$empty_bins[binned_grid]), na.rm = TRUE)
  if (empty_bins > 0) {
    messages <- c(messages, sprintf(
      "%s stored %d bins without observed values for: %s.",
      label,
      as.integer(empty_bins),
      paste(as.character(metadata$dimension[empty_grid]), collapse = ", ")
    ))
  }

  .sbf_emit_diagnostic_warnings(messages)
}

# Warn about binning prediction diagnostics stored in list output or attributes.
.sbf_warn_binning_prediction_diagnostics <- function(prediction,
                                                    label = "Binned prediction",
                                                    warn_diagnostics = TRUE) {
  if (!isTRUE(warn_diagnostics)) {
    return(invisible(character()))
  }
  if (is.list(prediction)) {
    diagnostics <- prediction$prediction_diagnostics
  } else {
    diagnostics <- attr(prediction, "prediction_diagnostics")
  }

  denom_adjusted <- .sbf_diag_sum(diagnostics, "denom_adjusted_count")
  support_zero <- .sbf_diag_sum(diagnostics, "kernel_support_zero_count")
  fallback <- .sbf_diag_sum(diagnostics, "bin_step_fallback_count")
  clamp_count <- if (is.list(prediction) && !is.null(prediction$clamp_count)) {
    .sbf_first_finite(prediction$clamp_count)
  } else {
    0
  }

  messages <- character()
  if (denom_adjusted > 0) {
    messages <- c(messages, sprintf(
      "%s used denominator flooring for %d query entries; see prediction diagnostics.",
      label,
      as.integer(denom_adjusted)
    ))
  }
  if (support_zero > 0) {
    messages <- c(messages, sprintf(
      "%s had empty kernel support for %d query points; see prediction diagnostics.",
      label,
      as.integer(support_zero)
    ))
  }
  if (fallback > 0) {
    messages <- c(messages, sprintf(
      "%s used bin-step fallback for %d query points; see prediction diagnostics.",
      label,
      as.integer(fallback)
    ))
  }
  if (is.finite(clamp_count) && clamp_count > 0) {
    messages <- c(messages, sprintf(
      "%s clamped %d hazard values; see prediction diagnostics.",
      label,
      as.integer(clamp_count)
    ))
  }

  .sbf_emit_diagnostic_warnings(messages)
}

# Warn about SBF prediction extrapolation, non-finite values, and clamping.
.sbf_warn_fit_prediction_diagnostics <- function(prediction,
                                                  label = "SBF formula prediction",
                                                  warn_diagnostics = TRUE) {
  if (!isTRUE(warn_diagnostics)) {
    return(invisible(character()))
  }

  if (is.list(prediction)) {
    diagnostics <- prediction$prediction_diagnostics
  } else {
    diagnostics <- attr(prediction, "prediction_diagnostics")
  }

  nonfinite_count <- .sbf_first_finite(diagnostics$nonfinite_count)
  extrapolated_count <- .sbf_first_finite(diagnostics$extrapolated_count)
  clamp_count <- .sbf_first_finite(diagnostics$clamp_count)

  messages <- character()
  if (is.finite(extrapolated_count) && extrapolated_count > 0) {
    messages <- c(messages, sprintf(
      "%s evaluated %d points outside effective support.",
      label,
      as.integer(extrapolated_count)
    ))
  }
  if (is.finite(nonfinite_count) && nonfinite_count > 0) {
    messages <- c(messages, sprintf(
      "%s returned %d non-finite values.",
      label,
      as.integer(nonfinite_count)
    ))
  }
  if (is.finite(clamp_count) && clamp_count > 0) {
    messages <- c(messages, sprintf(
      "%s clamped %d values.",
      label,
      as.integer(clamp_count)
    ))
  }

  .sbf_emit_diagnostic_warnings(messages)
}
