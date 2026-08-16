# Shared diagnostics helpers for active binning experiments.
#
# Denominator diagnostics use one schema for fit and prediction.  LC and sparse
# fits report value denominators only; LL also supplies slope denominators.

# Summarize raw denominators before flooring; adjusted_count can be supplied by
# callers that already know how many entries were floored in their guard step.
.sbf_binning_denominator_stats <- function(raw,
                                           floor_value = NA_real_,
                                           adjusted_count = NULL) {
  raw <- as.numeric(unlist(raw, use.names = FALSE))
  finite <- raw[is.finite(raw)]
  adjusted <- if (is.null(adjusted_count)) {
    if (is.finite(floor_value)) {
      as.integer(sum(!is.finite(raw) | raw < floor_value))
    } else {
      NA_integer_
    }
  } else {
    as.integer(adjusted_count)
  }

  list(
    min = if (length(finite) > 0L) min(finite) else NA_real_,
    p05 = if (length(finite) > 0L) as.numeric(stats::quantile(finite, probs = 0.05, names = FALSE, type = 8)) else NA_real_,
    floor = as.numeric(floor_value),
    zero_count = as.integer(sum(is.finite(raw) & raw <= 0)),
    adjusted_count = adjusted
  )
}

# Build the public diagnostics row.  denom_* summarizes all supplied
# denominators, while value_denom_* and slope_denom_* retain role-specific
# summaries for downstream experiment CSVs.
.sbf_binning_denominator_diagnostics <- function(value_raw,
                                                 slope_raw = NULL,
                                                 floor_value = NA_real_,
                                                 value_adjusted_count = NULL,
                                                 slope_adjusted_count = NULL) {
  value_stats <- .sbf_binning_denominator_stats(
    raw = value_raw,
    floor_value = floor_value,
    adjusted_count = value_adjusted_count
  )
  has_slope <- !is.null(slope_raw)
  slope_stats <- if (has_slope) {
    .sbf_binning_denominator_stats(
      raw = slope_raw,
      floor_value = floor_value,
      adjusted_count = slope_adjusted_count
    )
  } else {
    list(
      min = NA_real_,
      p05 = NA_real_,
      floor = NA_real_,
      zero_count = NA_integer_,
      adjusted_count = NA_integer_
    )
  }

  combined_adjusted <- if (has_slope) {
    adjusted_counts <- c(value_stats$adjusted_count, slope_stats$adjusted_count)
    adjusted_counts <- adjusted_counts[is.finite(adjusted_counts)]
    if (length(adjusted_counts) > 0L) {
      as.integer(sum(adjusted_counts))
    } else {
      NA_integer_
    }
  } else {
    value_stats$adjusted_count
  }
  combined_raw <- if (has_slope) {
    c(
      as.numeric(unlist(value_raw, use.names = FALSE)),
      as.numeric(unlist(slope_raw, use.names = FALSE))
    )
  } else {
    value_raw
  }
  denom_stats <- .sbf_binning_denominator_stats(
    raw = combined_raw,
    floor_value = floor_value,
    adjusted_count = combined_adjusted
  )

  data.frame(
    denom_min = denom_stats$min,
    denom_p05 = denom_stats$p05,
    denom_floor = denom_stats$floor,
    denom_zero_count = denom_stats$zero_count,
    denom_adjusted_count = denom_stats$adjusted_count,
    value_denom_min = value_stats$min,
    value_denom_p05 = value_stats$p05,
    value_denom_floor = value_stats$floor,
    value_denom_zero_count = value_stats$zero_count,
    value_denom_adjusted_count = value_stats$adjusted_count,
    slope_denom_min = slope_stats$min,
    slope_denom_p05 = slope_stats$p05,
    slope_denom_floor = slope_stats$floor,
    slope_denom_zero_count = slope_stats$zero_count,
    slope_denom_adjusted_count = slope_stats$adjusted_count,
    stringsAsFactors = FALSE
  )
}

# Convert component-level binning prediction diagnostics into the table stored
# in hazard/survival prediction output.
.sbf_binning_prediction_diagnostics_table <- function(diagnostics,
                                                     result) {
  component_info <- result$component_info
  if (length(diagnostics) == 0L) {
    n_components <- nrow(component_info)
    return(data.frame(
      component = as.integer(component_info$index),
      role = as.character(component_info$role),
      denom_min = rep(NA_real_, n_components),
      denom_p05 = rep(NA_real_, n_components),
      denom_floor = rep(NA_real_, n_components),
      denom_zero_count = rep(0, n_components),
      denom_adjusted_count = rep(0, n_components),
      value_denom_min = rep(NA_real_, n_components),
      value_denom_p05 = rep(NA_real_, n_components),
      value_denom_floor = rep(NA_real_, n_components),
      value_denom_zero_count = rep(0, n_components),
      value_denom_adjusted_count = rep(0, n_components),
      slope_denom_min = rep(NA_real_, n_components),
      slope_denom_p05 = rep(NA_real_, n_components),
      slope_denom_floor = rep(NA_real_, n_components),
      slope_denom_zero_count = rep(0, n_components),
      slope_denom_adjusted_count = rep(0, n_components),
      kernel_support_zero_count = rep(0, n_components),
      bin_step_fallback_count = rep(0, n_components),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(seq_len(nrow(component_info)), function(i) {
    item <- diagnostics[[i]]
    data.frame(
      component = as.integer(component_info$index[i]),
      role = as.character(component_info$role[i]),
      denom_min = as.numeric(item$denom_min[1L]),
      denom_p05 = as.numeric(item$denom_p05[1L]),
      denom_floor = as.numeric(item$denom_floor[1L]),
      denom_zero_count = as.numeric(item$denom_zero_count[1L]),
      denom_adjusted_count = as.numeric(item$denom_adjusted_count[1L]),
      value_denom_min = as.numeric(item$value_denom_min[1L]),
      value_denom_p05 = as.numeric(item$value_denom_p05[1L]),
      value_denom_floor = as.numeric(item$value_denom_floor[1L]),
      value_denom_zero_count = as.numeric(item$value_denom_zero_count[1L]),
      value_denom_adjusted_count = as.numeric(item$value_denom_adjusted_count[1L]),
      slope_denom_min = as.numeric(item$slope_denom_min[1L]),
      slope_denom_p05 = as.numeric(item$slope_denom_p05[1L]),
      slope_denom_floor = as.numeric(item$slope_denom_floor[1L]),
      slope_denom_zero_count = as.numeric(item$slope_denom_zero_count[1L]),
      slope_denom_adjusted_count = as.numeric(item$slope_denom_adjusted_count[1L]),
      kernel_support_zero_count = as.numeric(item$kernel_support_zero_count[1L]),
      bin_step_fallback_count = as.numeric(item$bin_step_fallback_count[1L]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
