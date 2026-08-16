# Shared binning helpers.
# They normalize binning requests and turn one numeric vector into grid values,
# bin indices, breaks, and support metadata.

# Return the canonical spec consumed by active fit precompute builders. 
.sbf_binning_prepare_spec <- function(time_bins = NULL,
                                      covariate_bins = NULL,
                                      time_binning_method = "quantile",
                                      covariate_binning_method = "quantile",
                                      representative = "midpoint") {
  list(
    time_method = time_binning_method,
    covariate_method = covariate_binning_method,
    time_bins = as.integer(time_bins),
    covariate_bins = as.integer(covariate_bins),
    representative = representative
  )
}

# Return the common binning result when requested bins reach the number of
# observed distinct values.
.sbf_binning_unique_grid_result <- function(x,
                                            bins,
                                            rule,
                                            representative) {
  rule <- as.character(rule)[1L]
  representative <- as.character(representative)[1L]

  reps <- sort(unique(x))
  if (length(reps) < 2L) {
    stop(
      "Binning exposure calculations require at least two distinct values for each component.",
      call. = FALSE
    )
  }
  bins_requested <- as.integer(bins[1L])
  mids <- 0.5 * (reps[-length(reps)] + reps[-1L])
  breaks <- c(min(x), mids, max(x))

  list(
    values = x,
    breaks = breaks,
    widths = diff(breaks),
    representatives = reps,
    effective_bins = length(reps),
    bin_index = match(x, reps),
    support_min = min(x),
    support_max = max(x),
    effective_support_min = min(reps),
    effective_support_max = max(reps),
    bins_requested = bins_requested,
    empty_bins = 0L,
    rule = rule,
    representative = representative,
    grid_mode = "unique"
  )
}

# Bin one finite numeric vector and return the common binning result.
.sbf_binning_bin_numeric_vector <- function(x,
                                           bins,
                                           rule,
                                           representative) {
  if (length(unique(x)) <= bins) {
    return(.sbf_binning_unique_grid_result(
      x = x,
      bins = bins,
      rule = rule,
      representative = representative
    ))
  }

  # equal_width_log: bins are equal width on the log1p scale (intended for the
  # time component), then mapped back to the original scale. log_u keeps the
  # log-scale breakpoints so a midpoint representative is the log-space midpoint.
  log_u <- NULL
  breaks <- if (identical(rule, "quantile")) {
    unique(as.numeric(stats::quantile(
      x,
      probs = seq(0, 1, length.out = bins + 1L),
      names = FALSE,
      type = 8
    )))
  } else if (identical(rule, "equal_width")) {
    unique(as.numeric(seq(min(x), max(x), length.out = bins + 1L)))
  } else if (identical(rule, "equal_width_log")) {
    if (any(x < 0)) {
      stop("equal_width_log binning requires nonnegative values.", call. = FALSE)
    }
    log_u <- seq(log1p(min(x)), log1p(max(x)), length.out = bins + 1L)
    edges <- as.numeric(expm1(log_u))
    edges[1L] <- min(x)
    edges[length(edges)] <- max(x)
    edges
  } else {
    stop("Internal error: unknown binning rule.", call. = FALSE)
  }
  if (length(unique(breaks)) < 2L) {
    return(.sbf_binning_unique_grid_result(
      x = x,
      bins = bins,
      rule = rule,
      representative = representative
    ))
  }

  bin_index <- cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  # Bin midpoints: log1p-scale for equal_width_log, original scale otherwise.
  # Shared by the midpoint representative and the empty-bin fallback of the mean
  # representative so a bin's center is defined consistently for each rule.
  midpoints <- if (!is.null(log_u)) {
    as.numeric(expm1(0.5 * (log_u[-1L] + log_u[-length(log_u)])))
  } else {
    0.5 * (breaks[-1L] + breaks[-length(breaks)])
  }

  reps <- if (identical(representative, "midpoint")) {
    midpoints
  } else if (identical(representative, "mean")) {
    vapply(
      seq_len(length(breaks) - 1L),
      function(k) {
        vals <- x[bin_index == k]
        if (length(vals) == 0L) {
          midpoints[k]
        } else {
          mean(vals)
        }
      },
      numeric(1)
    )
  } else {
    stop("Internal error: unknown bin representative.", call. = FALSE)
  }

  counts <- tabulate(bin_index, nbins = length(reps))
  used_bins <- which(counts > 0L)

  list(
    values = reps[bin_index],
    breaks = breaks,
    widths = diff(breaks),
    representatives = reps,
    effective_bins = length(reps),
    bin_index = bin_index,
    support_min = min(x),
    support_max = max(x),
    effective_support_min = min(reps[used_bins]),
    effective_support_max = max(reps[used_bins]),
    bins_requested = as.integer(bins),
    empty_bins = sum(counts == 0L),
    rule = rule,
    representative = representative,
    grid_mode = "binned"
  )
}
