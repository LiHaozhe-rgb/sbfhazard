graphics.off()

suppressPackageStartupMessages(library(survival))

source(file.path("experiments", "script_utils.R"))
project.root <- sbf_project_root()
sbf_load_package(project.root)

default_binning_layout_diagnostics_config <- function(n_train = 360L,
                                                      d = 4L,
                                                      fit_it = 50L,
                                                      model = 1L,
                                                      violate_cox = TRUE,
                                                      rho = 0.5,
                                                      additive_bandwidth = NULL,
                                                      multiplicative_bandwidth = NULL) {
  d <- as.integer(d)
  if (is.null(additive_bandwidth)) {
    additive_bandwidth <- c(0.4, rep(0.3, d - 1L))
  }
  if (is.null(multiplicative_bandwidth)) {
    multiplicative_bandwidth <- c(0.4, rep(0.2, d - 1L))
  }
  list(
    n_train = as.integer(n_train),
    d = d,
    fit_it = as.integer(fit_it),
    model = as.integer(model),
    violate_cox = isTRUE(violate_cox),
    rho = as.numeric(rho),
    additive_bandwidth = sbf_expand_bandwidth(additive_bandwidth, d),
    multiplicative_bandwidth = sbf_expand_bandwidth(multiplicative_bandwidth, d)
  )
}

default_binning_layout_diagnostics_specs <- function() {
  list(
    additive = list(
      list(
        label = "additive_equal_width",
        time_bins = 24L,
	        covariate_bins = 16L,
	        time_method = "equal_width",
	        covariate_method = "equal_width",
	        representative = "midpoint",
	        local_constant = TRUE
	      ),
      list(
        label = "additive_quantile",
        time_bins = 24L,
	        covariate_bins = 16L,
	        time_method = "quantile",
	        covariate_method = "quantile",
	        representative = "midpoint",
	        local_constant = TRUE
	      )
    ),
    multiplicative = list(
      list(
        label = "multiplicative_equal_width",
        time_bins = 24L,
        covariate_bins = 16L,
        time_method = "equal_width",
        covariate_method = "equal_width",
        representative = "midpoint"
      ),
      list(
        label = "multiplicative_quantile",
        time_bins = 24L,
        covariate_bins = 16L,
        time_method = "quantile",
        covariate_method = "quantile",
        representative = "midpoint"
      )
    )
  )
}

binning_layout_data_summary <- function(data, family, seed) {
  variables <- names(data)
  do.call(rbind, lapply(variables, function(variable) {
    values <- as.numeric(data[[variable]])
    finite <- values[is.finite(values)]
    role <- if (identical(variable, "time")) {
      "time"
    } else if (identical(variable, "status")) {
      "status"
    } else {
      "covariate"
    }
    data.frame(
      family = family,
      seed = seed,
      variable = variable,
      role = role,
      n = length(values),
      mean = mean(finite),
      sd = stats::sd(finite),
      min = min(finite),
      q05 = as.numeric(stats::quantile(finite, 0.05, names = FALSE, type = 8)),
      median = stats::median(finite),
      q95 = as.numeric(stats::quantile(finite, 0.95, names = FALSE, type = 8)),
      max = max(finite),
      stringsAsFactors = FALSE
    )
  }))
}

binning_layout_fit_row <- function(family, spec_label, fit, runtime_sec, error = NA_character_) {
  if (is.null(fit)) {
    return(data.frame(
      family = family,
      spec_label = spec_label,
      status = "error",
      runtime_sec = runtime_sec,
      converged = FALSE,
      iterations_used = NA_real_,
      final_delta = NA_real_,
      n_fit_rows = NA_real_,
      pairwise_total_cells = NA_real_,
      pairwise_nonzero_cells = NA_real_,
      pairwise_nonzero_share = NA_real_,
      pairwise_memory_mb = NA_real_,
      sparse_joint_cells = NA_real_,
      joint_exposure_nonzero_share = NA_real_,
      joint_exposure_memory_mb = NA_real_,
      error = error,
      stringsAsFactors = FALSE
    ))
  }

  binning <- fit$binning
  data.frame(
    family = family,
    spec_label = spec_label,
    status = "ok",
    runtime_sec = runtime_sec,
    converged = isTRUE(fit$converged),
    iterations_used = as.numeric(fit$iterations_used),
    final_delta = as.numeric(fit$final_delta),
    n_fit_rows = as.numeric(binning$n_fit_rows),
    pairwise_total_cells = if (!is.null(binning$pairwise_total_cells)) as.numeric(binning$pairwise_total_cells) else NA_real_,
    pairwise_nonzero_cells = if (!is.null(binning$pairwise_nonzero_cells)) as.numeric(binning$pairwise_nonzero_cells) else NA_real_,
    pairwise_nonzero_share = if (!is.null(binning$pairwise_nonzero_share)) as.numeric(binning$pairwise_nonzero_share) else NA_real_,
    pairwise_memory_mb = if (!is.null(binning$pairwise_memory_mb)) as.numeric(binning$pairwise_memory_mb) else NA_real_,
    sparse_joint_cells = if (!is.null(binning$sparse_joint_cells)) as.numeric(binning$sparse_joint_cells) else NA_real_,
    joint_exposure_nonzero_share = if (!is.null(binning$joint_exposure_nonzero_share)) as.numeric(binning$joint_exposure_nonzero_share) else NA_real_,
    joint_exposure_memory_mb = if (!is.null(binning$joint_exposure_memory_mb)) as.numeric(binning$joint_exposure_memory_mb) else NA_real_,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

binning_layout_fit_one <- function(data, family, bandwidth, iterations, spec) {
  start <- proc.time()[[3L]]
  fit <- tryCatch({
    if (identical(family, "additive")) {
      sbf_fit_binning(
        data = data,
        bandwidth = bandwidth,
        family = "additive",
        time_bins = spec$time_bins,
        covariate_bins = spec$covariate_bins,
	        time_binning_method = spec$time_method,
	        covariate_binning_method = spec$covariate_method,
	        representative = spec$representative,
	        local_constant = spec$local_constant,
        iterations = iterations,
        warn_nonconvergence = FALSE,
        warn_diagnostics = FALSE
      )
    } else {
      sbf_fit_binning(
        data = data,
        bandwidth = bandwidth,
        family = "multiplicative",
        time_bins = spec$time_bins,
        covariate_bins = spec$covariate_bins,
        time_binning_method = spec$time_method,
        covariate_binning_method = spec$covariate_method,
        representative = spec$representative,
        iterations = iterations,
        warn_nonconvergence = FALSE,
        warn_diagnostics = FALSE
      )
    }
  }, error = function(e) e)
  runtime <- as.numeric(proc.time()[[3L]] - start)

  if (inherits(fit, "error")) {
    return(list(
      fit = NULL,
      metadata = data.frame(),
      summary = binning_layout_fit_row(family, spec$label, NULL, runtime, conditionMessage(fit))
    ))
  }

  metadata <- fit$binning$metadata
  metadata$family <- family
  metadata$spec_label <- spec$label
  metadata <- metadata[, c("family", "spec_label", setdiff(names(metadata), c("family", "spec_label"))), drop = FALSE]
  list(
    fit = fit,
    metadata = metadata,
    summary = binning_layout_fit_row(family, spec$label, fit, runtime)
  )
}

binning_layout_plot_data <- function(data_sets, out_dir) {
  file <- file.path(out_dir, "binning_layout_data_distribution.png")
  grDevices::png(file, width = 1500, height = 1000, res = 150)
  old_par <- par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.5, 1))
  tryCatch({
    for (family in names(data_sets)) {
      data <- data_sets[[family]]
      features <- setdiff(names(data), c("time", "status"))
      hist(data$time, breaks = 24L, col = "#9ECAE1", border = "white",
           xlab = "time", main = sprintf("%s time", family))
      boxplot(data[, features, drop = FALSE], col = "#C7E9C0", outline = FALSE,
              ylab = "value", main = sprintf("%s covariates", family))
    }
  }, finally = {
    par(old_par)
    grDevices::dev.off()
  })
  normalizePath(file, winslash = "/", mustWork = TRUE)
}

binning_layout_plot_metadata <- function(metadata, out_dir) {
  file <- file.path(out_dir, "binning_layout_metadata.png")
  grDevices::png(file, width = 1600, height = 900, res = 150)
  old_par <- par(mfrow = c(1, 2), mar = c(7, 4.2, 3, 1))
  tryCatch({
    labels <- paste(metadata$family, metadata$spec_label, metadata$dimension, sep = "\n")
    barplot(metadata$effective_bins, names.arg = labels, las = 2, cex.names = 0.65,
            col = "#6BAED6", ylab = "stored bins", main = "Stored bins")
    barplot(metadata$empty_bins, names.arg = labels, las = 2, cex.names = 0.65,
            col = "#FDAE6B", ylab = "empty stored bins", main = "Empty stored bins")
  }, finally = {
    par(old_par)
    grDevices::dev.off()
  })
  normalizePath(file, winslash = "/", mustWork = TRUE)
}

binning_layout_plot_fit_summary <- function(summary, out_dir) {
  file <- file.path(out_dir, "binning_layout_cell_summary.png")
  grDevices::png(file, width = 1400, height = 800, res = 150)
  old_par <- par(mfrow = c(1, 2), mar = c(6, 4.2, 3, 1))
  tryCatch({
    labels <- paste(summary$family, summary$spec_label, sep = "\n")
    cell_count <- ifelse(
      is.finite(summary$pairwise_nonzero_cells),
      summary$pairwise_nonzero_cells,
      summary$sparse_joint_cells
    )
    nonzero_share <- ifelse(
      is.finite(summary$pairwise_nonzero_share),
      summary$pairwise_nonzero_share,
      summary$joint_exposure_nonzero_share
    )
    barplot(pmax(cell_count, 1), names.arg = labels, las = 2, cex.names = 0.75,
            log = "y", col = "#74C476", ylab = "cells, log scale", main = "Nonzero/observed cells")
    barplot(nonzero_share, names.arg = labels, las = 2, cex.names = 0.75,
            ylim = c(0, 1), col = "#9E9AC8", ylab = "share", main = "Nonzero share")
  }, finally = {
    par(old_par)
    grDevices::dev.off()
  })
  normalizePath(file, winslash = "/", mustWork = TRUE)
}

binning_layout_write_outputs <- function(data_summary, metadata, fit_summary, data_sets, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data_summary_file <- file.path(out_dir, "binning_layout_data_summary.csv")
  metadata_file <- file.path(out_dir, "binning_layout_metadata.csv")
  fit_summary_file <- file.path(out_dir, "binning_layout_fit_summary.csv")
  sbf_write_csv(data_summary, data_summary_file)
  sbf_write_csv(metadata, metadata_file)
  sbf_write_csv(fit_summary, fit_summary_file)
  plot_files <- c(
    binning_layout_plot_data(data_sets, out_dir),
    binning_layout_plot_metadata(metadata, out_dir),
    binning_layout_plot_fit_summary(fit_summary, out_dir)
  )
  list(
    data_summary = data_summary,
    metadata = metadata,
    fit_summary = fit_summary,
    paths = list(
      data_summary_csv = normalizePath(data_summary_file, winslash = "/", mustWork = TRUE),
      metadata_csv = normalizePath(metadata_file, winslash = "/", mustWork = TRUE),
      fit_summary_csv = normalizePath(fit_summary_file, winslash = "/", mustWork = TRUE),
      plots = plot_files
    )
  )
}

run_binning_layout_diagnostics <- function(base_seed = 2026060101L,
                                           config = default_binning_layout_diagnostics_config(),
                                           specs = default_binning_layout_diagnostics_specs(),
                                           out_dir = NULL) {
  if (is.null(out_dir)) {
    out_dir <- sbf_results_dir(
      project.root,
      "runs",
      "experiments",
      "diagnostics",
      "binning_layout_diagnostics"
    )
  }
  data_sets <- list(
    additive = sbf_simulate_data(
      n = config$n_train,
      d = config$d,
      rho = config$rho,
      family = "additive",
      model = config$model,
      violate_cox = config$violate_cox,
      seed = as.integer(base_seed)
    ),
    multiplicative = sbf_simulate_data(
      n = config$n_train,
      d = config$d,
      rho = config$rho,
      family = "multiplicative",
      model = config$model,
      violate_cox = config$violate_cox,
      seed = as.integer(base_seed) + 1L
    )
  )

  data_summary <- do.call(rbind, list(
    binning_layout_data_summary(data_sets$additive, "additive", as.integer(base_seed)),
    binning_layout_data_summary(data_sets$multiplicative, "multiplicative", as.integer(base_seed) + 1L)
  ))

  metadata_rows <- list()
  fit_rows <- list()
  idx <- 1L
  for (family in names(data_sets)) {
    bandwidth <- if (identical(family, "additive")) {
      config$additive_bandwidth
    } else {
      config$multiplicative_bandwidth
    }
    for (spec in specs[[family]]) {
      result <- binning_layout_fit_one(
        data = data_sets[[family]],
        family = family,
        bandwidth = bandwidth,
        iterations = config$fit_it,
        spec = spec
      )
      if (nrow(result$metadata) > 0L) {
        metadata_rows[[idx]] <- result$metadata
      }
      fit_rows[[idx]] <- result$summary
      idx <- idx + 1L
    }
  }

  binning_layout_write_outputs(
    data_summary = data_summary,
    metadata = do.call(rbind, metadata_rows),
    fit_summary = do.call(rbind, fit_rows),
    data_sets = data_sets,
    out_dir = out_dir
  )
}

if (sys.nframe() == 0L) {
  outputs <- run_binning_layout_diagnostics()
  cat("Binning layout diagnostics complete.\n")
  cat(sprintf("Data summary: %s\n", outputs$paths$data_summary_csv))
  cat(sprintf("Metadata: %s\n", outputs$paths$metadata_csv))
  cat(sprintf("Fit summary: %s\n", outputs$paths$fit_summary_csv))
}
