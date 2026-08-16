benchmark_runs_from_pointwise <- function(results) {
  pointwise <- do.call(rbind, lapply(results, function(result) result$pointwise))
  run_meta <- do.call(rbind, lapply(results, function(result) result$run_meta))
  insert_after <- match("converged_sbf", names(run_meta))
  run_meta <- cbind(
    run_meta[, seq_len(insert_after), drop = FALSE],
    paired_converged = run_meta$converged_origin & run_meta$converged_sbf,
    run_meta[, seq.int(insert_after + 1L, ncol(run_meta)), drop = FALSE]
  )
  metrics <- lapply(seq_len(nrow(run_meta)), function(i) {
    row <- run_meta[i, , drop = FALSE]
    rows <- pointwise[
      pointwise$algorithm == row$algorithm &
        pointwise$run == row$run,
      ,
      drop = FALSE
    ]
    effect <- rows[rows$curve_type == "effect", , drop = FALSE]
    survival <- rows[rows$curve_type == "survival", , drop = FALSE]
    data.frame(
      effect_mse_origin = sbf_mean_or_na((effect$origin_minus_truth)^2),
      effect_mse_sbf = sbf_mean_or_na((effect$sbf_minus_truth)^2),
      effect_mse_sbf_vs_origin = sbf_mean_or_na((effect$sbf_minus_origin)^2),
      effect_mean_abs_sbf_minus_origin = sbf_mean_or_na(effect$abs_sbf_minus_origin),
      survival_mse_origin = sbf_mean_or_na((survival$origin_minus_truth)^2),
      survival_mse_sbf = sbf_mean_or_na((survival$sbf_minus_truth)^2),
      survival_mse_sbf_vs_origin = sbf_mean_or_na((survival$sbf_minus_origin)^2),
      survival_mean_abs_sbf_minus_origin = sbf_mean_or_na(survival$abs_sbf_minus_origin),
      stringsAsFactors = FALSE
    )
  })
  cbind(run_meta, do.call(rbind, metrics))
}

benchmark_summary_from_runs <- function(results, runs) {
  do.call(rbind, lapply(results, function(result) {
    all_rows <- runs[runs$algorithm == result$algorithm, , drop = FALSE]
    rows <- all_rows[all_rows$paired_converged %in% TRUE, , drop = FALSE]
    time_origin_mean <- sbf_mean_or_na(rows$time_origin)
    time_sbf_mean <- sbf_mean_or_na(rows$time_sbf)
    speedup <- if (is.finite(time_origin_mean) && is.finite(time_sbf_mean) && time_sbf_mean != 0) {
      time_origin_mean / time_sbf_mean
    } else {
      NA_real_
    }
    data.frame(
      algorithm = result$algorithm,
      model_family = result$model_family,
      fit_variant = result$fit_variant,
      n = result$n,
      d = result$d,
      n_rep = result$n_rep,
      n_origin_converged = sum(all_rows$converged_origin, na.rm = TRUE),
      n_sbf_converged = sum(all_rows$converged_sbf, na.rm = TRUE),
      n_paired_valid = nrow(rows),
      time_origin_mean = time_origin_mean,
      time_sbf_mean = time_sbf_mean,
      speedup_origin_over_sbf = speedup,
      effect_mse_origin_mean = sbf_mean_or_na(rows$effect_mse_origin),
      effect_mse_sbf_mean = sbf_mean_or_na(rows$effect_mse_sbf),
      effect_mse_sbf_vs_origin_mean = sbf_mean_or_na(rows$effect_mse_sbf_vs_origin),
      effect_mean_abs_sbf_minus_origin = sbf_mean_or_na(rows$effect_mean_abs_sbf_minus_origin),
      survival_mse_origin_mean = sbf_mean_or_na(rows$survival_mse_origin),
      survival_mse_sbf_mean = sbf_mean_or_na(rows$survival_mse_sbf),
      survival_mse_sbf_vs_origin_mean = sbf_mean_or_na(rows$survival_mse_sbf_vs_origin),
      survival_mean_abs_sbf_minus_origin = sbf_mean_or_na(rows$survival_mean_abs_sbf_minus_origin),
      stringsAsFactors = FALSE
    )
  }))
}

benchmark_curve_mean <- function(pointwise, algorithm, curve_type, column) {
  rows <- pointwise[
    pointwise$algorithm == algorithm &
      pointwise$curve_type == curve_type,
    ,
    drop = FALSE
  ]
  grids <- sort(unique(rows$grid_value))
  values <- vapply(grids, function(g) {
    sbf_mean_or_na(rows[rows$grid_value == g, column])
  }, numeric(1))
  list(grid = grids, value = values)
}

benchmark_curve_runs <- function(pointwise, algorithm, curve_type, column) {
  rows <- pointwise[
    pointwise$algorithm == algorithm &
      pointwise$curve_type == curve_type,
    ,
    drop = FALSE
  ]
  runs <- sort(unique(rows$run))
  curves <- vector("list", length(runs))
  for (i in seq_along(runs)) {
    run_rows <- rows[rows$run == runs[[i]], , drop = FALSE]
    run_rows <- run_rows[order(run_rows$grid_value), , drop = FALSE]
    curves[[i]] <- list(grid = run_rows$grid_value, value = run_rows[[column]])
  }
  curves
}

benchmark_plot_result <- function(result, pointwise, out_dir) {
  valid_pointwise <- pointwise[pointwise$paired_converged %in% TRUE, , drop = FALSE]
  file <- file.path(out_dir, result$plot_file)
  grDevices::png(file, width = 1900, height = 760, res = 160)
  device_open <- TRUE
  old_par <- par(
    mfrow = c(1, 2),
    mar = c(4.3, 4.3, 3, 1.2),
    oma = c(0, 0, 2.2, 0),
    cex.axis = 1.15,
    cex.lab = 1.25,
    cex.main = 1.2
  )
  on.exit({
    if (device_open) {
      par(old_par)
      grDevices::dev.off()
    }
  }, add = TRUE)

  line_cols <- c(truth = "black", origin = "#D55E00", sbf = "#0072B2")
  cloud_cols <- c(
    origin = grDevices::adjustcolor(line_cols["origin"], alpha.f = 0.13),
    sbf = grDevices::adjustcolor(line_cols["sbf"], alpha.f = 0.13)
  )
  line_lty <- c(truth = 1, origin = 2, sbf = 3)
  line_lwd <- 2.6

  effect_truth <- benchmark_curve_mean(pointwise, result$algorithm, "effect", "truth")
  effect_origin <- benchmark_curve_mean(valid_pointwise, result$algorithm, "effect", "origin")
  effect_sbf <- benchmark_curve_mean(valid_pointwise, result$algorithm, "effect", "sbf")
  effect_origin_runs <- benchmark_curve_runs(valid_pointwise, result$algorithm, "effect", "origin")
  effect_sbf_runs <- benchmark_curve_runs(valid_pointwise, result$algorithm, "effect", "sbf")
  plot(effect_truth$grid, effect_truth$value, type = "l", col = line_cols["truth"], lwd = line_lwd,
       ylim = c(-2, 2), xlab = sprintf("z%d", result$target_k), ylab = "Effect",
       main = "Component effect")
  for (curve in effect_origin_runs) {
    lines(curve$grid, curve$value, col = cloud_cols["origin"], lwd = 1)
  }
  for (curve in effect_sbf_runs) {
    lines(curve$grid, curve$value, col = cloud_cols["sbf"], lwd = 1)
  }
  lines(effect_truth$grid, effect_truth$value, col = line_cols["truth"], lwd = line_lwd)
  lines(effect_origin$grid, effect_origin$value, col = line_cols["origin"], lwd = line_lwd, lty = line_lty["origin"])
  lines(effect_sbf$grid, effect_sbf$value, col = line_cols["sbf"], lwd = line_lwd, lty = line_lty["sbf"])

  survival_truth <- benchmark_curve_mean(pointwise, result$algorithm, "survival", "truth")
  survival_origin <- benchmark_curve_mean(valid_pointwise, result$algorithm, "survival", "origin")
  survival_sbf <- benchmark_curve_mean(valid_pointwise, result$algorithm, "survival", "sbf")
  survival_origin_runs <- benchmark_curve_runs(valid_pointwise, result$algorithm, "survival", "origin")
  survival_sbf_runs <- benchmark_curve_runs(valid_pointwise, result$algorithm, "survival", "sbf")
  plot(survival_truth$grid, survival_truth$value, type = "l", col = line_cols["truth"], lwd = line_lwd,
       ylim = c(0, 1), xlab = "time", ylab = "Survival",
       main = "Survival curve")
  for (curve in survival_origin_runs) {
    lines(curve$grid, curve$value, col = cloud_cols["origin"], lwd = 1)
  }
  for (curve in survival_sbf_runs) {
    lines(curve$grid, curve$value, col = cloud_cols["sbf"], lwd = 1)
  }
  lines(survival_truth$grid, survival_truth$value, col = line_cols["truth"], lwd = line_lwd)
  lines(survival_origin$grid, survival_origin$value, col = line_cols["origin"], lwd = line_lwd, lty = line_lty["origin"])
  lines(survival_sbf$grid, survival_sbf$value, col = line_cols["sbf"], lwd = line_lwd, lty = line_lty["sbf"])
  legend("topright", legend = c("truth", "origin mean", "SBF mean"),
         col = line_cols,
         lwd = line_lwd,
         lty = line_lty,
         cex = 1.0,
         seg.len = 2.2, x.intersp = 0.8, y.intersp = 0.9,
         bty = "n")

  mtext(result$plot_title, outer = TRUE, font = 2, cex = 1.25)
  par(old_par)
  grDevices::dev.off()
  device_open <- FALSE
  normalizePath(file, winslash = "/", mustWork = TRUE)
}

benchmark_write_outputs <- function(results, out_dir) {
  pointwise <- do.call(rbind, lapply(results, function(result) result$pointwise))
  runs <- benchmark_runs_from_pointwise(results)
  run_key <- paste(runs$algorithm, runs$run)
  pointwise$paired_converged <- runs$paired_converged[
    match(paste(pointwise$algorithm, pointwise$run), run_key)
  ]
  summary <- benchmark_summary_from_runs(results, runs)

  pointwise_file <- file.path(out_dir, "benchmark_sbf_vs_origin_pointwise.csv")
  runs_file <- file.path(out_dir, "benchmark_sbf_vs_origin_runs.csv")
  summary_file <- file.path(out_dir, "benchmark_sbf_vs_origin_summary.csv")
  pointwise_file <- sbf_write_csv(pointwise, pointwise_file)
  runs_file <- sbf_write_csv(runs, runs_file)
  summary_file <- sbf_write_csv(summary, summary_file)

  plot_files <- vapply(results, function(result) {
    benchmark_plot_result(result, pointwise, out_dir)
  }, character(1))

  list(
    pointwise = pointwise,
    runs = runs,
    summary = summary,
    paths = list(
      pointwise_csv = pointwise_file,
      runs_csv = runs_file,
      summary_csv = summary_file,
      plots = plot_files
    )
  )
}
