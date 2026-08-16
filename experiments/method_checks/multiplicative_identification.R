graphics.off()

source(file.path("experiments", "script_utils.R"))
project.root <- sbf_project_root()
sbf_load_package(project.root)

default_multiplicative_identification_method_check_config <- function(n_train = 600L,
                                                                       d = 5L,
                                                                       fit_it = 100L,
                                                                       target_k = 1L,
                                                                       model = 1L,
                                                                       violate_cox = TRUE,
                                                                       bandwidth = NULL,
                                                                       z_grid = seq(-1, 1, length.out = 300),
                                                                       t_grid = seq(0, 2, length.out = 300),
                                                                       rho = 0.5) {
  d <- as.integer(d)
  if (is.null(bandwidth)) {
    bandwidth <- c(0.2, rep(0.25, d - 1L))
  }
  list(
    n_train = as.integer(n_train),
    d = d,
    fit_it = as.integer(fit_it),
    target_k = as.integer(target_k),
    model = as.integer(model),
    rho = as.numeric(rho),
    violate_cox = isTRUE(violate_cox),
    bandwidth = sbf_expand_bandwidth(bandwidth, d),
    z_grid = as.numeric(z_grid),
    t_grid = as.numeric(t_grid)
  )
}

multiplicative_identification_truth_functions <- function(d, violate_cox = TRUE) {
  lapply(seq_len(as.integer(d) - 1L), function(k) {
    local({
      kk <- k
      function(z) sbf_true_multiplicative_phi(z, kk, violate_cox = violate_cox)
    })
  })
}

multiplicative_identification_predict <- function(fit, config) {
  z0 <- rep(0, config$d - 1L)
  component <- config$target_k + 1L
  effect <- sbf_predict(
    fit,
    type = "component",
    component = component,
    xout = config$z_grid,
    warn_diagnostics = FALSE
  )
  hazard <- sbf_predict(
    fit,
    type = "hazard",
    data = z0,
    times = config$t_grid,
    warn_diagnostics = FALSE
  )
  survival <- sbf_predict(
    fit,
    type = "survival",
    data = z0,
    times = config$t_grid,
    warn_diagnostics = FALSE
  )
  list(
    effect = log(pmax(1e-8, as.numeric(effect))),
    hazard = as.numeric(hazard$hazard[1L, ]),
    survival = as.numeric(survival$survival[1L, ])
  )
}

multiplicative_identification_curve_rows <- function(run, seed, identification, converged, config, truth, pred) {
  rbind(
    data.frame(
      run = run,
      seed = seed,
      identification = identification,
      converged = converged,
      curve_type = "effect",
      grid_value = config$z_grid,
      estimate = pred$effect,
      truth = truth$effect,
      stringsAsFactors = FALSE
    ),
    data.frame(
      run = run,
      seed = seed,
      identification = identification,
      converged = converged,
      curve_type = "survival",
      grid_value = config$t_grid,
      estimate = pred$survival,
      truth = truth$survival,
      stringsAsFactors = FALSE
    )
  )
}

multiplicative_identification_summary <- function(runs) {
  ids <- unique(runs$identification)
  do.call(rbind, lapply(ids, function(id) {
    all_rows <- runs[runs$identification == id, , drop = FALSE]
    rows <- all_rows[all_rows$converged %in% TRUE, , drop = FALSE]
    comparison_rows <- rows[rows$sample_mean_converged %in% TRUE, , drop = FALSE]
    data.frame(
      identification = id,
      n_runs = nrow(all_rows),
      n_converged = nrow(rows),
      mean_fit_runtime_sec = sbf_mean_or_na(rows$fit_runtime_sec),
      mean_iterations_used = sbf_mean_or_na(rows$iterations_used),
      mean_final_delta = sbf_mean_or_na(rows$final_delta),
      mean_effect_mse = sbf_mean_or_na(rows$effect_mse),
      mean_survival_mse = sbf_mean_or_na(rows$survival_mse),
      mean_hazard_max_abs_vs_sample_mean = sbf_mean_or_na(comparison_rows$hazard_max_abs_vs_sample_mean),
      mean_survival_max_abs_vs_sample_mean = sbf_mean_or_na(comparison_rows$survival_max_abs_vs_sample_mean),
      stringsAsFactors = FALSE
    )
  }))
}

multiplicative_identification_curve_mean <- function(curves, identification, curve_type, column) {
  rows <- curves[
    curves$identification == identification &
      curves$curve_type == curve_type,
    ,
    drop = FALSE
  ]
  grids <- sort(unique(rows$grid_value))
  values <- vapply(grids, function(g) {
    sbf_mean_or_na(rows[rows$grid_value == g, column])
  }, numeric(1))
  list(grid = grids, value = values)
}

multiplicative_identification_plot <- function(curves, out_dir) {
  file <- file.path(out_dir, "identification_curves.png")
  ids <- unique(curves$identification)
  curves <- curves[curves$converged %in% TRUE, , drop = FALSE]
  ids <- ids[ids %in% unique(curves$identification)]
  colors <- c(
    sample_mean = "#0072B2",
    integral = "#009E73",
    origin = "#CC79A7",
    jasa = "#E69F00"
  )
  display_labels <- c(
    sample_mean = "Sample mean",
    integral = "Integral",
    origin = "Single-point truth",
    jasa = "Two-point truth"
  )
  colors <- colors[ids]
  if (anyNA(colors)) {
    colors[is.na(colors)] <- grDevices::rainbow(sum(is.na(colors)))
  }

  grDevices::png(file, width = 1800, height = 760, res = 160)
  device_open <- TRUE
  old_par <- par(
    mfrow = c(1, 2),
    mar = c(4.3, 4.3, 3, 1.2),
    oma = c(0, 0, 2, 0),
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

  for (curve_type in c("effect", "survival")) {
    panel_title <- if (identical(curve_type, "effect")) "Component effect" else "Survival curve"
    if (length(ids) == 0L) {
      plot.new()
      box()
      title(main = panel_title)
      text(0.5, 0.5, "No converged fits")
      next
    }
    truth <- multiplicative_identification_curve_mean(curves, ids[[1L]], curve_type, "truth")
    ylim <- if (identical(curve_type, "survival")) c(0, 1) else range(truth$value, finite = TRUE)
    estimates <- lapply(ids, function(id) {
      multiplicative_identification_curve_mean(curves, id, curve_type, "estimate")
    })
    est_range <- range(unlist(lapply(estimates, `[[`, "value")), finite = TRUE)
    if (identical(curve_type, "effect")) {
      ylim <- range(c(ylim, est_range), finite = TRUE)
    }
    plot(
      truth$grid,
      truth$value,
      type = "l",
      col = "black",
      lwd = 2.5,
      ylim = ylim,
      xlab = if (identical(curve_type, "effect")) "Covariate value" else "Time",
      ylab = if (identical(curve_type, "effect")) "Log effect" else "Survival",
      main = panel_title
    )
    for (i in seq_along(ids)) {
      lines(estimates[[i]]$grid, estimates[[i]]$value, col = colors[[ids[[i]]]], lwd = 2, lty = i + 1L)
    }
    if (identical(curve_type, "survival")) {
      legend(
        "topright",
        legend = c("Truth", unname(display_labels[ids])),
        col = c("black", unname(colors)),
        lwd = c(2.5, rep(2, length(ids))),
        lty = c(1, seq_along(ids) + 1L),
        cex = 1.0,
        bty = "n"
      )
    }
  }
  mtext("Multiplicative Identification Rules", outer = TRUE, font = 2, cex = 1.2)
  par(old_par)
  grDevices::dev.off()
  device_open <- FALSE
  normalizePath(file, winslash = "/", mustWork = TRUE)
}

run_multiplicative_identification_method_check <- function(n_rep = 10L,
                                                           base_seed = 2026070201L,
                                                           identifications = c("sample_mean", "integral", "origin", "jasa"),
                                                           kernel = "epanechnikov",
                                                           config = NULL,
                                                           out_dir = NULL) {
  if (is.null(config)) {
    config <- default_multiplicative_identification_method_check_config()
  }
  identifications <- match.arg(identifications, several.ok = TRUE)
  if (!("sample_mean" %in% identifications)) {
    stop("identifications must include sample_mean for invariance comparisons.", call. = FALSE)
  }
  truth_functions <- multiplicative_identification_truth_functions(
    d = config$d,
    violate_cox = config$violate_cox
  )
  truth <- list(
    effect = sbf_true_multiplicative_phi(
      config$z_grid,
      config$target_k,
      violate_cox = config$violate_cox
    ),
    survival = sbf_true_multiplicative_survival(
      config$t_grid,
      rep(0, config$d - 1L),
      model = config$model,
      violate_cox = config$violate_cox
    )
  )
  if (is.null(out_dir)) {
    out_dir <- sbf_results_dir(project.root, "runs", "experiments", "method_checks", "multiplicative_identification")
  } else {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  }

  run_rows <- list()
  curve_rows <- list()
  for (run in seq_len(as.integer(n_rep))) {
    seed <- as.integer(base_seed) + run
    data <- sbf_simulate_data(
      n = config$n_train,
      d = config$d,
      rho = config$rho,
      family = "multiplicative",
      model = config$model,
      violate_cox = config$violate_cox,
      seed = seed
    )

    preds <- list()
    fits <- list()
    runtimes <- numeric(length(identifications))
    names(runtimes) <- identifications
    for (id in identifications) {
      args <- list(
        data = data,
        bandwidth = config$bandwidth,
        family = "multiplicative",
        iterations = config$fit_it,
        kernel = kernel,
        identification = id,
        warn_nonconvergence = FALSE
      )
      if (id %in% c("origin", "jasa")) {
        args$truth_functions <- truth_functions
      }
      start <- proc.time()[[3L]]
      fit <- do.call(sbf_fit, args)
      runtimes[[id]] <- as.numeric(proc.time()[[3L]] - start)
      fits[[id]] <- fit
      preds[[id]] <- multiplicative_identification_predict(fit, config)
      curve_rows[[length(curve_rows) + 1L]] <- multiplicative_identification_curve_rows(
        run = run,
        seed = seed,
        identification = id,
        converged = isTRUE(fit$converged),
        config = config,
        truth = truth,
        pred = preds[[id]]
      )
    }

    reference <- preds[["sample_mean"]]
    sample_mean_converged <- isTRUE(fits[["sample_mean"]]$converged)
    for (id in identifications) {
      run_rows[[length(run_rows) + 1L]] <- data.frame(
        run = run,
        seed = seed,
        identification = id,
        fit_runtime_sec = runtimes[[id]],
        converged = isTRUE(fits[[id]]$converged),
        sample_mean_converged = sample_mean_converged,
        iterations_used = as.numeric(fits[[id]]$iterations_used),
        final_delta = as.numeric(fits[[id]]$final_delta),
        effect_mse = sbf_mean_or_na((preds[[id]]$effect - truth$effect)^2),
        survival_mse = sbf_mean_or_na((preds[[id]]$survival - truth$survival)^2),
        hazard_max_abs_vs_sample_mean = max(abs(preds[[id]]$hazard - reference$hazard), na.rm = TRUE),
        survival_max_abs_vs_sample_mean = max(abs(preds[[id]]$survival - reference$survival), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  curves <- do.call(rbind, curve_rows)
  runs <- do.call(rbind, run_rows)
  summary <- multiplicative_identification_summary(runs)
  curves_file <- sbf_write_csv(curves, file.path(out_dir, "curves.csv"))
  runs_file <- sbf_write_csv(runs, file.path(out_dir, "runs.csv"))
  summary_file <- sbf_write_csv(summary, file.path(out_dir, "summary.csv"))
  plot_file <- multiplicative_identification_plot(curves, out_dir)

  list(
    curves = curves,
    runs = runs,
    summary = summary,
    curves_file = curves_file,
    runs_file = runs_file,
    summary_file = summary_file,
    plot_file = plot_file
  )
}
