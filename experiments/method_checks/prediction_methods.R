graphics.off()

source(file.path("experiments", "script_utils.R"))
project.root <- sbf_project_root()
sbf_load_package(project.root)

default_prediction_methods_method_check_config <- function(n_train = 600L,
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

prediction_method_truth_functions <- function(d, violate_cox = TRUE) {
  lapply(seq_len(as.integer(d) - 1L), function(k) {
    local({
      kk <- k
      function(z) sbf_true_multiplicative_phi(z, kk, violate_cox = violate_cox)
    })
  })
}

prediction_method_specs <- function() {
  data.frame(
    algorithm = c("additive_local_linear", "additive_local_constant", "multiplicative_local_constant"),
    family = c("additive", "additive", "multiplicative"),
    local_constant = c(FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
}

prediction_method_truth <- function(spec, config) {
  z0 <- rep(0, config$d - 1L)
  if (identical(spec$family, "additive")) {
    return(list(
      effect = sbf_true_additive_phi(config$z_grid, config$target_k, config$d - 1L, config$violate_cox),
      survival = sbf_true_additive_survival(config$t_grid, z0, config$model, config$violate_cox)
    ))
  }
  list(
    effect = sbf_true_multiplicative_phi(config$z_grid, config$target_k, config$violate_cox),
    survival = sbf_true_multiplicative_survival(config$t_grid, z0, config$model, config$violate_cox)
  )
}

prediction_method_fit <- function(data, spec, config, kernel, truth_functions) {
  args <- list(
    data = data,
    bandwidth = config$bandwidth,
    family = spec$family,
    iterations = config$fit_it,
    kernel = kernel,
    warn_nonconvergence = FALSE
  )
  if (identical(spec$family, "additive")) {
    args$local_constant <- isTRUE(spec$local_constant)
  } else {
    args$identification <- "origin"
    args$truth_functions <- truth_functions
  }
  start <- proc.time()[[3L]]
  fit <- do.call(sbf_fit, args)
  list(
    fit = fit,
    runtime_sec = as.numeric(proc.time()[[3L]] - start)
  )
}

prediction_method_predict <- function(fit, spec, config, prediction) {
  z0 <- rep(0, config$d - 1L)
  component <- config$target_k + 1L
  effect <- sbf_predict(
    fit,
    type = "component",
    component = component,
    xout = config$z_grid,
    prediction = prediction,
    warn_diagnostics = FALSE
  )
  survival <- sbf_predict(
    fit,
    type = "survival",
    data = z0,
    times = config$t_grid,
    prediction = prediction,
    warn_diagnostics = FALSE
  )
  effect <- as.numeric(effect)
  if (identical(spec$family, "multiplicative")) {
    effect <- log(pmax(1e-8, effect))
  }
  list(
    effect = effect,
    survival = as.numeric(survival$survival[1L, ])
  )
}

prediction_method_curve_rows <- function(algorithm, run, seed, prediction, converged, config, truth, pred) {
  rbind(
    data.frame(
      algorithm = algorithm,
      run = run,
      seed = seed,
      prediction = prediction,
      converged = converged,
      curve_type = "effect",
      grid_value = config$z_grid,
      estimate = pred$effect,
      truth = truth$effect,
      stringsAsFactors = FALSE
    ),
    data.frame(
      algorithm = algorithm,
      run = run,
      seed = seed,
      prediction = prediction,
      converged = converged,
      curve_type = "survival",
      grid_value = config$t_grid,
      estimate = pred$survival,
      truth = truth$survival,
      stringsAsFactors = FALSE
    )
  )
}

prediction_methods_summary <- function(runs) {
  groups <- unique(runs[, c("algorithm", "prediction"), drop = FALSE])
  do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    group <- groups[i, , drop = FALSE]
    all_rows <- runs[
      runs$algorithm == group$algorithm &
        runs$prediction == group$prediction,
      ,
      drop = FALSE
    ]
    rows <- all_rows[all_rows$converged %in% TRUE, , drop = FALSE]
    data.frame(
      algorithm = group$algorithm,
      prediction = group$prediction,
      n_runs = nrow(all_rows),
      n_converged = nrow(rows),
      mean_fit_runtime_sec = sbf_mean_or_na(rows$fit_runtime_sec),
      mean_iterations_used = sbf_mean_or_na(rows$iterations_used),
      mean_final_delta = sbf_mean_or_na(rows$final_delta),
      mean_effect_mse = sbf_mean_or_na(rows$effect_mse),
      mean_survival_mse = sbf_mean_or_na(rows$survival_mse),
      mean_effect_max_abs_formula_minus_interpolation = sbf_mean_or_na(rows$effect_max_abs_formula_minus_interpolation),
      mean_survival_max_abs_formula_minus_interpolation = sbf_mean_or_na(rows$survival_max_abs_formula_minus_interpolation),
      stringsAsFactors = FALSE
    )
  }))
}

prediction_method_curve_mean <- function(curves, algorithm, prediction, curve_type, column) {
  rows <- curves[
    curves$algorithm == algorithm &
      curves$prediction == prediction &
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

prediction_methods_plot <- function(curves, out_dir) {
  algorithm_labels <- c(
    additive_local_constant = "Additive local constant",
    additive_local_linear = "Additive local linear",
    multiplicative_local_constant = "Multiplicative local constant"
  )
  plot_names <- c(
    additive_local_constant = "prediction_additive_local_constant.png",
    additive_local_linear = "prediction_additive_local_linear.png",
    multiplicative_local_constant = "prediction_multiplicative_local_constant.png"
  )
  algorithms <- intersect(names(algorithm_labels), unique(curves$algorithm))
  curves <- curves[curves$converged %in% TRUE, , drop = FALSE]
  predictions <- c("formula", "interpolation")
  colors <- c(formula = "#0072B2", interpolation = "#D55E00")
  line_types <- c(truth = 1, formula = 2, interpolation = 3)
  line_width <- 2.6
  plot_files <- character(length(algorithms))

  for (a in seq_along(algorithms)) {
    algorithm <- algorithms[[a]]
    file <- file.path(out_dir, plot_names[[algorithm]])
    grDevices::png(file, width = 1900, height = 760, res = 160)
    old_par <- par(
      mfrow = c(1, 2),
      mar = c(4.3, 4.3, 3, 1.2),
      oma = c(0, 0, 2.2, 0),
      cex.axis = 1.15,
      cex.lab = 1.25,
      cex.main = 1.2
    )

    for (curve_type in c("effect", "survival")) {
      panel_title <- if (identical(curve_type, "effect")) "Component effect" else "Survival curve"
      algorithm_rows <- curves[curves$algorithm == algorithm, , drop = FALSE]
      if (nrow(algorithm_rows) == 0L) {
        plot.new()
        box()
        title(main = panel_title)
        text(0.5, 0.5, "No converged fits")
        next
      }
      truth <- prediction_method_curve_mean(curves, algorithm, predictions[[1L]], curve_type, "truth")
      estimates <- lapply(predictions, function(prediction) {
        prediction_method_curve_mean(curves, algorithm, prediction, curve_type, "estimate")
      })
      ylim <- if (identical(curve_type, "survival")) c(0, 1) else {
        range(c(truth$value, unlist(lapply(estimates, `[[`, "value"))), finite = TRUE)
      }
      plot(
        truth$grid,
        truth$value,
        type = "l",
        col = "black",
        lwd = line_width,
        lty = line_types[["truth"]],
        ylim = ylim,
        xlab = if (identical(curve_type, "effect")) "z1" else "time",
        ylab = if (identical(curve_type, "effect")) "Effect" else "Survival",
        main = panel_title
      )
      for (i in seq_along(predictions)) {
        prediction <- predictions[[i]]
        lines(
          estimates[[i]]$grid,
          estimates[[i]]$value,
          col = colors[[prediction]],
          lwd = line_width,
          lty = line_types[[prediction]]
        )
      }
      if (identical(curve_type, "survival")) {
        legend(
          "topright",
          legend = c("truth", predictions),
          col = c("black", unname(colors)),
          lwd = line_width,
          lty = unname(line_types),
          cex = 1.0,
          seg.len = 2.2,
          x.intersp = 0.8,
          y.intersp = 0.9,
          bty = "n"
        )
      }
    }
    mtext(algorithm_labels[[algorithm]], outer = TRUE, font = 2, cex = 1.25)
    par(old_par)
    grDevices::dev.off()
    plot_files[[a]] <- normalizePath(file, winslash = "/", mustWork = TRUE)
  }

  plot_files
}

run_prediction_methods_method_check <- function(n_rep = 10L,
                                                base_seed = 2026070301L,
                                                kernel = "epanechnikov",
                                                config = NULL,
                                                out_dir = NULL) {
  if (is.null(config)) {
    config <- default_prediction_methods_method_check_config()
  }
  if (is.null(out_dir)) {
    out_dir <- sbf_results_dir(project.root, "runs", "experiments", "method_checks", "prediction_methods")
  } else {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  }
  specs <- prediction_method_specs()
  truth_functions <- prediction_method_truth_functions(config$d, config$violate_cox)
  predictions <- c("formula", "interpolation")

  run_rows <- list()
  curve_rows <- list()
  for (run in seq_len(as.integer(n_rep))) {
    seed <- as.integer(base_seed) + run
    for (s in seq_len(nrow(specs))) {
      spec <- specs[s, , drop = FALSE]
      data <- sbf_simulate_data(
        n = config$n_train,
        d = config$d,
        rho = config$rho,
        family = spec$family,
        model = config$model,
        violate_cox = config$violate_cox,
        seed = seed
      )
      fit_obj <- prediction_method_fit(
        data = data,
        spec = spec,
        config = config,
        kernel = kernel,
        truth_functions = truth_functions
      )
      truth <- prediction_method_truth(spec, config)
      pred_values <- list()
      for (prediction in predictions) {
        pred_values[[prediction]] <- prediction_method_predict(
          fit = fit_obj$fit,
          spec = spec,
          config = config,
          prediction = prediction
        )
        curve_rows[[length(curve_rows) + 1L]] <- prediction_method_curve_rows(
          algorithm = spec$algorithm,
          run = run,
          seed = seed,
          prediction = prediction,
          converged = isTRUE(fit_obj$fit$converged),
          config = config,
          truth = truth,
          pred = pred_values[[prediction]]
        )
      }

      effect_diff <- max(abs(pred_values$formula$effect - pred_values$interpolation$effect), na.rm = TRUE)
      survival_diff <- max(abs(pred_values$formula$survival - pred_values$interpolation$survival), na.rm = TRUE)
      for (prediction in predictions) {
        run_rows[[length(run_rows) + 1L]] <- data.frame(
          algorithm = spec$algorithm,
          run = run,
          seed = seed,
          prediction = prediction,
          fit_runtime_sec = fit_obj$runtime_sec,
          converged = isTRUE(fit_obj$fit$converged),
          iterations_used = as.numeric(fit_obj$fit$iterations_used),
          final_delta = as.numeric(fit_obj$fit$final_delta),
          effect_mse = sbf_mean_or_na((pred_values[[prediction]]$effect - truth$effect)^2),
          survival_mse = sbf_mean_or_na((pred_values[[prediction]]$survival - truth$survival)^2),
          effect_max_abs_formula_minus_interpolation = effect_diff,
          survival_max_abs_formula_minus_interpolation = survival_diff,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  curves <- do.call(rbind, curve_rows)
  runs <- do.call(rbind, run_rows)
  summary <- prediction_methods_summary(runs)
  curves_file <- sbf_write_csv(curves, file.path(out_dir, "curves.csv"))
  runs_file <- sbf_write_csv(runs, file.path(out_dir, "runs.csv"))
  summary_file <- sbf_write_csv(summary, file.path(out_dir, "summary.csv"))
  plot_files <- prediction_methods_plot(curves, out_dir)

  list(
    curves = curves,
    runs = runs,
    summary = summary,
    curves_file = curves_file,
    runs_file = runs_file,
    summary_file = summary_file,
    plot_files = plot_files
  )
}
