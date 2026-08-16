benchmark_api_curves <- function(fit, family, component, z_grid, z0, t_grid) {
  effect <- sbf_predict(
    fit,
    type = "component",
    component = component,
    xout = z_grid,
    warn_diagnostics = FALSE
  )
  effect <- as.numeric(effect)
  if (identical(family, "multiplicative")) {
    effect <- log(effect)
  }

  survival <- sbf_predict(
    fit,
    type = "survival",
    data = z0,
    times = t_grid,
    warn_diagnostics = FALSE
  )
  list(effect = effect, survival = as.numeric(survival$survival[1L, ]))
}

benchmark_make_pointwise <- function(algorithm,
                                     model_family,
                                     fit_variant,
                                     run,
                                     seed,
                                     effect_grid,
                                     survival_grid,
                                     truth_effect,
                                     truth_survival,
                                     origin,
                                     sbf) {
  make_curve <- function(curve_type, grid, truth, origin_curve, sbf_curve) {
    origin_minus_truth <- as.numeric(origin_curve) - as.numeric(truth)
    sbf_minus_truth <- as.numeric(sbf_curve) - as.numeric(truth)
    sbf_minus_origin <- as.numeric(sbf_curve) - as.numeric(origin_curve)
    data.frame(
      algorithm = algorithm,
      model_family = model_family,
      fit_variant = fit_variant,
      run = run,
      seed = seed,
      curve_type = curve_type,
      grid_value = as.numeric(grid),
      truth = as.numeric(truth),
      origin = as.numeric(origin_curve),
      sbf = as.numeric(sbf_curve),
      origin_minus_truth = origin_minus_truth,
      sbf_minus_truth = sbf_minus_truth,
      sbf_minus_origin = sbf_minus_origin,
      abs_sbf_minus_origin = abs(sbf_minus_origin),
      stringsAsFactors = FALSE
    )
  }
  rbind(
    make_curve("effect", effect_grid, truth_effect, origin$effect, sbf$effect),
    make_curve("survival", survival_grid, truth_survival, origin$survival, sbf$survival)
  )
}

benchmark_make_run_meta <- function(algorithm,
                                    model_family,
                                    fit_variant,
                                    run,
                                    seed,
                                    origin,
                                    sbf) {
  data.frame(
    algorithm = algorithm,
    model_family = model_family,
    fit_variant = fit_variant,
    run = run,
    seed = seed,
    time_origin = origin$time,
    time_sbf = sbf$time,
    converged_origin = origin$converged,
    converged_sbf = sbf$converged,
    iter_origin = origin$iterations_used,
    iter_sbf = sbf$iterations_used,
    final_delta_origin = origin$final_delta,
    final_delta_sbf = sbf$final_delta,
    stringsAsFactors = FALSE
  )
}

run_additive_local_linear <- function(add_origin_fun,
                                      n,
                                      d,
                                      rho,
                                      model,
                                      violate_cox,
                                      iterations,
                                      n_rep,
                                      seed_base,
                                      target_k,
                                      bandwidth,
                                      z_grid,
                                      t_grid) {
  z0 <- rep(0, d - 1L)
  component <- target_k + 1L
  truth_effect <- sbf_true_additive_phi(z_grid, target_k, d - 1L, violate_cox)
  truth_survival <- sbf_true_additive_survival(t_grid, z0, model, violate_cox)

  pointwise <- vector("list", n_rep)
  run_meta <- vector("list", n_rep)
  seeds <- seed_base + seq_len(n_rep) - 1L

  for (i in seq_len(n_rep)) {
    data <- sbf_simulate_data(
      n = n,
      d = d,
      rho = rho,
      family = "additive",
      model = model,
      violate_cox = violate_cox,
      seed = seeds[[i]]
    )

    start <- proc.time()[[3L]]
    origin_fit <- sbf_fit_origin(
      data = data,
      bandwidth = bandwidth,
      family = "additive",
      origin_engine = add_origin_fun,
      integral_approx = "midd",
      iterations = iterations,
      kernel = "epanechnikov",
      local_constant = FALSE
    )
    origin_time <- as.numeric(proc.time()[[3L]] - start)
    origin_curves <- benchmark_api_curves(origin_fit, "additive", component, z_grid, z0, t_grid)
    origin <- list(
      time = origin_time,
      iterations_used = origin_fit$iterations_used,
      final_delta = origin_fit$final_delta,
      converged = origin_fit$converged,
      effect = origin_curves$effect,
      survival = origin_curves$survival
    )

    start <- proc.time()[[3L]]
    sbf_fit_obj <- sbf_fit(
      data = data,
      bandwidth = bandwidth,
      family = "additive",
      iterations = iterations,
      integral_approx = "midd",
      kernel = "epanechnikov",
      local_constant = FALSE,
      warn_nonconvergence = FALSE
    )
    sbf_time <- as.numeric(proc.time()[[3L]] - start)
    sbf_curves <- benchmark_api_curves(sbf_fit_obj, "additive", component, z_grid, z0, t_grid)
    sbf <- list(
      time = sbf_time,
      iterations_used = sbf_fit_obj$iterations_used,
      final_delta = sbf_fit_obj$final_delta,
      converged = sbf_fit_obj$converged,
      effect = sbf_curves$effect,
      survival = sbf_curves$survival
    )

    pointwise[[i]] <- benchmark_make_pointwise(
      "additive_local_linear", "additive", "local_linear", i, seeds[[i]],
      z_grid, t_grid, truth_effect, truth_survival, origin, sbf
    )
    run_meta[[i]] <- benchmark_make_run_meta(
      "additive_local_linear", "additive", "local_linear", i, seeds[[i]], origin, sbf
    )
  }

  list(
    algorithm = "additive_local_linear",
    model_family = "additive",
    fit_variant = "local_linear",
    n = n,
    d = d,
    n_rep = n_rep,
    target_k = target_k,
    plot_title = "Additive local linear",
    plot_file = "benchmark_sbf_vs_origin_additive_local_linear.png",
    pointwise = do.call(rbind, pointwise),
    run_meta = do.call(rbind, run_meta)
  )
}

run_additive_local_constant <- function(add_origin_fun,
                                        n,
                                        d,
                                        rho,
                                        model,
                                        violate_cox,
                                        iterations,
                                        n_rep,
                                        seed_base,
                                        target_k,
                                        bandwidth,
                                        z_grid,
                                        t_grid) {
  z0 <- rep(0, d - 1L)
  component <- target_k + 1L
  truth_effect <- sbf_true_additive_phi(z_grid, target_k, d - 1L, violate_cox)
  truth_survival <- sbf_true_additive_survival(t_grid, z0, model, violate_cox)

  pointwise <- vector("list", n_rep)
  run_meta <- vector("list", n_rep)
  seeds <- seed_base + seq_len(n_rep) - 1L

  for (i in seq_len(n_rep)) {
    data <- sbf_simulate_data(
      n = n,
      d = d,
      rho = rho,
      family = "additive",
      model = model,
      violate_cox = violate_cox,
      seed = seeds[[i]]
    )

    start <- proc.time()[[3L]]
    origin_fit <- sbf_fit_origin(
      data = data,
      bandwidth = bandwidth,
      family = "additive",
      origin_engine = add_origin_fun,
      integral_approx = "midd",
      iterations = iterations,
      kernel = "epanechnikov",
      local_constant = TRUE
    )
    origin_time <- as.numeric(proc.time()[[3L]] - start)
    origin_curves <- benchmark_api_curves(origin_fit, "additive", component, z_grid, z0, t_grid)
    origin <- list(
      time = origin_time,
      iterations_used = origin_fit$iterations_used,
      final_delta = origin_fit$final_delta,
      converged = origin_fit$converged,
      effect = origin_curves$effect,
      survival = origin_curves$survival
    )

    start <- proc.time()[[3L]]
    sbf_fit_obj <- sbf_fit(
      data = data,
      bandwidth = bandwidth,
      family = "additive",
      iterations = iterations,
      integral_approx = "midd",
      kernel = "epanechnikov",
      local_constant = TRUE,
      warn_nonconvergence = FALSE
    )
    sbf_time <- as.numeric(proc.time()[[3L]] - start)
    sbf_curves <- benchmark_api_curves(sbf_fit_obj, "additive", component, z_grid, z0, t_grid)
    sbf <- list(
      time = sbf_time,
      iterations_used = sbf_fit_obj$iterations_used,
      final_delta = sbf_fit_obj$final_delta,
      converged = sbf_fit_obj$converged,
      effect = sbf_curves$effect,
      survival = sbf_curves$survival
    )

    pointwise[[i]] <- benchmark_make_pointwise(
      "additive_local_constant", "additive", "local_constant", i, seeds[[i]],
      z_grid, t_grid, truth_effect, truth_survival, origin, sbf
    )
    run_meta[[i]] <- benchmark_make_run_meta(
      "additive_local_constant", "additive", "local_constant", i, seeds[[i]], origin, sbf
    )
  }

  list(
    algorithm = "additive_local_constant",
    model_family = "additive",
    fit_variant = "local_constant",
    n = n,
    d = d,
    n_rep = n_rep,
    target_k = target_k,
    plot_title = "Additive local constant",
    plot_file = "benchmark_sbf_vs_origin_additive_local_constant.png",
    pointwise = do.call(rbind, pointwise),
    run_meta = do.call(rbind, run_meta)
  )
}

run_multiplicative_local_constant <- function(mul_origin_fun,
                                              n,
                                              d,
                                              rho,
                                              model,
                                              iterations,
                                              n_rep,
                                              seed_base,
                                              target_k,
                                              bandwidth,
                                              z_grid,
                                              t_grid) {
  z0 <- rep(0, d - 1L)
  component <- target_k + 1L
  truth_effect <- sbf_true_multiplicative_phi(z_grid, target_k, violate_cox = TRUE)
  truth_survival <- sbf_true_multiplicative_survival(t_grid, z0, model, violate_cox = TRUE)
  truth_functions <- lapply(seq_len(d - 1L), function(k) {
    local({
      kk <- k
      function(z) sbf_true_multiplicative_phi(z, kk, violate_cox = TRUE)
    })
  })

  pointwise <- vector("list", n_rep)
  run_meta <- vector("list", n_rep)
  seeds <- seed_base + seq_len(n_rep) - 1L

  for (i in seq_len(n_rep)) {
    data <- sbf_simulate_data(
      n = n,
      d = d,
      rho = rho,
      family = "multiplicative",
      model = model,
      violate_cox = TRUE,
      seed = seeds[[i]]
    )

    start <- proc.time()[[3L]]
    origin_fit <- sbf_fit_origin(
      data = data,
      bandwidth = bandwidth,
      family = "multiplicative",
      origin_engine = mul_origin_fun,
      integral_approx = "midd",
      iterations = iterations,
      kernel = "epanechnikov",
      truth_functions = truth_functions
    )
    origin_time <- as.numeric(proc.time()[[3L]] - start)
    origin_curves <- benchmark_api_curves(origin_fit, "multiplicative", component, z_grid, z0, t_grid)
    origin <- list(
      time = origin_time,
      iterations_used = origin_fit$iterations_used,
      final_delta = origin_fit$final_delta,
      converged = origin_fit$converged,
      effect = origin_curves$effect,
      survival = origin_curves$survival
    )

    start <- proc.time()[[3L]]
    sbf_fit_obj <- sbf_fit(
      data = data,
      bandwidth = bandwidth,
      family = "multiplicative",
      iterations = iterations,
      integral_approx = "midd",
      kernel = "epanechnikov",
      identification = "origin",
      truth_functions = truth_functions,
      warn_nonconvergence = FALSE
    )
    sbf_time <- as.numeric(proc.time()[[3L]] - start)
    sbf_curves <- benchmark_api_curves(sbf_fit_obj, "multiplicative", component, z_grid, z0, t_grid)
    sbf <- list(
      time = sbf_time,
      iterations_used = sbf_fit_obj$iterations_used,
      final_delta = sbf_fit_obj$final_delta,
      converged = sbf_fit_obj$converged,
      effect = sbf_curves$effect,
      survival = sbf_curves$survival
    )

    pointwise[[i]] <- benchmark_make_pointwise(
      "multiplicative_local_constant", "multiplicative", "local_constant", i, seeds[[i]],
      z_grid, t_grid, truth_effect, truth_survival, origin, sbf
    )
    run_meta[[i]] <- benchmark_make_run_meta(
      "multiplicative_local_constant", "multiplicative", "local_constant", i, seeds[[i]], origin, sbf
    )
  }

  list(
    algorithm = "multiplicative_local_constant",
    model_family = "multiplicative",
    fit_variant = "local_constant",
    n = n,
    d = d,
    n_rep = n_rep,
    target_k = target_k,
    plot_title = "Multiplicative local constant",
    plot_file = "benchmark_sbf_vs_origin_multiplicative_local_constant.png",
    pointwise = do.call(rbind, pointwise),
    run_meta = do.call(rbind, run_meta)
  )
}
