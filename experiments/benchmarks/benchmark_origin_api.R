sbf_fit_origin <- function(data,
                           bandwidth,
                           family = c("additive", "multiplicative"),
                           origin_engine,
                           integral_approx = c("midd", "left", "right"),
                           iterations = 100L,
                           kernel = "epanechnikov",
                           convergence_tol = 0.001,
                           local_constant = FALSE,
                           truth_functions = NULL) {
  family <- match.arg(family)
  integral_approx <- match.arg(integral_approx)

  prepared <- .sbf_prepare_fit_inputs(data = data, formula = NULL)
  fit_data <- prepared$data
  d <- ncol(fit_data) - 1L
  bandwidth <- .sbf_assert_numeric_vector(
    bandwidth,
    "bandwidth",
    expected_length = d,
    positive = TRUE
  )
  iterations <- .sbf_assert_scalar_whole_number(
    iterations,
    "iterations",
    min_value = 1L
  )
  origin_iterations <- iterations + 1L
  convergence_tol <- .sbf_assert_scalar_numeric(
    convergence_tol,
    "convergence_tol",
    positive = TRUE
  )
  .sbf_assert_flag(local_constant, "local_constant")

  kernel_spec <- .sbf_kernel_spec(kernel, arg_name = "kernel")
  feature_names <- prepared$feature_names
  origin_formula <- stats::as.formula(sprintf(
    "survival::Surv(time, status) ~ %s",
    paste(feature_names, collapse = " + ")
  ))

  if (identical(family, "additive")) {
    fit <- origin_engine(
      formula = origin_formula,
      data = fit_data,
      bandwidth = bandwidth,
      weight = "sw",
      x.grid = NULL,
      n.grid.additional = 0L,
      integral.approx = integral_approx,
      it = origin_iterations,
      kern = kernel_spec$kernel,
      initial = NULL,
      convergence_tol = convergence_tol,
      kcorr = TRUE,
      LC = local_constant,
      wrong = FALSE,
      classic.backfit = FALSE,
      print = FALSE
    )

    fit$model_type <- "additive"
    fit$fit_settings <- list(
      bandwidth = bandwidth,
      integral_approx = integral_approx,
      iterations = iterations,
      convergence_tol = convergence_tol,
      kernel = kernel_spec$kernel,
      kernel_name = kernel_spec$name,
      kernel_correction = TRUE,
      local_constant = local_constant,
      wrong_taylor = FALSE,
      classic_backfit = FALSE,
      print_trace = FALSE
    )
    class(fit) <- unique(c("sbf_additive_fit", class(fit)))
  } else {
    if (!identical(local_constant, FALSE)) {
      stop("local_constant must be FALSE when family = \"multiplicative\".", call. = FALSE)
    }
    .sbf_assert_multiplicative_truth_functions(truth_functions, p = d - 1L)

    fit <- origin_engine(
      data = fit_data,
      b.grid = bandwidth,
      x.grid = NULL,
      n.grid.additional = 0L,
      integral.approx = integral_approx,
      it = origin_iterations,
      kern = kernel_spec$kernel,
      initial = NULL,
      convergence_tol = convergence_tol,
      verbose = FALSE,
      phi = truth_functions
    )

    fit$model_type <- "multiplicative"
    fit$fit_settings <- list(
      bandwidth = bandwidth,
      integral_approx = integral_approx,
      iterations = iterations,
      convergence_tol = convergence_tol,
      kernel = kernel_spec$kernel,
      kernel_name = kernel_spec$name
    )
    class(fit) <- unique(c("sbf_multiplicative_fit", class(fit)))
  }

  # Origin engines label initialization as state 1 and update in states 2:it.
  fit$origin_state_label <- as.integer(fit$l)
  fit$l <- NULL
  fit$iterations_used <- fit$origin_state_label - 1L
  fit$iterations_max <- iterations
  fit$converged <- isTRUE(fit$converged)
  fit$data <- fit_data
  fit["formula"] <- list(NULL)
  fit$formula_design <- NULL
  fit$feature_names <- feature_names

  fit$component_info <- .sbf_component_info_from_grid(
    x_grid = fit$x.grid,
    component_names = c("time", feature_names),
    roles = c("time", rep("covariate", length(feature_names))),
    support_min = NULL,
    support_max = NULL
  )

  fit
}
