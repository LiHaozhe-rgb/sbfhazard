# Shared fit input preparation for additive and multiplicative estimators.
#
# Formula handling converts raw data to fitted-scale numeric features. Prediction
# can rebuild that scale from raw RHS variables, otherwise it uses numeric input.

.sbf_as_numeric_matrix <- function(x,
                                  ncol_expected = NULL,
                                  arg_name = "data") {
  if (!(is.data.frame(x) || is.matrix(x) || is.atomic(x))) {
    .sbf_stop_bad_arg(arg_name, "a numeric vector, matrix, or data.frame", actual = sprintf("class=%s", paste(class(x), collapse = "/")))
  }
  mat <- as.matrix(x)
  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    .sbf_stop_bad_arg(arg_name, "a non-empty matrix-like object")
  }
  if (!is.null(ncol_expected) && ncol(mat) != ncol_expected) {
    .sbf_stop_bad_arg(arg_name, sprintf("an object with %d columns", ncol_expected), actual = sprintf("%d columns", ncol(mat)))
  }
  if (!is.numeric(mat)) {
    .sbf_stop_bad_arg(arg_name, "numeric values with no NA/NaN/Inf entries", actual = sprintf("storage mode=%s", storage.mode(mat)))
  }
  storage.mode(mat) <- "double"
  if (any(!is.finite(mat))) {
    .sbf_stop_bad_arg(arg_name, "numeric values with no NA/NaN/Inf entries")
  }
  mat
}

# Non-formula fit input must already be on fitted feature scale.
.sbf_normalize_fit_data <- function(data) {
  if (!is.data.frame(data)) {
    .sbf_stop_bad_arg("data", "a data.frame with time and status columns", actual = sprintf("class=%s", paste(class(data), collapse = "/")))
  }
  if (nrow(data) == 0L) {
    .sbf_stop_bad_arg("data", "a non-empty data.frame")
  }
  df <- data
  if (!all(c("time", "status") %in% names(df))) {
    .sbf_stop_bad_arg("data", "a data.frame with time and status columns")
  }
  feature_names <- setdiff(names(df), c("time", "status"))
  if (length(feature_names) == 0L) {
    .sbf_stop_bad_arg("data", "at least one feature column in addition to time and status")
  }
  df <- df[, c("time", feature_names, "status"), drop = FALSE]

  if (!all(vapply(df, is.numeric, logical(1)))) {
    .sbf_stop_bad_arg("data", "numeric time/status/feature columns")
  }
  if (any(!is.finite(as.matrix(df)))) {
    .sbf_stop_bad_arg("data", "numeric time/status/feature columns with no NA/NaN/Inf entries")
  }
  if (any(df$time < 0)) {
    .sbf_stop_bad_arg("data", "non-negative event times")
  }
  if (any(!(df$status %in% c(0, 1)))) {
    .sbf_stop_bad_arg("data", "status coded as 0/1")
  }
  df$status <- as.integer(df$status)
  if (!any(df$status == 1L)) {
    .sbf_stop_bad_arg("data", "at least one observed event with status = 1")
  }
  df
}

# Formula fit converts raw variables into the numeric design used by engines.
.sbf_formula_fit_data <- function(formula, data) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    .sbf_stop_bad_arg("formula", "a two-sided formula like survival::Surv(time, status) ~ V1 + V2")
  }
  if (!is.data.frame(data)) {
    .sbf_stop_bad_arg("data", "a data.frame for formula fit", actual = sprintf("class=%s", paste(class(data), collapse = "/")))
  }
  if (nrow(data) == 0L) {
    .sbf_stop_bad_arg("data", "a non-empty data.frame")
  }

  Terms <- stats::terms(x = formula, data = data)
  if (length(attr(Terms, "term.labels")) == 0L) {
    .sbf_stop_bad_arg("formula", "at least one right-hand-side feature")
  }

  mm <- stats::na.omit(stats::get_all_vars(stats::formula(Terms), data = data))
  if (NROW(mm) == 0L) {
    stop("No (non-missing) observations.", call. = FALSE)
  }

  response_values <- stats::model.response(stats::model.frame(stats::update(formula, ".~1"), data = mm))
  if (!is.matrix(response_values) || !all(c("time", "status") %in% colnames(response_values))) {
    .sbf_stop_bad_arg("formula", "a survival response with time and status columns")
  }

  design_matrix <- prodlim::model.design(
    Terms,
    data = mm,
    maxOrder = 1,
    dropIntercept = TRUE
  )[[1]]
  if (ncol(design_matrix) == 0L) {
    .sbf_stop_bad_arg("data", "formula design with at least one feature column")
  }
  if (is.null(colnames(design_matrix))) {
    colnames(design_matrix) <- sprintf("feature_%d", seq_len(ncol(design_matrix)))
  }

  training_levels <- attr(design_matrix, "levels")
  design_terms <- stats::delete.response(Terms)
  attr(design_terms, "intercept") <- 1L
  design_frame <- stats::model.frame(
    design_terms,
    data = mm,
    na.action = stats::na.pass
  )
  for (name in names(training_levels)) {
    values <- design_frame[[name]]
    design_frame[[name]] <- factor(
      values,
      levels = training_levels[[name]],
      ordered = is.ordered(values)
    )
  }
  training_matrix <- stats::model.matrix(design_terms, data = design_frame)
  formula_design <- list(
    terms = design_terms,
    xlevels = training_levels,
    contrasts = attr(training_matrix, "contrasts")
  )

  feature_names <- colnames(design_matrix)
  fit_data <- data.frame(
    time = as.numeric(response_values[, "time"]),
    as.data.frame(design_matrix, check.names = FALSE, stringsAsFactors = FALSE),
    status = as.numeric(response_values[, "status"]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  fit_data <- .sbf_normalize_fit_data(fit_data)
  list(
    formula = formula,
    formula_design = formula_design,
    data = fit_data,
    feature_names = feature_names
  )
}

# Public fit APIs share this branch between named data.frame input and formula input.
.sbf_prepare_fit_inputs <- function(data,
                                   formula = NULL) {
  if (is.null(formula)) {
    fit_data <- .sbf_normalize_fit_data(data)
    feature_names <- setdiff(names(fit_data), c("time", "status"))
    return(list(
      formula = NULL,
      formula_design = NULL,
      data = fit_data,
      feature_names = feature_names
    ))
  }

  .sbf_formula_fit_data(formula = formula, data = data)
}

# Try raw formula prediction first; return NULL to fall back to fitted-scale input.
.sbf_formula_prediction_matrix <- function(result, data, p) {
  if (!is.data.frame(data) || is.null(result$formula)) {
    return(NULL)
  }

  Terms <- stats::delete.response(stats::terms(result$formula))
  rhs_vars <- all.vars(stats::formula(Terms))
  if (length(rhs_vars) == 0L || !all(rhs_vars %in% names(data))) {
    return(NULL)
  }

  formula_design <- result$formula_design
  if (is.null(formula_design)) {
    .sbf_stop_bad_arg("result", "formula design metadata for raw-data prediction")
  }
  design_frame <- stats::model.frame(
    formula_design$terms,
    data = data,
    xlev = formula_design$xlevels,
    na.action = stats::na.pass
  )
  design_matrix <- stats::model.matrix(
    formula_design$terms,
    data = design_frame,
    contrasts.arg = formula_design$contrasts
  )
  design_matrix <- design_matrix[, -1L, drop = FALSE]
  if (ncol(design_matrix) != p) {
    .sbf_stop_bad_arg("data", sprintf("formula prediction data with %d fitted feature columns", p))
  }
  if (is.null(colnames(design_matrix))) {
    colnames(design_matrix) <- sprintf("feature_%d", seq_len(ncol(design_matrix)))
  }

  expected_names <- as.character(result$feature_names)
  if (length(expected_names) == p && !identical(colnames(design_matrix), expected_names)) {
    .sbf_stop_bad_arg("data", "formula prediction design matching fitted feature names")
  }
  if (any(!is.finite(design_matrix))) {
    .sbf_stop_bad_arg("data", "formula prediction design with no NA/NaN/Inf entries")
  }

  design_matrix
}

# Prediction accepts raw formula data when possible, otherwise fitted-scale covariates.
.sbf_prepare_prediction_feature_matrix <- function(result,
                                                  data,
                                                  p) {
  p <- .sbf_assert_scalar_whole_number(p, "p", min_value = 0L)
  formula_matrix <- .sbf_formula_prediction_matrix(result = result, data = data, p = p)
  if (!is.null(formula_matrix)) { # formula prediction succeeded, use that design matrix
    return(formula_matrix)
  }

  if (is.data.frame(data) && !is.null(result$feature_names) && ncol(data) == p) { # data has same number of columns as fitted features, check names if available
    data_names <- names(data)
    expected_names <- as.character(result$feature_names)
    if (!identical(data_names, expected_names)) {
      .sbf_stop_bad_arg(
        "data",
        sprintf("fitted-scale feature columns named %s", paste(expected_names, collapse = ", ")),
        actual = sprintf("columns: %s", paste(data_names, collapse = ", "))
      )
    }
  }

  if (is.atomic(data) && !is.matrix(data) && !is.data.frame(data)) { # single numeric vector, treat as one row of features
    data <- matrix(data, nrow = 1L)
  }
  .sbf_as_numeric_matrix(data, ncol_expected = p, arg_name = "data")
}
