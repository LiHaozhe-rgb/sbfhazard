# Common argument validation helpers.

# Build one consistent package-boundary argument error.
.sbf_stop_bad_arg <- function(arg_name, expected, actual = NULL, hint = NULL) {
  msg <- sprintf("%s is invalid. Expected %s.", arg_name, expected)
  if (!is.null(actual) && nzchar(actual)) {
    msg <- sprintf("%s Actual: %s.", msg, actual)
  }
  if (!is.null(hint) && nzchar(hint)) {
    msg <- sprintf("%s %s", msg, hint)
  }
  stop(msg, call. = FALSE)
}

# TRUE/FALSE options such as local_constant and warn_diagnostics.
.sbf_assert_flag <- function(x, arg_name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .sbf_stop_bad_arg(arg_name, "a single TRUE/FALSE value")
  }
  invisible(x)
}

# Function-valued inputs such as a custom kernel
.sbf_assert_function <- function(x, arg_name) {
  if (!is.function(x)) {
    .sbf_stop_bad_arg(arg_name, "a function", actual = sprintf("class=%s", paste(class(x), collapse = "/")))
  }
  invisible(x)
}

# Integer-like scalar arguments such as iterations, bins, seed, and component.
.sbf_assert_scalar_whole_number <- function(x,
                                           arg_name,
                                           min_value = NULL) {
  if (is.null(x)) {
    .sbf_stop_bad_arg(arg_name, "a single whole number")
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != round(x)) {
    .sbf_stop_bad_arg(arg_name, "a single whole number", actual = paste(capture.output(str(x)), collapse = " "))
  }
  if (!is.null(min_value) && x < min_value) {
    .sbf_stop_bad_arg(arg_name, sprintf("a whole number >= %s", min_value), actual = as.character(x))
  }
  invisible(as.integer(x))
}

# Real-valued scalar arguments such as convergence_tol, rho, and min_hazard.
.sbf_assert_scalar_numeric <- function(x,
                                       arg_name,
                                       positive = FALSE) {
  if (is.null(x)) {
    .sbf_stop_bad_arg(arg_name, "a single numeric value")
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    .sbf_stop_bad_arg(arg_name, "a single finite numeric value", actual = paste(capture.output(str(x)), collapse = " "))
  }
  if (positive && x <= 0) {
    .sbf_stop_bad_arg(arg_name, "a positive numeric value", actual = as.character(x))
  }
  invisible(as.numeric(x))
}

# Numeric vector arguments such as bandwidth, times, and xout.
.sbf_assert_numeric_vector <- function(x,
                                       arg_name,
                                       min_length = 1L,
                                       expected_length = NULL,
                                       positive = FALSE) {
  if (is.null(x)) {
    .sbf_stop_bad_arg(arg_name, "a numeric vector")
  }
  if (!is.numeric(x)) {
    .sbf_stop_bad_arg(arg_name, "a numeric vector", actual = sprintf("class=%s", paste(class(x), collapse = "/")))
  }
  if (length(x) < min_length) {
    .sbf_stop_bad_arg(arg_name, sprintf("a numeric vector with length >= %d", min_length), actual = sprintf("length=%d", length(x)))
  }
  if (!is.null(expected_length) && length(x) != expected_length) {
    .sbf_stop_bad_arg(arg_name, sprintf("a numeric vector with length %d", expected_length), actual = sprintf("length=%d", length(x)))
  }
  if (any(!is.finite(x))) {
    .sbf_stop_bad_arg(arg_name, "a numeric vector with only finite values")
  }
  if (positive && any(x <= 0)) {
    .sbf_stop_bad_arg(arg_name, "a numeric vector with only positive values")
  }
  invisible(as.numeric(x))
}
