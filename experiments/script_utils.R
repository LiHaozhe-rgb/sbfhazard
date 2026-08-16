# Shared helpers for package-facing research scripts.

# Truth-curve definitions (model phi/survival) live in their own file.
source(file.path("experiments", "truth_curves.R"))

sbf_project_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(root, "sbfhazard", "DESCRIPTION"))) {
    stop("Run this script from the repository root.", call. = FALSE)
  }
  root
}

sbf_load_package <- function(root = sbf_project_root()) {
  devtools::load_all(file.path(root, "sbfhazard"), export_all = TRUE, quiet = TRUE)
  invisible(TRUE)
}

sbf_results_dir <- function(root, tier, ...) {
  path <- file.path(root, "test_results", tier, ...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

sbf_mean_or_na <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  mean(x)
}

sbf_expand_bandwidth <- function(v, d) {
  v <- as.numeric(v)
  d <- as.integer(d)
  if (length(v) == 1L) {
    return(rep(v, d))
  }
  if (length(v) == d) {
    return(v)
  }
  stop("bandwidth must be scalar or length d.", call. = FALSE)
}

sbf_write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
