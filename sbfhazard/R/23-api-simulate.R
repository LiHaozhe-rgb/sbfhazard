#' Simulate example hazard data
#'
#' Generate small additive or multiplicative hazard simulation datasets for
#' examples and tests.
#'
#' @param n Number of rows.
#' @param d Total number of model components, including time.
#' @param rho Equicorrelation value for simulated covariates.
#' @param family Simulation family: `"additive"` or `"multiplicative"`.
#' @param model Simulation model identifier. Model 1 uses exponential event and
#'   censoring times. Model 2 uses a Makeham setting and requires package
#'   VGAM.
#' @param violate_cox Whether to use the nonlinear non-Cox-style effect setting.
#' @param seed Optional random seed.
#'
#' @return A data frame with `time`, simulated covariates, and `status`.
#'
#' @details
#' The simulation utilities require package MASS. Makeham model 2 also requires
#' package VGAM. The `family` argument selects the data-generating process; it
#' does not make the additive and multiplicative mechanisms identical.
#'
#' @examples
#' dat <- sbf_simulate_data(n = 20, d = 3, family = "additive", seed = 1)
#' head(dat)
#'
#' mdat <- sbf_simulate_data(n = 20, d = 3, family = "multiplicative", seed = 2)
#' head(mdat)
#'
#' @export

sbf_simulate_data <- function(n = 200L,
                              d = 5L,
                              rho = 0,
                              family = c("additive", "multiplicative"),
                              model = 1L,
                              violate_cox = TRUE,
                              seed = NULL) {
  family <- match.arg(family)
  n <- .sbf_assert_scalar_whole_number(n, "n", min_value = 1L)
  d <- .sbf_assert_scalar_whole_number(d, "d", min_value = 2L)
  rho <- .sbf_assert_scalar_numeric(rho, "rho")
  model <- .sbf_assert_scalar_whole_number(model, "model", min_value = 1L)
  .sbf_assert_flag(violate_cox, "violate_cox")
  if (!is.null(seed)) {
    seed <- .sbf_assert_scalar_whole_number(seed, "seed", min_value = 0L)
  }

  if (identical(family, "additive")) {
    return(.sbf_simulate_additive_data(
      n = n,
      d = d,
      rho = rho,
      model = model,
      violate_cox = violate_cox,
      seed = seed
    ))
  }

  .sbf_simulate_multiplicative_data(
    n = n,
    d = d,
    rho = rho,
    model = model,
    violate_cox = violate_cox,
    seed = seed
  )
}
