
\name{sbf_predict}
\alias{sbf_predict}
\title{Predict from smooth backfitting hazard models}
\usage{
sbf_predict(
  result,
  type = c("component", "hazard", "survival"),
  component = NULL,
  xout = NULL,
  data = NULL,
  times = NULL,
  clamp_nonnegative = TRUE,
  min_hazard = 1e-08,
  prediction = c("formula", "interpolation"),
  warn_diagnostics = TRUE
)
}
\arguments{
\item{result}{A fit object returned by \code{sbf_fit()} or
\code{sbf_fit_binning()}.}

\item{type}{Prediction type: \code{"component"}, \code{"hazard"}, or \code{"survival"}.}

\item{component}{Component index to evaluate when \code{type = "component"}.}

\item{xout}{Evaluation points when \code{type = "component"}. These are on the
fitted feature scale for formula fits.}

\item{data}{New covariate data when predicting hazards or survival. For
formula fits, a data frame containing the raw right-hand-side variables is
transformed automatically using the factor levels and contrasts actually
used by the training design; unseen factor levels are rejected.
Matrix/vector inputs are treated as fitted-scale covariates.}

\item{times}{Optional evaluation times for hazard or survival prediction.}

\item{clamp_nonnegative}{Whether to replace non-finite hazards and hazards
below \code{min_hazard} by \code{min_hazard}. If \code{FALSE}, finite hazards are left
unchanged and non-finite hazards raise an error.}

\item{min_hazard}{Positive lower bound used when clamping hazards.}

\item{prediction}{SBF-fit prediction method. \code{"formula"} evaluates the
fitted estimating equations at new points. \code{"interpolation"} linearly
interpolates fitted SBF-grid component values and uses the nearest fitted
value outside the grid. Binning fits only support \code{"formula"}.}

\item{warn_diagnostics}{Whether to emit diagnostic warnings.}
}
\value{
A component vector, hazard prediction object, or survival prediction
object.
}
\description{
Predict components, hazards, or survival curves from fitted smooth
backfitting hazard models.
}
\details{
\code{component = 1} is the baseline/time component. Components \code{2} and above
refer to the fitted feature components in \code{component_info}. Component
prediction \code{xout} is always on the fitted feature scale.

Formula and interpolation survival prediction use an integration grid formed
from the fitted support start, the baseline grid, and the requested times,
then apply trapezoidal integration. Interpolation continues with the
nearest-boundary hazard after the last fitted baseline point. At or before
the support start, survival is one.
}
\examples{
dat <- sbf_simulate_data(n = 50, d = 3, family = "additive", seed = 3)
fit <- sbf_fit(
  dat,
  bandwidth = c(0.35, 0.35, 0.35),
  family = "additive",
  formula = survival::Surv(time, status) ~ V1 + V2,
  local_constant = TRUE,
  iterations = 5,
  warn_nonconvergence = FALSE
)

sbf_predict(fit, type = "component", component = 2, xout = c(-0.2, 0, 0.2))
sbf_predict(
  fit,
  type = "hazard",
  data = data.frame(V1 = 0, V2 = 0),
  times = c(0.25, 0.5),
  warn_diagnostics = FALSE
)
sbf_predict(
  fit,
  type = "survival",
  data = data.frame(V1 = 0, V2 = 0),
  times = c(0.25, 0.5),
  warn_diagnostics = FALSE
)

}
