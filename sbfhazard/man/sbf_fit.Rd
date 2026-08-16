
\name{sbf_fit}
\alias{sbf_fit}
\alias{sbf_fit_binning}
\title{Fit smooth backfitting hazard models}
\usage{
sbf_fit(
  data,
  bandwidth,
  family = c("additive", "multiplicative"),
  formula = NULL,
  integral_approx = c("midd", "left", "right"),
  iterations = 100L,
  kernel = "epanechnikov",
  convergence_tol = 0.001,
  local_constant = FALSE,
  identification = c("sample_mean", "integral", "origin", "jasa"),
  truth_functions = NULL,
  warn_nonconvergence = TRUE
)

sbf_fit_binning(
  data,
  bandwidth,
  time_bins,
  covariate_bins,
  family = c("additive", "multiplicative"),
  time_binning_method = c("quantile", "equal_width"),
  covariate_binning_method = c("quantile", "equal_width"),
  representative = c("midpoint", "mean"),
  formula = NULL,
  local_constant = NULL,
  iterations = 100L,
  kernel = "epanechnikov",
  convergence_tol = 0.001,
  identification = c("sample_mean", "integral", "origin", "jasa"),
  truth_functions = NULL,
  warn_nonconvergence = TRUE,
  warn_diagnostics = TRUE
)
}
\arguments{
\item{data}{A numeric data frame with \code{time} and \code{status} columns, or a data
frame used with \code{formula}. Non-formula fits use all columns other than
\code{time} and \code{status} as fitted-scale features. The complete fitting data
must contain at least one observed event with \code{status = 1}.}

\item{bandwidth}{Numeric bandwidth vector with one value per model component,
including time.}

\item{family}{Model family: \code{"additive"} or \code{"multiplicative"}.}

\item{formula}{Optional formula of the form
\code{survival::Surv(time, status) ~ features}.
Right-hand side terms are converted to fitted-scale numeric features with
\code{prodlim::model.design(maxOrder = 1, dropIntercept = TRUE)}.}

\item{integral_approx}{SBF-grid integration rule.}

\item{iterations}{Maximum number of completed backfitting iterations. One
iteration updates every model component once.}

\item{kernel}{Kernel name or custom kernel function.}

\item{convergence_tol}{Positive convergence tolerance for backfitting.}

\item{local_constant}{For additive fits, whether to use local constant
updates. Multiplicative SBF fits require \code{FALSE}.}

\item{identification}{Multiplicative component identification rule.
The default \code{"sample_mean"} centers the log component over the training
sample.}

\item{truth_functions}{Truth log-effect functions for truth-based
multiplicative identification; one function per covariate component.}

\item{warn_nonconvergence}{Whether to warn when the fit does not reach the
convergence tolerance within \code{iterations}.}

\item{time_bins, covariate_bins}{Binning resolution.}

\item{time_binning_method, covariate_binning_method}{Binning method for time
and covariate components.}

\item{representative}{Bin representative: \code{"midpoint"} or \code{"mean"}.}

\item{warn_diagnostics}{Whether to emit diagnostic warnings for numerical or
binning issues.}
}
\value{
A fitted additive or multiplicative smooth backfitting object.
}
\description{
Fit additive or multiplicative smooth backfitting hazard models, or their
binning-based approximations.
}
\details{
\code{sbf_fit()} evaluates the SBF estimators on observation-value grids.
\code{sbf_fit_binning()} returns speed-oriented binned approximations with bin
layout metadata in the fitted object.
Common fit fields include \code{alpha_backfit}, \code{x.grid}, \code{component_info},
\code{fit_settings}, \code{fit_diagnostics}, and \code{feature_names}; binned fits also
include \code{binning}.
Multiplicative diagnostics include \code{component_floor_count}, the total number
of baseline and covariate raw-update values replaced by the fixed positive
floor across completed backfitting iterations.
Multiplicative \code{fit_settings} record this fixed value as \code{min_component}.
\code{iterations_used} records completed iterations and \code{iterations_max} records
the requested iteration limit for both SBF and binned fits.

The multiplicative identification rules are: Sample-mean identification
(\code{"sample_mean"}) for training-sample log-mean centering, Integral
identification (\code{"integral"}) for fitted-grid integral mean one,
Single-point truth-based identification (\code{"origin"}) for midpoint
truth-anchor alignment, and Two-point truth-based identification (\code{"jasa"})
for the source-script convention at -0.5 and 0.5. The \code{"jasa"} rule evaluates
the fitted component at the nearest grid points and evaluates truth at the
corresponding absolute distances, following the source script. The
two truth-based rules require \code{truth_functions} and are mainly for simulation
benchmarks.
}
\examples{
dat <- sbf_simulate_data(n = 50, d = 3, family = "additive", seed = 1)
fit <- sbf_fit(
  dat,
  bandwidth = c(0.35, 0.35, 0.35),
  family = "additive",
  local_constant = TRUE,
  iterations = 5,
  warn_nonconvergence = FALSE
)

bin_fit <- sbf_fit_binning(
  dat,
  bandwidth = c(0.35, 0.35, 0.35),
  family = "additive",
  time_bins = 5,
  covariate_bins = 5,
  iterations = 5,
  warn_nonconvergence = FALSE,
  warn_diagnostics = FALSE
)

mdat <- sbf_simulate_data(n = 50, d = 3, family = "multiplicative", seed = 2)
mfit <- sbf_fit(
  mdat,
  bandwidth = c(0.35, 0.35, 0.35),
  family = "multiplicative",
  iterations = 5,
  warn_nonconvergence = FALSE
)

mbin_fit <- sbf_fit_binning(
  mdat,
  bandwidth = c(0.35, 0.35, 0.35),
  family = "multiplicative",
  time_bins = 5,
  covariate_bins = 5,
  iterations = 5,
  warn_nonconvergence = FALSE,
  warn_diagnostics = FALSE
)

}
