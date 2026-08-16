
\name{sbf_simulate_data}
\alias{sbf_simulate_data}
\title{Simulate example hazard data}
\usage{
sbf_simulate_data(
  n = 200L,
  d = 5L,
  rho = 0,
  family = c("additive", "multiplicative"),
  model = 1L,
  violate_cox = TRUE,
  seed = NULL
)
}
\arguments{
\item{n}{Number of rows.}

\item{d}{Total number of model components, including time.}

\item{rho}{Equicorrelation value for simulated covariates.}

\item{family}{Simulation family: \code{"additive"} or \code{"multiplicative"}.}

\item{model}{Simulation model identifier. Model 1 uses exponential event and
censoring times. Model 2 uses a Makeham setting and requires package
VGAM.}

\item{violate_cox}{Whether to use the nonlinear non-Cox-style effect setting.}

\item{seed}{Optional random seed.}
}
\value{
A data frame with \code{time}, simulated covariates, and \code{status}.
}
\description{
Generate small additive or multiplicative hazard simulation datasets for
examples and tests.
}
\details{
The simulation utilities require package MASS. Makeham model 2 also requires
package VGAM. The \code{family} argument selects the data-generating process; it
does not make the additive and multiplicative mechanisms identical.
}
\examples{
dat <- sbf_simulate_data(n = 20, d = 3, family = "additive", seed = 1)
head(dat)

mdat <- sbf_simulate_data(n = 20, d = 3, family = "multiplicative", seed = 2)
head(mdat)

}
