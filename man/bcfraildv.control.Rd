% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcfraildv.control.R
\name{bcfraildv.control}
\alias{bcfraildv.control}
\title{Arguments for controlling bcfraildv fits.}
\usage{
bcfraildv.control(
  max.iter = 5000,
  max.iter2 = 350,
  tol = 1e-06,
  toll = 1e-05,
  lower = c(0, 0, 0),
  upper = c(Inf, Inf, 1),
  fastfit = TRUE,
  obt.se = TRUE,
  fscale = 1,
  print.level = 0,
  ndigit = 12,
  steptol = 1e-06,
  iterlim = 100,
  gradtol = 1e-08,
  check.analyticals = FALSE,
  nlminb_control = list()
)
}
\arguments{
\item{max.iter}{Maximum number of iterations allowed in gamma frailty fit. The default is 5000.}

\item{max.iter2}{Maximum number of iterations allowed in lognormal frailty fit. The default is 350.}

\item{tol}{A tolerance for convergence in gamma frailty fit i.e the maximum absolute differences between succssive iterations.The default is 1e-06.}

\item{toll}{tolerance for convergence in lognormal frailty fit.The default is 1e-05.}

\item{lower}{vectors of lower  bounds of the frailty parameters.}

\item{upper}{vectors of  upper bounds of the frailty parameters.}

\item{fastfit}{if true, an algorithm that make lognormal frailty fit more faster will be used. We sugest to leave it as it is.}

\item{obt.se}{Logical. If TRUE, for gamma fit, standard errors will be obtaind using the proposed method else observed information matrix will be used.. We sugest to leave it as it is.}

\item{fscale}{argument used to control \link{nlm} fits used.}

\item{print.level}{argument used to control \link{nlm} fits used.}

\item{ndigit}{argument used to control \link{nlm} fits used.}

\item{steptol}{argument used to control \link{nlm} fits used.}

\item{iterlim}{argument used to control \link{nlm} fits used.}

\item{gradtol}{argument used to control \link{nlm} fits used.}

\item{check.analyticals}{arguments used to control \link{nlm} fits used.}

\item{nlminb_control}{Arguments used to control \link{nlminb} fits used.}
}
\value{
The above control parameters in a list.
}
\description{
This is used to set various numeric parameters controlling a bcfraildv model fit as a single list.
}
\seealso{
\code{\link{bcfraildv}}
}