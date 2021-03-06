% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.bcgfrailev.R
\name{summary.bcgfrailev}
\alias{summary.bcgfrailev}
\title{Print bcgfrailev}
\usage{
\method{summary}{bcgfrailev}(object, ...)
}
\arguments{
\item{object}{A class \code{bcgfrailev} object.}

\item{...}{ignored}
}
\description{
Generics to print the S3 class bcgfrailev.
}
\details{
Calls \code{summary.bcgfrailev()}.
}
\examples{
set.seed(24)
simdata<-simbcfrailph(p.size=300, c.rate= c(0.3),fraildistrn=c("gamma"),frail.par=c(0.5,0.5),
bhaz.arg=list(distrn = c("weibull"),shape =c(5), scale = c(0.1)),
covar.arg=list(coefs=c(2),types = c("B"),size=c(1),prob=c(0.5)))
dataa<-simdata$data

fitbcgfrailev=bcgfrailev(Surv(time,censor)~ X1+frailty(PID) ,data=dataa)
fitbcgfrailev
summary(fitbcgfrailev)

}
