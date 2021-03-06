% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.coxphgmfr.R
\name{summary.coxphgmfr}
\alias{summary.coxphgmfr}
\title{Print coxphgmfr}
\usage{
\method{summary}{coxphgmfr}(object, ...)
}
\arguments{
\item{object}{A class \code{coxphgmfr} object.}

\item{...}{ignored}
}
\value{
An object of \code{summary.coxphgmfr}, with some more human-readable results from \code{coxphgmfr} object.
}
\description{
Generics to print the S3 class coxphgmfr.
}
\details{
Calls \code{summary.coxphgmfr()}.
}
\note{
The summary function is currently identical to the print function.
}
\examples{
set.seed(2)
n1=600;IID=array(1:n1)
X1<-runif(n1,  min=0, max=1)
z=rgamma(n1,shape=2,scale=0.5)
u1<-runif(n1,  min=0, max=1)
time<- 1/0.1*log(1-0.1*log(u1)/(0.0001*exp(3*X1)*z))
censor=rep(1,n1)
dataa <- data.frame(time=time, X1=X1,censor=censor,IID=IID)

fitcoxfr=coxphgmfr(Surv(time,censor)~ X1+frailty(IID) ,data=dataa)
fitcoxfr
summary(fitcoxfr)
names(fitcoxfr)

}
\seealso{
\code{\link{coxphgmfr}}
}
