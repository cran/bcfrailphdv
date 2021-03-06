% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcfraildv.R
\name{bcfraildv}
\alias{bcfraildv}
\title{Bivariate correlated frailty model with Proportional hazard structure when frailty variances of the two subjects are allowed to be different.}
\usage{
bcfraildv(
  formula,
  data,
  frail_distrn = c("gamma", "lognormal"),
  initfrailp = NULL,
  control = bcfraildv.control(),
  ...
)
}
\arguments{
\item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}

\item{data}{A dataframe contain survival time, censor, covariate etc with data in columns.}

\item{frail_distrn}{A type of frailty distribution to be used in fit. Either gamma or lognormal. The default is gamma.}

\item{initfrailp}{Initial estimates for the frailty parameters. The default is c(0.5,0.5,0.5).}

\item{control}{Arguments to control the fit. The default is \code{\link{bcfraildv.control}}.}

\item{...}{further arguments}
}
\value{
An object of  that contains  the following components.
\itemize{
  \item \code{coefficients} - {A vector of estimated Covariate coefficients.}
  \item \code{frailparest} - {A vector of estimated Frailty/random effect parameters i.e. frailty variances and correlation.}
  \item \code{vcov}- {Variance Covariance matrix of the Estimated Covariate coefficients obtained from the observed information matrix.}
  \item \code{vcovth}-{Variance Covariance matrix of the Estimated Frailty/random effect parameters obtained from the observed information matrix.}
  \item \code{stderr}-{A vector containing the Standard errors of the Estimated parameters both covariate coefficients and frailty parameters.}
  \item \code{loglilk0}- {Log likelihood of without frailty model.}
  \item \code{loglilk}-{Log likelihood of Cox PH model with frailty.}
  \item \code{Iloglilk}- {Log likelihood obtained after integrating out the frailty term in gamma fit.}
  \item \code{cbashaz}- {array containing Cummulative baseline hazard.}
  \item \code{X}-{Matrix of observed covariates.}
  \item \code{time}-{the observed survival time.}
  \item \code{censor}-{censoring indicator.}
  \item \code{resid}-{the martingale residuals.}
  \item \code{lin.prid}-{the vector of linear predictors.}
  \item \code{frail}-{estimated Frailty values.}
  \item \code{iteration}-{Number of outer iterations.}
  \item \code{e.time}-{the vector of unique event times.}
  \item \code{n.event}- {the number of events at each of the unique event times.}
  \item \code{converg}-  {TRUE if converge, FALSE otherwise.}
  }
}
\description{
Fit a semiparametric Bivariate correlated frailty models with Proportional Hazard structure where frailty/random effect variances of the two subject are allowed to be varied.
}
\note{
Parameters of Bivariate correlated gamma frailty model was estimated using EM approch. The EM approach used here is basically simmilar to the
EM-approach proposed by Iachine, I. A. (1995)with modifications. In the  proposed EM approach, the expectation part is modified and in maximization part, frailty parameters are estimated
from the marginal log likelihood function. Standard errors of the estimated covariate coefficients and the frailty parameters are obtained from the second derivatives of the marginal log likelihood function.

Parameters of Bivariate correlated lognormal frailty model was estimated in similar manner as
the penalized likelihood approach used by Ripatti and Palmgren (2000).
}
\examples{
set.seed(2)
simdata<-simbcfraildv(p.size=300, c.rate= c(0.3),
fraildistrn=c("lognormal"),frail.par=c(0.7,0.5,0.4),
bhaz.arg=list(distrn = c("weibull"),shape =c(4), scale = c(0.1)),
covar.arg=list(coefs=c(2),types = c("B"),size=c(1),prob=c(0.5)))
dataa<-simdata$data

fitbcfraildv=bcfraildv(Surv(time,censor)~ X1+frailty(PID) ,
data=dataa,frail_distrn=c("lognormal"))
fitbcfraildv

\dontshow{
set.seed(2)
simdata<-simbcfraildv(p.size=300, c.rate= c(0.3),
fraildistrn=c("lognormal"),frail.par=c(0.7,0.5,0.4),
bhaz.arg=list(distrn = c("weibull"),shape =c(4), scale = c(0.1)),
covar.arg=list(coefs=c(2),types = c("B"),size=c(1),prob=c(0.5)))
dataa<-simdata$data

fitbcfraildv=bcfraildv(Surv(time,censor)~ X1+frailty(PID),
data=dataa,frail_distrn=c("lognormal"))
fitbcfraildv

}

\donttest{
# for gamma fit

# simulate the data set
#Weibull baseline hazard with parameters shape= 4 and scale=0.1.
#a dataset with 500 pairs.
#gamma frailty distribution with frailty parameters are taken to
#be variance1=0.8, variance2=0.5 and rho=0.4.
#One binomial B(1,0.5) with regression coefficient 2.
#Each observed covariate for the two subjects in a
#pair is taken to be independent and 20 percent of the observations are censored.

set.seed(1)
simdata<-simbcfraildv(p.size=500, c.rate= c(0.2),fraildistrn=c("gamma"),
frail.par=c(0.8,0.5,0.4),
bhaz.arg=list(distrn = c("weibull"),shape =c(4), scale = c(0.1)),
covar.arg=list(coefs=c(2),types = c("B"),size=c(1),prob=c(0.5)))
dataa<-simdata$data

#fit
fitbcfraildv=bcfraildv(Surv(time,censor)~ X1+frailty(PID),
data=dataa,frail_distrn = c("gamma"))
fitbcfraildv

# or simply
fitbcfraildv=bcfraildv(Surv(time,censor)~ X1+frailty(PID),data=dataa)
fitbcfraildv

# the output looks like
# Call:
# bcfraildv(formula = Surv(time, censor) ~ X1 + frailty(PID), data = dataa,
#     frail_distrn = c("gamma"))
#
# n=  1000 and number of events= 805
#
# Regression Coefficients:
#    Estimate   StdErr      se2 z.value   p.value
# X1 1.961681 0.177475 0.084801  11.053 < 2.2e-16 ***
# ---
#
# Frailty Distribution:Bivariate Correlated  gamma
# Frailty variance 1 = 0.8232635 ( 0.1720332 )
# Frailty variance 2 = 0.4935778 ( 0.1516554 )
# Correlation Estimate = 0.4250812 ( 0.1523974 )
# Log likelihood = -4513.342

## one can set the initial parameter for the frailty parameters
##the default is initfrailp = c(0.5,0.5,0.5)
fitbcfraildv=bcfraildv(Surv(time,censor)~ X1+frailty(PID),
data=dataa,initfrailp = c(0.4,0.6,0.1))
fitbcfraildv

# Not run

#if covariates are not included
fitmoe=bcfraildv(Surv(time,censor)~0,data=dataa,
frail_distrn=c("lognormal"))
fitmoe
fitmoe=bcfraildv(Surv(time,censor)~1,data=dataa,
frail_distrn=c("lognormal"))
fitmoe

#if fraility id is not specified correctly
#or if it is not specified in a way that it indicates pairs.
ID=array(1:nrow(dataa))# this is not pair id rather it is individual id.
fitmoe=bcfraildv(Surv(time,censor)~ X1+frailty(ID),
data=dataa,frail_distrn=c("lognormal"))
fitmoe

# if control is not specified correctly.
# if one needs to change only max.iter to be 100,

fitmoe=bcfraildv(Surv(time,censor)~ X1+frailty(PID),data=dataa,
control=c(max.iter=100))
fitmoe

#the correct way is
fitmoe=bcfraildv(Surv(time,censor)~ X1+frailty(PID),data=dataa,
control=bcfraildv.control(max.iter=100))
fitmoe

#if initial frailty parameters are in the boundary of parameter space
fitmoe=bcfraildv(Surv(time,censor)~ X1,data=dataa,initfrailp=c(0.2,0.3,1))
fitmoe
fitmoe=bcfraildv(Surv(time,censor)~ X1,data=dataa,initfrailp=c(0,0.1,0.1))
fitmoe
#if gamma model correlation parameter space violation
fitmoe=bcfraildv(Surv(time,censor)~ X1,data=dataa,initfrailp=c(0.9,0.3,0.6))
fitmoe

#if a frailty distribution other than gamma and lognormal are specified

fitmoe=bcfraildv(Surv(time,censor)~ X1,data=dataa,frail_distrn=c("exp"))
fitmoe
# End Not run
}

}
\references{
Duchateau, L., Janssen, P. (2008) The Frailty Model. Springer, New York.

Iachine, I. A. (1995). Correlated frailty concept in the analysis of bivariate survival data. Bachelor project, Odense University, Department of Mathematics and Computer Science, Denmark.

Rippatti, S. and Palmgren, J (2000). Estimation of multivariate frailty models using penalized partial likelihood. Biometrics, 56: 1016-1022.

Wienke, A. (2011). Frailty Models in Survival Analysis. Chapman & Hall/CRC,Taylor & Francis Group
}
\seealso{
\code{\link{bcgfrailev}},\code{\link{simbcfraildv}}
}
