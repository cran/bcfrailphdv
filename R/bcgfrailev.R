#'
#' @name bcgfrailev
#' @title Bivariate correlated gamma frailty modeling with Proportional hazard.
#' @description Fit a semiparametric Bivariate correlated gamma frailty model with Proportional Hazard structure.
#'
#' @param formula A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.
#' @param data A dataframe contain survival time, censor, covariate etc with data in columns.
#' @param initfrailp Initial estimates for the frailty parameters. The default is c(0.5,0.5).
#' @param control Arguments to control the fit. The default is \code{bcfrailph.control}.
#' @param ... further arguments
#'
#' @return An object of  that contains  the following components.
#' \itemize{
#'   \item \code{coefficients} - {A vector of estimated Covariate coefficients.}
#'   \item \code{frailparest} - {A vector of estimated Frailty parameters i.e. frailty variance and correlation.}
#'   \item \code{vcov2}- {Variance Covariance matrix of the Estimated Covariate coefficients obtained from the observed information matrix.}
#'   \item \code{vcovth2}-Variance Covariance matrix of the Estimated Frailty parameters obtained from the observed information matrix of the mariginal likelihood.
#'   \item \code{stderr}-{A vector containing the Standard error of the Estimated parameters both covariate coefficients and  frailty parameters.}
#'   \item \code{loglilk0}- Log likelihood of without frailty model.
#'   \item \code{loglilk}-Log likelihood of Cox PH model with frailty.
#'   \item \code{Iloglilk}- Log likelihood of with frailty model after integrating out the frailty term.
#'   \item \code{cbashaz}- array containing Cummulative baseline hazard.
#'   \item \code{X}-{Matrix of observed covariates.}
#'   \item \code{time}-{the observed survival time.}
#'   \item \code{censor}-{censoring indicator.}
#'   \item \code{resid}-{the martingale residuals.}
#'   \item \code{lin.prid}-{the vector of linear predictors.}
#'   \item \code{frail}-{estimated Frailty values.}
#'   \item \code{iteration}-{Number of outer iterations.}
#'   \item \code{e.time}-{the vector of unique event times.}
#'   \item \code{n.event}- {the number of events at each of the unique event times.}
#'   \item \code{converg}-  {TRUE if converge, FALSE otherwise.}
#'   }
#'
#' @export bcgfrailev
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats nlminb
#' @importFrom stats nlm
#' @importFrom stats terms
#' @importFrom survival Surv
#' @importFrom survival coxph
#' @importFrom bcfrailph bcfrailph
#'
#' @note Parameters of Bivariate correlated gamma frailty model was estimated basically using the EM-approach proposed by Iachine, I. A. (1995) with modifications.
#' The main modification that made on the original EM-approach was similar to the modification made on EM approach for univariate gamma frailty model by Duchateau and Janssen (2008).
#' This means following more or less similar procedure as Duchateau and Janssen (2008), frailty parameters are estimated from the marginal log likelihood function.
#' The results of both EM- approach and the modified EM- approach are similar. The difference is that the modified one is much faster. The standard error of the estimated 
#' covariate coefficients and frailty parameters are based on the inverse of the observed information matrix constructed by taking partial
#' derivatives of minus the observable likelihood (the Log likelihood obtained after integrating out the frailty term).
#' 
#' @seealso \code{\link{bcfraildv}}
#'
#' @references
#'
#' Duchateau, L., Janssen, P. (2008) The Frailty Model. Springer, New York.
#'
#' Iachine, I. A. (1995). Correlated frailty concept in the analysis of bivariate survival data. Bachelor project, Odense University, Department of Mathematics and Computer Science, Denmark.
#'
#' Klein, J. P., and Moeschberger, M. L. (2003), Survival analysis: techniques for censored and truncated data, New York: Springer.
#'
#' @examples
#' set.seed(24)
#' simdata<-simbcfrailph(p.size=300, c.rate= c(0.3),fraildistrn=c("gamma"),frail.par=c(0.5,0.5),
#' bhaz.arg=list(distrn = c("weibull"),shape =c(5), scale = c(0.1)),
#' covar.arg=list(coefs=c(2),types = c("B"),size=c(1),prob=c(0.5)))
#' dataa<-simdata$data
#'
#' fitbcgfrailev=bcgfrailev(Surv(time,censor)~ X1+frailty(PID) ,data=dataa)
#' fitbcgfrailev
#'
#'
#' \dontshow{
#' set.seed(24)
#' simdata<-simbcfrailph(p.size=300, c.rate= c(0.3),fraildistrn=c("gamma"),frail.par=c(0.5,0.5),
#' bhaz.arg=list(distrn = c("weibull"),shape =c(5), scale = c(0.1)),
#' covar.arg=list(coefs=c(2),types = c("B"),size=c(1),prob=c(0.5)))
#' dataa<-simdata$data
#'
#' fitbcgfrailev=bcgfrailev(Surv(time,censor)~ X1+frailty(PID) ,data=dataa)
#' fitbcgfrailev
#' }
#'
#' \donttest{
#' # gamma fit in uncensored data
#'
#' # simulate the data set
#' set.seed(3)
#' simdata<-simbcfrailph(p.size=300, c.rate= c(0),fraildistrn=c("gamma"),frail.par=c(0.5,0.6),
#' bhaz.arg=list(distrn = c("weibull"),shape =c(5), scale = c(0.1)),
#' covar.arg=list(coefs=c(1.5),types = c("B"),size=c(1),prob=c(0.5)))
#' dataa<-simdata$data ## the simulated data set
#'
#' #fit
#' fitbcgfrailev=bcgfrailev(Surv(time,censor)~ X1+cluster(PID) ,data=dataa)
#' fitbcgfrailev
#'
#' ## one can set the initial parameter for the frailty parameters
#' ## the default is initfrailp = c(0.5,0.5)
#' fitbcgfrailev=bcgfrailev(Surv(time,censor)~ X1+X2+frailty(PID),data=dataa,initfrailp = c(0.1,0.5))
#' fitbcgfrailev
#'
#' # Not run
#'
#' #if covariates are not included
#' fitmoe=bcgfrailev(Surv(time,censor)~0,data=dataa)
#' fitmoe
#' fitmoe=bcgfrailev(Surv(time,censor)~1,data=dataa)
#' fitmoe
#'
#' #if fraility id is not specified correctly
#' #or if it is not specified in a way that it indicates pairs.
#' ID=array(1:nrow(dataa))# this is not pair id rather it is individual id.
#' fitmoe=bcgfrailev(Surv(time,censor)~ X1+frailty(ID),data=dataa)
#' fitmoe
#' fitmoe=bcgfrailev(Surv(time,censor)~ X1+cluster(ID),data=dataa)
#' fitmoe
#'
#' # if control is not specified correctly.
#' # if one needs to change only max.iter to be 100,
#'
#' fitmoe=bcgfrailev(Surv(time,censor)~ X1+frailty(PID),data=dataa,control=c(max.iter=100))
#' fitmoe
#'
#' #the correct way is
#' fitmoe=bcgfrailev(Surv(time,censor)~ X1+frailty(PID),data=dataa,
#' control=bcfrailph.control(max.iter=100))
#' fitmoe
#'
#' #if initial frailty parameters are in the boundary of parameter space
#' fitmoe=bcgfrailev(Surv(time,censor)~ X1,data=dataa,initfrailp=c(0.2,1))
#' fitmoe
#' fitmoe=bcgfrailev(Surv(time,censor)~ X1,data=dataa,initfrailp=c(0,0.1))
#' fitmoe
#'
#' # End Not run
#' }
#'
bcgfrailev<- function(formula, data,initfrailp = NULL,control,...)
{
Call <- match.call()
if(any(is.na(charmatch(names(Call), c("formula","data", "initfrailp", "control"),
nomatch = NA_integer_)))){
stop("There is/are unused argument(s)")}
if (missing(formula)) {stop("formula must be supplied")}
if (!inherits(formula,"formula")) stop("First argument must be a formula")
if (missing(data)) {stop("Data must be supplied")}
if (!is.data.frame(data) && !is.matrix(data)) stop("Data must be data frame or matrix.")
frail_distrn<-c("gamma")
if (missing(control)){control=bcfrailph::bcfrailph.control()}
fit1=bcfrailph::bcfrailph(formula=formula, data=data,initfrailp = initfrailp ,frail_distrn=frail_distrn,control=control)
fit<-adj.SE(fit1=fit1)
fit$call <- Call
fit$formula<- formula
fit$frail_distrn<-frail_distrn
class(fit) <- c("bcgfrailev",class(fit1)) 
fit
}



##
