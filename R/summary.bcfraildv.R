#'
#' @name summary.bcfraildv
#' @title Print bcfraildv
#' @description Generics to print the S3 class bcfraildv.
#' @details Calls \code{summary.bcfraildv()}.
#'
#' @param object A class \code{bcfraildv} object.
#' @param ... ignored
#'
#' @return An object of \code{print.bcfraildv}, with some more human-readable results from \code{bcfraildv} object.
#' @export
#  (deprecated) @S3method summary bcfraildv
#' @importFrom stats printCoefmat
#' @importFrom stats pnorm
#'
#' @note The summary function is currently identical to the print function.
#' @seealso \code{\link{bcfraildv}}
#'
#' @examples
#' set.seed(24)
#' simdata<-simbcfraildv(p.size=300, c.rate= c(0.3),fraildistrn=c("gamma"),frail.par=c(0.5,0.5,0.5),
#' bhaz.arg=list(distrn = c("weibull"),shape =c(5), scale = c(0.1)),
#' covar.arg=list(coefs=c(2),types = c("B"),size=c(1),prob=c(0.5)))
#' dataa<-simdata$data
#'
#' fitbcfraildv=bcfraildv(Surv(time,censor)~ X1+frailty(PID) ,data=dataa,frail_distrn=c("gamma"))
#' fitbcfraildv
#' summary(fitbcfraildv)
#' names(fitbcfraildv)
#'
summary.bcfraildv<- function(object, ...)
{
if(!inherits(object, "bcfraildv")){
stop("Argument must be the result of bcfraildv")}
cat("Call:\n")
print(object$call)
cat("\nn= ",length(object$censor),"and number of events=",sum(object$censor)," \n")
cat("\nRegression Coefficients:\n")
zval <- object$coefficients/object$stderr[1:length(object$coefficients)]
se=object$stderr
RTAB <- cbind(Estimate =object$coefficients,StdErr =se[1:length(object$coefficients)],
se2= c(sqrt(diag(object$vcov))),z.value =  zval,p.value = 2*(1-pnorm(abs(zval), mean=0,sd=1)))
printCoefmat(RTAB, P.values=TRUE, has.Pvalue=TRUE)
cat("\nFrailty Distribution:Bivariate Correlated ",object$frail_distrn,"\n")
if(object$frail_distrn==c("gamma")){
cat("Frailty variance 1 =",object$frailparest[1],"(",se[length(se)-2],")\n")
cat("Frailty variance 2 =",object$frailparest[2],"(",se[length(se)-1],")\n")
cat("Correlation Estimate =",object$frailparest[3],"(",se[length(se)],")\n")}
if(object$frail_distrn==c("lognormal")){
cat("Variance of random effect 1 =",object$frailparest[1],"(",se[length(se)-2],")\n")
cat("Variance of random effect 2 =",object$frailparest[2],"(",se[length(se)-1],")\n")
cat("Correlation Estimate of random effects =",object$frailparest[3],"(",se[length(se)],")\n")}
cat("Log likelihood =",object$Iloglilk,"\n")
}

