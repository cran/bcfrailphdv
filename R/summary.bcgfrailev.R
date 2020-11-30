#'
#' @name summary.bcgfrailev
#' @title Print bcgfrailev
#' @description Generics to print the S3 class bcgfrailev.
#' @details Calls \code{summary.bcgfrailev()}.
#'
#' @param object A class \code{bcgfrailev} object.
#' @param ... ignored
#'
#' @export
#  (deprecated) @S3method summary bcgfrailev
#' @importFrom stats printCoefmat
#' @importFrom stats pnorm
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
#' summary(fitbcgfrailev)
#'
summary.bcgfrailev<- function(object, ...)
  {
    if(!inherits(object, "bcgfrailev")){
      stop("Argument must be the result of bcgfrailev")}
    cat("Call:\n")
    print(object$call)
    cat("\nn= ",length(object$censor),"and number of events=",sum(object$censor)," \n")
    cat("\nRegression Coefficients:\n")
    zval <- object$coefficients/object$stderr[1:length(object$coefficients)]
    se=object$stderr
    RTAB <- cbind(Estimate =object$coefficients,StdErr =se[1:length(object$coefficients)],
                  se2= c(sqrt(diag(object$vcov2))),z.value =  zval,p.value = 2*(1-pnorm(abs(zval), mean=0,sd=1)))
    printCoefmat(RTAB, P.values=TRUE, has.Pvalue=TRUE)
    cat("\nFrailty Distribution:Bivariate Correlated ",object$frail_distrn,"\n")
    if(object$frail_distrn==c("gamma")){
      cat("Frailty variance =",object$frailparest[1],"(",se[length(se)-1],")\n")
      cat("Correlation Estimate =",object$frailparest[2],"(",se[length(se)],")\n")}
    if(object$frail_distrn==c("lognormal")){
      cat("Variance of random effect =",object$frailparest[1],"(",se[length(se)-1],")\n")
      cat("Correlation Estimate of random effects =",object$frailparest[2],"(",se[length(se)],")\n")}
    cat("Log likelihood with frailty =",object$Iloglilk,"\n")
    cat("Log likelihood without frailty=",object$loglilk0,"\n")
}
