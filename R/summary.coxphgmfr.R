#'
#' @name summary.coxphgmfr
#' @title Print coxphgmfr
#' @description Generics to print the S3 class coxphgmfr.
#' @details Calls \code{summary.coxphgmfr()}.
#'
#' @param object A class \code{coxphgmfr} object.
#' @param ... ignored
#'
#' @return An object of \code{summary.coxphgmfr}, with some more human-readable results from \code{coxphgmfr} object.
#' @export
#  (deprecated) @S3method summary coxphgmfr
#' @importFrom stats printCoefmat
#' @importFrom stats pnorm
#'
#' @note The summary function is currently identical to the print function.
#' @seealso \code{\link{coxphgmfr}}
#'
#' @examples
#' set.seed(2)
#' n1=600;IID=array(1:n1)
#' X1<-runif(n1,  min=0, max=1)
#' z=rgamma(n1,shape=2,scale=0.5)
#' u1<-runif(n1,  min=0, max=1)
#' time<- 1/0.1*log(1-0.1*log(u1)/(0.0001*exp(3*X1)*z))
#' censor=rep(1,n1)
#' dataa <- data.frame(time=time, X1=X1,censor=censor,IID=IID)
#'
#' fitcoxfr=coxphgmfr(Surv(time,censor)~ X1+frailty(IID) ,data=dataa)
#' fitcoxfr
#' summary(fitcoxfr)
#' names(fitcoxfr)
#'
summary.coxphgmfr<- function(object, ...)
{
if(!inherits(object, "coxphgmfr")){
stop("Argument must be the result of coxphgmfr")}
cat("Call:\n")
print(object$call)
cat("\nn= ",length(object$censor),"and number of events=",sum(object$censor)," \n")
cat("\nRegression Coefficients:\n")
zval <- object$coefficients/object$stderr[1:length(object$coefficients)]
se=object$stderr
RTAB <- cbind(Estimate =object$coefficients,StdErr =se[1:length(object$coefficients)],
se2= c(sqrt(diag(object$vcov2))),z.value =  zval,p.value = 2*(1-pnorm(abs(zval), mean=0,sd=1)))
printCoefmat(RTAB, P.values=TRUE, has.Pvalue=TRUE)
cat("\nFrailty Distribution:gamma\n")
cat("Frailty variance =",object$frailparest,"(",se[length(se)],")\n")
cat("Log likelihood with frailty =",object$Iloglilk,"\n")
cat("Log likelihood without frailty=",object$loglilk0,"\n")
}

