#'
#' @name print.coxphgmfr
#' @title Print coxphgmfr
#' @description Generics to print the S3 class coxphgmfr.
#' @details Calls \code{print.coxphgmfr()}.
#'
#' @param x A class \code{coxphgmfr} object.
#' @param ... ignored
#' @return An object of \code{print.coxphgmfr}, with some more human-readable results from \code{bcfrailph} object.
#' @export
#  (deprecated) @S3method print coxphgmfr
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
#'
print.coxphgmfr<- function(x, ...)
{
if(!inherits(x, "coxphgmfr")){
stop("Argument must be the result of coxphgmfr")}
cat("Call:\n")
print(x$call)
cat("\nn= ",length(x$censor),"and number of events=",sum(x$censor)," \n")
cat("\nRegression Coefficients:\n")
zval <- x$coefficients/x$stderr[1:length(x$coefficients)]
se=x$stderr
RTAB <- cbind(Estimate =x$coefficients,StdErr =se[1:length(x$coefficients)],
se2= c(sqrt(diag(x$vcov2))),z.value =  zval,p.value = 2*(1-pnorm(abs(zval), mean=0,sd=1)))
printCoefmat(RTAB, P.values=TRUE, has.Pvalue=TRUE)
cat("\nFrailty Distribution:gamma\n")
cat("Frailty variance =",x$frailparest,"(",se[length(se)],")\n")
cat("Log likelihood with frailty =",x$Iloglilk,"\n")
cat("Log likelihood without frailty=",x$loglilk0,"\n")
}

