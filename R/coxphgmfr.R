#'
#' @name coxphgmfr
#' @title Cox PH model with univariate and bivariate shared gamma frailty model.
#' @description Fit Cox PH model with univariate and bivariate shared gamma frailty model.
#'
#' @param formula A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.
#' @param data A dataframe contain survival time, censor, covariate etc with data in columns.
#' @param initfrailp Initial estimates for the frailty parameters. The default is c(0.5).
#' @param control Arguments to control the fit. The default is \code{\link{bcfraildv.control}}.
#' @param ... further arguments
#'
#' @return An object of  that contains  the following components.
#' \itemize{
#'   \item \code{coefficients} - {A vector of estimated Covariate coefficients.}
#'   \item \code{frailparest} - {A vector of estimated Frailty parameters i.e. frailty variance and correlation.}
#'   \item \code{vcov}- {Variance Covariance matrix of the Estimated Covariate coefficients obtained from the observed information matrix.}
#'   \item \code{vcov2}-Variance Covariance matrix of the Estimated Frailty parameters obtained from the observed information matrix of the mariginal likelihood.
#'   \item \code{stderr}-{A vector containing the Standard error of the Estimated parameters both covariate coefficients and  frailty parameter.}
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
#' @export coxphgmfr
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats nlminb
#' @importFrom stats terms
#' @importFrom survival Surv
#' @importFrom survival coxph
#'
#' @note This is just a \code{\link{coxph}} model with gamma frailty and the differences between
#' \code{\link{coxph}} with gamma frailty fit and \code{\link{coxphgmfr}} fit is on the standard errors of the
#' covariates cofficients. Here, the standard errors of the estimated covariate coefficients and the frailty variance parameter are obtained using
#' the standard errors estimation approach given in Klein and Moeschberger (2003).
#'
#' @seealso \code{\link{bcgfrailev}}
#'
#' @references
#'
#' Duchateau, L., Janssen, P. (2008) The Frailty Model. Springer, New York.
#'
#' Klein, J. P., and Moeschberger, M. L. (2003), Survival analysis: techniques for censored and truncated data, New York: Springer.
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
#' \donttest{
#' ### for SHARED bivariate gamma frailty fit
#' # simulate data
#' # 800 pairs,uniform covariate with coefficient 3.
#' # gamma frailty variable with parameters:variance =0.5 and mean 1
#' set.seed(3)
#' n=800; n1=n*2   ### 800 pairs
#' indic1=2*array(1:n)-1;indic2=2*array(1:n)
#' PID=1;e1=array(1:n);PID[indic1]=e1;PID[indic2]=e1  ### cluster
#' X1<-runif(n1,  min=0, max=1)
#' z=rgamma(n,shape=2,scale=0.5)
#' Z=1;Z[indic1]=z;Z[indic2]=z
#' u1<-runif(n1,  min=0, max=1)
#' time<- 1/0.1*log(1-0.1*log(u1)/(0.0001*exp(3*X1)*Z))
#' censor=rep(1,n1)
#' dataa <- data.frame(time=time, X1=X1,censor=censor,PID=PID)
#'
#' # fit
#' fitcoxfr=coxphgmfr(Surv(time,censor)~ X1+frailty(PID) ,data=dataa)
#' fitcoxfr
#' # Compare with coxph fit
#' cphfit <- coxph(Surv(time, censor, type = "right") ~ X1+frailty(ID),data =  dataa )
#' cphfit
#' # see the differences on the standard errors of the covariate coefficients
#'
#' # Not run
#'
#' #if data is not supplied
#' fitcoxfr=coxphgmfr(Surv(time,censor)~ X1+frailty(PID))
#' fitcoxfr
#'
#' #if covariates are not included
#' fitcoxfr=coxphgmfr(Surv(time,censor)~ 1+frailty(PID) ,data=dataa)
#' fitcoxfr
#'
#' # End Not run
#' }
#'
coxphgmfr<- function(formula, data,initfrailp = NULL,control=bcfraildv.control(),...)
{
Call <- match.call()
if (missing(formula)) {stop("formula must be supplied")}
if (!inherits(formula,"formula")) stop("First argument must be a formula")
if (missing(data)) {stop("Data must be supplied")}
if (!is.data.frame(data) && !is.matrix(data)) stop("Data must be data frame or matrix.")
special <- c("cluster","frailty","strata")
Terms <- terms(formula,special ,data)
mf <- model.frame(Terms,data)
mm <- model.matrix(Terms,mf)
pos_strata_mm <- grep(c("strata"), colnames(mm))
if(length(pos_strata_mm)>0){ stop(" strata is invalid in coxgamfr")}
pos_cluster_mm <- grep(c("cluster"), colnames(mm))
pos_frailty_mm <- grep("frailty", colnames(mm))
if( (length(pos_cluster_mm)>0)& (length(pos_frailty_mm)>0) ){
stop(" Simultaneous model fit using both frailty and cluster is invalid in coxphgmfr")}
pos_special_mm<-c(pos_cluster_mm,pos_frailty_mm)
fittype<-NULL
if(length(pos_special_mm)==0){fittype=c("Univ")}
if(length(pos_special_mm)>0){
uniq_id<-unique(mm[,pos_special_mm])
if(length(uniq_id)<(length(mm[,pos_special_mm])/2)){
stop("coxgamfr only support univariate and bivariate shared gamma frailty fits")}
if(length(uniq_id)==(length(mm[,pos_special_mm]))){fittype=c("Univ")}
if( (length(uniq_id)<length(mm[,pos_special_mm])) & (length(uniq_id)>(length(mm[,pos_special_mm])/2))){
stop("coxphgmfr only support univariate fit and shared gamma frailty fit for bivariate data")}
if(length(uniq_id)==(length(mm[,pos_special_mm])/2)){fittype=c("Shared")}
if(fittype==c("Shared")){
indic1<-2*array(1:(nrow(mm)/2))-1
indic2<-2*array(1:(nrow(mm)/2))
if(sum(mm[indic1,pos_special_mm]-mm[indic2,pos_special_mm])!=0){
order=sort(mm[,pos_special_mm], decreasing = FALSE, index.return = TRUE)
subject_indx=order$ix
mm<-mm[subject_indx,];mf<-mf[subject_indx,]}}}
X <- mm[, -c(1, pos_special_mm), drop = FALSE]
if(ncol(X)<1){stop("covariates must be included")}
mmatt <- attributes(mm)
attr(X, "assign") <- mmatt$assign[-c(1, pos_cluster_mm,pos_frailty_mm)]
attr(X, "contrasts") <- mmatt$contrasts
Y <- mf[[1]]
if (!inherits(Y, "Surv")){stop("Response is not a survival object")}
type <- attr(Y, "type")
if (type != "right"){stop(paste("coxphgmfr doesn't support \"", type,
"\" survival data", sep = ""))}
if (ncol(Y) != 3) {Y <- Surv(rep(0, nrow(Y)), Y[, 1], Y[, 2])}
control$lower<-control$lower[1]
control$upper<-control$upper[1]
if((control$lower<0)|(control$upper<0)){stop("Invalid lower or upper bound for frailty parameter")}
if(fittype==c("Univ")){
fit<-.univgamfit(X,Y,initfrailp,control)}
if(fittype==c("Shared")){
fit<-.sharedgamfit(X,Y,initfrailp,control)}
fit$call <- Call
fit$formula<- formula
class(fit) <- c("coxphgmfr")
fit
}



##
