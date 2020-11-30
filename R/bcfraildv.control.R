#'
#' @name bcfraildv.control
#' @title Arguments for controlling bcfraildv fits.
#' @description This is used to set various numeric parameters controlling a bcfraildv model fit as a single list.
#'
#' @param max.iter Maximum number of iterations allowed in gamma frailty fit. The default is 5000.
#' @param max.iter2 Maximum number of iterations allowed in lognormal frailty fit. The default is 350.
#' @param tol A tolerance for convergence in gamma frailty fit i.e the maximum absolute differences between succssive iterations.The default is 1e-06.
#' @param toll tolerance for convergence in lognormal frailty fit.The default is 1e-05.
#' @param lower vectors of lower  bounds of the frailty parameters.
#' @param upper vectors of  upper bounds of the frailty parameters.
#' @param fastfit if true, an algorithm that make lognormal frailty fit more faster will be used. We sugest to leave it as it is.
#' @param obt.se Logical. If TRUE, for gamma fit, standard errors will be obtaind using the proposed method else observed information matrix will be used.. We sugest to leave it as it is.
#' @param fscale argument used to control \link{nlm} fits used.
#' @param print.level argument used to control \link{nlm} fits used.
#' @param ndigit argument used to control \link{nlm} fits used.
#' @param steptol argument used to control \link{nlm} fits used.
#' @param iterlim argument used to control \link{nlm} fits used.
#' @param gradtol argument used to control \link{nlm} fits used.
#' @param check.analyticals arguments used to control \link{nlm} fits used.
#' @param nlminb_control Arguments used to control \link{nlminb} fits used.
#' @return The above control parameters in a list.
#'
#' @export bcfraildv.control
#'
#' @seealso \code{\link{bcfraildv}}
#'
bcfraildv.control<-function (max.iter=5000,max.iter2=350,tol=1e-06,toll=1e-05,
lower = c(0,0,0), upper = c(Inf,Inf,1),fastfit=TRUE,obt.se=TRUE,
fscale = 1, print.level = 0, ndigit = 12,steptol = 1e-06,
iterlim = 100,gradtol =1e-8,check.analyticals=FALSE,nlminb_control =list()){
if(!missing(nlminb_control)){nlminb_control=nlminb_control}
else{rel.tol=1e-10
nlminb_control=list(eval.max=300,iter.max=300,trace=0,abs.tol=1e-20,rel.tol= 1e-10,x.tol=1.5e-8,
xf.tol= 2.2e-14,step.min=1,step.max=1,sing.tol=rel.tol)}
res<-list(max.iter=max.iter,max.iter2=max.iter2,tol=tol,toll=toll,lower =lower,upper =upper,
fastfit=fastfit,obt.se=obt.se,gradtol=gradtol,fscale =fscale,print.level = print.level,ndigit=ndigit,steptol=steptol,
iterlim =iterlim,check.analyticals=check.analyticals,
nlminb_control=list(eval.max=nlminb_control$eval.max,iter.max=nlminb_control$iter.max,
trace=nlminb_control$trace,abs.tol=nlminb_control$abs.tol,rel.tol=nlminb_control$rel.tol,
x.tol=nlminb_control$x.tol,xf.tol=nlminb_control$xf.tol,step.min=nlminb_control$step.min,
step.max=nlminb_control$step.max,sing.tol=nlminb_control$sing.tol))
class(res) <- c("bcfraildv.control")
res
}
##
