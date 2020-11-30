#'
#' @name simbcfraildv
#' @title Simulate data from bivariate correlated frailty models.
#' @description Simulate data from bivariate correlated gamma and log-normal frailty models with observed covariates. It allows inclussion of one or two covariates. In addition, frailty variances of the two artificial subjects can be equal or different.
#'
#' @param p.size pair size.
#' @param c.rate censored rate. The default is zero..
#' @param fraildistrn A type of frailty distribution to be used. Either gamma or log-normal.
#' @param frail.par vector of frailty parameters, variance and correlation respectively. The default is c(0.5,0.5,0.25) meaning frailty variance of subject1 is 0.5, frailty variance of subject2 is 0.5, and correlation 0.25.
#' @param bhaz.arg is a \code{\link{list}} i.e,\code{list(distrn = c("weibull"),shape =c(2.5), scale = c(0.01), rate = c(0.5))}.\code{distrn} is the type of baseline hazard to be used. weibull, gompertz and exponential is possible. \code{shape}, \code{scale} and \code{rate} are the parameters of the coresponding baseline hazard parameters. rate needs to be specified if the baseline is exponential.
#' @param covar.arg is a \code{\link{list}} i.e,\code{coefs} is covariate coefficients, \code{types} is the type of covariate to be used. One of the following can be specified \code{c("BU","BB","UU","B","U")}. Use \code{"BU"} if one binomial and one uniform covariates are desired.\code{"BB"} is for two binomial covariates. \code{"B"} is for binomial and \code{"U"} is for uniform.\code{size} and \code{prob} needs to be specified if binomial covariate(s) are used and if uniform covariate(s) are used, then \code{min} and \code{max} should be specified.
#'
#' @return
#' An object of class \code{simbcfraildv} that contain the following:
#'
#' \itemize{
#'
#' \item{\code{data}}  {A data frame i.e, the simulated data set. IID is individual Id, PID is pair ID, time is the simulated survival time, censor is censoring indicator and X1 or X2 denote the simulated covariate.}
#'
#' \item{\code{numberofpair}}  {The specified number of pairs.}
#'
#' \item{\code{censoredrate} } {The specified censored rate.}
#'
#' \item{\code{fraildist} } {The specified frailty distribution.}
#'
#' \item{\code{frailpar}} {The specified frailty parameters.}
#' }
#'
#'
#' @export simbcfraildv
#'
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats rbinom
#' @importFrom stats quantile
#'
#' @seealso \code{\link{bcfraildv}}
#'
#' @examples
#' simdata<-simbcfraildv(p.size=1000, c.rate= c(0.2),fraildistrn=c("gamma"),frail.par=c(0.5,0.5,0.5),
#' bhaz.arg=list(distrn = c("gompertz"),shape =c(3), scale = c(0.1)),
#' covar.arg=list(coefs=c(1),types = c("U"),min=0,max=1))
#' simdata
#'
#' \donttest{
#' #Let us simulate a data set with the following parameters
#' #weibull baseline hazard with parameters shape= 2.5 and scale=0.01.
#' #a dataset with 1000 pairs. Frailty distribution is gamma
#' # and the frailty parameters are taken to
#' #be variance1=0.6,variance2=0.4 and correlation =0.2.
#' #One binomial covariate i.e,
#' #(Binomial (1,0.5)) with regression coefficient 0.5.
#' #Each observed covariate for the two individuals in a
#' #pair is taken to be independent and 20 percent
#' #observations are censored.
#'
#' #simulate the data set
#'
#' set.seed(1)
#' simdata<-simbcfraildv(p.size=1000, c.rate= c(0.2),
#' fraildistrn=c("gamma"),frail.par=c(0.6,0.4,0.2),
#' bhaz.arg=list(distrn = c("weibull"),shape =c(2.5), scale = c(0.01)),
#' covar.arg=list(coefs=c(0.5),types = c("B"),size=1,prob=0.5))
#'
#' #to extract the simulated data set
#' dataa<-simdata$data ## the simulated data set
#' dataa[1:4,] # the first four rows looks like
#'
#' #  IID PID     time censor X1
#' #1   1   1 2.704927      1  1
#' #2   2   1 5.418071      0  1
#' #3   3   2 4.602736      1  0
#' #4   4   2 6.205303      1  0
#'
#' # IID is individual indicator
#' # PID is pair indicator
#' # time is the simulated survival time
#' # censor is the simulated censoring indicator
#' # X1 is the simulated covariate
#'
#' # if log-normal frailty is desired
#'
#' simdata<-simbcfraildv(p.size=1000, c.rate= c(0.2),
#' fraildistrn=c("lognormal"),frail.par=c(0.6,0.4,0.2),
#' bhaz.arg=list(distrn = c("weibull"),shape =c(2.5), scale = c(0.01)),
#' covar.arg=list(coefs=c(0.5),types = c("B"),size=1,prob=0.5))
#'
#' dataa<-simdata$data # the simulated data set
#'
#' # if log-normal frailty with two covariates
#' #i.e., binomial (Binomial (1,0.5)) and uniform U[0,1] is desired
#'
#' simdata<-simbcfraildv(p.size=1000, c.rate= c(0.2),
#' fraildistrn=c("lognormal"),frail.par=c(0.6,0.4,0.2),
#' bhaz.arg=list(distrn = c("weibull"),shape =c(2.5), scale = c(0.01)),
#' covar.arg=list(coefs=c(0.5),types = c("BU"),
#' size=1,prob=0.5,min=c(0),max=c(1)))
#'
#' dataa<-simdata$data ## the simulated data set
#'
#' # Not run
#' # if p.size, pair size missed
#' simdata<-simbcfraildv( c.rate= c(0.2),fraildistrn=c("gamma"),
#' frail.par=c(0.6,0.4,0.2),
#' covar.arg=list(coefs=c(0.5),types = c("B"),size=1,prob=0.5))
#'
#' # if frailty distribution other than gamma and lognormal specified
#'
#' simdata<-simbcfraildv(p.size=100, c.rate= c(0.2),fraildistrn=c("exp"),frail.par=c(0.6,0.4,0.6),
#' covar.arg=list(coefs=c(0.5),types = c("B"),size=1,prob=0.5))
#' # End Not run
#' }
#'
simbcfraildv<-function(p.size, c.rate= c(0),fraildistrn,frail.par=c(0.5,0.5,0.25),
bhaz.arg=list(distrn = c("weibull"),shape =c(2.5), scale = c(0.01), rate = c(0.5)),
covar.arg=list(coefs=c(0.5,0),types = c("BU","BB","UU","B","U"),
size=c(1,1),prob=c(0.5,0.6),min=c(0,0),max=c(1,3))){
Call <- match.call()
if (missing(p.size)) {stop("number of pairs must be supplied")}
p.size<-round(p.size,digits=0)
p.size=p.size[1];c.rate=c.rate[1]
if (p.size<=0) {stop(" number of pairs must be positive")}
if ((c.rate<0) |(c.rate>=1)) {stop(" Censoring rate must be between 0 and 1")}
if (missing(fraildistrn)) {fraildistrn<-c("gamma")}
if(length(fraildistrn)>1){stop("simbcfraildv doesn't support more than one frailty distributions")}
if(!(c(fraildistrn)%in%c("gamma","lognormal"))){
stop("simbcfraildv only support gamma or lognormal frailty distributions")}
if (length(frail.par)<2) {stop(" frailty parameters i.e. variance and correlation must be supplied")}
if (length(frail.par)==2) {frail.par<-c(frail.par[1],frail.par[1],frail.par[2])}
if (length(frail.par)>3) {frail.par<-frail.par[1:3]}
if(any(frail.par<=0)){stop("At least one invalid frailty parameter specified")}
if (frail.par[3]>=1) {stop(" invalid frailty correlation parameter")}
if(fraildistrn==c("gamma")){
parbound=min(c((sqrt(frail.par[1])/sqrt(frail.par[2])),(sqrt(frail.par[2])/sqrt(frail.par[1]))))
if(frail.par[3]>=(parbound-0.0001)){ stop(" Invalid parameter specification or in boundary of parameter space for gamma model when frailty variances are different")}}
bhaz.arg$distrn=bhaz.arg$distrn[1];bhaz.arg$shape=bhaz.arg$shape[1]
bhaz.arg$scale=bhaz.arg$scale[1];bhaz.arg$rate=bhaz.arg$rate[1]
if(!(c(bhaz.arg$distrn)%in%c("weibull","gompertz","exponential"))){
stop("simbcfraildv only support weibull gompertz or exponential baseline hazards")}
if((pmatch(c(bhaz.arg$distrn),c("weibull","gompertz","exponential"))<3)){
if(any(c(bhaz.arg$shape,bhaz.arg$scale)<=0)){stop(" invalid baseline hazard parameters")}}
if((pmatch(c(bhaz.arg$distrn),c("weibull","gompertz","exponential"))==3)){
if(bhaz.arg$rate<=0){stop(" invalid rate parameter for exponential baseline hazard")}}
if (length(covar.arg$types)>1) {covar.arg$types<-covar.arg$types[1]}
if(!(c(covar.arg$types)%in%c("BU","BB","UU","B","U"))){
stop("incorrect covariates type and types shuld be either of BU, BB, UU, B or U")}
if(pmatch(c(covar.arg$types),c("BU","BB","UU","B","U"))<=3){
if(length(covar.arg$coefs)<1){stop(" provide covariate coefficients")}
if(length(covar.arg$coefs)==1){covar.arg$coefs<-c(covar.arg$coefs[1],covar.arg$coefs[1])}
if(length(covar.arg$coefs)>=2){covar.arg$coefs<-covar.arg$coefs[1:2]}
if(pmatch(c(covar.arg$types),c("BU","BB","UU","B","U"))==1){
if( (length(covar.arg$prob)<1)|(length(covar.arg$size)<1)){stop(" provide size and prob parameter for the binomial covariate")}
if((length(covar.arg$min)<1)|(length(covar.arg$max)<1)){stop(" provide min and max parameters for the uniform covariate")}
covar.arg$prob<-covar.arg$prob[1];covar.arg$size<-covar.arg$size[1]
covar.arg$min<-covar.arg$min[1];covar.arg$max<-covar.arg$max[1]
if(c(covar.arg$max-covar.arg$min)<=0){stop(" invalid min and max parameters for uniform covariate")}
if((covar.arg$prob<=0)|(covar.arg$prob>=1)){stop(" invalid prob parameter for the binomial covariate")}
if(covar.arg$size<=0){stop(" invalid size parameter for the binomial covariate")}}
if(pmatch(c(covar.arg$types),c("BU","BB","UU","B","U"))==2){
if((length(covar.arg$prob)<1)|(length(covar.arg$size)<1)){stop(" provide size and prob parameter for the binomial covariate")}
if(length(covar.arg$size)<2){covar.arg$size<-c(covar.arg$size[1],covar.arg$size[1])}
if(length(covar.arg$prob)<2){covar.arg$prob<-c(covar.arg$prob[1],covar.arg$prob[1])}
if(any(covar.arg$prob<=0)|any(covar.arg$prob>=1)){stop(" invalid prob parameter for the binomial covariate")}
if(any(covar.arg$size<=0)){stop(" invalid size parameter for the binomial covariate")}}
if(pmatch(c(covar.arg$types),c("BU","BB","UU","B","U"))==3){
if((length(covar.arg$min)<1)|(length(covar.arg$max)<1)){stop(" provide min and max parameters for the uniform covariate")}
if(length(covar.arg$min)<2){covar.arg$min<-c(covar.arg$min[1],covar.arg$min[1])}
if(length(covar.arg$max)<2){covar.arg$max<-c(covar.arg$max[1],covar.arg$max[1])}
if(length(covar.arg$max)>2){covar.arg$max<-covar.arg$max[1:2]}
if(length(covar.arg$min)>2){covar.arg$min<-covar.arg$min[1:2]}
if(any(c(covar.arg$max-covar.arg$min)<=0)){stop(" invalid min and max parameters for uniform covariate")}}}
if(pmatch(c(covar.arg$types),c("BU","BB","UU","B","U"))>3){
if(length(covar.arg$coefs)<1){stop(" provide covariate coefficients")}
if(length(covar.arg$coefs)>=2){covar.arg$coefs<-covar.arg$coefs[1]}
if(pmatch(c(covar.arg$types),c("BU","BB","UU","B","U"))==4){
if((length(covar.arg$prob)<1)|(length(covar.arg$size)<1)){stop(" provide size and prob parameter for the binomial covariate")}
if(length(covar.arg$size)>1){covar.arg$size<-c(covar.arg$size[1])}
if(length(covar.arg$prob)>1){covar.arg$prob<-c(covar.arg$prob[1])}
if((covar.arg$prob<=0)|(covar.arg$prob>=1)){stop(" invalid prob parameter for the binomial covariate")}
if(covar.arg$size<=0){stop(" invalid size parameter for the binomial covariate")}}
if(pmatch(c(covar.arg$types),c("BU","BB","UU","B","U"))==5){
if((length(covar.arg$min)<1)|(length(covar.arg$max)<1)){stop(" provide min and max parameters for the uniform covariate")}
if(length(covar.arg$min)>1){covar.arg$min<-c(covar.arg$min[1])}
if(length(covar.arg$max)>1){covar.arg$max<-c(covar.arg$max[1])}
if(c(covar.arg$max-covar.arg$min)<=0){stop(" invalid min and max parameters for uniform covariate")}}
}
dataa<-gener.datadv(p.size,c.rate,fraildistrn,frail.par,bhaz.arg,covar.arg)
gendat <-list(numberofpair=p.size,censoredrate= c.rate,
fraildist=fraildistrn,frailpar=frail.par,basehazdis=bhaz.arg$distrn,
covartype=covar.arg$types,covarcoef=covar.arg$coefs)
gendat$data<-dataa$data
gendat$call <- match.call()
class(gendat) <- c("simbcfraildv",class(gendat))
gendat}
###

