
########################################################################
############Bivariate Correlated gamma fit functions###########
.expfrailfun1 <- function(x,k0,k1,k2){
q1<-q2<-q3<-q4<-p1<-p2<-p3<-p4<-su<-NULL
B0<-B1<-B2<-frai<-NULL
if((x[1]==0) & (x[2]==0)){
B0<-(k0/x[5])
B1<-(k1/x[3]);B2<-(k2/x[4])
frai<-c(B0,B1,B2)}
if((x[1]==1) & (x[2]==0)){
q1=(k0/x[5]); q2=(k1/x[3])
p1=(q1/(q1+q2)); p2=(q2/(q1+q2))
B0=((p1*(k0+1)+p2*k0)/x[5])
B2=(k2/x[4])
B1=((p1*k1+p2*(k1+1))/x[3])
frai<-c(B0,B1,B2)}
if((x[1]==0) & (x[2]==1)){
q1=(k0/x[5]); q2=(k2/x[4])
p1=(q1/(q1+q2)); p2=(q2/(q1+q2))
B0=((p1*(k0+1)+p2*k0)/x[5])
B1=(k1/x[3])
B2=((p1*k2+p2*(k2+1))/x[4])
frai<-c(B0,B1,B2)}
if((x[1]==1) & (x[2]==1)){
q1=((k0*(k0+1))/x[5]^2); q2=((k0*k1)/(x[5]*x[3]))
q3=((k0*k2)/(x[5]*x[4])); q4=((k1*k2)/(x[3]*x[4]))
su=(q1+q2+q3+q4)
p1=(q1/su); p2=(q2/su); p3=(q3/su); p4=(q4/su)
B0=((p1*(k0+2)+(p2+p3)*(k0+1)+p4*k0)/(x[5]))
B1=(((p1+p3)*k1+(p2+p4)*(k1+1))/(x[3]))
B2=(((p1+p2)*k2+(p3+p4)*(k2+1))/(x[4]))
frai<-c(B0,B1,B2)}
frai
}

############ Expected values of frailty variables####
.Expefrail1=function(newtht,cen1,cen2,H1,H2){
a1<-newtht[1];a2<-newtht[2];R<-newtht[3]
sa1<-sqrt(a1);sa2<-sqrt(a2)
k0<-(R/(sa1*sa2));k1<-(1-(sa1/sa2)*R)/a1;k2<-(1-(sa2/sa1)*R)/a2
x<-cbind(cen1,cen2,(1/a1+H1),(1/a2+H2),(1/(sa1*sa2)+(sa1/sa2)*H1+(sa2/sa1)*H2))
frailcomp <- apply(x,1,.expfrailfun1,k0=k0,k1=k1,k2=k2)
frailcomp =t(as.matrix(frailcomp))
z1=(sa1/sa2)*c(frailcomp[,1])+c(frailcomp[,2])
z2=(sa2/sa1)*c(frailcomp[,1])+c(frailcomp[,3])
list(frail1=log(z1),frail2=log(z2))
}
####

############ The mariginal Log likelihood function and its first derivative w.r.t. frailty parameters ####

.Llikgammadv = function (theta, cen1,cen2,H1,H2){
a1 = abs(theta[1]);a2 = abs(theta[2]);R = abs(theta[3])
logli=sum((-((1/a1)-(R/(a1^(1/2)*a2^(1/2))))*log(1+a1*H1)-((1/a2)-(R/(a1^(1/2)*a2^(1/2))))*log(1+a2*H2)-
(R/(a1^(1/2)*a2^(1/2)))*log(1+a1*H1+a2*H2))+
cen1*cen2*(-log(1+a1*H1)-log(1+a2*H2)-2*log(1+a1*H1+a2*H2)+
log(((R+(a1^(1/2)*a2^(1/2)))*R*(1+a1*H1)*(1+a2*H2)+(((a1^(1/2)/a2^(1/2))*R-R^2)*(1+a1*H1)+
((a2^(1/2)/a1^(1/2))*R-R^2)*(1+a2*H2))*(1+a1*H1+a2*H2)+
(1-(a1^(1/2)/a2^(1/2))*R)*(1-(a2^(1/2)/a1^(1/2))*R)*((1+a1*H1+a2*H2))^2)))+
cen1*(1-cen2)*(-log(1+a1*H1)-log(1+a1*H1+a2*H2)+
log(((1+a1*H1)+(1-(a1^(1/2)/a2^(1/2))*R)*(1+a2*H2)-(1-(a1^(1/2)/a2^(1/2))*R))))+
(1-cen1)*cen2*(-log(1+a2*H2)-log(1+a1*H1+a2*H2)+
log(((1+a2*H2)+(1-(a2^(1/2)/a1^(1/2))*R)*(1+a1*H1)-(1-(a2^(1/2)/a1^(1/2))*R)))))
-logli
}


.fdLlikgammadv =function (theta, cen1,cen2,H1,H2){
a1 = abs(theta[1]);a2 = abs(theta[2]);R = abs(theta[3])
da1=sum(((1/a1^2 - R * (a1^((1/2) - 1) * (1/2) * a2^(1/2))/(a1^(1/2) * a2^(1/2))^2) * log(1 + a1 * H1) - ((1/a1) - (R/(a1^(1/2) *
a2^(1/2)))) * (H1/(1 + a1 * H1)) - R * (a1^((1/2) - 1) * (1/2) * a2^(1/2))/(a1^(1/2) * a2^(1/2))^2 * log(1 + a2 *
H2) - ((R/(a1^(1/2) * a2^(1/2))) * (H1/(1 + a1 * H1 + a2 * H2)) - R * (a1^((1/2) - 1) * (1/2) * a2^(1/2))/(a1^(1/2) *
a2^(1/2))^2 * log(1 + a1 * H1 + a2 * H2)) + cen1 * cen2 * (((a1^((1/2) - 1) * (1/2) * a2^(1/2) * R * (1 + a1 * H1) +
(R + (a1^(1/2) * a2^(1/2))) * R * H1) * (1 + a2 * H2) + ((a1^((1/2) - 1) * (1/2)/a2^(1/2) * R * (1 + a1 * H1) +
((a1^(1/2)/a2^(1/2)) * R - R^2) * H1 - a2^(1/2) * (a1^((1/2) - 1) * (1/2))/(a1^(1/2))^2 * R * (1 +
a2 * H2)) * (1 + a1 * H1 + a2 * H2) + (((a1^(1/2)/a2^(1/2)) * R - R^2) * (1 + a1 * H1) + ((a2^(1/2)/a1^(1/2)) *
R - R^2) * (1 + a2 * H2)) * H1) + (((1 - (a1^(1/2)/a2^(1/2)) * R) * (a2^(1/2) * (a1^((1/2) - 1) * (1/2))/(a1^(1/2))^2 *
R) - a1^((1/2) - 1) * (1/2)/a2^(1/2) * R * (1 - (a2^(1/2)/a1^(1/2)) * R)) * ((1 + a1 * H1 + a2 * H2))^2 + (1 - (a1^(1/2)/a2^(1/2)) *
R) * (1 - (a2^(1/2)/a1^(1/2)) * R) * (2 * (H1 * ((1 + a1 * H1 + a2 * H2))))))/((R + (a1^(1/2) * a2^(1/2))) *
R * (1 + a1 * H1) * (1 + a2 * H2) + (((a1^(1/2)/a2^(1/2)) * R - R^2) * (1 + a1 * H1) + ((a2^(1/2)/a1^(1/2)) * R -
R^2) * (1 + a2 * H2)) * (1 + a1 * H1 + a2 * H2) + (1 - (a1^(1/2)/a2^(1/2)) * R) * (1 - (a2^(1/2)/a1^(1/2)) *
R) * ((1 + a1 * H1 + a2 * H2))^2) - (H1/(1 + a1 * H1) + 2 * (H1/(1 + a1 * H1 + a2 * H2)))) + cen1 * (1 - cen2) *
((H1 - a1^((1/2) - 1) * (1/2)/a2^(1/2) * R * (1 + a2 * H2) + a1^((1/2) - 1) * (1/2)/a2^(1/2) * R)/((1 + a1 * H1) +
(1 - (a1^(1/2)/a2^(1/2)) * R) * (1 + a2 * H2) - (1 - (a1^(1/2)/a2^(1/2)) * R)) - (H1/(1 + a1 * H1) + H1/(1 +
a1 * H1 + a2 * H2))) + (1 - cen1) * cen2 * ((a2^(1/2) * (a1^((1/2) - 1) * (1/2))/(a1^(1/2))^2 * R * (1 + a1 * H1) +
(1 - (a2^(1/2)/a1^(1/2)) * R) * H1 - a2^(1/2) * (a1^((1/2) - 1) * (1/2))/(a1^(1/2))^2 * R)/((1 + a2 * H2) + (1 - (a2^(1/2)/a1^(1/2)) *
R) * (1 + a1 * H1) - (1 - (a2^(1/2)/a1^(1/2)) * R)) - H1/(1 + a1 * H1 + a2 * H2))))
da2=sum((cen1 * cen2 * ((a1^(1/2) * (a2^((1/2) - 1) * (1/2)) * R * (1 + a1 * H1) * (1 + a2 * H2) + (R + (a1^(1/2) * a2^(1/2))) *
R * (1 + a1 * H1) * H2 + ((a2^((1/2) - 1) * (1/2)/a1^(1/2) * R * (1 + a2 * H2) + ((a2^(1/2)/a1^(1/2)) * R - R^2) * H2 -
a1^(1/2) * (a2^((1/2) - 1) * (1/2))/(a2^(1/2))^2 * R * (1 + a1 * H1)) * (1 + a1 * H1 + a2 * H2) + (((a1^(1/2)/a2^(1/2)) *
R - R^2) * (1 + a1 * H1) + ((a2^(1/2)/a1^(1/2)) * R - R^2) * (1 + a2 * H2)) * H2) + ((a1^(1/2) * (a2^((1/2) - 1) * (1/2))/(a2^(1/2))^2 *
R * (1 - (a2^(1/2)/a1^(1/2)) * R) - (1 - (a1^(1/2)/a2^(1/2)) * R) * (a2^((1/2) - 1) * (1/2)/a1^(1/2) * R)) * ((1 + a1 *
H1 + a2 * H2))^2 + (1 - (a1^(1/2)/a2^(1/2)) * R) * (1 - (a2^(1/2)/a1^(1/2)) * R) * (2 * (H2 * ((1 + a1 * H1 + a2 * H2))))))/((R + (a1^(1/2) *
a2^(1/2))) * R * (1 + a1 * H1) * (1 + a2 * H2) + (((a1^(1/2)/a2^(1/2)) * R - R^2) * (1 + a1 * H1) + ((a2^(1/2)/a1^(1/2)) * R - R^2) *
(1 + a2 * H2)) * (1 + a1 * H1 + a2 * H2) + (1 - (a1^(1/2)/a2^(1/2)) * R) * (1 - (a2^(1/2)/a1^(1/2)) * R) * ((1 + a1 * H1 + a2 *
H2))^2) - (H2/(1 + a2 * H2) + 2 * (H2/(1 + a1 * H1 + a2 * H2)))) - (R * (a1^(1/2) * (a2^((1/2) - 1) * (1/2)))/(a1^(1/2) *
a2^(1/2))^2 * log(1 + a1 * H1) + (((1/a2) - (R/(a1^(1/2) * a2^(1/2)))) * (H2/(1 + a2 * H2)) - (1/a2^2 - R * (a1^(1/2) *
(a2^((1/2) - 1) * (1/2)))/(a1^(1/2) * a2^(1/2))^2) * log(1 + a2 * H2)) + ((R/(a1^(1/2) * a2^(1/2))) * (H2/(1 + a1 * H1 +
a2 * H2)) - R * (a1^(1/2) * (a2^((1/2) - 1) * (1/2)))/(a1^(1/2) * a2^(1/2))^2 * log(1 + a1 * H1 + a2 * H2))) + cen1 * (1 -
cen2) * ((a1^(1/2) * (a2^((1/2) - 1) * (1/2))/(a2^(1/2))^2 * R * (1 + a2 * H2) + (1 - (a1^(1/2)/a2^(1/2)) * R) * H2 -
a1^(1/2) * (a2^((1/2) - 1) * (1/2))/(a2^(1/2))^2 * R)/((1 + a1 * H1) + (1 - (a1^(1/2)/a2^(1/2)) * R) * (1 + a2 * H2) -
(1 - (a1^(1/2)/a2^(1/2)) * R)) - H2/(1 + a1 * H1 + a2 * H2)) + (1 - cen1) * cen2 * ((H2 - a2^((1/2) - 1) * (1/2)/a1^(1/2) *
R * (1 + a1 * H1) + a2^((1/2) - 1) * (1/2)/a1^(1/2) * R)/((1 + a2 * H2) + (1 - (a2^(1/2)/a1^(1/2)) * R) * (1 +
a1 * H1) - (1 - (a2^(1/2)/a1^(1/2)) * R)) - (H2/(1 + a2 * H2) + H2/(1 + a1 * H1 + a2 * H2)))))
dR=sum((1/(a1^(1/2) * a2^(1/2)) * log(1 + a1 * H1) + 1/(a1^(1/2) * a2^(1/2)) * log(1 + a2 * H2) - 1/(a1^(1/2) * a2^(1/2)) * log(1 + a1 *
H1 + a2 * H2) + cen1 * cen2 * (((R + (R + (a1^(1/2) * a2^(1/2)))) * (1 + a1 * H1) * (1 + a2 * H2) + (((a1^(1/2)/a2^(1/2)) - 2 *
R) * (1 + a1 * H1) + ((a2^(1/2)/a1^(1/2)) - 2 * R) * (1 + a2 * H2)) * (1 + a1 * H1 + a2 * H2) - ((1 - (a1^(1/2)/a2^(1/2)) *
R) * (a2^(1/2)/a1^(1/2)) + (a1^(1/2)/a2^(1/2)) * (1 - (a2^(1/2)/a1^(1/2)) * R)) * ((1 + a1 * H1 + a2 * H2))^2)/((R + (a1^(1/2) * a2^(1/2))) *
R * (1 + a1 * H1) * (1 + a2 * H2) + (((a1^(1/2)/a2^(1/2)) * R - R^2) * (1 + a1 * H1) + ((a2^(1/2)/a1^(1/2)) * R - R^2) *
(1 + a2 * H2)) * (1 + a1 * H1 + a2 * H2) + (1 - (a1^(1/2)/a2^(1/2)) *R) * (1 - (a2^(1/2)/a1^(1/2)) * R) * ((1 + a1 * H1 + a2 *
H2))^2)) - cen1 * (1 - cen2) * (((a1^(1/2)/a2^(1/2)) * (1 + a2 * H2) - (a1^(1/2)/a2^(1/2)))/((1 + a1 * H1) + (1 - (a1^(1/2)/a2^(1/2)) *
R) * (1 + a2 * H2) - (1 - (a1^(1/2)/a2^(1/2)) * R))) - (1 - cen1) * cen2 * (((a2^(1/2)/a1^(1/2)) * (1 + a1 * H1) - (a2^(1/2)/a1^(1/2)))/((1 +
a2 * H2) + (1 - (a2^(1/2)/a1^(1/2)) * R) * (1 + a1 * H1) - (1 - (a2^(1/2)/a1^(1/2)) * R)))))
gd=c(da1,da2,dR)
-gd
}

####

### parameter space restriction####
Uperfun<-function(theta){
rs=c(Inf,Inf,min(c((sqrt(theta[1])/sqrt(theta[2])),(sqrt(theta[2])/sqrt(theta[1])))))
return(rs)}
####

.bcgamfitdv<-function(X,Y,initfrailp,control){
time=Y[, 2];censor=Y[, 3]
order=sort(time, decreasing = FALSE, index.return = TRUE);indx=order$ix;timeo=order$x
uniq_tim<-unique(sort(time))
ind.haz=match(time,uniq_tim)
if(any(is.na(ind.haz))){
tord_diff<-as.array(diff(c(0,timeo)))
id.zero_tord <- which(apply(((tord_diff<0.0000001)&(tord_diff>0)),1, all))
if(length(id.zero_tord)>0){time[indx[id.zero_tord]]<- time[indx[id.zero_tord-1]]
order=sort(time, decreasing = FALSE, index.return = TRUE);indx=order$ix;timeo=order$x}
uniq_tim<-unique(sort(time))
ind.haz=match(time,uniq_tim)}
data.n1 <- length(time);data.n <-data.n1/2#### data.n is the number of pairs
risskset <- function(uniq_tim,x) ifelse(uniq_tim[1]<=x,1,0)
RI <- apply(as.array(uniq_tim),1,risskset,x=time)
RI <-t(RI)
numbevent <- function(uniq_tim,time,censor){
indic0=match(time,uniq_tim[1])*array(1:data.n1)
sum(censor[indic0[!is.na(indic0)]])}
n_eve<-apply(as.array(uniq_tim),1,numbevent,time,censor)
indic1<-2*array(1:data.n)-1;indic2<-2*array(1:data.n)
cph0 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset = NULL, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
bet<-cph0$coefficients
x_bet<-X%*%bet
svexp_bet_xo=as.vector(RI%*%(exp(x_bet)))
H0<-cumsum(c(n_eve/svexp_bet_xo))###obtain initial parameter for H0(t)
H_bet_x=c(H0[ind.haz]*exp(c(x_bet)))
cen1=censor[indic1];cen2=censor[indic2]
if(length(initfrailp)==0){newtht<-c(0.5,0.5,0.5)}
if(length(initfrailp)>0){newtht<-initfrailp}
estim0=c(newtht)
W<-NULL;new.diff=newdiff=1
upb=itero=newthto=0
e_frail<-.Expefrail1(newtht=newtht,cen1=cen1,cen2=cen2,
H1=H_bet_x[indic1],H2=H_bet_x[indic2])
W[indic1]<-e_frail$frail1;W[indic2]<-e_frail$frail2
if(any(is.na(W))){W[is.na(W)]<-mean(W[!is.na(W)])}
cph1 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset =W, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
bet<-cph1$coefficients
x_bet<-X%*%bet
svexp_bet_xo=as.vector(RI%*%(exp(x_bet+W)))
H0<-cumsum(c(n_eve/svexp_bet_xo))###obtain initial parameter for H0(t)
H_bet_x=c(H0[ind.haz]*exp(c(x_bet)))
iter=0
repeat{
iter=iter+1
e_frail<-.Expefrail1(newtht=newtht,cen1=cen1,cen2=cen2,
H1=H_bet_x[indic1],H2=H_bet_x[indic2])
W[indic1]<-e_frail$frail1;W[indic2]<-e_frail$frail2
if(any(is.na(W))){W[is.na(W)]<-mean(W[!is.na(W)])}
cph1 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset =W, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
bet<-cph1$coefficients
x_bet<-X%*%bet
svexp_bet_xo=as.vector(RI%*%(exp(x_bet+W)))
H0<-cumsum(c(n_eve/svexp_bet_xo))###obtain initial parameter for H0(t)
H_bet_x=c(H0[ind.haz]*exp(c(x_bet)))
fittr=do.call(nlminb, args=c(list(start=newtht, objective=.Llikgammadv,
gradient = .fdLlikgammadv,cen1=cen1,cen2=cen2,H1=H_bet_x[indic1],
H2=H_bet_x[indic2],lower = control$lower, upper = Uperfun(theta=newtht),
control=control$nlminb_control)))
newtht=fittr$par##### obtain new estimates of sigma^2 and row
upb=Uperfun(theta=newtht)
if(newtht[3]>upb[3]){
itero=0
repeat{
itero=itero+1
newthto=newtht
newthto[3]=abs(upb[3]-0.00000000001)
fittr=do.call(nlminb, args=c(list(start=newthto, objective=.Llikgammadv,
gradient = .fdLlikgammadv,cen1=cen1,cen2=cen2,H1=H_bet_x[indic1],
H2=H_bet_x[indic2],lower = control$lower, upper = Uperfun(theta=newtht),
control=control$nlminb_control)))
newtht=fittr$par
upb=Uperfun(theta=newtht)
if((itero>= 20) | (newtht[3]<upb[3])) break}}
new.diff=max(abs(abs(c(newtht))-abs(estim0)))
estim0=c(newtht)
if((new.diff < control$tol)  |  (iter >= control$max.iter)) break}
if (iter > control$max.iter){warning("Ran out of iterations and did not converge")}
h0=diff(c(0,H0));nonzero_h0=h0[h0>0]
H=H0[ind.haz];h=h0[ind.haz];exp_bet_x=exp(x_bet)
lik=(sum(censor*(c(x_bet)))+sum(log(nonzero_h0))-fittr$objective+sum(censor))
vcov=cph1$var
colnames(vcov) <- rownames(vcov) <- colnames(X)
if(control$obt.se){
adj_se=try(.SEbcfrailgamdv(bet=bet,newtht=newtht,n_eve=n_eve,etime=uniq_tim,h0=h0,
censor=censor,time=time,X=X,H=c(H_bet_x)),silent=TRUE)
vcovth<-solve(adj_se$vco)
colnames(vcovth) <- rownames(vcovth) <- c("theta1","theta2","Row")
if(length(attr(adj_se,"class"))==0){adjse=c(adj_se$se)
if(any(is.na(adjse))){adjse=c(sqrt(diag(vcov)),sqrt(diag(vcovth)))
warning("Check the approperiatness of the model.")
warning("frailty parameters might be at the boundary of parameter space.")
warning("In obtaining SE Na is produced")}}
if(length(attr(adj_se,"class"))!=0){
adjse=c(sqrt(diag(vcov)),sqrt(diag(vcovth)))
warning("Check the approperiatness of the model.")
warning("frailty parameters might be at the boundary of parameter space.")
warning("In obtaining SE Na is produced")}}
if(!control$obt.se){adjse=c(sqrt(diag(vcov)),sqrt(diag(vcovth)))}
res <-list(coefficients=bet,frailparest= c(theta1=newtht[1],theta2=newtht[2],Row=newtht[3]),
vcov = vcov,vcovth= vcovth,stderr=adjse,loglilk0=cph0$loglik[1],loglilk=cph1$loglik[2],
Iloglilk=lik,cbasehaz=cbind(haz=H0,time=uniq_tim),X=X,time=time,censor=censor,
resid=cph1$residuals,lin.prid=cph1$linear.predictors,
frail=exp(W),iteration=iter,inh=ind.haz,indx=indx,e.time=uniq_tim,n.event=n_eve,
converg = ifelse(new.diff<control$tol,TRUE,FALSE))
res$call <- match.call()
class(res) <- c(".bcgamfitdv")
res
}


############ End Bivariate Correlated gamma fit functions###
########################################################################


########################################################################
############Bivariate Correlated lognormal fit functions###

.hesfunc1 = function(par,bet,censoro,Xo,RI,INC){
n=nrow(Xo)
i1=INC[,1];i2=INC[,2];i3=INC[,3]
W=par
W1<-W[i1]
W2<-W[i2]
x_bet<-c(Xo%*%bet+W)
expx_bet<-exp(c(Xo%*%bet+W))
svexp_bet_xo=as.vector(RI%*%(exp(x_bet)))
as0=cumsum(c(censoro/svexp_bet_xo))
as1=((cumsum(c(censoro/(svexp_bet_xo^2))))*c(expx_bet^2)-(as0*c(expx_bet)))
as2=cumsum(c(censoro/svexp_bet_xo^2))
as3=as2[i3]*(expx_bet[i1])*(expx_bet[i2])
xx<-matrix(c(W1,W2,as1[i1],as1[i2],as3),(n/2),5)
xx
}

.llpenn1= function (par,newtht,censoro,sumx,Xo,RI,INC){
n4=ncol(Xo)
n=nrow(Xo)
i1=INC[,1];i2=INC[,2];i3=INC[,3]
a1=newtht[1];a2=newtht[2];roh=newtht[3]
sa1=sqrt(a1);sa2=sqrt(a2)
bet=par[array(1:n4)]
W=par[(1+n4):(n+n4)]
x_bet<-c(Xo%*%bet+W)
svexp_bet_xo=as.vector(RI%*%(exp(x_bet)))
.llpenn1=(-sum(censoro*(x_bet-log(svexp_bet_xo)))+
(1/(2*a1*a2*(1-roh^2)))*sum(a2*W[i1]^2+a1*W[i2]^2-2*roh*sa1*sa2*W[i1]*W[i2]))
Dw<-(-censoro+exp(c(x_bet))*cumsum(censoro/svexp_bet_xo))
Dw[i1]<-Dw[i1]+(a2*W[i1]-roh*sa1*sa2*W[i2])*(1/(a1*a2*(1-roh^2)))
Dw[i2]<-Dw[i2]+(a1*W[i2]-roh*sa1*sa2*W[i1])*(1/(a1*a2*(1-roh^2)))
attr(.llpenn1,"gradient")<-c(as.vector(-sumx+colSums( c(censoro/svexp_bet_xo)*(RI%*%(Xo*c(exp(x_bet)))))),c(Dw))
.llpenn1
}

.Lliklogndv = function (thetaa, xx){
xx=as.matrix(xx)
W1=xx[,1];W2=xx[,2];S1=xx[,3];S2=xx[,4];S3=xx[,5]
a1=abs(thetaa[1]);a2=abs(thetaa[2]);e=thetaa[3]
sa1=sqrt(a1);sa2=sqrt(a2)
n = nrow(xx)
P11<-(S1*S2-(1/(a2*(1-e^2)))*S1-(1/(a1*(1-e^2)))*S2+
(1/(a1*a2*(1-e^2)^2))-S3^2-2*(e/(sa1*sa2*(1-e^2)))*S3-(e^2/(a1*a2*(1-e^2)^2)))
logli=(-(1/(2*a1*a2*(1-e^2)))*sum(a2*W1^2+a1*W2^2-2*e*sa1*sa2*W1*W2)-
(n/2)*log(a1)-(n/2)*log(a2)-(n/2)*log(1-e^2)-(1/2)*sum(log(P11)))
-logli
}


.fdLliklogndv = function (thetaa, xx){
xx=as.matrix(xx)
W1=xx[,1];W2=xx[,2];S1=xx[,3];S2=xx[,4];S3=xx[,5]
a1=thetaa[1];a2=thetaa[2];e=thetaa[3]
sa1=sqrt(a1);sa2=sqrt(a2)
n = nrow(xx)
P11<-(S1*S2-(1/(a2*(1-e^2)))*S1-(1/(a1*(1-e^2)))*S2+(1/(a1*a2*(1-e^2)^2))-
S3^2-2*(e/(sa1*sa2*(1-e^2)))*S3-(e^2/(a1*a2*(1-e^2)^2)))
da1P11<-((1/(a1^2*(1-e^2)))*S2-(1/(a1^2*a2*((1-e^2)^2)))+
((e/(a1^(3/2)*sa2*(1-e^2))))*S3+(e^2/(a1^2*a2*((1-e^2)^2))))
da2P11<-((1/(a2^2*(1-e^2)))*S1-(1/(a2^2*a1*((1-e^2)^2)))+
((e/(a2^(3/2)*sa1*(1-e^2))))*S3+(e^2/(a2^2*a1*((1-e^2)^2))))
deP11<-(-((2*e*a2)/(a2*(1-e^2))^2)*S1-((2*e*a1)/(a1*(1-e^2))^2)*S2+
(4*(a1*a2*(1-e^2)*e)/(a1*a2*(1-e^2)^2)^2)-(2/(sa1*sa2*(1-e^2)))*S3-
((4*e^2*sa1*sa2)/(sa1*sa2*(1-e^2))^2)*S3-
(2*e/(a1*a2*(1-e^2)^2)+(e^3*4*a1*a2*(1-e^2))/(a1*a2*(1-e^2)^2)^2))
grr= c(((1/(2*a1^2*a2*(1-e^2)))*sum(a2*W1^2+a1*W2^2-2*e*sa1*sa2*W1*W2)-
(1/(2*a1*a2*(1-e^2)))*sum(W2^2-e*(sa2/sa1)*W1*W2)-(n/(2*a1))-(1/2)*sum(da1P11/P11)),
((1/(2*a2^2*a1*(1-e^2)))*sum(a2*W1^2+a1*W2^2-2*e*sa1*sa2*W1*W2)-
(1/(2*a1*a2*(1-e^2)))*sum(W1^2-e*(sa1/sa2)*W1*W2)-(n/(2*a2))-(1/2)*sum(da2P11/P11)),
( -((a1*a2*4*e)/(2*a1*a2*(1-e^2))^2)*sum(a2*W1^2+a1*W2^2-2*e*sa1*sa2*W1*W2)+
(1/(a1*a2*(1-e^2)))*sum(sa1*sa2*W1*W2)+n*(e/(1-e^2))-(1/2)*sum(deP11/P11)))
-grr
}


.sdLliklogndv= function (thetaa, xx){
xx=as.matrix(xx)
W1=xx[,1];W2=xx[,2];S1=xx[,3];S2=xx[,4];S3=xx[,5]
a1=thetaa[1];a2=thetaa[2];e=thetaa[3]
sa1=sqrt(a1);sa2=sqrt(a2)
n = nrow(xx)
P11<-(S1*S2-(1/(a2*(1-e^2)))*S1-(1/(a1*(1-e^2)))*S2+(1/(a1*a2*(1-e^2)^2))-
S3^2-2*(e/(sa1*sa2*(1-e^2)))*S3-(e^2/(a1*a2*(1-e^2)^2)))
da1P11<-((1/(a1^2*(1-e^2)))*S2-(1/(a1^2*a2*((1-e^2)^2)))+
((e/(a1^(3/2)*sa2*(1-e^2))))*S3+(e^2/(a1^2*a2*((1-e^2)^2))))
da2P11<-((1/(a2^2*(1-e^2)))*S1-(1/(a2^2*a1*((1-e^2)^2)))+
((e/(a2^(3/2)*sa1*(1-e^2))))*S3+(e^2/(a2^2*a1*((1-e^2)^2))))
deP11<-(-((2*e*a2)/(a2*(1-e^2))^2)*S1-((2*e*a1)/(a1*(1-e^2))^2)*S2+
(4*(a1*a2*(1-e^2)*e)/(a1*a2*(1-e^2)^2)^2)-(2/(sa1*sa2*(1-e^2)))*S3-
((4*e^2*sa1*sa2)/(sa1*sa2*(1-e^2))^2)*S3-
(2*e/(a1*a2*(1-e^2)^2)+(e^3*4*a1*a2*(1-e^2))/(a1*a2*(1-e^2)^2)^2))
d2a1P11<-( (-(2 * a1 * (1 - e^2)/(a1^2 * (1 - e^2))^2))*S2-
(-(2 * a1 * a2 * ((1 - e^2)^2)/(a1^2 * a2 * ((1 - e^2)^2))^2))+
(-(e * (a1^((3/2) - 1) * (3/2) * a2^(1/2) * (1 - e^2))/(a1^(3/2) *
    a2^(1/2) * (1 - e^2))^2))*S3+
(-(e^2 * (2 * a1 * a2 * ((1 - e^2)^2))/(a1^2 * a2 * ((1 - e^2)^2))^2)))
da1a2P11<-( -(-(a1^2 * ((1 - e^2)^2)/(a1^2 * a2 * ((1 - e^2)^2))^2))+
(-(e * (a1^(3/2) * (a2^((1/2) - 1) * (1/2)) * (1 - e^2))/(a1^(3/2) *
    a2^(1/2) * (1 - e^2))^2))*S3+
(-(e^2 * (a1^2 * ((1 - e^2)^2))/(a1^2 * a2 * ((1 - e^2)^2))^2)))
da1eP11<-( (a1^2 * (2 * e)/(a1^2 * (1 - e^2))^2)*S2-
(a1^2 * a2 * (2 * (2 * e * (1 - e^2)))/(a1^2 * a2 * ((1 - e^2)^2))^2)+
(1/(a1^(3/2) * a2^(1/2) * (1 - e^2)) + e * (a1^(3/2) * a2^(1/2) *
    (2 * e))/(a1^(3/2) * a2^(1/2) * (1 - e^2))^2)*S3+
(2 * e/(a1^2 * a2 * ((1 - e^2)^2)) + e^2 * (a1^2 * a2 * (2 * (2 *
    e * (1 - e^2))))/(a1^2 * a2 * ((1 - e^2)^2))^2))
d2a2P11<-((-(2 * a2 * (1 - e^2)/(a2^2 * (1 - e^2))^2))*S1-
(-(2 * a2 * a1 * ((1 - e^2)^2)/(a2^2 * a1 * ((1 - e^2)^2))^2))+
(-(e * (a2^((3/2) - 1) * (3/2) * a1^(1/2) * (1 - e^2))/(a2^(3/2) *
    a1^(1/2) * (1 - e^2))^2))*S3+
(-(e^2 * (2 * a2 * a1 * ((1 - e^2)^2))/(a2^2 * a1 * ((1 - e^2)^2))^2)))
da2eP11<-((a2^2 * (2 * e)/(a2^2 * (1 - e^2))^2)*S1-
(a2^2 * a1 * (2 * (2 * e * (1 - e^2)))/(a2^2 * a1 * ((1 - e^2)^2))^2)+
(1/(a2^(3/2) * a1^(1/2) * (1 - e^2)) + e * (a2^(3/2) * a1^(1/2) *
    (2 * e))/(a2^(3/2) * a1^(1/2) * (1 - e^2))^2)*S3+
(2 * e/(a2^2 * a1 * ((1 - e^2)^2)) + e^2 * (a2^2 * a1 * (2 * (2 *
    e * (1 - e^2))))/(a2^2 * a1 * ((1 - e^2)^2))^2))
d2eP11<-(-(2 * a2/(a2 * (1 - e^2))^2 + (2 * e * a2) * (2 * (a2 * (2 * e) *
    (a2 * (1 - e^2))))/((a2 * (1 - e^2))^2)^2)*S1-
(2 * a1/(a1 * (1 - e^2))^2 + (2 * e * a1) * (2 * (a1 * (2 * e) *
    (a1 * (1 - e^2))))/((a1 * (1 - e^2))^2)^2)*S2+
(4 * (a1 * a2 * (1 - e^2) - a1 * a2 * (2 * e) * e)/(a1 * a2 *
    (1 - e^2)^2)^2 + 4 * (a1 * a2 * (1 - e^2) * e) * (2 * (a1 *
    a2 * (2 * (2 * e * (1 - e^2))) * (a1 * a2 * (1 - e^2)^2)))/((a1 *
    a2 * (1 - e^2)^2)^2)^2)-(2 * (sa1 * sa2 * (2 * e))/(sa1 * sa2 * (1 - e^2))^2)*S3-
(4 * (2 * e) * sa1 * sa2/(sa1 * sa2 * (1 - e^2))^2 + (4 * e^2 *
    sa1 * sa2) * (2 * (sa1 * sa2 * (2 * e) * (sa1 * sa2 * (1 -
    e^2))))/((sa1 * sa2 * (1 - e^2))^2)^2)*S3-
(2/(a1 * a2 * (1 - e^2)^2) + 2 * e * (a1 * a2 * (2 * (2 * e *
    (1 - e^2))))/(a1 * a2 * (1 - e^2)^2)^2 + ((3 * e^2 * 4 *
    a1 * a2 * (1 - e^2) - e^3 * 4 * a1 * a2 * (2 * e))/(a1 *
    a2 * (1 - e^2)^2)^2 + (e^3 * 4 * a1 * a2 * (1 - e^2)) * (2 *
    (a1 * a2 * (2 * (2 * e * (1 - e^2))) * (a1 * a2 * (1 - e^2)^2)))/((a1 *
    a2 * (1 - e^2)^2)^2)^2)))
D2a1=(((1/(2 * a1^2 * a2 * (1 - e^2))) * sum(W2^2 - 2 * e * (a1^((1/2) -
    1) * (1/2)) * a2^(1/2) * W1 * W2) - 2 * (2 * a1) * a2 * (1 -
    e^2)/(2 * a1^2 * a2 * (1 - e^2))^2 * sum(a2 * W1^2 + a1 * W2^2 -
    2 * e * a1^(1/2) * a2^(1/2) * W1 * W2))-
((1/(2 * a1 * a2 * (1 - e^2))) * (e * (a2^(1/2) * (a1^((1/2) -
    1) * (1/2))/(a1^(1/2))^2) * sum(W1 * W2)) - 2 * a2 * (1 - e^2)/(2 *
    a1 * a2 * (1 - e^2))^2 * sum(W2^2 - e * (a2^(1/2)/a1^(1/2)) *
    W1 * W2))+n/(2*a1^2)-(1/2)*sum(d2a1P11/P11-da1P11^2/P11^2))
D2a2=(((1/(2 * a2^2 * a1 * (1 - e^2))) * sum(W1^2 - 2 * e * a1^(1/2) *
    (a2^((1/2) - 1) * (1/2)) * W1 * W2) - 2 * (2 * a2) * a1 *
    (1 - e^2)/(2 * a2^2 * a1 * (1 - e^2))^2 * sum(a2 * W1^2 + a1 *
    W2^2 - 2 * e * a1^(1/2) * a2^(1/2) * W1 * W2))-
((1/(2 * a1 * a2 * (1 - e^2))) * (e * (a1^(1/2) * (a2^((1/2) -
    1) * (1/2))/(a2^(1/2))^2) * sum(W1 * W2)) - 2 * a1 * (1 - e^2)/(2 *
    a1 * a2 * (1 - e^2))^2 * sum(W1^2 - e * (a1^(1/2)/a2^(1/2)) *
    W1 * W2))+n/(2*a2^2)-(1/2)*sum(d2a2P11/P11-da2P11^2/P11^2))
D2a1a2=(((1/(2 * a1^2 * a2 * (1 - e^2))) * sum(W1^2 - 2 * e * a1^(1/2) *
    (a2^((1/2) - 1) * (1/2)) * W1 * W2) - 2 * a1^2 * (1 - e^2)/(2 *
    a1^2 * a2 * (1 - e^2))^2 * sum(a2 * W1^2 + a1 * W2^2 - 2 * e *
    a1^(1/2) * a2^(1/2) * W1 * W2))-
(-((1/(2 * a1 * a2 * (1 - e^2))) * (e * (a2^((1/2) - 1) * (1/2)/a1^(1/2)) *
    sum(W1 * W2)) + 2 * a1 * (1 - e^2)/(2 * a1 * a2 * (1 - e^2))^2 *
    sum(W2^2 - e * (a2^(1/2)/a1^(1/2)) * W1 * W2)))-
(1/2)*sum(da1a2P11/P11-(da2P11*da1P11)/P11^2))
Da1R=((2 *a1^2*a2*(2*e)/(2* a1^2 *a2 *(1-e^2))^2*sum(a2 *W1^2+a1*W2^2-2*e*a1^(1/2)*a2^(1/2)*W1*W2)-
(1/(2 * a1^2 * a2 * (1 - e^2))) * sum(2 * a1^(1/2) * a2^(1/2) * W1 * W2))-
(2 * a1 * a2 * (2 * e)/(2 * a1 * a2 * (1 - e^2))^2 * sum(W2^2 - e *
    (a2^(1/2)/a1^(1/2)) * W1 * W2) - (1/(2 * a1 * a2 * (1 - e^2))) *
    sum((a2^(1/2)/a1^(1/2)) * W1 * W2))-
(1/2)*sum(da1eP11/P11-(da1P11*deP11)/P11^2))
Da2R=((2 * a2^2 * a1 * (2 * e)/(2 * a2^2 * a1 * (1 - e^2))^2 * sum(a2 *
    W1^2 + a1 * W2^2 - 2 * e * a1^(1/2) * a2^(1/2) * W1 * W2) -
    (1/(2 * a2^2 * a1 * (1 - e^2))) * (2 * a1^(1/2) * a2^(1/2) *
       sum(W1 * W2)))-
(2 * a1 * a2 * (2 * e)/(2 * a1 * a2 * (1 - e^2))^2 * sum(W1^2 - e *
    (a1^(1/2)/a2^(1/2)) * W1 * W2) - (1/(2 * a1 * a2 * (1 - e^2))) *
    ((a1^(1/2)/a2^(1/2)) * sum(W1 * W2)))-
(1/2)*sum(da2eP11/P11-(da2P11*deP11)/P11^2))
D2R=((((a1 * a2 * 4 * e)/(2 * a1 * a2 * (1 - e^2))^2) * (2 * a1^(1/2) *
    a2^(1/2) * sum(W1 * W2)) - (a1 * a2 * 4/(2 * a1 * a2 * (1 - e^2))^2 +
    (a1 * a2 * 4 * e) * (2 * (2 * a1 * a2 * (2 * e) * (2 * a1 *
        a2 * (1 - e^2))))/((2 * a1 * a2 * (1 - e^2))^2)^2) *
    sum(a2 * W1^2 + a1 * W2^2 - 2 * e * a1^(1/2) * a2^(1/2) * W1 * W2))+
(a1 * a2 * (2 * e)/(a1 * a2 * (1 - e^2))^2 * a1^(1/2) * a2^(1/2) * sum(W1 * W2))+
 n*( 1/(1-e^2)+2*e^2/(1-e^2)^2)-(1/2)*sum(d2eP11/P11-deP11^2/P11^2))
hesi=matrix(c(D2a1,D2a1a2,Da1R,D2a1a2,D2a2,Da2R,Da1R,Da2R,D2R),3,3)
-hesi
}

.SEpenal1= function(par,newtht,censoro,Xo,RI,INC){
n4=ncol(Xo)
n=nrow(Xo)
i1=INC[,1];i2=INC[,2];i3=INC[,3]
a1=newtht[1];a2=newtht[2];roh=newtht[3]
sa1=sqrt(a1);sa2=sqrt(a2)
bet=par[1:n4]
W=par[(1+n4):(n+n4)]
x_bet<-c(Xo%*%bet+W)
vexpx_bet<-c(exp(x_bet))
svexp_bet_xo=as.vector(RI%*%(exp(x_bet)))
dbexp_bet_xo=(RI%*%(vexpx_bet*Xo))
interacmat <- function(Xo,u){Xo*u}
resinteracmat<-apply(as.matrix(Xo),2,interacmat,u=Xo)
Xoxo<-matrix(resinteracmat,n,n4^2)
as0=cumsum(c(censoro/svexp_bet_xo))*vexpx_bet
ma0<-matrix(0,n,n)
ma0[lower.tri(ma0, diag = TRUE)]<-1
ma1<-t((ma0*vexpx_bet))
ma2<-ma1*(-cumsum(c(censoro/(svexp_bet_xo^2))))*vexpx_bet
digg<-diag(ma2)+as0
ma3<-t(ma2)
D2W<-ma3+ma2
digg[i1]<-digg[i1]+(1/(a1*(1-roh^2)))
digg[i2]<-digg[i2]+(1/(a2*(1-roh^2)))
diag(D2W)<-digg
D2W[cbind(i1,i2)]<-D2W[cbind(i1,i2)]-((roh*sa1*sa2)/(a1*a2*(1-roh^2)))
D2W[cbind(i2,i1)]<-D2W[cbind(i2,i1)]-((roh*sa1*sa2)/(a1*a2*(1-roh^2)))
dwbpart20<-(dbexp_bet_xo)*(censoro/(svexp_bet_xo^2))
cumsummat <- function(dwbpart20){cumsum(dwbpart20)}
dwbpart21<-apply(as.matrix(dwbpart20),2,cumsummat)
DWB<-(as0*Xo-dwbpart21*(vexpx_bet))
d2bpart1=matrix(c(colSums((censoro/svexp_bet_xo)*(RI%*%(vexpx_bet*Xoxo)))),n4,n4)
d2bpart2=t((censoro/svexp_bet_xo)*dbexp_bet_xo)%*%(dbexp_bet_xo*(censoro/svexp_bet_xo))
D2B<-d2bpart1-d2bpart2
MA=rbind(D2B,DWB);MB=rbind(t(DWB),D2W)
HES=cbind(MA,MB)
INVE<-solve(HES)
se=sqrt(diag(INVE))
se[1:n4]
}

.llppl1= function (par,newtht,censoro,Xo,RI){
n4=ncol(Xo)
n=nrow(Xo)
a1=newtht[1];a2=newtht[2];roh=newtht[3]
sa1=sqrt(a1);sa2=sqrt(a2)
bet=par[array(1:n4)]
W=par[(1+n4):(n+n4)]
x_bet<-c(Xo%*%bet+W)
svexp_bet_xo=as.vector(RI%*%(exp(x_bet)))
logl=sum(censoro*(x_bet-log(svexp_bet_xo)))
logl
}

.machinelern = function (newtht, newtht0,direct,lo,G){
curdirect=estimdiff=c(newtht-newtht0)
curdirect[curdirect<0]<-c(-1);curdirect[curdirect==0]<-0;curdirect[curdirect>0]<-1
prevdirect=direct
newtht0<-newtht
if(G==0){
for(j in 1:length(newtht)){
if(prevdirect[j]==1){
if(curdirect[j]==1){
lo[j]=lo[j]+2;newtht0[j]<-newtht[j]+lo[j]*estimdiff[j]}
if(curdirect[j]!=1){lo[j]=0.4*lo[j];newtht0[j]<-newtht[j]+lo[j]*estimdiff[j]}}
if(prevdirect[j]==-1){
if(curdirect[j]==-1){
lo[j]=lo[j]+2;newtht0[j]<-newtht[j]+lo[j]*estimdiff[j]}
if(curdirect[j]!=-1){lo[j]=0.4*lo[j];newtht0[j]<-newtht[j]+lo[j]*estimdiff[j]}}
if(prevdirect[j]==0){newtht0[j]<-newtht[j];lo[j]=0.5*lo[j]}}
if(any(newtht0<0.00000001)){
newtht0[newtht0<0.00000001]<-newtht[newtht0<0.00000001]
curdirect[newtht0<0.00000001]<-0
lo[newtht0<0.00000001]<-0}
if(newtht0[length(newtht0)]>0.9999999){
newtht0[length(newtht0)]<-newtht[length(newtht0)]
curdirect[length(newtht0)]<-0;lo[length(newtht0)]<-0}}
if(G==2){newtht0<-newtht}
list(newt=newtht0,direct=curdirect,lo=lo)}



.bclognfitdv<-function(X,Y,initfrailp,control){
time=Y[, 2];censor=Y[, 3]
order=sort(time, decreasing = FALSE, index.return = TRUE);indx=order$ix;timeo=order$x
uniq_tim<-unique(sort(time))
ind.haz=match(time,uniq_tim)
if(any(is.na(ind.haz))){
tord_diff<-as.array(diff(c(0,timeo)))
id.zero_tord <- which(apply(((tord_diff<0.0000001)&(tord_diff>0)),1, all))
if(length(id.zero_tord)>0){time[indx[id.zero_tord]]<- time[indx[id.zero_tord-1]]
order=sort(time, decreasing = FALSE, index.return = TRUE);indx=order$ix;timeo=order$x}
uniq_tim<-unique(sort(time))
ind.haz=match(time,uniq_tim)}
data.n1 <- nrow(X);data.n <-data.n1/2#### data.n is the number of pairs
indic1<-2*array(1:data.n)-1;indic2<-2*array(1:data.n)
PPID<-array(1:data.n1)
pid=match(PPID,indx);i1<-pid[indic1];i2<-pid[indic2]
i3 <-pmin(i1,i2)
INC<-matrix(c(i1,i2,i3),data.n,3)
risskset <- function(timeo,x) ifelse(timeo[1]>=x,1,0)
RI <- apply(as.array(timeo),1,risskset,x=timeo)
numbevent <- function(uniq_tim,time,censor){
indic0=match(time,uniq_tim[1])*array(1:data.n1)
sum(censor[indic0[!is.na(indic0)]])}
n_eve<-apply(as.array(uniq_tim),1,numbevent,time,censor)
Xo=X[indx,]
Xo=as.matrix(Xo)
censoro=censor[indx]
cph0 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset = NULL, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
bet<-cph0$coefficients
ncovar_coef=length(bet)
W<-par<-rep(0,(data.n1+ncovar_coef))
par[1:ncovar_coef]<-bet
par0<-par
sumx=as.vector(colSums(censoro*Xo))
newtht0=newtht<-NULL
if(length(initfrailp)==0){
newtht0<-newtht<-c(0.5,0.5,0.5)
newthtinit=newtht0}
if(length(initfrailp)>0){
newtht0<-newtht<-initfrailp
newthtinit=newtht0}
fittrpen=do.call(nlm, args=c(list(f=.llpenn1,p=par,newtht=newtht0,
censoro=censoro,sumx=sumx,Xo=Xo,RI=RI,INC=INC,fscale = control$fscale,
print.level = control$print.level, ndigit = control$ndigit,steptol = control$steptol,
iterlim= control$iterlim,gradtol = control$gradtol,check.analyticals = control$check.analyticals)))
par=fittrpen$estimate
xx<-.hesfunc1(par=par[(1+ncovar_coef):(data.n1+ncovar_coef)],bet=par[1:ncovar_coef],
censoro=censoro,Xo=Xo,RI=RI,INC=INC)
fittr=do.call(nlminb, args=c(list(start=newtht0, objective=.Lliklogndv,
gradient = .fdLliklogndv,xx=xx,lower = control$lower, upper = control$upper,
control=control$nlminb_control)))
newtht=fittr$par##### obtain new
direct<-newtht-newtht0
direct[direct<0]<--1
direct[direct>0]<-1
direct[direct==0]<-0
estim0=newtht0<-c(newtht)
new.diff=1
oldlik<-fittr$objective
lo=rep(5,length(newtht))
G=0
if(control$fastfit){
iter=0
repeat{
iter=iter+1
fittrpen=do.call(nlm, args=c(list(f=.llpenn1,p=par,newtht=newtht0,
censoro=censoro,sumx=sumx,Xo=Xo,RI=RI,INC=INC,fscale = control$fscale,
print.level = control$print.level, ndigit = control$ndigit,steptol =control$steptol ,
iterlim= control$iterlim,gradtol = control$gradtol,check.analyticals = control$check.analyticals)))
par=fittrpen$estimate
xx<-.hesfunc1(par=par[(1+ncovar_coef):(data.n1+ncovar_coef)],bet=par[1:ncovar_coef],
censoro=censoro,Xo=Xo,RI=RI,INC=INC)
fittr=do.call(nlminb, args=c(list(start=newtht0, objective=.Lliklogndv,
gradient = .fdLliklogndv,xx=xx,lower = control$lower, upper = control$upper,
control=control$nlminb_control)))
newtht=fittr$par##### obtain new estimates of sigma^2 and row
newdif=max(abs(newtht-estim0))
if(iter>200){G=2}
newtht00<-.machinelern(newtht,newtht0,direct,lo,G)
newtht0<-newtht00$newt;direct<-newtht00$direct;lo<-newtht00$lo
estim0=newtht0
if((newdif < control$toll)  |  (iter >=control$max.iter2)) break}}
if(!control$fastfit){
iter=0
repeat{
iter=iter+1
fittrpen=do.call(nlm, args=c(list(f=.llpenn1,p=par,newtht=newtht0,
censoro=censoro,sumx=sumx,Xo=Xo,RI=RI,INC=INC,fscale = control$fscale,
print.level = control$print.level, ndigit = control$ndigit,steptol =control$steptol ,
iterlim= control$iterlim,gradtol = control$gradtol,check.analyticals = control$check.analyticals)))
par=fittrpen$estimate
llppl<-.llppl1(par=par,newtht=newtht0,censoro=censoro,Xo=Xo,RI=RI)
xx<-.hesfunc1(par=par[(1+ncovar_coef):(data.n1+ncovar_coef)],bet=par[1:ncovar_coef],
censoro=censoro,Xo=Xo,RI=RI,INC=INC)
fittr=do.call(nlminb, args=c(list(start=newtht0, objective=.Lliklogndv,
gradient = .fdLliklogndv,xx=xx,lower = control$lower, upper = control$upper,
control=control$nlminb_control)))
newtht=fittr$par##### obtain new estimates of sigma^2 and row
newdif=max(abs(newtht-estim0))
newtht0<-newtht
estim0<-newtht0
if((newdif < control$toll)  |  (iter >=control$max.iter2)) break}}
if (iter >= control$max.iter2){warning("Ran out of iterations and did not converge")}
W<-1
W[indx]<-par[(1+ncovar_coef):(data.n1+ncovar_coef)]
cph1 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset = W, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
bet<-par[1:ncovar_coef]
llppl<-.llppl1(par=par,newtht=newtht,censoro=censoro,Xo=Xo,RI=RI)
lik=(-fittr$objective+llppl)
fishinfmat<-.sdLliklogndv(thetaa=newtht, xx=xx)
vcovth<-solve(fishinfmat)
if(any(is.nan(sqrt(diag(vcovth))))){
warning("Check the approperiatness of the model.")
warning("frailty parameters might be at the boundary of parameter space.")
warning("In obtaining SE Na is produced")}
colnames(vcovth) <- rownames(vcovth) <- c("theta1","theta2","Row")
vcov=cph1$var
colnames(vcov) <- rownames(vcov) <- colnames(X)
adj_se=.SEpenal1(par=par,newtht=newtht,censoro=censoro,Xo=Xo,RI=RI,INC=INC)
adjse=c(adj_se,sqrt(diag(vcovth)))
res <-list(coefficients=bet,frailparest= c(theta1=newtht[1],theta2=newtht[2],Row=newtht[3]),
vcov = vcov,vcovth= vcovth,stderr=adjse,loglilk0=cph0$loglik[1],loglilk=cph1$loglik,
penloglilk=-fittrpen$minimum,Iloglilk=lik,X=X,time=time,censor=censor,
resid=cph1$residuals,lin.prid=cph1$linear.predictors,
frail=exp(W),newthtinit=newthtinit,
iteration=iter,inh=ind.haz,indx=indx,e.time=uniq_tim,n.event=n_eve,
converg = ifelse(new.diff<control$tol,TRUE,FALSE))
res$call <- match.call()
class(res) <- c(".bclognfitdv")
res
}


############ End Bivariate Correlated lognormal fit functions###
########################################################################



########################################################################
############Univariate and bivariate shared gamma frailty fit functions###
######

############Univariate  gamma frailty fit functions###

llpengamuniv = function (par,censoro,Xo,RI0,thtfrail){
Xo=as.matrix(Xo);RI0=RI0;RI0=as.matrix(RI0)
n1=nrow(Xo);n4=ncol(Xo)
a=thtfrail;bet=matrix(c(par[array(1:n4)]),n4,1)
w0=matrix(c(par[array((1+n4):(n1+n4))]),n1,1)
d0=Xo%*%bet;g0=exp(w0+d0);G0=RI0%*%g0
logli=(c(colSums(c(censoro)*(d0+w0-log(G0))))+sum( c(1/a)*(c(w0)-exp(c(w0)))))
return(-logli)}

grdpengamuniv = function (par,censoro,Xo,RI0,thtfrail){
Xo=as.matrix(Xo);RI0=RI0;RI0=as.matrix(RI0)
n1=nrow(Xo);n4=ncol(Xo)
a=thtfrail;bet=matrix(c(par[array(1:n4)]),n4,1)
w0=matrix(c(par[array((1+n4):(n1+n4))]),n1,1)
d0=Xo%*%bet;g0=exp(w0+d0);G0=RI0%*%g0
Xog=Xo*c(g0);DG0=(RI0%*%Xog)*c(1/c(G0))
Dbeta=matrix(matrix(c(colSums((Xo*c(censoro)))),n4,1)-matrix(c(colSums((DG0*c(censoro)))),n4,1),n4,1)
cs=c(censoro/c(G0))
Dwo=c(c(censoro)-c(c(g0*cumsum(cs))))+c(1/a)*(1-exp(c(w0)))
Dlogli=c(c(Dbeta),Dwo)
return(-Dlogli)}


.Llikgamuniv= function (thetaa, xx){
xx=as.matrix(xx)
H=xx[,1];cen=xx[,2]
a = abs(thetaa)
logli= sum(-(1/a+cen)*log(1+a*H))
-logli}

.fdLlikgamuniv = function (thetaa, xx){
xx=as.matrix(xx)
H=xx[,1];cen=xx[,2]
a = abs(thetaa)
grr=c(sum((1/a^2)*log(1+a*H)-(1/a+cen)*(H/(1+a*H))))
-grr
}


.SEunivgam=function(bet,newtht,n_eve,etime,H,h0,censor,time,X){
a<-as.vector(newtht);n_eve<-as.vector(n_eve);etime<-as.vector(etime)
h0<-as.vector(h0);time<-as.vector(time);HH=H;HH<-as.vector(HH)
X <- as.matrix(X)
n_cov_coef= ncol(X);data.n1= nrow(X)
g0<-c(exp(X%*%bet))
AA=(1+a*HH);AA1=(1/a+HH)
n_eve0=as.numeric(n_eve>0)
trevntimein=n_eve0*array(1:length(n_eve0))
trevntimein1=trevntimein[trevntimein>0]
trevntime=etime[trevntimein1]
nonzero_h0<-h0[trevntimein1]
nev1<-n_eve[trevntimein1]
rissksetdr <- function(time,x) ifelse(time[1]>=x,1,0)
RIE <- apply(as.array(time),1,rissksetdr,x=trevntime)
RIE00 <-(t(RIE)*c(g0))
RIE0 <-(RIE00*c(1/AA1))
RIE1 <-RIE0*c(-(1/a+censor))
D2h <-(t(RIE1)%*%RIE0)
diag(D2h)<-diag(D2h)+(nev1/nonzero_h0^2)
RG2=t(D2h)
D2h[lower.tri(D2h, diag = FALSE)]<-RG2[lower.tri(RG2, diag = FALSE)]
DhDa <- c(colSums(c(c(1/a+censor)*c((1/a^2)/AA1^2))*RIE00-RIE0*c(1/a^2)))
interacmatza <- function(X,u){X*u}
resinteracmatx<-apply(as.matrix(X),2,interacmatza,u=RIE00)
IMX<-matrix(resinteracmatx,data.n1,n_cov_coef*ncol(RIE00))
DN0=((c(HH)*X)*c(-(1/a+censor)/AA1^2))
DbDh=matrix(c(colSums(IMX*c(c(1/a+censor)*c(1/AA1)))),ncol(RIE00),n_cov_coef)+ t(RIE00)%*%DN0
interacmatt <- function(X,u){X*u}
resinteracmat0<-apply(as.matrix(X),2,interacmatt,u=X)
IM<-matrix(resinteracmat0,data.n1,n_cov_coef^2)
DM=((c(HH)*IM)*c((1/a+censor)/AA1))
DN=(c(HH)*X)
D2b=matrix(c(colSums(DM)),n_cov_coef,n_cov_coef)+((t(DN))%*%DN0)
D2a=sum((2/a^3)*log(AA)-( (1/a+censor)*(HH^2/AA^2)+(2/a^2)*(HH/AA)))
DbDa=colSums(c(-(1/a^2)/AA1)*(c(HH)*X) +c(1/a+censor)*(c((1/a^2)/AA1^2)*(c(HH)*X)))
MA=rbind(D2b,DbDh,DbDa);MB=rbind(t(DbDh),D2h,DhDa);MC=c(c(DbDa),c(DhDa),c(D2a))
HES=cbind(MA,MB,MC)
INVE<-solve(HES)
seofthet=sqrt(diag(INVE))
SEofbetandsig=c(c(seofthet[1:n_cov_coef]),seofthet[length(seofthet)])
ii=cbind(rep(c(1:n_cov_coef),each=n_cov_coef),rep(c(1:n_cov_coef),n_cov_coef))
D2b0=D2b
D2b0[ii]<-INVE[ii]
list(se=SEofbetandsig,vco=D2b0)
}




.univgamfit<-function(X,Y,initfrailp,control){
time=Y[, 2];censor=Y[, 3]
order=sort(time, decreasing = FALSE, index.return = TRUE);indx=order$ix;timeo=order$x
uniq_tim<-unique(sort(time))
ind.haz=match(time,uniq_tim)
if(any(is.na(ind.haz))){
tord_diff<-as.array(diff(c(0,timeo)))
id.zero_tord <- which(apply(((tord_diff<0.0000001)&(tord_diff>0)),1, all))
if(length(id.zero_tord)>0){time[indx[id.zero_tord]]<- time[indx[id.zero_tord-1]]
order=sort(time, decreasing = FALSE, index.return = TRUE);indx=order$ix;timeo=order$x}
uniq_tim<-unique(sort(time))
ind.haz=match(time,uniq_tim)}
data.n1 <- nrow(X)
rissksetff <- function(timeo,x) ifelse(timeo[1]>=x,1,0)
RI0 <- apply(as.array(timeo),1,rissksetff,x=timeo)
numbevent <- function(uniq_tim,time,censor){
indic0=match(time,uniq_tim[1])*array(1:data.n1)
sum(censor[indic0[!is.na(indic0)]])}
n_eve<-apply(as.array(uniq_tim),1,numbevent,time,censor)
cph0 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset = NULL, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
bet<-cph0$coefficients
n4=ncovar_coef=length(bet)
par<-rep(0,(data.n1+ncovar_coef))
par[1:ncovar_coef]<-bet
newtht0=newtht<-NULL
if(length(initfrailp)==0){newtht0<-newtht<-c(0.5)}
if(length(initfrailp)>0){newtht0<-newtht<-initfrailp
newthtinit=newtht0}
Xo=X
Xo[(1:data.n1),]=X[indx,]
rissksety <- function(uniq_tim,x) ifelse(uniq_tim[1]<=x,1,0)
RII <- apply(as.array(uniq_tim),1,rissksety,x=timeo)
RII <-t(RII)
censoro=censor[indx]
new.diff=estim0=1
iter=0
repeat{
iter=iter+1
fittr=do.call(nlminb, args=c(list(start=par, objective=llpengamuniv,
gradient = grdpengamuniv,censoro=censoro,Xo=Xo,RI0=RI0,thtfrail=newtht,
control=control$nlminb_control)))
par=fittr$par;bet<-fittr$par[1:n4];xo_bet<-Xo%*%bet
w0=c(par[array((1+n4):(data.n1+n4))])
z=exp(w0);mean(z)
svexp_bet_xo=as.vector(RII%*%(exp(xo_bet+w0)))
H0<-cumsum(c(n_eve/svexp_bet_xo))###obtain initial parameter for H0(t)
H_bet_x=c(H0[ind.haz]*exp(c(X%*%bet)))
xx=matrix(c(H_bet_x,censor),data.n1,2)
fittr=do.call(nlminb, args=c(list(start=newtht, objective=.Llikgamuniv,
gradient = .fdLlikgamuniv,xx=xx,lower = c(0), upper = Inf,
control=control$nlminb_control)))
newtht;newtht=fittr$par ##### obtain new estimates of sigma^2 and row
new.diff=max(abs(abs(c(newtht))-abs(estim0)))
estim0=c(newtht)
if((new.diff < control$tol)  |  (iter >= control$max.iter)) break}
if (iter > control$max.iter){warning("Ran out of iterations and did not converge")}
h0=diff(c(0,H0));nonzero_h0=h0[h0>0]
H=H0[ind.haz];h=h0[ind.haz];x_bet=X%*%bet
H_bet_x=c(H0[ind.haz]*exp(c(X%*%bet)))
lik=(sum(censor*(c(x_bet)))+sum(log(nonzero_h0))-fittr$objective+sum(censor))
adjj_se=.SEunivgam(bet=bet,newtht=newtht,n_eve=n_eve,etime=uniq_tim,
H=c(H_bet_x),h0=h0,censor=censor,time=time,X=X)
adjse=c(adjj_se$se)
vcov=adjj_se$vco
colnames(vcov) <- rownames(vcov) <- colnames(X)
W=w0
W[indx]=w0
cph1 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset =W, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
vcov2=cph1$var
colnames(vcov2) <- rownames(vcov2) <- colnames(X)
res <-list(coefficients=bet,frailparest=newtht,
vcov = vcov,vcov2 = vcov2,stderr=adjse,loglilk0=cph0$loglik[1],loglilk=cph1$loglik[2],
Iloglilk=lik,cbasehaz=cbind(haz=H0,time=uniq_tim),X=X,time=time,censor=censor,
resid=cph1$residuals,lin.prid=cph1$linear.predictors,
frail=exp(W),iteration=iter,e.time=uniq_tim,n.event=n_eve,
converg = ifelse(new.diff< control$tol,TRUE,FALSE))
res$call <- match.call()
class(res) <- c(".univgamfit")
res
}




############bivariate shared gamma frailty fit functions###

llpengamshared = function (par,censor,X,RI,INC,indx,thtfrail){
X=as.matrix(X);RI=RI;RI=as.matrix(RI);INC=INC;INC=as.matrix(INC)
i1=INC[,1];i2=INC[,2];indic1=INC[,3];indic2=INC[,4]
n1=nrow(X);n4=ncol(X);n=n1/2
W=c(par[array((1+n4):(n+n4))])
WW=W;WW[indic1]<-W;WW[indic2]<-W
a=thtfrail;bet=matrix(c(par[array(1:n4)]),n4,1)
w=matrix(WW,n1,1)
d0=X%*%bet;g0=exp(w+d0);G0=c(RI%*%g0)
d00=c(d0)
logli=(sum(c(censor[indic1])*(d00[indic1]+W-log(G0[indic1]))+
c(censor[indic2])*(d00[indic2]+W-log(G0[indic2])))+sum(c(1/a)*(c(W)-exp(c(W)))))
return(-logli)}

grdpengamshared = function (par,censor,X,RI,INC,indx,thtfrail){
X=as.matrix(X);RI=RI;RI=as.matrix(RI);INC=INC;INC=as.matrix(INC)
i1=INC[,1];i2=INC[,2];indic1=INC[,3];indic2=INC[,4]
n1=nrow(X);n4=ncol(X);n=n1/2
W=c(par[array((1+n4):(n+n4))])
WW=W;WW[indic1]<-W;WW[indic2]<-W
a=thtfrail;bet=matrix(c(par[array(1:n4)]),n4,1)
w=matrix(WW,n1,1)
d0=X%*%bet;g0=exp(w+d0);G0=c(RI%*%g0)
g00=c(g0)
Xg=X*c(g0);DG0=(RI%*%Xg)*c(1/c(G0))
Dbeta=matrix(matrix(c(colSums((X*c(censor)))),n4,1)-matrix(c(colSums((DG0*c(censor)))),n4,1),n4,1)
cs=c(c(censor[indx])/c(G0[indx]))
cs0<-cumsum(cs)
Dwo=c((censor[indic1]+censor[indic2])- (c(g00[indic1]*cs0[i1])+c(g00[indic2]*cs0[i2])))+(c(1/a)*(1-exp(c(W))))
Dlogli=c(c(Dbeta),Dwo)
return(-Dlogli)}



.Llikgamshared= function (thetaa, xx){
xx=as.matrix(xx)
H1=xx[,1];H2=xx[,2];cen1=xx[,3];cen2=xx[,4]
di=cen1+cen2
mdi=cen1*cen2
a = abs(thetaa)
logli= sum(mdi*log(1+a))+sum(-(1/a+di)*log(1+a*(H1+H2)))
-logli}


.fdLlikgamshared = function (thetaa, xx){
xx=as.matrix(xx)
H1=xx[,1];H2=xx[,2];cen1=xx[,3];cen2=xx[,4]
di=cen1+cen2
mdi=cen1*cen2
a = abs(thetaa)
grr=c( sum(mdi*(1/(1+a)))+ sum((1/a^2)*log(1+a*(H1+H2))-(1/a+di)*((H1+H2)/(1+a*(H1+H2)))))
-grr
}


.SEsharedgam=function(bet,newtht,n_eve,etime,H,h0,censor,time,X){
a<-as.vector(newtht);n_eve<-as.vector(n_eve);etime<-as.vector(etime)
h0<-as.vector(h0);time<-as.vector(time);HH=H;HH<-as.vector(HH)
X <- as.matrix(X)
n_cov_coef=n4= ncol(X);data.n1=n1= nrow(X);data.n=n=n1/2
indic1<-2*array(1:data.n)-1;indic2<-2*array(1:data.n)
g0<-c(exp(X%*%bet))
HH1=HH[indic1];HH2=HH[indic2]
g01=g0[indic1];g02=g0[indic2]
cen1=censor[indic1];cen2=censor[indic2]
AA=(1+a*(HH1+HH2));AA1=(1/a+(HH1+HH2))
di=(cen1+cen2)
X1=X[indic1,];X2=X[indic2,]
X1<- as.matrix(X1)
X2<- as.matrix(X2)
n_eve0=as.numeric(n_eve>0)
trevntimein=n_eve0*array(1:length(n_eve0))
trevntimein1=trevntimein[trevntimein>0]
trevntime=etime[trevntimein1]
nonzero_h0<-h0[trevntimein1]
nev1<-n_eve[trevntimein1]
rissksetj <- function(time,x) ifelse(time[1]>=x,1,0)
RIU <- apply(as.array(time),1,rissksetj,x=trevntime)
RIU <- as.matrix(RIU)
RIU0<-t(((t(RIU))*c(g0)))
funbf <- function(RIU,u){
re0=RIU[u]
re0}
RRW1 <- apply(RIU,1,funbf,u=indic1)
RRW2 <- apply(RIU,1,funbf,u=indic2)
RRW1 <-RRW1*c(g01)
RRW2 <-RRW2*c(g02)
risskf <- function(RIU0,indic1,indic2){
re0=RIU0[indic1]+RIU0[indic2]
re0}
RG <- apply(RIU0,1,risskf,indic1=indic1,indic2=indic2)
RG <- as.matrix(RG)
RG0<- (RG*c(c(-1/AA1^2)*c((1/a+di))))
D2h <-((t(RG))%*%RG0)
diag(D2h)<-diag(D2h)+(nev1/nonzero_h0^2)
RG2=t(D2h)
D2h[lower.tri(D2h, diag = FALSE)]<-RG2[lower.tri(RG2, diag = FALSE)]
interacmat1 <- function(X1,u){X1*u}
resinteracmat1<-apply(as.matrix(X1),2,interacmat1,u=X1)
IM1<-matrix(resinteracmat1,data.n,n_cov_coef^2)
interacmat2 <- function(X2,u){X2*u}
resinteracmat2<-apply(as.matrix(X2),2,interacmat2,u=X2)
IM2<-matrix(resinteracmat2,data.n,n_cov_coef^2)
resinteracmat2<-apply(as.matrix(X2),2,interacmat2,u=X1)
IM21<-matrix(resinteracmat2,data.n,n_cov_coef^2)
resinteracmat1<-apply(as.matrix(X1),2,interacmat1,u=X2)
IM12<-matrix(resinteracmat1,data.n,n_cov_coef^2)
DM=(c(a)*(c(HH1)*IM1+c(HH2)*IM2))*c((1/a+di)/AA)
DN=c(a)*(c(HH1)*X1+c(HH2)*X2)
DN0=(c(a)*(c(HH1)*X1+c(HH2)*X2))*c(-(1/a+di)/AA^2)
D2b=matrix(c(colSums(DM)),n_cov_coef,n_cov_coef)+((t(DN))%*%DN0)
D2a=sum(cen1*cen2*(1/(1+a)^2))+sum((2/a^3)*log(AA))- sum((1/a+di)*((HH1+HH2)^2/AA^2)+(2/a^2)*((HH1+HH2)/AA))
DbDa=colSums(c(-(1/a^2)/AA1)*(c(HH1)*X1+c(HH2)*X2) +c(1/a+di)*(c((1/a^2)/AA1^2)*(c(HH1)*X1+c(HH2)*X2)))
interacmatF1 <- function(X1,u){X1*u}
resinteracmatx1<-apply(as.matrix(X1),2,interacmatF1,u=RRW1)
IMx1<-matrix(resinteracmatx1,data.n,n_cov_coef*ncol(RRW1))
interacmatF2 <- function(X2,u){X2*u}
resinteracmatx2<-apply(as.matrix(X2),2,interacmatF2,u=RRW2)
IMx2<-matrix(resinteracmatx2,data.n,n_cov_coef*ncol(RRW2))
DbDh=matrix(c(colSums((IMx1+IMx2)*c(c(a/AA)*c((1/a+di))))),ncol(RRW2),n_cov_coef)+((t(c(a)*RG))%*%DN0)
DhDa=colSums((-1/a^2)*((RG*c(1/AA1)))+(RG*c(c((1/a^2)/AA1^2)*c((1/a+di)))))
MA=rbind(D2b,DbDh,DbDa);MB=rbind(t(DbDh),D2h,DhDa);MC=c(c(DbDa),c(DhDa),c(D2a))
HES=cbind(MA,MB,MC)
INVE<-solve(HES)
seofthet=sqrt(diag(INVE))
SEofbetandsig=c(c(seofthet[1:n_cov_coef]),seofthet[length(seofthet)])
SEofbetandsig
ii=cbind(rep(c(1:n_cov_coef),each=n_cov_coef),rep(c(1:n_cov_coef),n_cov_coef))
D2b0=D2b
D2b0[ii]<-INVE[ii]
list(se=SEofbetandsig,vco=D2b0)
}


.sharedgamfit<-function(X,Y,initfrailp,control){
time=Y[, 2];censor=Y[, 3]
order=sort(time, decreasing = FALSE, index.return = TRUE);indx=order$ix;timeo=order$x
uniq_tim<-unique(sort(time))
ind.haz=match(time,uniq_tim)
if(any(is.na(ind.haz))){
tord_diff<-as.array(diff(c(0,timeo)))
id.zero_tord <- which(apply(((tord_diff<0.0000001)&(tord_diff>0)),1, all))
if(length(id.zero_tord)>0){time[indx[id.zero_tord]]<- time[indx[id.zero_tord-1]]
order=sort(time, decreasing = FALSE, index.return = TRUE);indx=order$ix;timeo=order$x}
uniq_tim<-unique(sort(time))
ind.haz=match(time,uniq_tim)}
data.n1 <- nrow(X);data.n <-data.n1/2#### data.n is the number of pairs
rissksetff <- function(time,x) ifelse(time[1]>=x,1,0)
RI <- apply(as.array(time),1,rissksetff,x=time)
numbevent <- function(uniq_tim,time,censor){
indic0=match(time,uniq_tim[1])*array(1:data.n1)
sum(censor[indic0[!is.na(indic0)]])}
n_eve<-apply(as.array(uniq_tim),1,numbevent,time,censor)
indic1<-2*array(1:data.n)-1;indic2<-2*array(1:data.n)
cph0 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset = NULL, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
bet<-cph0$coefficients
n4=ncovar_coef=length(bet)
par<-rep(0,(data.n+ncovar_coef))
par[1:ncovar_coef]<-bet
newtht0=newtht<-NULL
if(length(initfrailp)==0){newtht0<-newtht<-c(0.5)}
if(length(initfrailp)>0){newtht0<-newtht<-initfrailp
newthtinit=newtht0}
Xo=X
Xo[(1:data.n1),]=X[indx,]
rissksety <- function(uniq_tim,x) ifelse(uniq_tim[1]<=x,1,0)
RII <- apply(as.array(uniq_tim),1,rissksety,x=timeo)
RII <-t(RII)
n=data.n
order=sort(indx, decreasing = FALSE, index.return = TRUE)
indx0=order$ix
cen1=censor[indic1];cen2=censor[indic2]
INC=matrix(c(indx0[indic1],indx0[indic2],indic1,indic2),data.n,4)
new.diff=estim0=1
iter=0
repeat{
iter=iter+1
fittr=do.call(nlminb, args=c(list(start=par, objective=llpengamshared,
gradient = grdpengamshared,censor=censor,X=X,RI=RI,INC=INC,indx=indx,thtfrail=newtht,
control=control$nlminb_control)))
bet;par=fittr$par;bet<-fittr$par[1:ncovar_coef];bet;xo_bet<-Xo%*%bet
W=c(par[array((1+ncovar_coef):(data.n+ncovar_coef))])
WW=W;WW[indic1]<-W;WW[indic2]<-W
w0=WW
w0=WW[indx]
z=exp(w0);mean(z)
svexp_bet_xo=as.vector(RII%*%(exp(xo_bet+w0)))
H0<-cumsum(c(n_eve/svexp_bet_xo))###obtain initial parameter for H0(t)
H_bet_x=c(H0[ind.haz]*exp(c(X%*%bet)))
xx=matrix(c(H_bet_x[indic1],H_bet_x[indic2],cen1,cen2),data.n,4)
fittr=do.call(nlminb, args=c(list(start=newtht, objective=.Llikgamshared,
gradient = .fdLlikgamshared,xx=xx,lower = control$lower, upper = control$upper,
control=control$nlminb_control)))
newtht;newtht=fittr$par ##### obtain new estimates of sigma^2 and row
new.diff=max(abs(abs(c(newtht))-abs(estim0)))
estim0=c(newtht)
if((new.diff < control$tol)  |  (iter >= control$max.iter)) break}
if (iter > control$max.iter){warning("Ran out of iterations and did not converge")}
h0=diff(c(0,H0));nonzero_h0=h0[h0>0]
H=H0[ind.haz];h=h0[ind.haz];x_bet=X%*%bet
lik=(sum(censor*(c(x_bet)))+sum(log(nonzero_h0))-fittr$objective+sum(censor))
adjj_se=.SEsharedgam(bet=bet,newtht=newtht,n_eve=n_eve,etime=uniq_tim,
H=c(H_bet_x),h0=h0,censor=censor,time=time,X=X)
adjse=c(adjj_se$se)
vcov=adjj_se$vco
colnames(vcov) <- rownames(vcov) <- colnames(X)
W=w0
W[indx]=w0
cph1 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset =W, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
vcov2=cph1$var
colnames(vcov2) <- rownames(vcov2) <- colnames(X)
res <-list(coefficients=bet,frailparest=newtht,
vcov = vcov,vcov2 = vcov2,stderr=adjse,loglilk0=cph0$loglik[1],loglilk=cph1$loglik[2],
Iloglilk=lik,cbasehaz=cbind(haz=H0,time=uniq_tim),X=X,time=time,censor=censor,
resid=cph1$residuals,lin.prid=cph1$linear.predictors,
frail=exp(W),iteration=iter,e.time=uniq_tim,n.event=n_eve,
converg = ifelse(new.diff< control$tol,TRUE,FALSE))
res$call <- match.call()
class(res) <- c(".sharedgamfit")
res
}

############ END Univariate and bivariate shared gamma frailty fit functions###
########################################################################

##################################################
##########SE for bcfraildv gamma fit###########

.SEbcfrailgamdv=function(bet,newtht,n_eve,etime,h0,censor,time,X,H){
n_eve<-as.vector(n_eve);etime<-as.vector(etime);censor<-as.vector(censor)
h0<-as.vector(h0);time<-as.vector(time)
newtht<-as.vector(newtht);a1 = (newtht[1]);a2 = (newtht[2]);R = (newtht[3])
X <- as.matrix(X);g0<-c(exp(X%*%bet))
HH=H;HH<-as.vector(HH)
n_cov_coef=n4= ncol(X);data.n1=n1= nrow(X);data.n=n=n1/2
indic1<-2*array(1:data.n)-1;indic2<-2*array(1:data.n)
k0<-(R*a1^(-1/2)*a2^(-1/2));k1<-(1/a1-k0);k2<-(1/a2-k0)
c1<-(R^2+R*a1^(1/2)*a2^(1/2));c2<-(R*(a2^(1/2)/a1^(1/2))-R^2)
c3<-(R*(a1^(1/2)/a2^(1/2))-R^2);c4<-(1-R*(a2^(1/2)/a1^(1/2))-R*(a1^(1/2)/a2^(1/2))+R^2)
da1k0<-(-(1/2)*(R*a1^(-3/2)*a2^(-1/2)));da2k0<-(-(1/2)*(R*a1^(-1/2)*a2^(-3/2)));dRk0<-(a1^(-1/2)*a2^(-1/2))
d2a1k0<-((3/4)*(R*a1^(-5/2)*a2^(-1/2)));da1a2k0<-((1/4)*(R*a1^(-3/2)*a2^(-3/2)));d2a2k0<-((3/4)*(R*a1^(-1/2)*a2^(-5/2)))
da1Rk0<-(-(1/2)*(a1^(-3/2)*a2^(-1/2)));da2Rk0<-(-(1/2)*(a1^(-1/2)*a2^(-3/2)))
da1k1<-(-(1/a1^2)-da1k0);da2k1<-(-da2k0);dRk1<-(-dRk0);d2a1k1<-((2/a1^3)-d2a1k0)
d2a2k1<-(-d2a2k0);da1a2k1<-(-da1a2k0);da1Rk1<-(-da1Rk0);da2Rk1<-(-da2Rk0)
da1k2<-(-da1k0);da2k2<-(-(1/a2^2)-da2k0);dRk2<-(-dRk0);d2a1k2<-(-d2a1k0);d2a2k2<-((2/a2^3)-d2a2k0)
da1a2k2<-(-da1a2k0);da1Rk2<-(-da1Rk0);da2Rk2<-(-da2Rk0)
da1c1=((1/2)*R*a1^(-1/2)*a2^(1/2));da2c1=((1/2)*R*a1^(1/2)*a2^(-1/2));dRc1<-(2*R+a1^(1/2)*a2^(1/2))
d2a1c1=(-(1/4)*R*a1^(-3/2)*a2^(1/2));da1a2c1=((1/4)*R*a1^(-1/2)*a2^(-1/2));d2a2c1=(-(1/4)*R*a1^(1/2)*a2^(-3/2))
d2Rc1<-(2);da1Rc1=((1/2)*a1^(-1/2)*a2^(1/2));da2Rc1=((1/2)*a1^(1/2)*a2^(-1/2))
da1c2<-((-1/2)*R*(a2^(1/2)*a1^(-3/2)));da2c2<-((1/2)*R*(a2^(-1/2)*a1^(-1/2)));dRc2<-((a2^(1/2)/a1^(1/2))-2*R)
d2a1c2<-((3/4)*R*(a2^(1/2)*a1^(-5/2)));d2a2c2<-(-(1/4)*R*(a2^(-3/2)*a1^(-1/2)));d2Rc2<-(-2)
da1a2c2<-((-1/4)*R*(a2^(-1/2)*a1^(-3/2)));da1Rc2<-((-1/2)*(a2^(1/2)*a1^(-3/2)));da2Rc2<-((1/2)*(a2^(-1/2)*a1^(-1/2)))
da1c3<-((1/2)*R*(a2^(-1/2)*a1^(-1/2)));da2c3<-((-1/2)*R*(a1^(1/2)*a2^(-3/2)));dRc3<-((a1^(1/2)/a2^(1/2))-2*R)
d2a1c3<-(-(1/4)*R*(a2^(-1/2)*a1^(-3/2)));d2a2c3<-((3/4)*R*(a1^(1/2)*a2^(-5/2)));d2Rc3<-(-2)
da1a2c3<-(-(1/4)*R*(a2^(-3/2)*a1^(-1/2)));da1Rc3<-((1/2)*(a2^(-1/2)*a1^(-1/2)));da2Rc3<-((-1/2)*(a1^(1/2)*a2^(-3/2)))
da1c4<-((1/2)*R*(a2^(1/2)*a1^(-3/2))-(1/2)*R*(a1^(-1/2)*a2^(-1/2)));da2c4<-(-(1/2)*R*(a2^(-1/2)*a1^(-1/2))+(1/2)*R*(a1^(1/2)*a2^(-3/2)))
dRc4<-(-(a2^(1/2)*a1^(-1/2))-(a1^(1/2)*a2^(-1/2))+2*R);d2Rc4<-(2);da1Rc4<-((1/2)*(a2^(1/2)*a1^(-3/2))-(1/2)*(a1^(-1/2)*a2^(-1/2)))
da2Rc4<-(-(1/2)*(a2^(-1/2)*a1^(-1/2))+(1/2)*(a1^(1/2)*a2^(-3/2)));d2a1c4<-(-(3/4)*R*(a2^(1/2)*a1^(-5/2))+(1/4)*R*(a1^(-3/2)*a2^(-1/2)))
d2a2c4<-((1/4)*R*(a2^(-3/2)*a1^(-1/2))-(3/4)*R*(a1^(1/2)*a2^(-5/2)));da1a2c4<-((1/4)*R*(a2^(-1/2)*a1^(-3/2))+(1/4)*R*(a1^(-1/2)*a2^(-3/2)))
H1=HH[indic1];H2=HH[indic2]
time1=time[indic1];time2=time[indic2]
g01=g0[indic1];g02=g0[indic2]
cen1=censor[indic1];cen2=censor[indic2]
di=(cen1+cen2)
X1=X[indic1,];X2=X[indic2,]
X1<- as.matrix(X1)
X2<- as.matrix(X2)
AA=(1+a1*H1+a2*H2);AA1=(1+a1*H1);AA2=(1+a2*H2)
n_eve0=as.numeric(n_eve>0)
trevntimein=n_eve0*array(1:length(n_eve0))
trevntimein1=trevntimein[trevntimein>0]
trevntime=etime[trevntimein1]
nonzero_h0<-h0[trevntimein1]
nev1<-n_eve[trevntimein1]
rissksetdr1 <- function(time1,x) ifelse(time1[1]>=x,1,0)
RIE1 <- apply(as.array(time1),1,rissksetdr1,x=trevntime)
RIE001 <-(t(RIE1)*c(g01))
rissksetdr2 <- function(time2,x) ifelse(time2[1]>=x,1,0)
RIE2 <- apply(as.array(time2),1,rissksetdr2,x=trevntime)
RIE002 <-(t(RIE2)*c(g02))
interacmatza <- function(X,u){X*u}
resinteracmatx<-apply(as.matrix(X1),2,interacmatza,u=RIE001)
IMX1<-matrix(resinteracmatx,data.n,n_cov_coef*ncol(RIE001))
resinteracmatx<-apply(as.matrix(X2),2,interacmatza,u=RIE001)
IMX12<-matrix(resinteracmatx,data.n,n_cov_coef*ncol(RIE001))
resinteracmatx<-apply(as.matrix(X2),2,interacmatza,u=RIE002)
IMX2<-matrix(resinteracmatx,data.n,n_cov_coef*ncol(RIE002))
resinteracmatx<-apply(as.matrix(X1),2,interacmatza,u=RIE002)
IMX21<-matrix(resinteracmatx,data.n,n_cov_coef*ncol(RIE002))
interacmat1 <- function(X1,u){X1*u}
resinteracmat1<-apply(as.matrix(X1),2,interacmat1,u=X1)
IM1<-matrix(resinteracmat1,data.n,n_cov_coef^2)
interacmat2 <- function(X2,u){X2*u}
resinteracmat2<-apply(as.matrix(X2),2,interacmat2,u=X2)
IM2<-matrix(resinteracmat2,data.n,n_cov_coef^2)
resinteracmat2<-apply(as.matrix(X2),2,interacmat2,u=X1)
IM21<-matrix(resinteracmat2,data.n,n_cov_coef^2)
resinteracmat1<-apply(as.matrix(X1),2,interacmat1,u=X2)
IM12<-matrix(resinteracmat1,data.n,n_cov_coef^2)
Dhp00i=c(-c(k1+cen1)*(1/AA1))*(c(a1)*RIE001)-c(c(k2+cen2)*(1/AA2))*(c(a2)*RIE002)-c(c(k0+di)*c(1/AA))*(c(a1)*RIE001+c(a2)*RIE002)
Dbp00i=(c(-c(k1+cen1)*(1/AA1))*(c(a1*H1)*X1)+c(-c(k2+cen2)*(1/AA2))*(c(a2*H2)*X2)+c(-c(k0+di)*(1/AA))*(c(a1*H1)*X1+c(a2*H2)*X2))
Da1p00i= -c(da1k1)*log(AA1)-c(k1+cen1)*(H1/AA1)-c(da1k2)*log(AA2)-c(da1k0)*log(AA)-c(k0+di)*(H1/AA)
Da2p00i= -c(da2k2)*log(AA2)-c(k2+cen2)*(H2/AA2)-c(da2k1)*log(AA1)-c(da2k0)*log(AA)-c(k0+di)*(H2/AA)
DRp00i=-c(dRk1)*log(AA1)-c(dRk2)*log(AA2)-c(dRk0)*log(AA)
D2hp00=(t((c(k1+cen1)*c(1/AA1^2)*c(a1)*RIE001))%*%(c(a1)*RIE001)+ t((c(k2+cen2)*c(1/AA2^2)*c(a2)*RIE002))%*%(c(a2)*RIE002)+
t((c(k0+di)*c(1/AA^2)*(c(a1)*RIE001+c(a2)*RIE002)))%*%(c(a1)*RIE001+c(a2)*RIE002))#####
Dhbp00i= (matrix(c(colSums(c(-c(k1+cen1)*c(1/AA1))*(c(a1)*IMX1))),ncol(RIE001),n_cov_coef)+ t(c(a1)*RIE001)%*%( c(c(k1+cen1)*c(1/AA1^2))*(c(a1*H1)*X1))+
matrix(c(colSums(c(-c(k2+cen2)*c(1/AA2))*(c(a2)*IMX2))),ncol(RIE002),n_cov_coef)+ t(c(a2)*RIE002)%*%(c(c(k2+cen2)*c(1/AA2^2))*(c(a2*H2)*X2))+
matrix(c(colSums(-c(c(k0+di)*c(1/AA))*(c(a1)*IMX1+c(a2)*IMX2))),ncol(RIE001),n_cov_coef)+
t(c(a1)*RIE001+c(a2)*RIE002 )%*%( c(c(k0+di)*c(1/AA^2))*(c(a1*H1)*X1+c(a2*H2)*X2)))####
Dha1p00i=colSums(c(-c(da1k1)*(1/AA1))*(c(a1)*RIE001)  -c(k1+cen1)*(c(1/AA1)*RIE001 -c(H1/AA1^2)*c(a1)*RIE001)-c(c(da1k2)*(1/AA2))*(c(a2)*RIE002)-
c(c(da1k0)*c(1/AA))*(c(a1)*RIE001+c(a2)*RIE002)-c(k0+di)*(c(1/AA)*RIE001 -c(H1/AA^2)*(c(a1)*RIE001+c(a2)*RIE002)))#####
Dha2p00i=colSums(c(-c(da2k2)*(1/AA2))*(c(a2)*RIE002)  -c(k2+cen2)*(c(1/AA2)*RIE002 -c(H2/AA2^2)*c(a2)*RIE002)-c(c(da2k1)*(1/AA1))*(c(a1)*RIE001)-
c(c(da2k0)*c(1/AA))*(c(a1)*RIE001+c(a2)*RIE002)-c(k0+di)*(c(1/AA)*RIE002 -c(H2/AA^2)*(c(a1)*RIE001+c(a2)*RIE002)))######
DhRp00i=colSums(c(-c(dRk1)*(1/AA1))*(c(a1)*RIE001)-c(c(dRk2)*(1/AA2))*(c(a2)*RIE002)-c(c(dRk0)*c(1/AA))*(c(a1)*RIE001+c(a2)*RIE002))#####
DN1=(c(a1*H1)*X1);DN2=(c(a2*H2)*X2);DN12=(c(a1*H1)*X1+c(a2*H2)*X2)
DN01=((c(a1*H1)*X1)*c(-(k1+cen1)/AA1^2));DN02=((c(a2*H2)*X2)*c(-(k2+cen2)/AA2^2))
DN012=((c(a1*H1)*X1+c(a2*H2)*X2)*c((k0+di)/AA^2))
DM1=((c(a1*H1)*IM1)*c(-(k1+cen1)/AA1));DM2=((c(a2*H2)*IM2)*c(-(k2+cen2)/AA2))
DM12=((c(a1*H1)*IM1+c(a2*H2)*IM2)*c(-(k0+di)/AA))
D2bp00i=(matrix(c(colSums(DM1)),n_cov_coef,n_cov_coef)-((t(DN1))%*%DN01)+matrix(c(colSums(DM2)),n_cov_coef,n_cov_coef)-((t(DN2))%*%DN02)+
matrix(c(colSums(DM12)),n_cov_coef,n_cov_coef)+((t(DN12))%*%DN012))######
Dba1p00i=colSums(c(-c(da1k1)*(1/AA1))*(c(a1*H1)*X1)-c(k1+cen1)*(c(1/AA1)*(c(H1)*X1)- c(1/AA1^2)*(c(a1*H1^2)*X1))+
c(-c(da1k2)*(1/AA2))*(c(a2*H2)*X2)-c(c(da1k0)*(1/AA))*(c(a1*H1)*X1+c(a2*H2)*X2)-c(k0+di)*(c(1/AA)*(c(H1)*X1)-c(H1/AA^2)*(c(a1*H1)*X1+c(a2*H2)*X2)))##
Dba2p00i=colSums(c(-c(da2k2)*(1/AA2))*(c(a2*H2)*X2)-c(k2+cen2)*(c(1/AA2)*(c(H2)*X2)- c(1/AA2^2)*(c(a2*H2^2)*X2))+
c(-c(da2k1)*(1/AA1))*(c(a1*H1)*X1)-c(c(da2k0)*(1/AA))*(c(a1*H1)*X1+c(a2*H2)*X2)-c(k0+di)*(c(1/AA)*(c(H2)*X2)-c(H2/AA^2)*(c(a1*H1)*X1+c(a2*H2)*X2)))##
DbRp00i=colSums(c(-c(dRk1)*(1/AA1))*(c(a1*H1)*X1)+c(-c(dRk2)*(1/AA2))*(c(a2*H2)*X2)+c(-c(dRk0)*(1/AA))*(c(a1*H1)*X1+c(a2*H2)*X2))##
D2a1p00i= sum(-c(d2a1k1)*log(AA1)-c(da1k1)*(2*H1/AA1)+c(k1+cen1)*(H1^2/AA1^2)-c(d2a1k2)*log(AA2)+
(-c(d2a1k0)*log(AA)-c(da1k0)*(2*H1/AA))+c(k0+di)*(H1^2/AA^2))
Da1a2p00i=sum(-c(da1a2k1)*log(AA1)-c(da2k1)*(H1/AA1)-c(da1a2k2)*log(AA2)-c(da1k2)*(H2/AA2)-
c(da1a2k0)*log(AA)-c(da1k0)*(H2/AA)-c(da2k0)*(H1/AA)+c(k0+di)*(H1*H2/AA^2))
Da1Rp00i= sum(-c(da1Rk1)*log(AA1)-c(dRk1)*(H1/AA1)-c(da1Rk2)*log(AA2)-c(da1Rk0)*log(AA)-c(dRk0)*(H1/AA))
D2a2p00i= sum(-c(d2a2k2)*log(AA2)-c(da2k2)*(2*H2/AA2)+c(k2+cen2)*(H2^2/AA2^2)-c(d2a2k1)*log(AA1)+
(-c(d2a2k0)*log(AA)-c(da2k0)*(2*H2/AA))+c(k0+di)*(H2^2/AA^2))
Da2Rp00i= sum(-c(da2Rk2)*log(AA2)-c(dRk2)*(H2/AA2)-c(da2Rk1)*log(AA1)-c(da2Rk0)*log(AA)-c(dRk0)*(H2/AA))
D2Rp00i=0
Dthtp00i=matrix(c(D2a1p00i,Da1a2p00i,Da1Rp00i,Da1a2p00i,D2a2p00i,Da2Rp00i,Da1Rp00i,Da2Rp00i,D2Rp00i),3,3)
P11=( c1*AA1*AA2+c2*AA2*AA+c3*AA1*AA+c4*AA^(2))
DhlogP11=((cen1*cen2)/P11)*( c1*( c(AA2)*(c(a1)*RIE001)+c(AA1)*(c(a2)*RIE002))+
c2*(c(AA2)*(c(a1)*RIE001+c(a2)*RIE002)+ c(AA)*(c(a2)*RIE002))+c3*(c(AA1)*(c(a1)*RIE001+c(a2)*RIE002)+ c(AA)*(c(a1)*RIE001))+
c4*2*AA*(c(a1)*RIE001+c(a2)*RIE002))
DblogP11=((cen1*cen2)/P11)*(c(c1*a1*H1*AA2)*X1+c(c1*AA1*a2*H2)*X2+
c(c2*a2*H2*AA)*X2+c(c2*AA2)*(a1*H1*X1+a2*H2*X2)+c(c3*a1*H1*AA)*X1+c(c3*AA1)*(a1*H1*X1+a2*H2*X2)+c4*2*AA*(a1*H1*X1+a2*H2*X2))
da1logP11=((cen1*cen2)/P11)*(da1c1*AA1*AA2+c1*H1*AA2+da1c2*AA2*AA+c2*H1*AA2+da1c3*AA1*AA+c3*(AA1*H1+H1*AA)+da1c4*AA^(2)+c4*2*AA*H1)
da2logP11=((cen1*cen2)/P11)*(da2c1*AA1*AA2+c1*H2*AA1+da2c2*AA2*AA+c2*(H2*AA+H2*AA2)+da2c3*AA1*AA+c3*AA1*H2+da2c4*AA^(2)+c4*2*AA*H2)
dRlogP11=((cen1*cen2)/P11)*(dRc1*AA1*AA2+dRc2*AA2*AA+dRc3*AA1*AA+dRc4*AA^(2))
D2hlogP110=(t(c(c(cen1*cen2/P11)*c(2)*c1)*(c(a1)*RIE001))%*%(c(a2)*RIE002)+
t(c(c(cen1*cen2/P11)*c(2)*c2)*(c(a2)*RIE002))%*%(c(a1)*RIE001+c(a2)*RIE002)+
t(c(c(cen1*cen2/P11)*c(2)*c3)*(c(a1)*RIE001))%*%(c(a1)*RIE001+c(a2)*RIE002)+
t(c(c(cen1*cen2/P11)*c(2)*c4)*(c(a1)*RIE001+c(a2)*RIE002))%*%(c(a1)*RIE001+c(a2)*RIE002))
D2hlogP11=D2hlogP110-t(DhlogP11)%*%DhlogP11###############
DhblogP110=((cen1*cen2)/P11)*(c1*((c(a1*a2*H2)*IMX12)+c(AA2)*(c(a1)*IMX1)+(c(a2*a1*H1)*IMX21)+c(AA1)*(c(a2)*IMX2))+
c2*((c(a1*a2*H2)*IMX12+c(a2*a2*H2)*IMX2)+ c(AA2)*(c(a1)*IMX1+c(a2)*IMX2)+(c(a1*a2*H1)*IMX21+c(a2*a2*H2)*IMX2)+c(AA)*(c(a2)*IMX2))+
c3*((c(a1*a1*H1)*IMX1+c(a1*a2*H1)*IMX21)+ c(AA1)*(c(a1)*IMX1+c(a2)*IMX2)+(c(a1*a2*H1)*IMX1+c(a2*a2*H2)*IMX12)+c(AA)*(c(a1)*IMX1))+
c4*2*((c(a1*a1*H1)*IMX1+c(a2*a1*H1)*IMX21+c(a1*a2*H2)*IMX12+c(a2*a2*H2)*IMX2)+AA*(c(a1)*IMX1+c(a2)*IMX2)))
DhblogP11=matrix(c(colSums(DhblogP110)),ncol(RIE001),n_cov_coef)-t(DhlogP11)%*%DblogP11#############
Dha1logP110=((cen1*cen2)/P11)*( da1c1*( c(AA2)*(c(a1)*RIE001)+c(AA1)*(c(a2)*RIE002))+c1*(c(AA2)*RIE001+c(H1)*(c(a2)*RIE002))+
da1c2*(c(AA2)*(c(a1)*RIE001+c(a2)*RIE002)+ c(AA)*(c(a2)*RIE002))+c2*(c(AA2)*RIE001+c(H1)*(c(a2)*RIE002))+
da1c3*(c(AA1)*(c(a1)*RIE001+c(a2)*RIE002)+ c(AA)*(c(a1)*RIE001))+c3*(c(H1)*(c(a1)*RIE001+c(a2)*RIE002)+c(AA1)*RIE001+c(H1)*(c(a1)*RIE001)+c(AA)*RIE001)+
da1c4*2*AA*(c(a1)*RIE001+c(a2)*RIE002)+c4*2*(H1*(c(a1)*RIE001+c(a2)*RIE002)+c(AA)*RIE001))
Dha1logP11=colSums((Dha1logP110)-c(da1logP11)*DhlogP11)######
Dha2logP110=((cen1*cen2)/P11)*( da2c1*( c(AA2)*(c(a1)*RIE001)+c(AA1)*(c(a2)*RIE002))+c1*(c(H2)*(c(a1)*RIE001)+c(AA1)*RIE002)+
da2c2*(c(AA2)*(c(a1)*RIE001+c(a2)*RIE002)+ c(AA)*(c(a2)*RIE002))+c2*(c(H2)*(c(a1)*RIE001+c(a2)*RIE002)+c(AA2)*RIE002+c(H2)*(c(a2)*RIE002)+c(AA)*RIE002)+
da2c3*(c(AA1)*(c(a1)*RIE001+c(a2)*RIE002)+ c(AA)*(c(a1)*RIE001))+c3*(c(AA1)*RIE002+ c(H2)*(c(a1)*RIE001))+
da2c4*2*AA*(c(a1)*RIE001+c(a2)*RIE002)+c4*2*(H2*(c(a1)*RIE001+c(a2)*RIE002)+c(AA)*RIE002))
Dha2logP11=colSums((Dha2logP110)-c(da2logP11)*DhlogP11)######
DhRlogP110=((cen1*cen2)/P11)*( dRc1*( c(AA2)*(c(a1)*RIE001)+c(AA1)*(c(a2)*RIE002))+
dRc2*(c(AA2)*(c(a1)*RIE001+c(a2)*RIE002)+ c(AA)*(c(a2)*RIE002))+dRc3*(c(AA1)*(c(a1)*RIE001+c(a2)*RIE002)+ c(AA)*(c(a1)*RIE001))+
dRc4*2*AA*(c(a1)*RIE001+c(a2)*RIE002))
DhRlogP11=colSums((DhRlogP110)-c(dRlogP11)*DhlogP11)######
D2blogP110=((cen1*cen2)/P11)*(c(c1*a1)*(c(H1*AA2)*IM1+c(H1*a2*H2)*IM21)+ c(c1*a2)*(c(a1*H1*H2)*IM12+c(AA1*H2)*IM2)+
c(c2*a2)*(c(H2*AA)*IM2+c(H2)*(a1*H1*IM12+a2*H2*IM2))+c(c2*a2*H2)*(a1*H1*IM21+a2*H2*IM2)+c(c2*AA2)*(a1*H1*IM1+a2*H2*IM2)+
c(c3*a1)*(c(H1*AA)*IM1+c(H1)*(a1*H1*IM1+a2*H2*IM21))+c(c3*a1*H1)*(a1*H1*IM1+a2*H2*IM12)+c(c3*AA1)*(a1*H1*IM1+a2*H2*IM2)+
c4*2*(c((a1*H1)^2)*IM1+c(a1*H1*a2*H2)*IM21+c(a1*H1*a2*H2)*IM12+c((a2*H2)^2)*IM2)+c4*2*AA*(a1*H1*IM1+a2*H2*IM2))
D2blogP11=matrix(c(colSums(D2blogP110)),n_cov_coef,n_cov_coef)-t(DblogP11)%*%DblogP11############
Dba1logP110=((cen1*cen2)/P11)*(c((da1c1*a1+c1)*H1*AA2)*X1+c((da1c1*AA1+c1*H1)*a2*H2)*X2+
c((da1c2*AA+c2*H1)*a2*H2)*X2+c(da1c2*AA2)*(a1*H1*X1+a2*H2*X2)+c(c2*AA2)*(H1*X1)+c((da1c3*a1+c3)*H1*AA+(c3*a1)*H1^2)*X1+
c(da1c3*AA1+c3*H1)*(a1*H1*X1+a2*H2*X2)+c(c3*AA1)*(H1*X1)+(da1c4*2*AA+c4*2*H1)*(a1*H1*X1+a2*H2*X2)+(c4*2*AA)*(H1*X1))
Dba1logP11=colSums((Dba1logP110)-c(da1logP11)*DblogP11)######
Dba2logP110=((cen1*cen2)/P11)*(c(a1*H1*(da2c1*AA2+c1*H2))*X1+c((da2c1*a2+c1)*AA1*H2)*X2+c((da2c2*a2+c2)*H2*AA+(c2*a2)*H2^2)*X2+
c(da2c2*AA2+c2*H2)*(a1*H1*X1+a2*H2*X2)+c(c2*AA2)*(H2*X2)+c((da2c3*AA+c3*H2)*a1*H1)*X1+c(da2c3*AA1)*(a1*H1*X1+a2*H2*X2)+
c(c3*AA1)*(H2*X2)+(da2c4*2*AA+c4*2*H2)*(a1*H1*X1+a2*H2*X2)+(c4*2*AA)*(H2*X2))
Dba2logP11=colSums((Dba2logP110)-c(da2logP11)*DblogP11)######
DbRlogP110=((cen1*cen2)/P11)*(c(dRc1*a1*H1*AA2)*X1+c(dRc1*AA1*a2*H2)*X2+
c(dRc2*a2*H2*AA)*X2+c(dRc2*AA2)*(a1*H1*X1+a2*H2*X2)+c(dRc3*a1*H1*AA)*X1+c(dRc3*AA1)*(a1*H1*X1+a2*H2*X2)+dRc4*2*AA*(a1*H1*X1+a2*H2*X2))
DbRlogP11=colSums((DbRlogP110)-c(dRlogP11)*DblogP11)######
d2a1logP110=((cen1*cen2)/P11)*(d2a1c1*AA1*AA2+2*da1c1*(H1*AA2)+d2a1c2*AA2*AA+2*da1c2*AA2*H1+
d2a1c3*AA1*AA+2*da1c3*(H1*AA+AA1*H1)+c3*2*H1^2+d2a1c4*AA^(2)+4*da1c4*AA*H1+c4*2*H1^2)
d2a1logP11=sum(d2a1logP110-da1logP11^2)
da1a2logP110=((cen1*cen2)/P11)*( da1a2c1*AA1*AA2+da1c1*AA1*H2+da2c1*H1*AA2+c1*H1*H2+da1a2c2*AA2*AA+da1c2*(H2*AA+AA2*H2)+da2c2*H1*AA2+c2*H1*H2+
da1a2c3*AA1*AA+da1c3*AA1*H2+da2c3*(AA1*H1+H1*AA)+c3*(H1*H2)+da1a2c4*AA^(2)+da1c4*2*AA*H2+da2c4*2*AA*H1+c4*2*H2*H1)
da1a2logP11=sum(da1a2logP110-(da1logP11*da2logP11))
da1RlogP110=((cen1*cen2)/P11)*( da1Rc1*AA1*AA2+dRc1*H1*AA2+da1Rc2*AA2*AA+dRc2*H1*AA2+da1Rc3*AA1*AA+dRc3*(AA1*H1+H1*AA)+da1Rc4*AA^(2)+dRc4*2*AA*H1)
da1RlogP11=sum(da1RlogP110-(da1logP11*dRlogP11))
d2a2logP110=((cen1*cen2)/P11)*(d2a2c1*AA1*AA2+2*da2c1*AA1*H2+d2a2c2*AA2*AA+2*da2c2*(H2*AA+H2*AA2)+c2*2*H2^2+
d2a2c3*AA1*AA+2*da2c3*AA1*H2+d2a2c4*AA^(2)+4*da2c4*AA*H2+ c4*2*H2^2)
d2a2logP11=sum(d2a2logP110-da2logP11^2)
da2RlogP110=((cen1*cen2)/P11)*(da2Rc1*AA1*AA2+dRc1*H2*AA1+da2Rc2*AA2*AA+dRc2*(H2*AA+H2*AA2)+da2Rc3*AA1*AA+dRc3*AA1*H2+da2Rc4*AA^(2)+dRc4*2*AA*H2)
da2RlogP11=sum(da2RlogP110-(da2logP11*dRlogP11))
d2RlogP110=((cen1*cen2)/P11)*(d2Rc1*AA1*AA2+d2Rc2*AA2*AA+d2Rc3*AA1*AA+d2Rc4*AA^(2))
d2RlogP11=sum(d2RlogP110-dRlogP11^2)
DthtP11=matrix(c(d2a1logP11,da1a2logP11,da1RlogP11,da1a2logP11,d2a2logP11,da2RlogP11,da1RlogP11,da2RlogP11,d2RlogP11),3,3)###
P10=(a1^(1/2)*a2^(-1/2))*R*AA1+(1-(a1^(1/2)*a2^(-1/2))*R)*AA
DhlogP10=(cen1*(1-cen2)/P10)*(c((a1^(1/2)*a2^(-1/2))*R)*(c(a1)*RIE001)+c((1-(a1^(1/2)*a2^(-1/2))*R))*(c(a1)*RIE001+c(a2)*RIE002))
DblogP10=(cen1*(1-cen2)/P10)*( c((a1^(3/2)*a2^(-1/2))*R*H1)*X1+c((1-(a1^(1/2)*a2^(-1/2))*R))*(c(a1*H1)*X1+c(a2*H2)*X2))
Da1logP10=(cen1*(1-cen2)/P10)*( (1/2)*(a1^(-1/2)*a2^(-1/2))*R*AA1+(a1^(1/2)*a2^(-1/2))*R*H1-(1/2)*(a1^(-1/2)*a2^(-1/2))*R*AA+
(1-(a1^(1/2)*a2^(-1/2))*R)*H1)
Da2logP10=(cen1*(1-cen2)/P10)*((-1/2)*(a1^(1/2)*a2^(-3/2))*R*AA1+(1/2)*(a1^(1/2)*a2^(-3/2))*R*AA +(1-(a1^(1/2)*a2^(-1/2))*R)*H2)
DRlogP10=(cen1*(1-cen2)/P10)*((a1^(1/2)*a2^(-1/2))*AA1-(a1^(1/2)*a2^(-1/2))*AA)
D2hlogP10=-t(DhlogP10)%*%DhlogP10####
DhblogP100=(cen1*(1-cen2)/P10)*(c((a1^(1/2)*a2^(-1/2))*R)*(c(a1)*IMX1)+c((1-(a1^(1/2)*a2^(-1/2))*R))*(c(a1)*IMX1+c(a2)*IMX2))
DhblogP10=matrix(c(colSums(DhblogP100)),ncol(RIE002),n_cov_coef)- t(DhlogP10)%*%DblogP10######
Dha1logP100=(cen1*(1-cen2)/P10)*(c((1/2)*(a1^(-1/2)*a2^(-1/2))*R)*(c(a1)*RIE001)+c((a1^(1/2)*a2^(-1/2))*R)*(RIE001)+
c(-(1/2)*(a1^(-1/2)*a2^(-1/2))*R)*(c(a1)*RIE001+c(a2)*RIE002)+c((1-(a1^(1/2)*a2^(-1/2))*R))*(RIE001))
Dha1logP10=colSums((Dha1logP100)-c(Da1logP10)*DhlogP10)######
Dha2logP100=(cen1*(1-cen2)/P10)*(c((-1/2)*(a1^(1/2)*a2^(-3/2))*R)*(c(a1)*RIE001)+c((1-(a1^(1/2)*a2^(-1/2))*R))*(RIE002)+
c((1/2)*(a1^(1/2)*a2^(-3/2))*R)*(c(a1)*RIE001+c(a2)*RIE002))
Dha2logP10=colSums((Dha2logP100)-c(Da2logP10)*DhlogP10)######
DhRlogP100=(cen1*(1-cen2)/P10)*(c((a1^(1/2)*a2^(-1/2)))*(c(a1)*RIE001)-c((a1^(1/2)*a2^(-1/2)))*(c(a1)*RIE001+c(a2)*RIE002))
DhRlogP10=colSums((DhRlogP100)-c(DRlogP10)*DhlogP10)######
D2blogP100=(cen1*(1-cen2)/P10)*(c((a1^(3/2)*a2^(-1/2))*R*H1)*IM1+c((1-(a1^(1/2)*a2^(-1/2))*R))*(c(a1*H1)*IM1+c(a2*H2)*IM2))
D2blogP10=matrix(c(colSums(D2blogP100)),n_cov_coef,n_cov_coef)-t(DblogP10)%*%DblogP10######
Dba1logP100=(cen1*(1-cen2)/P10)*( c((3/2)*(a1^(1/2)*a2^(-1/2))*R*H1)*X1+c(-(1/2)*(a1^(-1/2)*a2^(-1/2))*R)*(c(a1*H1)*X1+c(a2*H2)*X2)+
c((1-(a1^(1/2)*a2^(-1/2))*R))*(c(H1)*X1))
Dba1logP10=colSums(Dba1logP100-c(Da1logP10)*DblogP10)###
Dba2logP100=(cen1*(1-cen2)/P10)*( c((-1/2)*(a1^(3/2)*a2^(-3/2))*R*H1)*X1+c((1/2)*(a1^(1/2)*a2^(-3/2))*R)*(c(a1*H1)*X1+c(a2*H2)*X2)+
c((1-(a1^(1/2)*a2^(-1/2))*R))*(c(H2)*X2))
Dba2logP10=colSums(Dba2logP100-c(Da2logP10)*DblogP10)####
DbRlogP100=(cen1*(1-cen2)/P10)*( c((a1^(3/2)*a2^(-1/2))*H1)*X1+c(-(a1^(1/2)*a2^(-1/2)))*(c(a1*H1)*X1+c(a2*H2)*X2))
DbRlogP10=colSums(DbRlogP100-c(DRlogP10)*DblogP10)####
D2a1logP100=(cen1*(1-cen2)/P10)*( -(1/4)*(a1^(-3/2)*a2^(-1/2))*R*AA1+(1/4)*(a1^(-3/2)*a2^(-1/2))*R*AA)
D2a1logP10=sum(D2a1logP100-Da1logP10^2)
Da1a2logP100=(cen1*(1-cen2)/P10)*( -(1/4)*(a1^(-1/2)*a2^(-3/2))*R*AA1+(1/4)*(a1^(-1/2)*a2^(-3/2))*R*AA-(1/2)*(a1^(-1/2)*a2^(-1/2))*R*H2)
Da1a2logP10=sum(Da1a2logP100-(Da1logP10*Da2logP10))
Da1RlogP100=(cen1*(1-cen2)/P10)*( (1/2)*(a1^(-1/2)*a2^(-1/2))*AA1-(1/2)*(a1^(-1/2)*a2^(-1/2))*AA)
Da1RlogP10=sum(Da1RlogP100-(Da1logP10*DRlogP10))
D2a2logP100=(cen1*(1-cen2)/P10)*((3/4)*(a1^(1/2)*a2^(-5/2))*R*AA1-(3/4)*(a1^(1/2)*a2^(-5/2))*R*AA+(a1^(1/2)*a2^(-3/2))*R*H2)
D2a2logP10=sum(D2a2logP100-Da2logP10^2)
Da2RlogP100=(cen1*(1-cen2)/P10)*((-1/2)*(a1^(1/2)*a2^(-1/2))*AA1+(1/2)*(a1^(1/2)*a2^(-3/2))*AA-(a1^(1/2)*a2^(-1/2))*H2)
Da2RlogP10=sum(Da2RlogP100-(Da2logP10*DRlogP10))
D2RlogP10=-sum(DRlogP10^2)
DthtP10=matrix(c(D2a1logP10,Da1a2logP10,Da1RlogP10,Da1a2logP10,D2a2logP10,Da2RlogP10,Da1RlogP10,Da2RlogP10,D2RlogP10),3,3)###
P01=(a2^(1/2)*a1^(-1/2))*R*AA2+(1-(a2^(1/2)*a1^(-1/2))*R)*AA
DhlogP01=((1-cen1)*cen2/P01)*(c((a2^(1/2)*a1^(-1/2))*R)*(c(a2)*RIE002)+c((1-(a2^(1/2)*a1^(-1/2))*R))*(c(a1)*RIE001+c(a2)*RIE002))
DblogP01=((1-cen1)*cen2/P01)*( c((a2^(1/2)*a1^(-1/2))*R*a2*H2)*X2+c((1-(a2^(1/2)*a1^(-1/2))*R))*(c(a1*H1)*X1+c(a2*H2)*X2))
Da1logP01=(((1-cen1)*cen2)/P01)*( -(1/2)*(a1^(-3/2)*a2^(1/2))*R*AA2+(1-(a1^(-1/2)*a2^(1/2))*R)*H1 +(1/2)*(a1^(-3/2)*a2^(1/2))*R*AA)
Da2logP01=(((1-cen1)*cen2)/P01)*( (1/2)*(a1^(-1/2)*a2^(-1/2))*R*AA2+(a1^(-1/2)*a2^(1/2))*R*H2-(1/2)*(a1^(-1/2)*a2^(-1/2))*R*AA+(1-(a1^(-1/2)*a2^(1/2))*R)*H2)
DRlogP01=(((1-cen1)*cen2)/P01)*((a2^(1/2)*a1^(-1/2))*AA2-(a2^(1/2)*a1^(-1/2))*AA)
D2hlogP01=-t(DhlogP01)%*%DhlogP01####
DhblogP010=((1-cen1)*cen2/P01)*(c((a2^(1/2)*a1^(-1/2))*R)*(c(a2)*IMX2)+c((1-(a2^(1/2)*a1^(-1/2))*R))*(c(a1)*IMX1+c(a2)*IMX2))
DhblogP01=matrix(c(colSums(DhblogP010)),ncol(RIE002),n_cov_coef)- t(DhlogP01)%*%DblogP01######
Dha1logP010=((1-cen1)*cen2/P01)*(c( (-1/2)*(a2^(1/2)*a1^(-3/2))*R)*(c(a2)*RIE002)+
(1/2)*(a2^(1/2)*a1^(-3/2))*R*(c(a1)*RIE001+c(a2)*RIE002)+c((1-(a2^(1/2)*a1^(-1/2))*R))*(RIE001))
Dha1logP01=colSums((Dha1logP010)-c(Da1logP01)*DhlogP01)######
Dha2logP010=((1-cen1)*cen2/P01)*(c((1/2)*(a2^(-1/2)*a1^(-1/2))*R)*(c(a2)*RIE002)+c((a2^(1/2)*a1^(-1/2))*R)*(RIE002)+
-(1/2)*(a2^(-1/2)*a1^(-1/2))*R*(c(a1)*RIE001+c(a2)*RIE002)+c((1-(a2^(1/2)*a1^(-1/2))*R))*(RIE002))
Dha2logP01=colSums((Dha2logP010)-c(Da2logP01)*DhlogP01)######
DhRlogP010=((1-cen1)*cen2/P01)*(c((a2^(1/2)*a1^(-1/2)))*(c(a2)*RIE002)-(a2^(1/2)*a1^(-1/2))*(c(a1)*RIE001+c(a2)*RIE002))
DhRlogP01=colSums((DhRlogP010)-c(DRlogP01)*DhlogP01)#####
D2boP01=((cen2*(1-cen1))/P01)*( c((a2^(1/2)*a1^(-1/2))*R*a2*H2)*IM2+c(1-(a2^(1/2)*a1^(-1/2))*R)*(a1*H1*IM1+a2*H2*IM2))
D2blogP01=matrix(c(colSums(D2boP01)),n_cov_coef,n_cov_coef)-t(DblogP01)%*%DblogP01######
Dba1logP010=((1-cen1)*cen2/P01)*( c((-1/2)*(a2^(1/2)*a1^(-3/2))*R*a2*H2)*X2+
(1/2)*(a2^(1/2)*a1^(-3/2))*R*(c(a1*H1)*X1+c(a2*H2)*X2)+c((1-(a2^(1/2)*a1^(-1/2))*R))*(c(H1)*X1) )
Dba1logP01=colSums((Dba1logP010)-c(Da1logP01)*DblogP01)######
Dba2logP010=((1-cen1)*cen2/P01)*( c((3/2)*(a2^(1/2)*a1^(-1/2))*R*H2)*X2+
c(-(1/2)*(a2^(-1/2)*a1^(-1/2))*R)*(c(a1*H1)*X1+c(a2*H2)*X2)+c((1-(a2^(1/2)*a1^(-1/2))*R))*(c(H2)*X2))
Dba2logP01=colSums((Dba2logP010)-c(Da2logP01)*DblogP01)######
DbRlogP010=((1-cen1)*cen2/P01)*( c((a2^(1/2)*a1^(-1/2))*a2*H2)*X2-(a2^(1/2)*a1^(-1/2))*(c(a1*H1)*X1+c(a2*H2)*X2))
DbRlogP01=colSums((DbRlogP010)-c(DRlogP01)*DblogP01)####
D2a1logP010=(((1-cen1)*cen2)/P01)*((3/4)*(a1^(-5/2)*a2^(1/2))*R*AA2+(a1^(-3/2)*a2^(1/2)*R)*H1-(3/4)*(a1^(-5/2)*a2^(1/2))*R*AA)
D2a1logP01=sum(D2a1logP010-Da1logP01^2)
Da1a2logP010=(((1-cen1)*cen2)/P01)*( -(1/4)*(a1^(-3/2)*a2^(-1/2))*R*AA2-(1/2)*(a1^(-3/2)*a2^(1/2))*R*H2+
(-(1/2)*(a1^(-1/2)*a2^(-1/2))*R)*H1 +(1/4)*(a1^(-3/2)*a2^(-1/2))*R*AA+(1/2)*(a1^(-3/2)*a2^(1/2))*R*H2)
Da1a2logP01=sum(Da1a2logP010-(Da1logP01*Da2logP01))
D2a2logP010=(((1-cen1)*cen2)/P01)*( -(1/4)*(a1^(-1/2)*a2^(-3/2))*R*AA2+(1/4)*(a1^(-1/2)*a2^(-3/2))*R*AA)
D2a2logP01=sum(D2a2logP010-Da2logP01^2)
Da1RlogP010=(((1-cen1)*cen2)/P01)*( -(1/2)*(a1^(-3/2)*a2^(1/2))*AA2-(a1^(-1/2)*a2^(1/2))*H1 +(1/2)*(a1^(-3/2)*a2^(1/2))*AA)
Da1RlogP01=sum(Da1RlogP010-(Da1logP01*DRlogP01))
Da2RlogP010=(((1-cen1)*cen2)/P01)*( (1/2)*(a1^(-1/2)*a2^(-1/2))*AA2-(1/2)*(a1^(-1/2)*a2^(-1/2))*AA)
Da2RlogP01=sum(Da2RlogP010-(Da2logP01*DRlogP01))
D2RlogP01=-sum(DRlogP01^2)
DthtP01=matrix(c(D2a1logP01,Da1a2logP01,Da1RlogP01,Da1a2logP01,D2a2logP01,Da2RlogP01,Da1RlogP01,Da2RlogP01,D2RlogP01),3,3)###
D2b=D2bp00i+D2blogP11+D2blogP10+D2blogP01
DbDh=Dhbp00i+DhblogP11+DhblogP10+DhblogP01
DbDa=rbind(c(Dba1p00i+Dba1logP11+Dba1logP10+Dba1logP01),c(Dba2p00i+Dba2logP11+Dba2logP10+Dba2logP01),c(DbRp00i+DbRlogP11+DbRlogP10+DbRlogP01))
D2h=(D2hp00+D2hlogP11+D2hlogP10+D2hlogP01)
diag(D2h)<-diag(D2h)-(nev1/nonzero_h0^2)
DhDa=cbind(c(Dha1p00i+Dha1logP11+Dha1logP10+Dha1logP01),c(Dha2p00i+Dha2logP11+Dha2logP10+Dha2logP01),c(DhRp00i+DhRlogP11+DhRlogP10+DhRlogP01))
D2a=Dthtp00i+DthtP11+DthtP10+DthtP01
MA=rbind(D2b,DbDh,DbDa);MB=rbind(t(DbDh),D2h,t(DhDa));MC=rbind(t(DbDa),DhDa,D2a)
HES=-cbind(MA,MB,MC)
RG2=t(HES)
HES[lower.tri(HES, diag = FALSE)]<-RG2[lower.tri(RG2, diag = FALSE)]
INVE<-solve(HES)
seofthet=sqrt(diag(INVE))
SEofbetandsig=c(seofthet[1:n_cov_coef],seofthet[(length(seofthet)-2):length(seofthet)])
list(se=SEofbetandsig,vco=-D2a)}
########## End SE for bcfraildv gamma fit###########
##################################################


########################################################################
############ Functions to obtain the adjusted SE for bcfrailev fit###########

.SEbcfrailphgamcv=function(bet,newtht,n_eve,etime,h0,censor,time,X,H){
beto=bet;beto<-as.vector(beto);n_cov_coef=length(beto);bet=matrix(c(beto),n_cov_coef,1)
n_eve=n_eve;n_eve<-as.vector(n_eve);etime=etime;etime<-as.vector(etime);censor=censor;censor<-as.vector(censor)
h0=h0;h0<-as.vector(h0);time=time;time<-as.vector(time);newtht=newtht;newtht<-as.vector(newtht)
data.n1= length(time);data.n=data.n1/2
a = (newtht[1]);R = (newtht[2])
X <- X;X <-matrix(X,data.n1,n_cov_coef);g0<-c(exp(X%*%bet))
HH=H;HH<-as.vector(HH)
indic1<-2*array(1:data.n)-1;indic2<-2*array(1:data.n)
k0<-c(R*a^(-1));k<-c((1-R)*a^(-1))
c1<-c(R^2+R*a);c2<-c(R-R^2);c3<-c2;c4<-c(1-2*R+R^2)
dak0<-c(-R*a^(-2));dRk0<-c(a^(-1));d2ak0<-c(2*R*a^(-3));daRk0<-c(-a^(-2))
dak<-c(-(1-R)*a^(-2));dRk<-c(-a^(-1));d2ak<-c(2*(1-R)*a^(-3));daRk<-c(a^(-2))
H1=HH[indic1];H2=HH[indic2];time1=time[indic1];time2=time[indic2]
g01=g0[indic1];g02=g0[indic2];cen1=censor[indic1];cen2=censor[indic2]
di=(cen1+cen2)
X1=X[indic1,];X2=X[indic2,]
X1<- matrix(X1,data.n,n_cov_coef);X2<- matrix(X2,data.n,n_cov_coef)
AA=(1+a*(H1+H2));AA1=(1+a*H1);AA2=(1+a*H2)
n_eve0=as.numeric(n_eve>0)
trevntimein=n_eve0*array(1:length(n_eve0))
trevntimein1=trevntimein[trevntimein>0]
trevntime=etime[trevntimein1]
nonzero_h0<-h0[trevntimein1]
nev1<-n_eve[trevntimein1]
rissksetdr1 <- function(time1,x) ifelse(time1[1]>=x,1,0)
RIE1 <- apply(as.array(time1),1,rissksetdr1,x=trevntime)
RIE001 <-(t(RIE1)*c(g01))
rissksetdr2 <- function(time2,x) ifelse(time2[1]>=x,1,0)
RIE2 <- apply(as.array(time2),1,rissksetdr2,x=trevntime)
RIE002 <-(t(RIE2)*c(g02))
interacmatza <- function(X,u){X*u}
resinteracmatx<-apply(as.matrix(X1),2,interacmatza,u=RIE001)
IMX1<-matrix(resinteracmatx,data.n,n_cov_coef*ncol(RIE001))
resinteracmatx<-apply(as.matrix(X2),2,interacmatza,u=RIE001)
IMX12<-matrix(resinteracmatx,data.n,n_cov_coef*ncol(RIE001))
resinteracmatx<-apply(as.matrix(X2),2,interacmatza,u=RIE002)
IMX2<-matrix(resinteracmatx,data.n,n_cov_coef*ncol(RIE002))
resinteracmatx<-apply(as.matrix(X1),2,interacmatza,u=RIE002)
IMX21<-matrix(resinteracmatx,data.n,n_cov_coef*ncol(RIE002))
interacmat1 <- function(X1,u){X1*u}
resinteracmat10<-apply(as.matrix(X1),2,interacmat1,u=X1)
IM1<-matrix(resinteracmat10,data.n,n_cov_coef^2)
interacmat2 <- function(X2,u){X2*u}
resinteracmat20<-apply(as.matrix(X2),2,interacmat2,u=X2)
IM2<-matrix(resinteracmat20,data.n,n_cov_coef^2)
resinteracmat21<-apply(as.matrix(X2),2,interacmat2,u=X1)
IM21<-matrix(resinteracmat21,data.n,n_cov_coef^2)
resinteracmat12<-apply(as.matrix(X1),2,interacmat1,u=X2)
IM12<-matrix(resinteracmat12,data.n,n_cov_coef^2)
Dhp00i=c(-c(k+cen1)*(1/AA1))*(c(a)*RIE001)-c(c(k+cen2)*(1/AA2))*(c(a)*RIE002)-c(c(k0+di)*c(1/AA))*(c(a)*RIE001+c(a)*RIE002)
Dbp00i=(c(-c(k+cen1)*(1/AA1))*(c(a*H1)*X1)+c(-c(k+cen2)*(1/AA2))*(c(a*H2)*X2)+c(-c(k0+di)*(1/AA))*(c(a*H1)*X1+c(a*H2)*X2))
Dap00i= -c(dak)*log(AA1)-c(k+cen1)*(H1/AA1)-c(dak)*log(AA2)-c(k+cen2)*(H2/AA2)-c(dak0)*log(AA)-c(k0+di)*((H1+H2)/AA)
DRp00i=-c(dRk)*log(AA1)-c(dRk)*log(AA2)-c(dRk0)*log(AA)
D2hp00=(t((c(k+cen1)*c(1/AA1^2)*c(a)*RIE001))%*%(c(a)*RIE001)+ t((c(k+cen2)*c(1/AA2^2)*c(a)*RIE002))%*%(c(a)*RIE002)+
t((c(k0+di)*c(1/AA^2)*(c(a)*RIE001+c(a)*RIE002)))%*%(c(a)*RIE001+c(a)*RIE002))#####
Dhbp00i= (matrix(c(colSums(c(-c(k+cen1)*c(1/AA1))*(c(a)*IMX1))),ncol(RIE001),n_cov_coef)+ t(c(a)*RIE001)%*%( c(c(k+cen1)*c(1/AA1^2))*(c(a*H1)*X1))+
matrix(c(colSums(c(-c(k+cen2)*c(1/AA2))*(c(a)*IMX2))),ncol(RIE002),n_cov_coef)+ t(c(a)*RIE002)%*%(c(c(k+cen2)*c(1/AA2^2))*(c(a*H2)*X2))+
matrix(c(colSums(-c(c(k0+di)*c(1/AA))*(c(a)*IMX1+c(a)*IMX2))),ncol(RIE001),n_cov_coef)+
t(c(a)*RIE001+c(a)*RIE002 )%*%( c(c(k0+di)*c(1/AA^2))*(c(a*H1)*X1+c(a*H2)*X2)))####
Dhap00i=colSums(c(-c(dak)*(1/AA1))*(c(a)*RIE001)-c(k+cen1)*((1/AA1)*(RIE001)-(H1/AA1^2)*(c(a)*RIE001))-
c(c(dak)*(1/AA2))*(c(a)*RIE002)-c(k+cen2)*((1/AA2)*(RIE002)-(H2/AA2^2)*(c(a)*RIE002))-
c(c(dak0)*c(1/AA))*(c(a)*RIE001+c(a)*RIE002)-c(k0+di)*(c(1/AA)*(RIE001+RIE002)-c((H1+H2)/AA^2)*(c(a)*RIE001+c(a)*RIE002)))
DhRp00i=colSums(c(-c(dRk)*(1/AA1))*(c(a)*RIE001)-c(c(dRk)*(1/AA2))*(c(a)*RIE002)-c(c(dRk0)*c(1/AA))*(c(a)*RIE001+c(a)*RIE002))#####
DN1=(c(a*H1)*X1);DN2=(c(a*H2)*X2);DN12=(c(a*H1)*X1+c(a*H2)*X2)
DN01=((c(a*H1)*X1)*c(-(k+cen1)/AA1^2));DN02=((c(a*H2)*X2)*c(-(k+cen2)/AA2^2))
DN012=((c(a*H1)*X1+c(a*H2)*X2)*c((k0+di)/AA^2))
DM1=((c(a*H1)*IM1)*c(-(k+cen1)/AA1));DM2=((c(a*H2)*IM2)*c(-(k+cen2)/AA2))
DM12=((c(a*H1)*IM1+c(a*H2)*IM2)*c(-(k0+di)/AA))
D2bp00i=(matrix(c(colSums(DM1)),n_cov_coef,n_cov_coef)-((t(DN1))%*%DN01)+matrix(c(colSums(DM2)),n_cov_coef,n_cov_coef)-((t(DN2))%*%DN02)+
matrix(c(colSums(DM12)),n_cov_coef,n_cov_coef)+((t(DN12))%*%DN012))######
Dbap00i=colSums((c(-c(dak)*(1/AA1))*(c(a*H1)*X1) -c(k+cen1)*( (1/AA1)*(c(H1)*X1)-(H1/AA1^2)*(c(a*H1)*X1)) +
c(-c(dak)*(1/AA2))*(c(a*H2)*X2)-c(k+cen2)*( (1/AA2)*(c(H2)*X2)-(H2/AA2^2)*(c(a*H2)*X2))-c(dak0)*(1/AA)*(c(a*H1)*X1+c(a*H2)*X2)-
c(k0+di)*( (1/AA)*(c(H1)*X1+c(H2)*X2)-((H1+H2)/AA^2)*(c(a*H1)*X1+c(a*H2)*X2))))
DbRp00i=colSums((c(-c(dRk)*(1/AA1))*(c(a*H1)*X1)+c(-c(dRk)*(1/AA2))*(c(a*H2)*X2)+
c(-c(dRk0)*(1/AA))*(c(a*H1)*X1+c(a*H2)*X2)))##
D2ap00i= sum( -c(d2ak)*(log(AA1)+log(AA2))-c(dak)*(2*H1/AA1+2*H2/AA2)+c(k+cen1)*(H1^2/AA1^2)+c(k+cen2)*(H2^2/AA2^2)-
c(d2ak0)*log(AA)-c(dak0)*(2*(H1+H2)/AA)+c(k0+di)*((H1+H2)^2/AA^2))
DaRp00i= sum(-c(daRk)*(log(AA1)+log(AA2))-c(dRk)*(H1/AA1+H2/AA2)-c(daRk0)*log(AA)-c(dRk0)*((H1+H2)/AA))
D2Rp00i=0
Dthtp00i=matrix(c(D2ap00i,DaRp00i,DaRp00i,D2Rp00i),2,2)
P11=( c1*AA1*AA2+c2*AA2*AA+c3*AA1*AA+c4*AA^(2))
DhlogP11=((cen1*cen2)/P11)*( c1*( c(AA2)*(c(a)*RIE001)+c(AA1)*(c(a)*RIE002))+
c2*(c(AA2)*(c(a)*RIE001+c(a)*RIE002)+ c(AA)*(c(a)*RIE002))+c3*(c(AA1)*(c(a)*RIE001+c(a)*RIE002)+ c(AA)*(c(a)*RIE001))+
c4*2*AA*(c(a)*RIE001+c(a)*RIE002))
DblogP11=((cen1*cen2)/P11)*(c(c1*a*H1*AA2)*X1+c(c1*AA1*a*H2)*X2+
c(c2*a*H2*AA)*X2+c(c2*AA2)*(a*H1*X1+a*H2*X2)+c(c3*a*H1*AA)*X1+c(c3*AA1)*(a*H1*X1+a*H2*X2)+c4*2*AA*(a*H1*X1+a*H2*X2))
dalogP11=((cen1*cen2)/P11)*( R*AA1*AA2+c1*(H1*AA2+AA1*H2)+c2*(H2*AA+AA2*(H1+H2))+c3*(H1*AA+AA1*(H1+H2))+2*c4*AA*(H1+H2))
dRlogP11=((cen1*cen2)/P11)*((2*R+a)*AA1*AA2+(1-2*R)*AA2*AA+(1-2*R)*AA1*AA+(-2+2*R)*AA^(2))
D2hlogP110=(t(c(c(cen1*cen2/P11)*c(2)*c1)*(c(a)*RIE001))%*%(c(a)*RIE002)+
t(c(c(cen1*cen2/P11)*c(2)*c2)*(c(a)*RIE002))%*%(c(a)*RIE001+c(a)*RIE002)+
t(c(c(cen1*cen2/P11)*c(2)*c3)*(c(a)*RIE001))%*%(c(a)*RIE001+c(a)*RIE002)+
t(c(c(cen1*cen2/P11)*c(2)*c4)*(c(a)*RIE001+c(a)*RIE002))%*%(c(a)*RIE001+c(a)*RIE002))
D2hlogP11=D2hlogP110-t(DhlogP11)%*%DhlogP11###############
DhblogP110=((cen1*cen2)/P11)*(c1*((c(a*a*H2)*IMX12)+c(AA2)*(c(a)*IMX1)+(c(a*a*H1)*IMX21)+c(AA1)*(c(a)*IMX2))+
c2*((c(a*a*H2)*IMX12+c(a*a*H2)*IMX2)+ c(AA2)*(c(a)*IMX1+c(a)*IMX2)+(c(a*a*H1)*IMX21+c(a*a*H2)*IMX2)+c(AA)*(c(a)*IMX2))+
c3*((c(a*a*H1)*IMX1+c(a*a*H1)*IMX21)+ c(AA1)*(c(a)*IMX1+c(a)*IMX2)+(c(a*a*H1)*IMX1+c(a*a*H2)*IMX12)+c(AA)*(c(a)*IMX1))+
c4*2*((c(a*a*H1)*IMX1+c(a*a*H1)*IMX21+c(a*a*H2)*IMX12+c(a*a*H2)*IMX2)+AA*(c(a)*IMX1+c(a)*IMX2)))
DhblogP11=matrix(c(colSums(DhblogP110)),ncol(RIE001),n_cov_coef)-t(DhlogP11)%*%DblogP11#############
DhalogP110=((cen1*cen2)/P11)*( R*( c(AA2)*(c(a)*RIE001)+c(AA1)*(c(a)*RIE002))+
c1*(c(AA2+a*H2)*RIE001+c(AA1+a*H1)*RIE002)+
c2*(c(AA2+a*H2)*(RIE001+RIE002)+c(AA+a*(H1+H2))*RIE002)+
c3*(c(AA1+a*H1)*(RIE001+RIE002)+c(AA+a*(H1+H2))*RIE001)+
c4*2*c(AA+a*(H1+H2))*(RIE001+RIE002))
DhalogP11=colSums((DhalogP110)-c(dalogP11)*DhlogP11)######
DhRlogP110=((cen1*cen2)/P11)*( (2*R+a)*( c(AA2)*(c(a)*RIE001)+c(AA1)*(c(a)*RIE002))+
(1-2*R)*(c(AA2)*(c(a)*RIE001+c(a)*RIE002)+ c(AA)*(c(a)*RIE002))+(1-2*R)*(c(AA1)*(c(a)*RIE001+c(a)*RIE002)+ c(AA)*(c(a)*RIE001))+
(-2+2*R)*2*AA*(c(a)*RIE001+c(a)*RIE002))
DhRlogP11=colSums((DhRlogP110)-c(dRlogP11)*DhlogP11)######
D2blogP110=((cen1*cen2)/P11)*(c(c1*a)*(c(H1*AA2)*IM1+c(H1*a*H2)*IM21)+ c(c1*a)*(c(a*H1*H2)*IM12+c(AA1*H2)*IM2)+
c(c2*a)*(c(H2*AA)*IM2+c(H2)*(a*H1*IM12+a*H2*IM2))+c(c2*a*H2)*(a*H1*IM21+a*H2*IM2)+c(c2*AA2)*(a*H1*IM1+a*H2*IM2)+
c(c3*a)*(c(H1*AA)*IM1+c(H1)*(a*H1*IM1+a*H2*IM21))+c(c3*a*H1)*(a*H1*IM1+a*H2*IM12)+c(c3*AA1)*(a*H1*IM1+a*H2*IM2)+
c4*2*(c((a*H1)^2)*IM1+c(a*H1*a*H2)*IM21+c(a*H1*a*H2)*IM12+c((a*H2)^2)*IM2)+c4*2*AA*(a*H1*IM1+a*H2*IM2))
D2blogP11=matrix(c(colSums(D2blogP110)),n_cov_coef,n_cov_coef)-t(DblogP11)%*%DblogP11############
DbalogP110=(((cen1*cen2)/P11)*( c(R*a*H1*AA2)*X1+c(R*AA1*a*H2)*X2+c(c1*H1*(AA2+a*H2))*X1+c(c1*(AA1+a*H1)*H2)*X2+
c(c2*H2*(AA+a*(H1+H2)))*X2+c(c2*(AA2+a*H2))*(H1*X1+H2*X2)+c(c3*H1*(AA+a*(H1+H2)))*X1+c(c3*(AA1+a*H1))*(H1*X1+H2*X2)+
c4*2*(AA+a*(H1+H2))*(H1*X1+H2*X2)))
DbalogP11=colSums((DbalogP110)-c(dalogP11)*DblogP11)######
DbRlogP110=((cen1*cen2)/P11)*((2*R+a)*(c(a*H1*AA2)*X1+c(AA1*a*H2)*X2)+
(1-2*R)*(c(a*H2*AA)*X2+c(AA2)*(a*H1*X1+a*H2*X2)+c(a*H1*AA)*X1+c(AA1)*(a*H1*X1+a*H2*X2)) +(-2+2*R)*2*AA*(a*H1*X1+a*H2*X2))
DbRlogP11=colSums((DbRlogP110)-c(dRlogP11)*DblogP11)######
d2alogP110=((cen1*cen2)/P11)*( 2*R*(H1*AA2+AA1*H2)+2*c1*(H1*H2)+2*c2*(H2*(H1+H2))+2*c3*(H1*(H1+H2))+2*c4*(H1+H2)^2)
d2alogP11=sum(d2alogP110-dalogP11^2)
daRlogP110=((cen1*cen2)/P11)*( AA1*AA2+(2*R+a)*(H1*AA2+AA1*H2)+(1-2*R)*((H2*AA+AA2*(H1+H2))+(H1*AA+AA1*(H1+H2)))+2*(-2+2*R)*AA*(H1+H2))
daRlogP11=sum(daRlogP110-(dalogP11*dRlogP11))
d2RlogP110=((cen1*cen2)/P11)*(2*AA1*AA2-2*AA2*AA-2*AA1*AA+2*AA^(2))
d2RlogP11=sum(d2RlogP110-dRlogP11^2)
DthtP11=matrix(c(d2alogP11,daRlogP11,daRlogP11,d2RlogP11),2,2)
P10=R*AA1+(1-R)*AA
DhlogP10=(cen1*(1-cen2)/P10)*(c(R)*(c(a)*RIE001)+c((1-R))*(c(a)*RIE001+c(a)*RIE002))
DblogP10=(cen1*(1-cen2)/P10)*( c(a*R*H1)*X1+c((1-R))*(c(a*H1)*X1+c(a*H2)*X2))
DalogP10=(cen1*(1-cen2)/P10)*(R*H1+(1-R)*(H1+H2))
DRlogP10=(cen1*(1-cen2)/P10)*(AA1-AA)
D2hlogP10=-t(DhlogP10)%*%DhlogP10####
DhblogP100=(cen1*(1-cen2)/P10)*(c(R)*(c(a)*IMX1)+c((1-R))*(c(a)*IMX1+c(a)*IMX2))
DhblogP10=matrix(c(colSums(DhblogP100)),ncol(RIE002),n_cov_coef)- t(DhlogP10)%*%DblogP10######
DhalogP100=(cen1*(1-cen2)/P10)*(c(R)*(RIE001)+c((1-R))*(RIE001+RIE002))
DhalogP10=colSums((DhalogP100)-c(DalogP10)*DhlogP10)######
DhRlogP100=(cen1*(1-cen2)/P10)*((c(a)*RIE001)-(c(a)*RIE001+c(a)*RIE002))
DhRlogP10=colSums((DhRlogP100)-c(DRlogP10)*DhlogP10)######
D2blogP100=(cen1*(1-cen2)/P10)*(c(a*R*H1)*IM1+c((1-R))*(c(a*H1)*IM1+c(a*H2)*IM2))
D2blogP10=matrix(c(colSums(D2blogP100)),n_cov_coef,n_cov_coef)-t(DblogP10)%*%DblogP10######
DbalogP100=(cen1*(1-cen2)/P10)*( c(R*H1)*X1+c((1-R))*(c(H1)*X1+c(H2)*X2))
DbalogP10=colSums(DbalogP100-c(DalogP10)*DblogP10)###
DbRlogP100=(cen1*(1-cen2)/P10)*( c(a*H1)*X1-(c(a*H1)*X1+c(a*H2)*X2))
DbRlogP10=colSums(DbRlogP100-c(DRlogP10)*DblogP10)####
D2alogP10=sum((cen1*(1-cen2)/P10^2)*(-(R*H1+(1-R)*(H1+H2))^2))
DaRlogP100=(cen1*(1-cen2)/P10)*(H1-(H1+H2))
DaRlogP10=sum(DaRlogP100-(DalogP10*DRlogP10))
D2RlogP10=sum((cen1*(1-cen2)/P10^2)*(-(AA1-AA)^2))
DthtP10=matrix(c(D2alogP10,DaRlogP10,DaRlogP10,D2RlogP10),2,2)
P01=R*AA2+(1-R)*AA
DhlogP01=((1-cen1)*cen2/P01)*(c(R)*(c(a)*RIE002)+c((1-R))*(c(a)*RIE001+c(a)*RIE002))
DblogP01=((1-cen1)*cen2/P01)*( c(a*R*H2)*X2+c((1-R))*(c(a*H1)*X1+c(a*H2)*X2))
DalogP01=((1-cen1)*cen2/P01)*(R*H2+(1-R)*(H1+H2))
DRlogP01=((1-cen1)*cen2/P01)*(AA2-AA)
D2hlogP01=-t(DhlogP01)%*%DhlogP01####
DhblogP010=((1-cen1)*cen2/P01)*(c(R)*(c(a)*IMX2)+c((1-R))*(c(a)*IMX1+c(a)*IMX2))
DhblogP01=matrix(c(colSums(DhblogP010)),ncol(RIE002),n_cov_coef)- t(DhlogP01)%*%DblogP01######
DhalogP010=((1-cen1)*cen2/P01)*(c(R)*(RIE002)+c((1-R))*(RIE001+RIE002))
DhalogP01=colSums((DhalogP010)-c(DalogP01)*DhlogP01)######
DhRlogP010=((1-cen1)*cen2/P01)*((c(a)*RIE002)-(c(a)*RIE001+c(a)*RIE002))
DhRlogP01=colSums((DhRlogP010)-c(DRlogP01)*DhlogP01)######
D2blogP010=((1-cen1)*cen2/P01)*(c(a*R*H2)*IM2+c((1-R))*(c(a*H1)*IM1+c(a*H2)*IM2))
D2blogP01=matrix(c(colSums(D2blogP010)),n_cov_coef,n_cov_coef)-t(DblogP01)%*%DblogP01######
DbalogP010=((1-cen1)*cen2/P01)*( c(R*H2)*X2+c((1-R))*(c(H1)*X1+c(H2)*X2))
DbalogP01=colSums(DbalogP010-c(DalogP01)*DblogP01)###
DbRlogP010=((1-cen1)*cen2/P01)*( c(a*H2)*X2-(c(a*H1)*X1+c(a*H2)*X2))
DbRlogP01=colSums(DbRlogP010-c(DRlogP01)*DblogP01)####
D2alogP01=sum(((1-cen1)*cen2/P01^2)*(-(R*H2+(1-R)*(H1+H2))^2))
DaRlogP01=sum(((1-cen1)*cen2/P01)*(H2-(H1+H2))-(DRlogP01*DalogP01))
D2RlogP01=sum(((1-cen1)*cen2/P01^2)*(-(AA2-AA)^2))
DthtP01=matrix(c(D2alogP01,DaRlogP01,DaRlogP01,D2RlogP01),2,2)
D2b=D2bp00i+D2blogP11+D2blogP10+D2blogP01
DbDh=Dhbp00i+DhblogP11+DhblogP10+DhblogP01
DbDa=rbind(c(Dbap00i+DbalogP11+DbalogP10+DbalogP01),c(DbRp00i+DbRlogP11+DbRlogP10+DbRlogP01))
D2h=(D2hp00+D2hlogP11+D2hlogP10+D2hlogP01)
diag(D2h)<-diag(D2h)-(nev1/nonzero_h0^2)
DhDa=cbind(c(Dhap00i+DhalogP11+DhalogP10+DhalogP01),c(DhRp00i+DhRlogP11+DhRlogP10+DhRlogP01))
D2a=Dthtp00i+DthtP11+DthtP10+DthtP01
MA=rbind(D2b,DbDa,DbDh);MB=rbind(t(DbDa),D2a,DhDa);MC=rbind(t(DbDh),t(DhDa),D2h)
HES=-cbind(MA,MB,MC)
RG2=t(HES)
HES[lower.tri(HES, diag = FALSE)]<-RG2[lower.tri(RG2, diag = FALSE)]
INVE<-solve(HES)
seofthet=sqrt(diag(INVE))
SEofbetandsig=c(seofthet[1:(n_cov_coef+2)])
list(se=SEofbetandsig,vco=-D2a)}


adj.SE<-function(fit1){
if(!inherits(fit1, "bcfrailph")){ stop("Argument must be the result of bcfrailph")}
bet<-fit1$coefficients;newtht<-fit1$frailparest;n_eve<-fit1$n.event
etime<-fit1$e.time;H0<-fit1$cbasehaz[,1];h0=diff(c(0,H0))
censor<-fit1$censor;X<-fit1$X;time<-fit1$time
indx=fit1$indx
timeo<-time[indx]
ind.haz=match(time,etime)
if(any(is.na(ind.haz))){
tord_diff<-as.array(diff(c(0,timeo)))
id.zero_tord <- which(apply(((tord_diff<0.0000001)&(tord_diff>0)),1, all))
if(length(id.zero_tord)>0){time[indx[id.zero_tord]]<- time[indx[id.zero_tord-1]]
order=sort(time, decreasing = FALSE, index.return = TRUE);indx=order$ix;timeo=order$x}
etime<-unique(sort(time))
ind.haz=match(time,etime)}
x_bet<-X%*%bet
H=c(H0[ind.haz]*exp(c(x_bet)))
adjj_se=.SEbcfrailphgamcv(bet=bet,newtht=newtht,n_eve=n_eve,etime=etime,h0=h0,
censor=censor,time=time,X=X,H=H)
adjse=c(adjj_se$se)
vcovth<-solve(adjj_se$vco)
if(any(is.na(sqrt(diag(vcovth))))){
warning("Check the approperiatness of the model.")
warning("frailty parameters might be at the boundary of parameter space.")
warning("In obtaining SE Na is produced")}
colnames(vcovth) <- rownames(vcovth) <- c("theta","Row")
fit1$stderr<-adjse
fit1$vcovth2<-vcovth
fit1
}
############ End functions to obtain the adjusted SE for bcfrailev fit###########


###################################################
####### function used in simbcfraildv###########

gener.datadv<-function(p.size,c.rate,fraildistrn,frail.par,bhaz.arg,covar.arg){
n<-p.size; n1=n*2
IID=array(1:n1);indic2=2*array(1:n);indic1=indic2-1
PID=1;e1=array(1:n);PID[indic1]=e1;PID[indic2]=e1
if(fraildistrn==c("gamma")){
lam0=1/(sqrt(frail.par[1])*sqrt(frail.par[2]))
k0=frail.par[3]/(sqrt(frail.par[1])*sqrt(frail.par[2]))
k1=(1/frail.par[1])-k0;k2=(1/frail.par[2])-k0
lam1=(1/frail.par[1]);lam2=(1/frail.par[2])
y0=rgamma(n,shape=k0,scale=1/lam0) ;y1=rgamma(n,shape=k1,scale=1/lam1)
y2=rgamma(n,shape=k2,scale=1/lam2)
z1=((lam0/lam1)*y0+y1);z2=((lam0/lam2)*y0+y2)}
if(fraildistrn==c("lognormal")){
k0=frail.par[3]*(sqrt(frail.par[1])*sqrt(frail.par[2]))
k1=frail.par[1]-k0
k2=frail.par[2]-k0
y0=rnorm(n,mean=0,sd=sqrt(k0))
y1=rnorm(n,mean=0,sd=sqrt(k1))
y2=rnorm(n,mean=0,sd=sqrt(k2))
z1=exp((y0+y1)); z2=exp((y0+y2))}
if(pmatch(c(covar.arg$types),c("BU","BB","UU","B","U"))==1){
x11<-rbinom(n,size=covar.arg$size[1],prob=covar.arg$prob[1])
x21<-runif(n,min=covar.arg$min[1],max=covar.arg$max[1])
sj1cov<-matrix(c(x11,x21),n,2)
x12<-rbinom(n,size=covar.arg$size[1],prob=covar.arg$prob[1])
x22<-runif(n,min=covar.arg$min[1],max=covar.arg$max[1])
sj2cov<-matrix(c(x12,x22),n,2)}
if(pmatch(c(covar.arg$types),c("BU","BB","UU","B","U"))==2){
x11<-rbinom(n,size=covar.arg$size[1],prob=covar.arg$prob[1])
x21<-rbinom(n,size=covar.arg$size[2],prob=covar.arg$prob[2])
sj1cov<-matrix(c(x11,x21),n,2)
x12<-rbinom(n,size=covar.arg$size[1],prob=covar.arg$prob[1])
x22<-rbinom(n,size=covar.arg$size[2],prob=covar.arg$prob[2])
sj2cov<-matrix(c(x12,x22),n,2)}
if(pmatch(c(covar.arg$types),c("BU","BB","UU","B","U"))==3){
x11<-runif(n,min=covar.arg$min[1],max=covar.arg$max[1])
x21<-runif(n,min=covar.arg$min[2],max=covar.arg$max[2])
sj1cov<-matrix(c(x11,x21),n,2)
x12<-runif(n,min=covar.arg$min[1],max=covar.arg$max[1])
x22<-runif(n,min=covar.arg$min[2],max=covar.arg$max[2])
sj2cov<-matrix(c(x12,x22),n,2)}
if(pmatch(c(covar.arg$types),c("BU","BB","UU","B","U"))==4){
x11<-rbinom(n,size=covar.arg$size[1],prob=covar.arg$prob[1])
sj1cov<-matrix(c(x11),n,1)
x12<-rbinom(n,size=covar.arg$size[1],prob=covar.arg$prob[1])
sj2cov<-matrix(c(x12),n,1)}
if(pmatch(c(covar.arg$types),c("BU","BB","UU","B","U"))==5){
x11<-runif(n,min=covar.arg$min[1],max=covar.arg$max[1])
sj1cov<-matrix(c(x11),n,1)
x12<-runif(n,min=covar.arg$min[1],max=covar.arg$max[1])
sj2cov<-matrix(c(x12),n,1)}
u1<-runif(n,  min=0, max=1) # u1 & u2 are random variables represent distribution functions (Cdf)
u2<-runif(n,  min=0, max=1)
if (bhaz.arg$distrn==c("weibull")){
T1 <- (-log(u1) / ((bhaz.arg$scale)*z1*exp(c(sj1cov%*%covar.arg$coefs))))^(1/(bhaz.arg$shape))
T2 <- (-log(u2) / ((bhaz.arg$scale)*z2*exp(c(sj2cov%*%covar.arg$coefs))))^(1/(bhaz.arg$shape))}
if (bhaz.arg$distrn==c("gompertz")){
T1 <- 1/(bhaz.arg$shape)*log(1-(bhaz.arg$shape)*log(u1)/((bhaz.arg$scale)*exp(c(sj1cov%*%covar.arg$coefs))*z1))
T2 <- 1/(bhaz.arg$shape)*log(1-(bhaz.arg$shape)*log(u2)/((bhaz.arg$scale)*exp(c(sj2cov%*%covar.arg$coefs))*z2))}
if (bhaz.arg$distrn==c("exponential")){
T1 <- ((-log(u1)) / ((bhaz.arg$rate)*z1*exp(c(sj1cov%*%covar.arg$coefs))))
T2 <- ((-log(u2)) / ((bhaz.arg$rate)*z2*exp(c(sj2cov%*%covar.arg$coefs))))}
cen1=cen2<-NULL
if(c.rate==0){t1<-T1;t2<-T2;cen1=cen2<-rep(1,n)}
if(c.rate>0){
order=sort(T1, decreasing = FALSE, index.return = TRUE);indx1=order$ix
order=sort(T2, decreasing = FALSE, index.return = TRUE);indx2=order$ix
cr1=round(c(0.05*n),digits=0);I5p=array(1:cr1)
q1005<-quantile(T1, probs = c(0.05));q1099<-quantile(T1, probs = c(0.9999))
q2005<-quantile(T2, probs = c(0.05));q2099<-quantile(T2, probs = c(0.9999))
cen1=rbinom(n, size=1, prob=c(1-c.rate/0.95));cen1[indx1[I5p]]<-1
cen2=rbinom(n, size=1, prob=c(1-c.rate/0.95));cen2[indx2[I5p]]<-1
CC01<-runif(n,  min=q1005, max=q1099);CC1<-c(((q1005+min(CC01))/2),CC01)
CC02<-runif(n,  min=q2005, max=q2099);CC2<-c(((q2005+min(CC02))/2),CC02)
C<-c(CC1,CC2);t1<-T1;t2<-T2;CL=1
for(j in 1:n){if(cen1[j]==0){CL=C[C<T1[j]]; if(length(CL)<=1){t1[j]<-0.95*T1[j]}
if(length(CL)>1){t1[j]<-sample(CL, size=1)}}
if(cen2[j]==0){CL=C[C<T2[j]];if(length(CL)<=1){t2[j]<-0.95*T2[j]}
if(length(CL)>1){t2[j]<-sample(CL, size=1)}}}}
X1<-X2<-time<-censor<-NULL
IID0=array(1:n)
time[indic1]=t1;time[indic2]=t2;censor[indic1]=cen1;censor[indic2]=cen2
if(ncol(sj1cov)==1){X1[indic1]=sj1cov[,1];X1[indic2]=sj2cov[,1]
data<- data.frame(IID=IID,PID=PID,time=time,censor=censor,X1=X1)
DAT<-list(data=data)}
if(ncol(sj1cov)==2){X1[indic1]=sj1cov[,1];X1[indic2]=sj2cov[,1];X2[indic1]=sj1cov[,2];X2[indic2]=sj2cov[,2]
data<- data.frame(IID=IID,PID=PID,time=time,censor=censor,X1=X1,X2=X2)
DAT<-list(data=data)}
DAT}


####### End function used in simbcfraildv###########
###################################################

