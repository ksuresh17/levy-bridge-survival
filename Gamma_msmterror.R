##########################
#Gamma bridge process (with measurement error) 
##########################
rm(list=ls())
#Pacakges
library(numDeriv) #grad function
library(fAsianOptions) #kummerM function 

#######################################
#mu_param: parameters in mean of process
#var_param: parameters in variance of process
#bet: coefficients in conditional cumulative baseline hazard H(t) [additional parameter if Gamma, Weibull]
#race: covariate 
#mod: model for H(t) ["exp","weib","gamma"]
#T_val: time of last observed measurement 
#v: true value of marker at T_val 
#z: observed value of marker at T_val (with measurement error)
#M=mu_par^2/var_par (reparameterization)
#g0: [gamma in paper] measurement error parameter
#######################################

###########################
#Simulation setup 
###########################
tolerance<-1e-10
k_max<-100
n_max<-10000

mu_param<-c(-1.1,0.3)
var_param<-c(-2.1,0.5)
bet_param<-c(-3.6,0.6)
mod<-"exp"

n<-200
cens_horiz<-20
g0<-20

seed<-1101

set.seed(seed)

#H(t): cumulative baseline hazard
cumh<-function(t,bet,mod,race)
{
  exp_param<-exp(bet[1]+race*bet[2])
  if(mod=="exp")
  {
    Ht<- exp_param*t #Exp
  } else if(mod=="weib") {
    Ht<- (t/exp_param)^(bet[3]) #Weibull
  } else if(mod=="gamma")
  {
    Ht<- -log(1-(gamma(bet[3])-gammainc(bet[3],exp_param*t))/gamma(bet[3])) #Gamma
  }
  return(Ht)
}

#Generate race 
R_sim<-rbinom(n,1,0.30)

#Simulate gamma process (Avramidis, et al. "Efficient simulation of Gamma and Variance-Gamma processes")
dat<-NULL 
k<--log(.01/(cens_horiz+1))/log(2) #2^k time values
h=2^(-k)*(cens_horiz+1) #h: length of time increment
dat_times<-seq(0,cens_horiz+1,by=h) #times at which to generate process 
for(j in 1:n)
{
  vec<-c(0,1:(2^k))
  mu_sim<-exp(mu_param[1]+R_sim[j]*mu_param[2])
  var_sim<-exp(var_param[1]+R_sim[j]*var_param[2])
  for(i in 2:((2^k)+1))
  {
    Q<-rgamma(n=1,shape=mu_sim^2*h/var_sim,scale=var_sim/mu_sim)
    vec[i]<-vec[i-1]+Q
  }
  dat<-rbind(dat,vec)
}

#Simulate survival times
froot<-function(u,t,Wdat,bet,race)
{
  Wtime_ind<-which(abs(dat_times-round(t,2))<10e-5)
  1-exp(-Wdat[Wtime_ind]*cumh(t,bet,mod,race))-u
}

u<-runif(n,0,1)
ETimes<-NULL
for(i in 1:n)
{
  rootval<-try(uniroot(froot,u=u[i],Wdat=dat[i,],bet=bet_param,race=R_sim[i],
                       interval=c(1e-14,cens_horiz+1),extendInt="yes")$root,TRUE)
  ETimes<-c(ETimes,ifelse(inherits(rootval, "try-error"),Inf,rootval))
}


CTimes<-runif(n,0,cens_horiz)
Time<-pmin(ETimes,CTimes)
event<-ifelse(Time==ETimes,1,0)

#Generate V_T (last value of Gamma process)
VT_ind<-unlist(sapply(Time,function(x) which(abs(dat_times-round(x,2))<10e-5)))
VT<-NULL
for(i in 1:n)
{
  VT<-c(VT,dat[i,VT_ind[i]])
}

VT[VT==0]<-10e-100

#Generate observed marker value (with measurement error)
ZT<-rgamma(n,shape=g0,scale=VT/g0)
ZT[ZT<10e-5]<-10e-5

#Full dataset
fulldat<-data.frame("id"=1:n,"Time"=Time,"event"=event,"z"=ZT,"vT"=VT,"race"=R_sim)

######################################################
#Unconditional survival and hazard functions
######################################################
#h(t): baseline hazard 
fh<-function(t,bet,mod,race)
{
  exp_param<-exp(bet[1]+race*bet[2])
  if(mod=="exp")
  {
    ht<- exp_param #Exp
  } else if(mod=="weib") {
    ht<- (bet[3]/exp_param)*(t/exp_param)^(bet[3]-1) #Weibull
  } else if(mod=="gamma") {
    ht<- exp_param*(exp_param*t)^(bet[3]-1)*exp(-exp_param*t)/(gamma(bet[3])-(gamma(bet[3])-gammainc(bet[3],exp_param*t))) #Gamma
  } 
  
  return(ht)
}

#S(t): unconditional survival function 
fS<-function(t,mu,var,bet,mod,race)
{
  mu_par<-exp(mu[1]+race*mu[2])
  var_par<-exp(var[1]+race*var[2])
  M<-mu_par^2/var_par
  kappa<-var_par/mu_par
  
  S_t<-(1+kappa*cumh(t,bet,mod,race))^(-M*t)
  return(S_t)
}

#dLambda(t): unconditional hazard 
fdLambda<-function(t,mu,var,bet,mod,race)
{
  mu_par<-exp(mu[1]+race*mu[2])
  var_par<-exp(var[1]+race*var[2])
  M<-mu_par^2/var_par
  kappa<-var_par/mu_par
  
  ret<-M*log(1+kappa*cumh(t,bet,mod,race))+kappa*M*t/(1+kappa*cumh(t,bet,mod,race))*fh(t,bet,mod,race)
  return(ret)
}

#Lambda(t): unconditional cumulative hazard
fLambda<-function(t,mu,var,bet,race)
{
  ret<--log(fS(t,mu,var,bet,race))
  return(ret)
}

######################################################
#Derivation of functions using closed-form solutions 
######################################################
#f(Z): distribution of Z
fZ<-function(z,T_val,mu,var,g0,race)
{
  z<-as.numeric(z)
  T_val<-as.numeric(T_val)
  race<-as.numeric(race)
  
  mu_par<-exp(mu[1]+race*mu[2])
  var_par<-exp(var[1]+race*var[2])
  M<-mu_par^2/var_par
  kappa<-var_par/mu_par
  
  if(2/exp(lgamma(g0))==0) {ret<-0} else {
    ret<-2/exp(lgamma(g0))*1/exp(lgamma(M*T_val))*1/z*(g0*z/kappa)^((g0+M*T_val)/2)*besselK(2*sqrt(g0*z/kappa),M*T_val-g0)
  }
  return(ret)
}

#E_V[S(t|V)*f(z|V)]: numerator of conditional survival function G(t|z)
numG<-function(z,t,T_val,mu,var,g0,bet,mod,race)
{
  z<-as.numeric(z)
  t<-as.numeric(t)
  T_val<-as.numeric(T_val)
  race<-as.numeric(race)
  
  mu_par<-exp(mu[1]+race*mu[2])
  var_par<-exp(var[1]+race*var[2])
  M<-mu_par^2/var_par
  kappa<-var_par/mu_par
  
  k<-seq(0,k_max)
  sum_temp<-exp(lgamma(M*t+k))/exp(lgamma(M*t))*exp(lgamma(M*T_val))/exp(lgamma(M*T_val+k))*1/exp(lgamma(k+1))*((-1)*cumh(t,bet,mod,race))^k*2/exp(lgamma(g0))*1/exp(lgamma(M*T_val))*1/z*(g0*z/kappa)^((g0+M*T_val)/2)*(g0*z*kappa)^(k/2)*besselK(2*sqrt(g0*z/kappa),M*T_val+k-g0)
  sum_temp<-sum_temp[!is.na(sum_temp)]
  count_it<-1
  check<-ifelse(length(sum_temp)>2,abs(sum(sum_temp)-sum(sum_temp[1:(length(sum_temp)-1)])),NA)
  while(!is.na(check)&check>tolerance&count_it<10)
  {
    count_it<-count_it+1
    ind_k<-length(k)
    k<-seq(ind_k+1,ind_k+100)
    sum_temp<-c(sum_temp,exp(lgamma(M*t+k))/exp(lgamma(M*t))*exp(lgamma(M*T_val))/exp(lgamma(M*T_val+k))*1/exp(lgamma(k+1))*((-1)*cumh(t,bet,mod,race))^k*2/exp(lgamma(g0))*1/exp(lgamma(M*T_val))*1/z*(g0*z/kappa)^((g0+M*T_val)/2)*(g0*z*kappa)^(k/2)*besselK(2*sqrt(g0*z/kappa),M*T_val+k-g0))
    sum_temp<-sum_temp[!is.na(sum_temp)]
    check<-ifelse(length(sum_temp)>2,abs(sum(sum_temp)-sum(sum_temp[1:(length(sum_temp)-1)])),NA)
  }
  ret<-sum(sum_temp)
  ret<-ifelse(check>tolerance|is.na(check),NA,ret)
  
  return(ret)
}

#G(t|Z): conditional survival 
fGcond<-function(z,t,T_val,mu,var,g0,bet,mod,race)
{
  z<-as.numeric(z)
  t<-as.numeric(t)
  T_val<-as.numeric(T_val)
  race<-as.numeric(race)
  
  mu_par<-exp(mu[1]+race*mu[2])
  var_par<-exp(var[1]+race*var[2])
  M<-mu_par^2/var_par
  kappa<-var_par/mu_par
  
  if(t==0&T_val==0) {
    ret<-1
  } else {
    k<-seq(1,k_max)
    sum_temp<-1/k*1/beta(M*t,k)*exp(lgamma(M*T_val))/exp(lgamma(M*T_val+k))*((-1)*cumh(t,bet,mod,race)*sqrt(g0*z*kappa))^k*besselK(2*sqrt(g0*z/kappa),M*T_val+k-g0)/besselK(2*sqrt(g0*z/kappa),M*T_val-g0)
    sum_temp<-sum_temp[!is.na(sum_temp)]
    count_it<-1
    check<-ifelse(length(sum_temp)>2,abs(sum(sum_temp)-sum(sum_temp[1:(length(sum_temp)-1)])),NA)
    while(!is.na(check)&check>tolerance&count_it<10)
    {
      count_it<-count_it+1
      ind_k<-length(k)
      k<-seq(ind_k+1,ind_k+100)
      sum_temp<-c(sum_temp,1/k*1/beta(M*t,k)*exp(lgamma(M*T_val))/exp(lgamma(M*T_val+k))*((-1)*cumh(t,bet,mod,race)*sqrt(g0*z*kappa))^k*besselK(2*sqrt(g0*z/kappa),M*T_val+k-g0)/besselK(2*sqrt(g0*z/kappa),M*T_val-g0))
      sum_temp<-sum_temp[!is.na(sum_temp)]
      check<-ifelse(length(sum_temp)>2,abs(sum(sum_temp)-sum(sum_temp[1:(length(sum_temp)-1)])),NA)
    }
    ret<-1+sum(sum_temp)
    ret<-ifelse(check>tolerance|is.na(check)|ret>1,0,ret)
  }
  return(ret)
}

#d/dt G(t|Z): derivative of conditional survival distribution 
fderiv_Gcond<-function(z,t,T_val,mu,var,g0,bet,mod,race)
{
  z<-as.numeric(z)
  t<-as.numeric(t)
  T_val<-as.numeric(T_val)
  race<-as.numeric(race)
  
  mu_par<-exp(mu[1]+race*mu[2])
  var_par<-exp(var[1]+race*var[2])
  M<-mu_par^2/var_par
  kappa<-var_par/mu_par
  
  k<-seq(1,k_max)
  temp1<-1/k*1/beta(M*t,k)*exp(lgamma(M*T_val))/exp(lgamma(M*T_val+k))*((-1)*cumh(t,bet,mod,race)*sqrt(g0*z*kappa))^k*besselK(2*sqrt(g0*z/kappa),M*T_val+k-g0)/besselK(2*sqrt(g0*z/kappa),M*T_val-g0)
  temp2<-M*(digamma(M*t+k)-digamma(M*t))+k*fh(t,bet,mod,race)/cumh(t,bet,mod,race)
  sum_temp<-temp1*temp2
  sum_temp<-sum_temp[!is.na(sum_temp)]
  count_it<-1
  check<-ifelse(length(sum_temp)>2,abs(sum(sum_temp)-sum(sum_temp[1:(length(sum_temp)-1)])),NA)
  while(!is.na(check)&check>tolerance&count_it<10)
  {
    count_it<-count_it+1
    ind_k<-length(k)
    k<-seq(ind_k+1,ind_k+100)
    temp1<-1/k*1/beta(M*t,k)*exp(lgamma(M*T_val))/exp(lgamma(M*T_val+k))*((-1)*cumh(t,bet,mod,race)*sqrt(g0*z*kappa))^k*besselK(2*sqrt(g0*z/kappa),M*T_val+k-g0)/besselK(2*sqrt(g0*z/kappa),M*T_val-g0)
    temp2<-M*(digamma(M*t+k)-digamma(M*t))+k*fh(t,bet,mod,race)/cumh(t,bet,mod,race)
    temp<-temp1*temp2
    sum_temp<-c(sum_temp,temp)
    sum_temp<-sum_temp[!is.na(sum_temp)]
    check<-ifelse(length(sum_temp)>2,abs(sum(sum_temp)-sum(sum_temp[1:(length(sum_temp)-1)])),NA)
  }
  ret<-sum(sum_temp)
  ret<-ifelse(check>tolerance|is.na(check),0,ret)
  return(ret)
}

#dLambda(t|z): conditional hazard T|Z
fdLambda_Gcond<-function(z,t,T_val,mu,var,g0,bet,mod,race)
{
  z<-as.numeric(z)
  t<-as.numeric(t)
  T_val<-as.numeric(T_val)
  race<-as.numeric(race)
  
  if(fderiv_Gcond(z,t,T_val,mu,var,g0,bet,mod,race)==0){
    ret<-0
  } else {
    ret<--1/fGcond(z,t,T_val,mu,var,g0,bet,mod,race)*fderiv_Gcond(z,t,T_val,mu,var,g0,bet,mod,race)
  }
  return(ret)
}

######################################################
#Derivation of functions using numerical integration/differentiation 
######################################################
# #Checked and gives same as Rcond 
fZ_int<-function(z,T_val,mu,var,g0,race)
{
  z<-as.numeric(z)
  T_val<-as.numeric(T_val)
  race<-as.numeric(race)
  
  mu_par<-exp(mu[1]+race*mu[2])
  var_par<-exp(var[1]+race*var[2])
  M<-mu_par^2/var_par
  kappa<-var_par/mu_par
  
  integrate(function(v,z,mu,var,g0,T_val,M,kappa) {
    dgamma(z,shape=g0,scale=v/g0)*dgamma(v,shape=M*T_val,scale=kappa)
  },z=z,g0=g0,mu=mu,var=var,T_val=T_val,M=M,kappa=kappa,lower=0,upper=Inf,stop.on.error = FALSE)$val
}

#E_V[S(t|V)*f(z|V)]: numerator of conditional survival function G(t|z)
numG_int<-function(z,t,T_val,mu,var,g0,bet,mod,race)
{
  z<-as.numeric(z)
  T_val<-as.numeric(T_val)
  race<-as.numeric(race)
  
  mu_par<-exp(mu[1]+race*mu[2])
  var_par<-exp(var[1]+race*var[2])
  M<-mu_par^2/var_par
  kappa<-var_par/mu_par
  
  integrate(function(v,z,t,T_val,mu,var,g0,bet,race,M,kappa,mod) {
    Re(kummerM(x=(-1)*v*cumh(t,bet,mod,race),a=M*t,b=M*T_val))*dgamma(z,shape=g0,scale=v/g0)*dgamma(v,shape=M*T_val,scale=kappa)
  },z=z,t=t,mu=mu,var=var,bet=bet,g0=g0,T_val=T_val,race=race,mod=mod,M=M,kappa=kappa,lower=0,upper=Inf,rel.tol = 1e-9,stop.on.error = FALSE)$val
}

#G(t|Z): conditional survival
fGcond_int<-function(z,t,T_val,mu,var,g0,bet,mod,race)
{
  numG_int(z,t,T_val,mu,var,g0,bet,mod,race)/fZ_int(z,T_val,mu,var,g0,race)
}

#Lambda(t|z): conditional cumulative hazard T|Z
fLambda_Gcond_int<-function(z,t,T_val,mu,var,g0,bet,mod,race)
{
  ret<--log(fGcond_int(z,t,T_val,mu,var,g0,bet,mod,race))
  return(ret)
}

#dLambda(t|z): conditional hazard T|Z
fdLambda_Gcond_int<-function(z,t,T_val,mu,var,g0,bet,mod,race)
{
  ret<-grad(fLambda_Gcond_int,t,T_val=T_val,z=z,mu=mu,var=var,bet=bet,g0=g0,mod=mod,race=race)
  return(ret)
}

##################################################
##################################################
#Estimation (parameters: mu, var, g0, bet)
######################################################
#log-likelihood
######################################################
#param: c(mu,var,g0,bet)
#df: dataset (with headings "Time","z","race")
#mod: model for baseline hazard 
loglik_full<-function(param,df,mod)
{
  dat_names<-names(df)
  t_ind<-which(dat_names=="Time")
  T_ind<-which(dat_names=="Time")
  z_ind<-which(dat_names=="z")
  race_ind<-which(dat_names=="race")
  
  df_event<-subset(df,df$event==1)
  df_noevent<-subset(df,df$event==0)
  
  #Restricting parameters to be positive 
  mu_param<-c(param[1],param[2])
  var_param<-c(param[3],param[4])
  g0_param<-exp(param[5])
  # g0_param<-param[5]*100
  if(mod=="exp") {
    bet_param<-param[c(6,7)] #Exp
  } else if (mod=="weib"|mod=="gamma") {
    bet_param<-c(param[6:7],exp(param[8]))
  }
  
  #Contribution for those who experience the event 
  dat_event<-apply(df_event,1,function(x) {
    c(fdLambda_Gcond(t=x[t_ind],T_val=x[T_ind],z=x[z_ind],mu=mu_param,var=var_param,bet=bet_param,g0=g0_param,mod=mod,race=x[race_ind]),
      fGcond(t=x[t_ind],T_val=x[T_ind],z=x[z_ind],mu=mu_param,var=var_param,bet=bet_param,g0=g0_param,mod=mod,race=x[race_ind]),
      fZ(z=x[z_ind],T_val=x[T_ind],mu=mu_param,var=var_param,g0=g0_param,race=x[race_ind]))
  })
  
  dLambda_cond_event<-dat_event[1,]
  Gcond_event<-dat_event[2,]
  fZ_event<-dat_event[3,]
  
  cont1<-sum(log(dLambda_cond_event)+log(Gcond_event)+log(fZ_event))
  
  #Contribution for those who do not experience the event 
  dat_noevent<-apply(df_noevent,1,function(x) {
    c(fGcond(t=x[t_ind],T_val=x[T_ind],z=x[z_ind],mu=mu_param,var=var_param,bet=bet_param,g0=g0_param,mod=mod,race=x[race_ind]),
      fZ(z=x[z_ind],T_val=x[T_ind],mu=mu_param,var=var_param,g0=g0_param,race=x[race_ind]))
  })
  
  Gcond_noevent<-dat_noevent[1,]
  fZ_noevent<-dat_noevent[2,]
  
  cont2<-sum(log(Gcond_noevent)+log(fZ_noevent))
  
  ret<-cont1+cont2
  
  return(-ret)
}

start.vals<-c(mu_param,var_param,log(g0),bet_param) #starting values

param<-start.vals

opt_est<-try(optim(start.vals,loglik_full,df=fulldat, hessian=FALSE,
                   mod=mod,control=list(maxit=20000,trace=FALSE)))

opt_est


