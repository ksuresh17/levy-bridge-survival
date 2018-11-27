##########################
#Gamma bridge process (no measurement error)
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
#v: observed value at T_val 
#M=mu_par^2/var_par (reparameterization)
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

seed<-1101

set.seed(seed)

#H(t): conditional cumulative baseline hazard
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

#Simulate censoring times
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

#Generate marker
ZT<-VT
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
#f(v): distribution of V (assumed from same distribution as process, gamma)
fV<-function(v,T_val,mu,var,race)
{
  v<-as.numeric(v)
  T_val<-as.numeric(T_val)
  race<-as.numeric(race)
  
  mu_par<-exp(mu[1]+race*mu[2])
  var_par<-exp(var[1]+race*var[2])
  M<-mu_par^2/var_par
  kappa<-var_par/mu_par
  dgamma(v,shape=M*T_val,scale=kappa)
}

#S(t|Z): conditional survival [Eq.(S2)]
fGcond<-function(v,t,T_val,mu,var,bet,mod,race)
{
  v<-as.numeric(v)
  t<-as.numeric(t)
  T_val<-as.numeric(T_val)
  race<-as.numeric(race)
  
  mu_par<-exp(mu[1]+race*mu[2])
  var_par<-exp(var[1]+race*var[2])
  M<-mu_par^2/var_par
  kappa<-var_par/mu_par
  if(t==0&T_val==0)
  {
    ret<-1 
  } else {
    ret<-Re(kummerM(-cumh(t,bet,mod,race)*v,M*t,M*T_val))
  }
  return(ret)
}

#d/dt S(t|v): derivative of conditional survival distribution 
fderiv_Gcond<-function(v,t,T_val,mu,var,bet,mod,race)
{
  mu_par<-exp(mu[1]+race*mu[2])
  var_par<-exp(var[1]+race*var[2])
  M<-mu_par^2/var_par
  kappa<-var_par/mu_par
  ret<-Re(kummerM(-cumh(t,bet,mod,race)*v,M*t+1,M*T_val+1))*(M*t)/(M*T_val)*(-v*fh(t,bet,mod,race))
  return(ret)
}

#dLambda(t|v): conditional hazard T|Z
fdLambda_Gcond<-function(v,t,T_val,mu,var,bet,mod,race)
{
  v<-as.numeric(v)
  t<-as.numeric(t)
  T_val<-as.numeric(T_val)
  race<-as.numeric(race)
  
  mu_par<-exp(mu[1]+race*mu[2])
  var_par<-exp(var[1]+race*var[2])
  M<-mu_par^2/var_par
  ret<--1/fGcond(v,t,T_val,mu,var,bet,mod,race)*fderiv_Gcond(v,t,T_val,mu,var,bet,mod,race)
  return(ret)
}

##################################################
##################################################
#Estimation (parameters: mu, var, bet)
######################################################
#log-likelihood
######################################################
#param: c(mu,var,bet)
#df: dataset (with headings "Time","z","race")
#mod: model for baseline hazard 
loglik_full<-function(param,df,mod)
{
  dat_names<-names(df)
  t_ind<-which(dat_names=="Time")
  T_ind<-which(dat_names=="Time")
  z_ind<-which(dat_names=="vT")
  race_ind<-which(dat_names=="race")
  
  #Restricting parameters to be positive 
  mu_param<-c(param[1],param[2])
  var_param<-c(param[3],param[4])
  if(mod=="exp") {
    bet_param<-param[c(5,6)] #Exp
  } else if (mod=="weib"|mod=="gamma") {
    bet_param<-c(param[5:6],exp(param[7]))
  }
  
  df_event<-subset(df,df$event==1)
  
  dLambda_cond_event<-NULL
  fz_event<-NULL
  
  dat_event<-apply(df_event,1,function(x) {
    c(fdLambda_Gcond(t=x[t_ind],T_val=x[T_ind],v=x[z_ind],race=x[race_ind],mu=mu_param,var=var_param,bet=bet_param,mod=mod),
      fGcond(t=x[t_ind],T_val=x[T_ind],v=x[z_ind],race=x[race_ind],mu=mu_param,var=var_param,bet=bet_param,mod=mod),
      fV(v=x[z_ind],T_val=x[T_ind],race=x[race_ind],mu=mu_param,var=var_param))
  })
  
  dLambda_cond_event<-dat_event[1,]
  Gcond_event<-dat_event[2,]
  fZ_event<-dat_event[3,]
  
  #Contribution for those who experience the event 
  cont1<-sum(log(dLambda_cond_event)+log(Gcond_event)+log(fZ_event))
  
  df_noevent<-subset(df,df$event==0)
  
  #Contribution for those who do not experience the event 
  dat_noevent<-apply(df_noevent,1,function(x) {
    c(fGcond(t=x[t_ind],T_val=x[T_ind],v=x[z_ind],race=x[race_ind],mu=mu_param,var=var_param,bet=bet_param,mod=mod),
      fV(v=x[z_ind],T_val=x[T_ind],race=x[race_ind],mu=mu_param,var=var_param))
  })
  
  Gcond_noevent<-dat_noevent[1,]
  fZ_noevent<-dat_noevent[2,]
  
  cont2<-sum(log(Gcond_noevent)+log(fZ_noevent))
  
  ret<-cont1+cont2
  
  return(-ret)
}

start.vals<-c(mu_param,var_param,bet_param)

param<-start.vals

opt_est<-try(optim(start.vals,loglik_full,df=fulldat, hessian=FALSE,
                   mod=mod,control=list(maxit=20000,trace=FALSE))) 


opt_est




