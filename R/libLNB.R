## Logit-normal distribution functions ####
# logit <- function(x) log(x/(1-x))

# dlogitnorm<-function(x,m,s) dnorm(logit(x),m,s)/(x*(1-x));
# 

## rlogitnorm from package logitnorm doesn't have parameter n
rlogitnorm<-function(n,m,s) plogis(rnorm(n,m,s))

# logist<-function(x) plogis(x) #1/(1+exp(-x))

# mean.logitnorm<-function(m,s) integrate(function(x,m,s){dlogitnorm(x,m,s)*x},0,1,m=m,s=s)$value
# 
# var.logitnorm<-function(m,s) {
#   mean<-mean.logitnorm(m,s)
#   f<-function(x,m,s){(mean-x)^2*dlogitnorm(x,m,s)}
#   s<-integrate(f,0,1,m=m,s=s)
#   s$value
# }

## Logit-normal binomial distribution functions ####

dlnbinom <- function(x, size, m, s){
  x = as.integer(x)
  size = as.integer(size)
  f=function(P,x,size,m,s) (1-P)^(size-x-1)*P^(x-1)*exp(-(logit(P)-m)^2/(2*s^2))
  F=function(x,size,m,s) {
    if(s>0) return(integrate(f,0,1,x,size,m,s,rel.tol=.Machine$double.eps^0.25,subdivisions=100)$value*choose(size,x)/(sqrt(2*pi*s^2)))
    if(s==0) return(dbinom(x,size,exp(m)))
    if(s<0) return(0)
  }
  #mapply(F,x,size,m,s)
  tryCatch(mapply(F,x,size,m,s), error=function(cond){return(dlnbinom.mc(x,size,m,s, nruns=10000))})
  #tryCatch(mapply(F,x,size,m,s), error=function(cond){warning("IF"); return(NaN)})
}

dlnbinom.mc <- function(x,size,m,s,nruns=100000){
  #  warning(paste("Solving with Monte-Carlo,",nruns,"runs"))
  f = function(x,size,m,s) freq(rlnbinom(nruns,size,m,s),x)/nruns
  return(mapply(f,x,size,m,s))
}

mean.lnbinom<-function(n,m,s) sum(dlnbinom(c(0:n),n,m,s)*c(0:n))

var.lnbinom<-function(n,m,s) {
  mean<-mean.lnbinom(n,m,s)
  #sum(((mean-c(0:n))*dlnbinom(c(0:n),n,m,s))^2)
  mean((mean-c(0:n))^2*dlnbinom(c(0:n),n,m,s)*n)
}


plnbinom<-function(x,size,m,s) sum(dlnbinom(c(0:x),size,m,s))
rlnbinom<-function(n,size,m,s) rbinom(n,size,rlogitnorm(n,m,s))

## Zero-runcated LNB ####
dlnbinom.zt<-function(x,size,m,s){
  d0<-(1-dlnbinom(0,size,m,s))
  (x>0)*dlnbinom(x,size,m,s)/d0
}

dlnbinom.zot<-function(x,size,m,s){
  d0<-(1-sum(dlnbinom(c(1,2),size,m,s)))
  (x>1)*dlnbinom(x,size,m,s)/d0
}

# Logit-normal geometric distribution functions ####

dlngeom<-function(k,m,s) {
  f<-function(P,k,m,s) dgeom(k,P)*dlogitnorm(P,m,s)
  F<-function(k,m,s) integrate(f,0,1,k,m,s)$value
  mapply(F,k,m,s)
}

plngeom<-function(k,m,s) {
  F<-function(k,m,s) 1-dlnbinom(0,k+1,m,s)
  mapply(F,k,m,s)
}

qlngeom<-function(q,m,s){
  f<-function(q,m,s){
    i<-0; 
    while(plngeom(i-1,m,s)<q) {i<-i+1};   
    return(i)
  }
  return(mapply(f,q,m,s))
}

## Hazard function ####

hlngeom<-function(k,m,s) dlngeom(k,m,s)/(1-plngeom(k,m,s))
hgeom<-function(k,p) dgeom(k,p)/(1-pgeom(k,p))

## Maximum Likelihood estimation ####
nloglik.lnbinom.zt<-function(p,n,K){ 
  K<-K[K>0]
  if(p[2]>=0){
    range<-c(c(1:n))
    frq<-freq(K,range)
    filter<-frq>0
    ml<-log(dlnbinom.zt(range[filter],n,p[1],p[2]))
    nloglik<-(-sum(ml*frq[filter]))
    return(nloglik)
  }else{return(-Inf)}
}

fitLNBzt<-function(ms,n,startval=c(-1,2)){
  ## takes a margin sum data set and the size (number of sessions) and fits by ML
  ## for extreme sparse or strongly dispersed data, the fit may fail.
  ## Experiment with the start values (startval) in that case.
  maxLik<-optim(fn=nloglik.lnbinom.zt, par=startval, 
                K=ms, n=n, method="BFGS")
  LNBfit<-list( nlogLik=maxLik$value, mu=maxLik$par[1], sd=maxLik$par[2],
                n=n, discovered=sum(ms>0), ms=ms)
  LNBfit$AIC=2*(-LNBfit$nlogLik)+2*2
  LNBfit <- structure(LNBfit, class = c("LNBfit"))
  return(LNBfit)
}


print.LNBfit <- function(object) {
    dhat <- Dhat(object)
    dnull <- Dnull(object)
    cat("Observations:\n")
    print(data.frame(Trials = object$n,
                   Discovered = object$discovered), 
          digits = 0,
          row.names = F)
    cat("\nCoefficients:\n")
    print(data.frame(mu = object$mu,
                   sd = object$sd), 
          digits = 3,
          row.names = F)
  
    cat("\nProcess:\n")
    print(data.frame(Completion = dhat,
                   Unobserved = dnull), 
          digits = 2,
          row.names = F)
    cat("\nnlogLik: ", object$nlogLik, "AIC: ", object$AIC,"\n")
    cat("\n")
}

## estimated process completion ####
Dhat<-function(LNBfit) 1-dlnbinom(0,LNBfit$n,LNBfit$mu, LNBfit$sd)

## number of undiscovered problems ####
Dnull<- function(LNBfit) LNBfit$discovered/Dhat(LNBfit) - LNBfit$discovered


## Bayesian Estimation ####

## TODO: make prior configurable
lnbzt_posterior <- function(param, K, MS){
  if(param[2]>=0){
    occ = unique(MS)
    frq = freq(MS,occ)
    ml = log(dlnbinom.zt(occ,K,param[1],param[2]))
    loglik = sum(ml*frq)
  }else{loglik = NaN}
#  loglik = -nloglik.lnbinom.zt.fast(param, K, MS)
#  mprior = dnorm(param[1], mean=-2, sd=5, log=T)
#  sprior = dgamma(param[2], 3, 1, log=T)
  mprior = dnorm(param[1], mean=-2, sd=5, log=T)
  sprior = dunif(param[2], 0, 10, log=T)
  loglik + mprior + sprior
}

## helper function
freq<-function(ms, occ = NULL){
  if(is.null(occ)) occ = unique(ms)
  sapply(occ, function(x) sum(ms == x))
}
