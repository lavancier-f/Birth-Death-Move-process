library(spatstat)

#####################
#### Functions for the distance between point patterns
#####################

#Hausdorff distance between two point patterns
Hausdorffdist=function(x,y){
  #x and y : two point patterns (class ppp of spatstat)
  res=crossdist(x,y)
  return(max(apply(res,1,min),apply(res,2,min)))
}

#Cardinality distance : indicator function for the equality of the cardinalities of two point patterns 
npointsdist=function(x,y){
  #x and y : two point patterns (class ppp of spatstat)
  return(1-as.numeric(x$n==y$n))
}

#Smoothed cardinality distance : absolute distance between the cardinalities of two point patterns
npointsdist2=function(x,y){
  return(abs(x$n-y$n))
}



#####################
#### Estimation of intensities for continuous-time observations
#####################

#Estimation of the jump rate for continuous-time observations
alphaest=function(x,deltaT,h,dist,K=dnorm){
  #x : we want to estimate alpha(x) where x is point pattern (class ppp of spatstat)
  #deltaT : vector of inter-jump times (the last value being the residual time T-T_{N_T})
  #h : bandwidth  
  #dist : vector of distances between x and the observed point patterns at the jump times (same length as deltaT)
  #K : kernel function
  
  kt=K(dist,0,h)
  hatT=sum(kt*deltaT)
  n=length(kt)
  return(sum(kt[-n])/hatT)
}

#Estimation of the birth rate for continuous-time observations
betaest=function(x,deltaTbirths,h,distbirths,K=dnorm){
  #deltaTbirths : vector of inter-birth times (the last value being the residual time)
  #distbirths : vector of distances between x and the observed point patterns at the birth times (same length as deltaTbirths)
  #x, h, K : same as for alphaest
  return(alphaest(x,deltaTbirths,h,distbirths,K))
}

#Estimation of the death rate for continuous-time observations
deltaest=function(x,deltaTdeaths,h,distdeaths,K=dnorm){
  #deltaTdeaths : vector of inter-death times (the last value being the residual time)
  #distdeaths : vector of distances between x and the observed point patterns at the death times (same length as deltaTdeaths)
  #x, h, K : same as for alphaest
  return(alphaest(x,deltaTdeaths,h,distdeaths,K))}



#log-likelihood used for the choice of h by cross-validation
loglik=function(h,XX,deltaT,distmat,K=dnorm){
  #XX : list of observed point patterns at the jump times (of length n, say)
  #distmat : matrice of distances between the observed point patterns at the jump times
  #h, deltaT, K : same entries as alphaest or betaest or deltaest
  resalphai=NULL 
  n=length(XX)
  for(i in 1:n){
    #estimation of alpha(XX[[i]]) without using XX[[i]]
    resalphai[i]=alphaest(XX[[i]],deltaT[-i],h,distmat[i,][-i],K=K)
  }
  return(-sum(log(resalphai))+sum(deltaT*resalphai))
}

#Choice of the bandwidth by likelihood cross-validation
hoptim=function(XX,deltaT,distmat,K=dnorm){
  #same entries as for loglik
  init=quantile(distmat,0.1)/2
  up=max(distmat)
  return(optim(par=init,fn=loglik,XX=XX,deltaT=deltaT,distmat=distmat,lower=0.0001,upper=up,method="Brent")$par)
}



#####################
#### Estimation of intensities for discrete-time observations
#####################

#Estimation of the jump rate for discrete-time observations
alphaestdisc=function(x,deltaN,deltat,h,dist,K=dnorm){
  #x : we want to estimate alpha(x) where x is point pattern (class ppp of spatstat)
  #deltaN : number of observed jumps between each observed point pattern 
  #deltat : vector of inter-times (of the same length as deltaN, but if they are all equal one value is enough)
  #h : bandwidth  
  #dist : vector of distances between x and the observed point patterns (length=length(deltat)+1)
  #K : kernel function
  
  kt=K(dist,0,h)
  n=length(kt)
  hatT=sum(kt[-n]*deltat)
  return(sum(deltaN*kt[-n])/hatT)
}

#Estimation of the birth rate for discrete-time observations
betaestdisc=function(x,deltaNbirths,deltat,h,dist,K=dnorm){
  #deltaNbirths : number of observed births between each observed point pattern
  #x, deltat, h, dist, K : same as for alphaestdisc
  return(alphaestdisc(x,deltaNbirths,deltat,h,dist,K))
}

#Estimation of the death rate for discrete-time observations
deltaestdisc=function(x,deltaNdeaths,deltat,h,dist,K=dnorm){
  #deltaNdeaths : number of observed deaths between each observed point pattern
  #x, deltat, h, dist, K : same as for alphaestdisc
  return(alphaestdisc(x,deltaNdeaths,deltat,h,dist,K))
}


#log-likelihood used for the choice of h by cross-validation
loglikdisc=function(h,XX,deltaN,deltat,distmat,K=dnorm){
  #XX : list of observed point patterns 
  #distmat : matrice of distances between the observed point patterns
  #h, deltaN, deltat, K : same entries as for alphaestdisc or betaestdisc or deltaestdisc
  n=length(XX)
  if(length(deltat)==1){deltat=rep(deltat,n-1)}
  resalphai=NULL 
  for(i in 1:(n-1)){
    #estimation of alpha(XX[[i]]) without using XX[[i]]
    resalphai[i]=alphaestdisc(XX[[i]],deltaN[-i],deltat[-i],h,distmat[i,][-i],K=dnorm)
  }
  return(-sum(deltaN*log(resalphai))+sum(deltat*resalphai))
}

#Choice of the bandwidth by likelihood cross-validation
hoptimdisc=function(XX,deltaN,deltat,distmat,K=dnorm){
  #same entries as for loglikdisc
  init=quantile(distmat,0.1)/2
  up=max(distmat)
  return(optim(par=init,fn=loglikdisc,XX=XX,deltaN=deltaN,deltat=deltat,distmat=distmat,lower=0.0001,upper=up,method="Brent")$par)
}


