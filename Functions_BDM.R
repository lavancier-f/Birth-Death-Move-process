library(spatstat)

#####################
#### Functions for the distance between point patterns
#####################


#Cardinality distance : indicator function for the equality of the cardinalities of two point patterns 
npointsdist=function(x,y){
  #x and y : two point patterns (class ppp of spatstat)
  return(1-as.numeric(x$n==y$n))
}

#Smoothed cardinality distance : absolute distance between the cardinalities of two point patterns
npointsdist2=function(x,y){
  #x and y : two point patterns (class ppp of spatstat)
  return(abs(x$n-y$n))
}

#Hausdorff distance between two point patterns
Hausdorffdist=function(x,y){
  #x and y : two point patterns (class ppp of spatstat)
  res=crossdist(x,y)
  return(max(apply(res,1,min),apply(res,2,min)))
}

#Optimal matching distance : see ?pppdist
optimalmatchingdist=function(x,y){
  #x and y : two point patterns defined on the same window (class ppp of spatstat)
  diam=diameter.owin(x)
  return(pppdist(x,y,cutoff=diam,matching=F))
}

  

#####################
#### Estimation of intensities for continuous-time observations 
#### This is for a pure birth-death process (or neglecting the motions between jumps)
#####################

#Estimation of the jump (or birth or death) intensity (i.e. alpha, beta or delta) of a pure BD process for continuous-time observations 
intensityest=function(deltaT,h,dist,index,K=dnorm){
  #We want to estimate alpha(x) or beta(x) or delta(x) where x is a point pattern 
  #deltaT : vector of inter-jump times (the last value being the residual time T-T_{N_T})
  #h : bandwidth  
  #dist : vector of distances between x and the observed point patterns at 0 and at the jump times (same length as deltaT)
  #index : boolean vector indicating for each jump if it must be taken into account (length=length(deltaT)-1), i.e.: 
  #       TRUE for all of them for the estimation of alpha; 
  #       TRUE only for the births for the estimation of beta;
  #       TRUE only for the deaths for the estimation of delta.
  #K : kernel function
  
  kt=K(dist,0,h)
  hatT=sum(kt*deltaT)
  n=length(kt)
  return(sum(kt[-n][index])/hatT)
}

#log-likelihood used for the choice of h by cross-validation
loglik=function(h,XX,deltaT,distmat,index,K=dnorm){
  #XX : list of observed point patterns at 0 and at the jump times
  #distmat : matrix of distances between the elements of XX
  #h, deltaT, index, K : same entries as for intensityest
  res=NULL 
  n=length(XX)
  for(i in 1:(n-1)){
    #estimation of beta(XX[[i]]) without using XX[[i]]
    res[i]=intensityest(deltaT[-i],h,distmat[i,][-i],index[-i],K=K)
  }
  res[n]=intensityest(deltaT,h,distmat[n,],index,K=K)
  return(-sum(log(res[-n][index]))+sum(deltaT*res))
}

#Choice of the bandwidth by likelihood cross-validation
hoptim=function(XX,deltaT,distmat,index,K=dnorm){
  #same entries as for loglik
  init=quantile(distmat,0.1)/2
  up=max(distmat)
  return(optim(par=init,fn=loglik,XX=XX,deltaT=deltaT,distmat=distmat,index=index,lower=0.0001,upper=up,method="Brent")$par)
}



#####################
#### Estimation of intensities for discrete-time observations
#####################

#Estimation of the jump rate for discrete-time observations
alphaestdisc=function(deltaN,deltat,h,dist,K=dnorm){
  #We want to estimate alpha(x) where x is point pattern 
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
betaestdisc=function(deltaNbirths,deltat,h,dist,K=dnorm){
  #We want to estimate beta(x) where x is point pattern 
  #deltaNbirths : number of observed births between each observed point pattern
  #deltat, h, dist, K : same as for alphaestdisc
  return(alphaestdisc(deltaNbirths,deltat,h,dist,K))
}

#Estimation of the death rate for discrete-time observations
deltaestdisc=function(deltaNdeaths,deltat,h,dist,K=dnorm){
  #We want to estimate delta(x) where x is point pattern 
  #deltaNdeaths : number of observed deaths between each observed point pattern
  #deltat, h, dist, K : same as for alphaestdisc
  return(alphaestdisc(deltaNdeaths,deltat,h,dist,K))
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
    resalphai[i]=alphaestdisc(deltaN[-i],deltat[-i],h,distmat[i,][-i],K=dnorm)
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


