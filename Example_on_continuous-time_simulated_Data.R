#####
## This is an example of how to estimate the intensities of a Birth-death-move process 
## from a continuous-time observed trajectory (simulated Data)
#####
load(file="Simulated_BDM.RData")
#This is the result from the simulation of a BDM process on [0,T] with T=1000 where
# the jump intensity only depends on the cardinalities and is alphan=function(n) exp(5*(n/100-1))
# the birth intensity and the death intensity are alphan/2, the motion of each point between jumps is Brownian.
# This is the simulated trajectory used in the article.
#The environment contains
# XX : the pre-jumps point patterns, in chronological order
# deltaT : vector of inter-jump times (the last value being the residual time T-T_{N_T})  
# dist1, dist2, dist3, dist4 : the matrices of distances between all point patterns in XX for
#   the cardinality distance, the smoothed cardinality distance, the Hausdorff distance, the optimal matching distance, respectively


#Load the functions
source("Functions_BDM.R")
n=length(deltaT)

# Construction of the distance matrices between all pre-jump configurations
# NOT RUN (already saved in Simulated_BDM.RData)
#
# dist1=matrix(0,n,n)
# dist2=matrix(0,n,n)
# for(i in 1:(n-1)){cat(i)
#   for(j in (i+1):n){
#     dist1[i,j]<-npointsdist(XX[[i]],XX[[j]])
#     dist2[i,j]<-npointsdist2(XX[[i]],XX[[j]])
#   }
# }
# dist1=t(dist1)+dist1
# dist2=t(dist2)+dist2
# 
# 
# dist3=matrix(0,n,n)
# for(i in 1:(n-1)){cat(i)
#   for(j in (i+1):n){
#     dist3[i,j]<-Hausdorffdist(XX[[i]],XX[[j]])
#   }
# }
# dist3=t(dist3)+dist3
# 
# dist4=matrix(0,n,n)
# diam=diameter.owin(XX[[1]])
# for(i in 1:(n-1)){cat(i)
#   for(j in (i+1):n){
#     dist4[i,j]<-optimalmatchingdist(XX[[i]],XX[[j]])
#   }
# }
# dist4=t(dist4)+dist4


#True intensity:
nn=unlist(lapply(XX,npoints))
alphan=function(n) exp(5*(n/100-1))
plot(sapply(nn,alphan),type='l')

index=rep(TRUE,n-1) #indicates which jump counts for the estimation (all of them here since we are estimating the total intensity) 

#Estimation based on the number of points (without smoothing)
resalpha1=NULL
for(i in 1:n){
  resalpha1[i]=intensityest(deltaT,h=0,dist1[i,],index,K=function(t,a,b) 1-t) #h is useless
}
plot(sapply(nn,alphan),type='l')
lines(resalpha1,col=2)

#Estimation based on the number of points, but through a smoothing kernel
h2=hoptim(XX,deltaT,dist2,index,K=dnorm) #choice of the bandwidth by likelihood CV
resalpha2=NULL
for(i in 1:n){
  resalpha2[i]=intensityest(deltaT,h=h2,dist2[i,],index)
}
plot(sapply(nn,alphan),type='l')
lines(resalpha2,col=2)


#Estimation using the Hausdorff distance
h3=hoptim(XX,deltaT,dist3,index,K=dnorm) #choice of the bandwidth by likelihood CV
resalpha3=NULL
for(i in 1:n){
  resalpha3[i]=intensityest(deltaT,h=h3,dist3[i,],index)
}
plot(sapply(nn,alphan),type='l')
lines(resalpha3,col=2)


#Estimation using the optimal matching distance
h4=hoptim(XX,deltaT,dist4,index,K=dnorm) #choice of the bandwidth by likelihood CV
resalpha4=NULL
for(i in 1:n){
  resalpha4[i]=intensityest(deltaT,h=h4,dist4[i,],index)
}
plot(sapply(nn,alphan),type='l')
lines(resalpha4,col=2)

#All plots together
plot(sapply(nn,alphan),type='l') #ground truth
lines(resalpha1,col=2) #Estimation from the cardinality distance
lines(resalpha2,col=3) #Estimation from the smoothed cardinality distance
lines(resalpha3,col=4) #Estimation from the Hausdorff distance
lines(resalpha4,col=5) #Estimation from the d_k distance

#Intensities wrt the number of points
nn=unlist(lapply(XX,npoints))
plot(nn,resalpha3,col=4) #from the Hausdorff distance
points(nn,resalpha1,col=2) #from the cardinality distance
points(nn,resalpha2,col=3) #from the smoothed cardinality distance
points(nn,resalpha4,col=5) #from the d_k distance
lines(nn,sapply(nn,alphan)) #ground truth
