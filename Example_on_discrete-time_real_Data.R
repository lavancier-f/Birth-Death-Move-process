#####
## Estimation of the birth rate and the death rate on real data.
## They correspond to the locations of two type of proteins (Langerin and Rab11 proteins) in a cell, observed in discrete time.
## The data for each type of proteins are recorded in a separate .RData file
## These are the data analysed in the article.
#####

#Load the functions
source("Functions_BDM.R")


### I Estimation of intensities for Langerin data
load(file="Utrack_Langerin.RData")
#This environment contains :
# XXdisc : the observed point patterns, in chronological order
# deltaNbirth : the observed number of births between each observed point pattern 
# deltaNdeath : the observed number of deaths between each observed point pattern 
# dist1, dist2, dist_dk : the matrices of distances between all observed point patterns for
#   the cardinality distance, the smoothed cardinality distance, the d_k distance, respectively

ndisc=length(XXdisc)

# I.1 Estimation of the birth rate
# I.1.1 Estimation of the birth rate based on the number of points (without smoothing)
resalphadiscLbirth1=NULL
for(i in 1:ndisc){
  resalphadiscLbirth1[i]=betaestdisc(XXdisc[[i]],deltaNbirth,deltat=0.14,dist=dist1[i,],K=function(t,a,b) 1-t)
}
plot(resalphadiscLbirth1,type='l',col=2)

# I.1.2 Estimation of the birth rate based on the number of points, but through a smoothing kernel
h2=hoptimdisc(XXdisc,deltaNbirth,0.14,dist2)
resalphadiscLbirth2=NULL
for(i in 1:ndisc){
  resalphadiscLbirth2[i]=betaestdisc(XXdisc[[i]],deltaNbirth,deltat=0.14,h=h2,dist=dist2[i,])
}
lines(resalphadiscLbirth2,col=3)

# I.1.3 Estimation of the birth rate based on the d_k distance 
h_dk=hoptimdisc(XXdisc,deltaNbirth,0.14,dist_dk)
resalphadiscLbirth_dk=NULL
for(i in 1:ndisc){
  resalphadiscLbirth_dk[i]=betaestdisc(XXdisc[[i]],deltaNbirth,deltat=0.14,h=h_dk,dist=dist_dk[i,])
}
lines(resalphadiscLbirth_dk,col=5)


# I.2 Estimation of the death rate
# I.2.1 Estimation of the death rate based on the number of points (without smoothing)
resalphadiscLdeath1=NULL
for(i in 1:ndisc){
  resalphadiscLdeath1[i]=deltaestdisc(XXdisc[[i]],deltaNdeath,deltat=0.14,dist=dist1[i,],K=function(t,a,b) 1-t)
}
plot(resalphadiscLdeath1,type='l',col=2)

# I.2.2 Estimation of the death rate based on the number of points, but through a smoothing kernel
h2=hoptimdisc(XXdisc,deltaNdeath,0.14,dist2)
resalphadiscLdeath2=NULL
for(i in 1:ndisc){
  resalphadiscLdeath2[i]=deltaestdisc(XXdisc[[i]],deltaNdeath,deltat=0.14,h=h2,dist=dist2[i,])
}
lines(resalphadiscLdeath2,col=3)

# I.2.3 Estimation of the death rate based on the d_k distance 
h_dk=hoptimdisc(XXdisc,deltaNdeath,0.14,dist_dk)
resalphadiscLdeath_dk=NULL
for(i in 1:ndisc){
  resalphadiscLdeath_dk[i]=deltaestdisc(XXdisc[[i]],deltaNdeath,deltat=0.14,h=h_dk,dist=dist_dk[i,])
}
lines(resalphadiscLdeath_dk,col=5)

#Death rates of Langerin wrt to the number of points
nndisc=unlist(lapply(XXdisc,npoints))
plot(nndisc,resalphadiscLdeath1,col=2)
points(nndisc,resalphadiscLdeath2,col=3)
points(nndisc,resalphadiscLdeath_dk,col=4)


### II Estimation of intensities for Rab11 data
load(file="Utrack_Rab11.RData")
#This environment 
# XXdisc : the observed point patterns, in chronological order
# deltaNbirth : the observed number of births between each observed point pattern 
# deltaNdeath : the observed number of deaths between each observed point pattern 
# dist1, dist2, dist_dk : the matrices of distances between all observed point patterns for
#   the cardinality distance, the smoothed cardinality distance, the d_k distance, respectively

ndisc=length(XXdisc)

# II.1 Estimation of the birth rate
# II.1.1 Estimation of the birth rate based on the number of points (without smoothing)
resalphadiscRbirth1=NULL
for(i in 1:ndisc){
  resalphadiscRbirth1[i]=betaestdisc(XXdisc[[i]],deltaNbirth,deltat=0.14,dist=dist1[i,],K=function(t,a,b) 1-t)
}
plot(resalphadiscRbirth1,type='l',col=2)

# II.1.2 Estimation of the birth rate based on the number of points, but through a smoothing kernel
h2=hoptimdisc(XXdisc,deltaNbirth,0.14,dist2)
resalphadiscRbirth2=NULL
for(i in 1:ndisc){
  resalphadiscRbirth2[i]=betaestdisc(XXdisc[[i]],deltaNbirth,deltat=0.14,h=h2,dist=dist2[i,])
}
lines(resalphadiscRbirth2,col=3)

# II.1.3 Estimation of the birth rate based on the d_k distance 
h_dk=hoptimdisc(XXdisc,deltaNbirth,0.14,dist_dk)
resalphadiscRbirth_dk=NULL
for(i in 1:ndisc){
  resalphadiscRbirth_dk[i]=betaestdisc(XXdisc[[i]],deltaNbirth,deltat=0.14,h=h_dk,dist=dist_dk[i,])
}
lines(resalphadiscRbirth_dk,col=5)


# II.2 Estimation of the death rate
# II.2.1 Estimation of the death rate based on the number of points (without smoothing)
resalphadiscRdeath1=NULL
for(i in 1:ndisc){
  resalphadiscRdeath1[i]=deltaestdisc(XXdisc[[i]],deltaNdeath,deltat=0.14,dist=dist1[i,],K=function(t,a,b) 1-t)
}
plot(resalphadiscRdeath1,type='l',col=2)

# II.2.2 Estimation of the death rate based on the number of points, but through a smoothing kernel
h2=hoptimdisc(XXdisc,deltaNdeath,0.14,dist2)
resalphadiscRdeath2=NULL
for(i in 1:ndisc){
  resalphadiscRdeath2[i]=deltaestdisc(XXdisc[[i]],deltaNdeath,deltat=0.14,h=h2,dist=dist2[i,])
}
lines(resalphadiscRdeath2,col=3)

# II.2.3 Estimation of the death rate based on the d_k distance 
h_dk=hoptimdisc(XXdisc,deltaNdeath,0.14,dist_dk)
resalphadiscRdeath_dk=NULL
for(i in 1:ndisc){
  resalphadiscRdeath_dk[i]=deltaestdisc(XXdisc[[i]],deltaNdeath,deltat=0.14,h=h_dk,dist=dist_dk[i,])
}
lines(resalphadiscRdeath_dk,col=5)

#Death rates of Rab11 wrt to the number of points
nndisc=unlist(lapply(XXdisc,npoints))
plot(nndisc,resalphadiscRdeath1,col=2)
points(nndisc,resalphadiscRdeath2,col=3)
points(nndisc,resalphadiscRdeath_dk,col=4)


# III Comparison between the estimated death rates of Langerin and Rab11
plot(resalphadiscLdeath2,type='l',ylim=range(c(resalphadiscLdeath2,resalphadiscRdeath2)))
lines(resalphadiscRdeath2,col=2)
legend("top",c('Langerin','Rab11'),col=1:2,lty=1)
title('Estimated death rates')

plot(ccf(resalphadiscLdeath2,resalphadiscRdeath2,lag.max=20))
abline(v=0,col=gray(0.5,0.5),lty=1)

