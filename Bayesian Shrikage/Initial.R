#Inverse Gaussian
library(statmod)
#This is for inverse gamma
library(invgamma)
#this is for polya-gamma
library(pgdraw)
# Generage Omega
library(purrr)
#multivariate normal
library(MASS)
#make positive
library(corpcor)
#bernoulli 
library(statip)
#test positive definite
library(matrixcalc)
#For force symmetric
library(Matrix)
#for conditional normal
library(condMVNorm)

#setwd("C:/Users/think/Desktop/R/MCMC code")
#x=read.csv("x.csv")
#y=read.csv("y.csv")
#x=as.matrix(x[,-1])
#y=as.matrix(y[,-1])


#train = c(1,2,sample(5:nrow(x), (nrow(x)-5)*0.3))
#xtest = x[-train,]
#ytest = y[-train,]


### 2.1 Initial data
n=nrow(y)
p=ncol(x)
n_iter=1000
n_burnin=500
### 2.2 Initial data
beta=matrix(rep(10,p),nrow=p)
sigma=1   #variance
a_sigma=1
b_sigma=1
tau=rep(1,p)
lambda=rep(1,p)
alpha=log(lambda)
mu=0
v=1
N=10000
chi=rep(1,p)
phi=rep(1,p)
#positive-definite symmetric matrix
set.seed(10) 


#A_sim, D_sim, S_sim is like the paper to get A_sim to be omega
A_sim=diag(p)
D_sim=rep(0,p)
for(j in c(1:p)){
  for(k in c(1:p)){
    D_sim[j]=D_sim[j]+G0[j,k]
  }
}
S_sim=matrix(rep(0,p*p),nrow = p)
for(k in c(2:p)){
  for(j in c(1:(k-1))){
    if(j<=q & k<=q){S_sim[j,k]=1}
    else {S_sim[j,k]=rbern(1, 0.5)}
    if(G0[j,k]==1){
      A_sim[j,k]= -S_sim[j,k]/(max(D_sim[j],D_sim[k])*1.1+0.1)
      A_sim[k,j]=A_sim[j,k]
    }
  }
}

# SIGMAX_sim=solve(A_sim)
# SIGMAX_sim=SIGMAX_sim%*%(diag(1/diag(SIGMAX_sim)))
OMEGA=A_sim
omega=matrix(rep(0,p*(p-1)),nrow=p-1)

# omega = diag(p)
# for(j in c(1:p)){
#   for(i in c(1:p)){
#     omega[i,j]=omega[j,i]}}
is.positive.definite(OMEGA)
# eta=1
# epsilon=1
eta=10
epsilon=10
SIGMA=rep(1,p)

### 3. calculation before loop
one=matrix(rep(1,p),nrow=p)

###4. result
finalbeta=matrix(rep(1,n_iter*p),ncol=p)
finalchi=matrix(rep(1,n_iter*p),ncol=p)
finalSIGMA=matrix(rep(1,n_iter*p),ncol=p)
finaltau=matrix(rep(1,n_iter*p),ncol=p)
finalphi=matrix(rep(1,n_iter*p),ncol=p)
finalalpha=matrix(rep(1,n_iter*p),ncol=p)
finalOMEGA=list()
finalsigma_2=array()
finalA=list()