#multivariate normal
library(MASS)
#library(corpcor)
#bernoulli 
library(statip)

#Step 1: x
#setwd("C:/Users/think/Desktop/R/MCMC code")
set.seed(563) 
#G0 is the true G
n=50
p=100
g=30 # g is the number of virtual pathways
q=5  #q=5 important variables
ngene=rnbinom(g, size=30,mu=30) # 
# g pathways
path=matrix(NA,nrow=g,ncol=max(ngene))
path[1,]=c(c(1:q),rep(0,max(ngene)-q))
for(j in c(2:g)){
  path[j,]=c(sample((q+1):p,size=ngene[j],replace=FALSE),rep(0,max(ngene)-ngene[j]))}

#Randomly choose two genes and insert an edge between the two genes (G0)
G0=matrix(rep(0,p*p),nrow = p,ncol=p)
for(i in c(1:g)){ 
  pathtry=path[i,][path[i,]>0]
  connected=sample(pathtry,2,replace=F)
  unconnected=pathtry[!pathtry %in% connected]
  G0[connected[1],connected[2]]=1
  G0[connected[2],connected[1]]=1
  #Choose a connected and an unconnected 
  for (j in c(1:length(pathtry))){
    rancon=sample(connected,size=1)
    ranuncon=sample(unconnected,size=1)
    G0[rancon,ranuncon]=1
    G0[ranuncon,rancon]=1
    connected=c(connected,ranuncon)
    unconnected=unconnected[unconnected!=ranuncon]
    if (length(unconnected)==1){
      rancon=sample(connected,size=1)
      ranuncon=unconnected
      G0[rancon,ranuncon]=1
      G0[ranuncon,rancon]=1
      connected=c(connected,ranuncon)
      unconnected=unconnected[unconnected!=ranuncon]
      break
    }
  }
}
#(d)
for(i in c(1:n)){
  for(j in c((i+1):n)){
    if(G0[i,j]==0){
      if(rbinom(1, 1, 0.05)==1){
        G0[i,j]==1
        G0[j,i]==1}
    }
  }
}
diag(G0)=1

#Get Sigma_x 
A_x=diag(p)
D=rep(0,p)
for(j in c(1:p)){
  for(k in c(1:p)){
    D[j]=D[j]+G0[j,k]
  }
}
S=matrix(rep(0,p*p),nrow = p)
q=5
for(k in c(2:p)){
  for(j in c(1:(k-1))){
    if(j<=q & k<=q){S[j,k]=1}
    else {S[j,k]=rbern(1, 0.5)}
    if(G0[j,k]==1){
      A_x[j,k]= -S[j,k]/(max(D[j],D[k])*1.1+0.1)
      A_x[k,j]=A_x[j,k]
    }
  }
}
SIGMAX=solve(A_x)
SIGMAX=SIGMAX%*%(diag(1/diag(SIGMAX)))
x=matrix(NA,nrow = n,ncol=p)
for(i in c(1:n)){
  x[i,]=mvrnorm(n=1,mu=rep(0,p),Sigma = SIGMAX)
}
#GT=matrix(as.numeric(rbernoulli(p*p,0.7)),nrow=p)
#for(j in c(1:p)){
#  for(i in c(1:p)){GT[i,j]=GT[j,i]}}
#for(j in c(1:p)){GT[j,j]=1}



#Step 2: sigma
#a_sigma=1
#b_sigma=1
#library(invgamma)
#sigma_2=rinvgamma(n=1,shape=a_sigma,scale =b_sigma)


#Step 3: beta
#p=ncol(x)
#lambda=c(1:p)
#beta=matrix(rep(1,p),nrow=p)
#library(rmutil)
#for(j in c(1:p)){
#beta[j]=rlaplace(1, m=0, s=sqrt(sigma_2)/lambda[j])}
beta=matrix(c(rep(1,5),rep(0,p-5)),nrow=p)

#Step 4 : y
sigma_epi=1
epi=rnorm(n,mean=0,sd=sqrt(sigma_epi))
y=x%*%beta+epi

#write.csv(x,"x.csv")
#write.csv(y,"y.csv")







