source("Generate.R")
source("Initial.R")

n_iter=10
timestart<-Sys.time()
#G0_index = matrix(NA,nrow = (p-1), ncol = p)
#G0_index for finding where G0 is 0
# for(j in c(1:p-1)){
# G0_index[,j]=G0[-j,j]
# }

for(t in c(1:n_iter)){
  D_tau=diag(tau)
  B=t(x)%*% x+D_tau
  beta=mvrnorm(n=1,mu=solve(B)%*%t(x)%*%y,sigma^2*solve(B))
  # rgamma's parameter is n, shape, rate=1/scale)
  sigma_2=rinvgamma(n=1,shape=(n+p-1)/2+a_sigma,scale =1/(0.5*t(y-x%*%beta)%*%(y-x%*%beta)+0.5%*%t(beta)%*%solve(D_tau)%*%beta+b_sigma))
  for (j in c(1:p)){
    chi[j]=log(2*N/tau[j])}
  for (j in c(1:p)){
    # rinvgauss's parameter is mu and lambda 
    tau[j]=1/(rinvgauss(n=1,m=lambda[j]/abs(beta[j]),s=lambda[j]^2))
    phi[j]=pgdraw(N,abs(2*alpha[j]-chi[j]))
    }
  alpha=mvrnorm(n=1,mu=solve(4*diag(phi)+OMEGA/v)%*%((2-N+mu/v)*one+2*diag(phi)%*%chi),Sigma=solve(4*diag(phi)+OMEGA/v))
  A=epsilon*diag(rep(1,p))+eta*matrix(rep(1,p*p),nrow=p)+(alpha-mu)%*%t(alpha-mu)/v
  
  #set \omega_j from \Omega
  for (j in c(1:p-1)){
    omega[,j]=OMEGA[-j,j]
  }
  for(j in c(1:p)){
  #try without considering G_0
    if(j==1){
      posi=round(OMEGA[-j,-j]*(A[j,j]^(-1)),6)
      print(is.positive.definite(posi))
      #omega[,j]=mvrnorm(n=1,mu=t(A[-j,j])%*%solve(A[j,j]*solve(OMEGA[-j,-j])),Sigma=posi)
      omega[,j]=round(mvrnorm(n=1,mu=-t(posi)%*%A[-j,j],Sigma=posi),6)
    }
    
    if(j>1){
      posi=round(OMEGA[-j,-j]*A[j,j]^(-1),6)
      print(is.positive.definite(posi))
      given_num=c(1:(j-1))
      given_x=omega[j-1,c(1:(j-1))]
      # omega[,j]=rcmvnorm(n=1,mean=t(A[-j,j])%*%solve(A[j,j]*solve(OMEGA[-j,-j])),sigma=posi,
      #                    dep=c(1:p-1),given=given_num,X=given_x, method = "eigen")
      omega[,j]=round(rcmvnorm(n=1,mean=-t(posi)%*%A[-j,j],sigma=posi,
                         dep=c(1:p-1),given=given_num,X=given_x, method = "svd"),6)
    }
    #SIGMA[j]=round((rgamma(n=1,shape=(1-epsilon-eta)/2,rate=A[j,j]/2))^{-1},6)  #scale=1/(2/A[j,j])
    SIGMA[j]=round(rinvgamma(n=1, shape=(epsilon+eta+3)/2, rate=A[j,j]/2),6)
    OMEGA[j,j]=round(1/SIGMA[j]+t(omega[,j])%*%solve(OMEGA[-j,-j])%*%omega[,j],6)
    
    #OMEGA[j,j]=1/SIGMA[j]+t(OMEGA[-j,j])%*%solve(OMEGA[-j,-j])%*%OMEGA[-j,j]
   }
  for(j in c(1:p)){
  OMEGA[-j,j]=omega[,j]}
  
  OMEGA <- matrix(OMEGA,nrow = p, ncol = p)
  #Save other parameters
  #OMEGA <- matrix(forceSymmetric(OMEGA),nrow = p, ncol = p)
  #OMEGA <- make.positive.definite(OMEGA, tol=1e-3)
  #OMEGA1 <- make.positive.definite(OMEGA, tol=1)
  # finalbeta[t,]=as.numeric(beta)
  # finalchi[t,]=as.numeric(chi)
  # finalSIGMA[t,]=as.numeric(SIGMA)
  # finaltau[t,]=as.numeric(tau)
  # finalphi[t,]=as.numeric(phi)
  # finalalpha[t,]=as.numeric(alpha)
  # finalOMEGA[[t]]=OMEGA
  # finalsigma_2[t]=sigma_2
  # finalA[[t]]=A
}

timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)


# sy=matrix(,ncol = p,nrow = p)
# for(i in c(1:p)){
#   for (j in c(1:p)){
#     sy[i,j]=(OMEGA[i,j]==OMEGA[j,i])
#   }
#   
# }
# 
# compare=matrix(,nrow=100,ncol=2)
# for(i in c(1:p)){compare[i,1]=OMEGA[i,1]}
# for(i in c(1:p)){compare[i,2]=OMEGA[1,i]}


