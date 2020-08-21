#supporting functions for GP regression:
library(fields); library(FastGP); library(MASS)

# For continuous function estimation:
h_j=function(x,my_knot,delta_N){
  N=length(my_knot)
  h=(1 - abs((x - my_knot)/delta_N))*(abs((x - my_knot)/delta_N) <= 1)
  return(h)
}

# For monotone function estimation:
psi_j=function(x,my_knot,delta_N)
{
  N=length(my_knot)
  k=rep(0,N)
  i=max(which(my_knot<=x))
  if(i==1)
  {
    k[1]=x-0.5*(x^2)/delta_N
    k[2]=x-my_knot[2]*x/delta_N+0.5*x^2/delta_N
  }
  if(i==2)
  {
    k[1]=delta_N/2
    k[2]=delta_N/2+(x-my_knot[2])*(1+my_knot[2]/delta_N)-0.5*(x^2-my_knot[2]^2)/delta_N
    k[3]=(x-my_knot[2])*(1-my_knot[3]/delta_N)+0.5*(x^2-my_knot[2]^2)/delta_N
  }
  if(i==N)
  {
    k[1]=delta_N/2
    k[2:(N-1)]=delta_N
    k[N]=delta_N/2
  }
  if(i!=1 && i!=2 && i!=N)
  {
    k[1]=delta_N/2
    k[2:(i-1)]=delta_N
    k[i]=delta_N/2+(x-my_knot[i])*(1+my_knot[i]/delta_N)-0.5*(x^2-my_knot[i]^2)/delta_N
    k[i+1]=(x-my_knot[i])*(1-my_knot[i+1]/delta_N)+0.5*(x^2-my_knot[i]^2)/delta_N
  }
  return(k)
}

# For convex function estimation:
phi_j=function(x,my_knot,delta_N)
{
  N=length(my_knot)
  k=rep(0,N)
  if(x>=my_knot[1] && x<my_knot[2])
  {
    k[1]=(x^2/2)-0.5*(x^3/3)/delta_N
    k[2]=0.5*(x^3/3)/delta_N
  }
  if(x>=my_knot[2] && x<my_knot[3])
  {
    k[1]=(my_knot[2]^2/2)-0.5*(my_knot[2]^3/3)/delta_N+0.5*delta_N*(x-my_knot[2])
    k[2]=my_knot[2]^3/(3*delta_N)+0.5*delta_N*(x-my_knot[2])+(x^2-my_knot[2]^2)-1.5*my_knot[2]*(x-my_knot[2])-0.5*(x^3/3)/delta_N
    k[3]=(1-my_knot[3]/delta_N)*(x^2/2-x*my_knot[2]+my_knot[2]^2/2)+0.5*(x^3/3-x*my_knot[2]^2+2*my_knot[2]^3/3)/delta_N
  }
  if(x>=my_knot[3] && x<my_knot[N])
  {
    k[1]=(my_knot[2]^2/2)-0.5*(my_knot[2]^3/3)/delta_N+0.5*delta_N*(x-my_knot[2])
    k[2]=my_knot[2]^3/(3*delta_N)+0.5*delta_N^2+my_knot[3]^2-2.5*my_knot[2]^2-0.5*(my_knot[3]^3/3)/delta_N+delta_N*(x-my_knot[3])
    k[3]=(1-my_knot[3]/delta_N)*(my_knot[3]^2/2-my_knot[3]*my_knot[2]+my_knot[2]^2/2)+0.5*(my_knot[3]^3/3-my_knot[3]*my_knot[2]^2+2*my_knot[2]^3/3)/delta_N+delta_N*(x-my_knot[3])
    
    if(N>4){
      for(j in 4:(N-1)){
        if(x<my_knot[j-1])
          k[j]=0
        if(x>=my_knot[j-1] && x<my_knot[j])
          k[j]=(1-my_knot[j]/delta_N)*(0.5*(x^2-my_knot[j-1]^2)-my_knot[j-1]*(x-my_knot[j-1]))+0.5*((x^3/3-my_knot[j-1]^3/3)-my_knot[j-1]^2*(x-my_knot[j-1]))/delta_N
        if(x>=my_knot[j] && x<my_knot[j+1])
          k[j]=(1-my_knot[j]/delta_N)*(0.5*(my_knot[j]^2-my_knot[j-1]^2)-my_knot[j-1]*(my_knot[j]-my_knot[j-1]))+0.5*((my_knot[j]^3/3-my_knot[j-1]^3/3)-my_knot[j-1]^2*(my_knot[j]-my_knot[j-1]))/delta_N+0.5*delta_N*(x-my_knot[j])+
            (1+my_knot[j]/delta_N)*(0.5*(x^2-my_knot[j]^2)-my_knot[j]*(x-my_knot[j]))-0.5*((x^3/3-my_knot[j]^3/3)-my_knot[j]^2*(x-my_knot[j]))/delta_N
        if(x>=my_knot[j+1])
          k[j]=(1-my_knot[j]/delta_N)*(0.5*(my_knot[j]^2-my_knot[j-1]^2)-my_knot[j-1]*(my_knot[j]-my_knot[j-1]))+0.5*((my_knot[j]^3/3-my_knot[j-1]^3/3)-my_knot[j-1]^2*(my_knot[j]-my_knot[j-1]))/delta_N+0.5*delta_N*(my_knot[j+1]-my_knot[j])+
            (1+my_knot[j]/delta_N)*(0.5*(my_knot[j+1]^2-my_knot[j]^2)-my_knot[j]*(my_knot[j+1]-my_knot[j]))-0.5*((my_knot[j+1]^3/3-my_knot[j]^3/3)-my_knot[j]^2*(my_knot[j+1]-my_knot[j]))/delta_N+delta_N*(x-my_knot[j+1])
      }
    }
    
    if(x>=my_knot[N-1] && x<my_knot[N])
      k[N]=(1-my_knot[N]/delta_N)*(0.5*(x^2-my_knot[N-1]^2)-my_knot[N-1]*(x-my_knot[N-1]))+0.5*((x^3/3-my_knot[N-1]^3/3)-my_knot[N-1]^2*(x-my_knot[N-1]))/delta_N
  }
  if(x>=my_knot[N])
  {
    k[1]=(my_knot[2]^2/2)-0.5*(my_knot[2]^3/3)/delta_N+0.5*delta_N*(x-my_knot[2])
    k[2]=my_knot[2]^3/(3*delta_N)+0.5*delta_N^2+my_knot[3]^2-2.5*my_knot[2]^2-0.5*(my_knot[3]^3/3)/delta_N+delta_N*(x-my_knot[3])
    k[3]=(1-my_knot[3]/delta_N)*(my_knot[3]^2/2-my_knot[3]*my_knot[2]+my_knot[2]^2/2)+0.5*(my_knot[3]^3/3-my_knot[3]*my_knot[2]^2+2*my_knot[2]^3/3)/delta_N+delta_N*(x-my_knot[3])
    for(j in 4:(N-1)){
      k[j]=(1-my_knot[j]/delta_N)*(0.5*(my_knot[j]^2-my_knot[j-1]^2)-my_knot[j-1]*(my_knot[j]-my_knot[j-1]))+0.5*((my_knot[j]^3/3-my_knot[j-1]^3/3)-my_knot[j-1]^2*(my_knot[j]-my_knot[j-1]))/delta_N+0.5*delta_N*(my_knot[j+1]-my_knot[j])+
        (1+my_knot[j]/delta_N)*(0.5*(my_knot[j+1]^2-my_knot[j]^2)-my_knot[j]*(my_knot[j+1]-my_knot[j]))-0.5*((my_knot[j+1]^3/3-my_knot[j]^3/3)-my_knot[j]^2*(my_knot[j+1]-my_knot[j]))/delta_N+delta_N*(x-my_knot[j+1])
    }
    k[N]=(1-my_knot[N]/delta_N)*(0.5*(my_knot[N]^2-my_knot[N-1]^2)-my_knot[N-1]*(my_knot[N]-my_knot[N-1]))+0.5*((my_knot[N]^3/3-my_knot[N-1]^3/3)-my_knot[N-1]^2*(my_knot[N]-my_knot[N-1]))/delta_N+0.5*delta_N*(x-my_knot[N])
  }
  return(k)
}

### Function to form design matrix:
des.mat1=function(x,my_knot,delta_N){
  # Function to form basis matrix for monotone constraint
  n=length(x)
  N=length(my_knot)-1
  # design matrix \Psi(n X N+1)
  X=matrix(0,n,N+1)
  for(l in 1:n){
    X[l,1:(N+1)]=psi_j(x[l],my_knot,delta_N)
  }
  return(X)
}

des.mat2=function(x,my_knot,delta_N){
  # Function to form basis matrix for convex constraint
  n=length(x)
  N=length(my_knot)-1
  # design matrix \Phi(n X N+1)
  X=matrix(0,n,N+1)
  for(l in 1:n){
    X[l,1:(N+1)]=phi_j(x[l],my_knot,delta_N)
  }
  return(X)
}


#Given a \nu (smoothness parameter of matern kernel) finding a value of 
# l (length-scale parameter) such that the correlation between the 
# maximum seperation is some small value, say 0.05

#Matern kernel with smoothness nu and length-scale l:
MK = function(x, y ,l, nu){
  ifelse(abs(x-y)>0, (sqrt(2*nu)*abs(x-y)/l)^nu/(2^(nu-1)*gamma(nu))*besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu), 1.0)
}

# function for uniroot:
fl=function(l,para){ 
  #para[1]=x, para[2]=y and para[3]=nu of MK : Matern kernel function;
  #para[4]=pre-specified value of the correlation
  a=MK(para[1],para[2],l,para[3])
  return(a-para[4])
}

# function for estimating l:
l_est=function(nu,range,val){
  # nu : smoothness; range : c(min, max) of the range of variable
  # val : pre-specified value of the correlation between the maximum seperation
  para=c(range[1],range[2],nu,val)
  rl=uniroot(f=fl,interval = c(0.000001,100000),para)
  return(rl$root)
}

# Covariance matrix
covmat=function(knot,nu,l){
  return(MK(rdist(knot),0,l,nu))
}

# Function that creates data to use in rstan()
data_stan <- function(y,d,tm,tmax,N,nu,l,val=0.05){
  # y : response
  # d : dose
  # tm : time/weeks/visits
  # tmax : maximum time value; extracted from 'tm', if missing
  # N : number of basis/knots
  # nu : smoothness parameter of Matern kernel
  # l : length-scale parameter of Matern kernel
  
  if(length(y) != length(tm)){
    stop("Weeks and response has different lengths")
  }
  n = length(y)
  if(missing(tmax)){
    tmax = max(tm)
  }
  x = tm/tmax
  if(missing(N)){
    N = ceiling(n/2) - 1 # can use (ceiling(n/4) - 1) or (n-1) instead.
  }
  #forming the basis matrix:
  int_length=1/N
  my_knots=seq(0,1,by=int_length)
  X=des.mat1(x,my_knots,int_length)
  cj=psi_j(1,my_knots,int_length)
  
  #forming prior matrix
  if(missing(nu)){
    nu = 0.75 # can use nu={0.5,1,1.5}
  } 
  if(missing(l)){
    l=l_est(nu,c(my_knots[1],my_knots[length(my_knots)]),val)
    #can use val={0.01,0.1}
    #or can use l={0.3,0.5,1}
  }
  # prior covariance K:
  K = covmat(my_knots,nu,l)
  
  if(!missing(d)){
    return(list("n"=n,"y"=y,"N"=N,"X"=X,"K"=K,"cj"=cj,"d"=d))
  }else{
    return(list("n"=n,"y"=y,"N"=N,"X"=X,"K"=K,"cj"=cj))
  }
}



