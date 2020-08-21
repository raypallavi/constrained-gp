### fitting nonparametric GP using rstan()
## code for longitudinal part only
source("supp_funs.R")

# data generating functions:
#function 1 ---> ITP
f1 <- function(x, xmax, b1){
  (1 - exp(-b1 * x))/(1 - exp(-b1 * xmax))
}
#function 2 ---> strictly monotone increasing
f2 <- function(x, xmax, b1=0.1){
  log(1+b1*x)/log(1+b1*xmax)
}


# generating the data:
set.seed(1234)
n <- 200
tmax <- 52 #maximum time T
tm <- sort(runif(n,0,tmax))
#beta1 <- 0.1
#y.tr <- f1(tm,xmax = tmax, b1 = beta1)
y.tr <- f2(tm,xmax = tmax)
#error standard deviation:
sig.tr <- 0.07
err <- rnorm(n,0,sig.tr)
y <- y.tr + err
plot(tm,y,type="p",pch=20,col="green3",ylim=range(c(y,y.tr)))
points(tm,y.tr,type="l",lwd=2)

# forming the data for rstan:
dat1 <- data_stan(y, tm=tm, tmax=tmax)

library(rstan)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = TRUE)

# defining model for rstan:
mod_code <- "
 data{
   int<lower=1> n;
   vector[n] y;
   int<lower=1,upper=n> N;
   matrix[n,N+1] X;
   matrix[N+1,N+1] K;
   vector[N+1] cj;
 }
 parameters{
   simplex[N+1] xitemp;
   real<lower=0> sig;
   real<lower=0> tau;
 }
 transformed parameters{
   vector[N+1] xi = xitemp ./ cj;
   vector[n] mu = X * xi;
 }
 model{
   y ~ multi_normal(mu, diag_matrix(rep_vector(square(sig),n)));
   xitemp ~ multi_normal(rep_vector(0,N+1), square(tau) * K);
   target += -2 * log(sig);
   target += -2 * log(tau);
 }
"

# drawing posterior samples using rstan:
system.time({
  fit1 <- stan(model_code = mod_code, data = dat1, chains = 3,
               iter = 3000, warmup = 1000, thin = 1, verbose = FALSE,
               seed = "12345", algorithm = "HMC")
})
#extracting posterior samples:
sam.out <- extract(fit1)
pm <- colMeans(sam.out$mu)
lcl <- apply(sam.out$mu,2,quantile,probs=0.025)
ucl <- apply(sam.out$mu,2,quantile,probs=0.975)

# plot for longitudinal profile:
plot(tm,y,type="p",pch=20,col="green3",ylim=range(c(y,y.tr,pm,lcl,ucl)))
points(tm,y.tr,type="l",lwd=2,col=2)
points(tm,pm,type="l",lwd=2,col="blue")
points(tm,lcl,type="l",lwd=2,lty=2,col="blue")
points(tm,ucl,type="l",lwd=2,lty=2,col="blue")

#error measures:
(rmse <- sqrt(mean((y.tr - pm)^2)))
(mae <- mean(abs(y.tr - pm)))

### diagnostic checks:
library(mcmcse); library(coda); library(mcmcplots)
efsam.xi = efsam.mu = integer(); for(i in 1:ncol(sam.out$xi)){
  efsam.xi[i] = mcmcse::ess(sam.out$xi[,i])
  efsam.mu[i] = mcmcse::ess(sam.out$mu[,i])
}
#checking effective sample size to MCMC sample size ratio:
min(efsam.mu/nrow(sam.out$mu))
min(efsam.xi/nrow(sam.out$xi))
mcmcse::ess(sam.out$sig)/length(sam.out$sig)
mcmcse::ess(sam.out$tau)/length(sam.out$tau)
#checking traceplots:
par(mfrow=c(3,3))
xi_mat=matrix(sam.out$xi[,1:9],ncol = 9,
              dimnames = list(NULL,c(paste("xi[",1:9,"]",sep = ""))))
traplot(as.mcmc(xi_mat))
mu_mat=matrix(sam.out$mu[,1:9],ncol = 9,
              dimnames = list(NULL,c(paste("mu[",1:9,"]",sep = ""))))
traplot(as.mcmc(mu_mat))
par(mfrow=c(1,1))
traplot(as.mcmc(as.vector(sam.out$sig)),main="Traceplot for sig")
traplot(as.mcmc(as.vector(sam.out$tau)),main="Traceplot for tau")
