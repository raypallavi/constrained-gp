### fitting dose-response and longitudinal model combined together
### fitting nonparametric GP using rstan() on the longitudinal part:
source("supp_funs.R")
library(dplyr, ggplot2)

#strictly monotone increasing
g0 <- function(x, xmax, b1=0.1){
  log(1+b1*x)/log(1+b1*xmax)
}

# data:
n_cohorts = c(10, 10, 10, 10) # number of subjects in each cohort
dose = c(.25, .5, .75, 1.5) # dose administered to each cohort
# linear dose-response model:
b0.tr = 0.5 # intercept
b1.tr = 1 # slope
tm = c(4, 12, 24, 52)
tmax = 52

### Forming the noisy data with patient ID information:
sig.tr = 0.1 # standard deviation
set.seed(1234)
err = rnorm(sum(n_cohorts)*length(tm), mean = 0, sd = sig.tr)
dl.dat = data.frame("PatientID"=rep(1:40, each=length(tm)),
                    "dose"=rep(dose, length(tm)*n_cohorts),
                    "time"=rep(tm, sum(n_cohorts))) %>% mutate(
                      resp.tr = b0.tr + (b1.tr * dose * g0(time,tmax)),
                      resp = resp.tr + err
                    )

#Plot for smaple means and standard errors:
msr <- dl.dat %>% group_by(dose,time) %>% 
  summarise("mean"=mean(resp),"sd"=sd(resp),"N"=length(resp))
msr$se <- msr$sd / sqrt(msr$N)
ggplot(msr, aes(x=time, y=mean, color=as.factor(dose))) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) + geom_line() +
  labs(color="Dose",x="Time",y="Mean response") + geom_point()
#Plot for the maximum time point:
msr %>% dplyr::filter(time == 52) %>% 
  ggplot(data = ., aes(x=dose, y=mean)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) + 
  geom_point() + geom_line() + labs(y="Mean response",x="Dose",title="Time = 52") +
  scale_x_continuous(breaks = dose)

#forming the GP part:
# forming the data as list() for rstan:
dat1 <- data_stan(y=dl.dat$resp,d=dl.dat$dose,tm=dl.dat$time,tmax=tmax)

library(rstan)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = TRUE)

# defining model for rstan:
#b0 and b1 are dose-model parameters
#rest are parameters associated with longitudinal model

mod_code <- "
 data{
   int<lower=1> n;
   vector[n] y;
   int<lower=1,upper=n> N;
   matrix[n,N+1] X;
   matrix[N+1,N+1] K;
   vector[N+1] cj;
   vector[n] d;
 }
 parameters{
   simplex[N+1] xitemp;
   real b0;
   real b1;
   real<lower=0> sig;
   real<lower=0> tau;
 }
 transformed parameters{
   vector[N+1] xi = xitemp ./ cj;
   vector[n] mu1 = b1 * d;
   vector[n] mu2 = X * xi;
   vector[n] mu = b0 + (mu1 .* mu2);
 }
 model{
   y ~ multi_normal(mu, diag_matrix(rep_vector(square(sig),n)));
   xitemp ~ multi_normal(rep_vector(0,N+1), square(tau) * K);
   b0 ~ normal(0, 1);
   b1 ~ normal(0, 1);
   target += -2 * log(sig);
   target += -2 * log(tau);
 }
"

# drawing posterior samples using rstan:
system.time({
  fit1 <- stan(model_code = mod_code, data = dat1, chains = 3,
               iter = 3000, warmup = 1000, thin = 2, verbose = FALSE,
               seed = "12345", algorithm = "NUTS")
})
print(fit1,pars=c("b0","b1","sig","tau"))
sam.out <- rstan::extract(fit1)
pm <- colMeans(sam.out$mu)
lcl <- apply(sam.out$mu,2,quantile,probs=0.025)
ucl <- apply(sam.out$mu,2,quantile,probs=0.975)

#error measures:
(rmse <- sqrt(mean((dl.dat$resp - pm)^2)))
(mae <- mean(abs(dl.dat$resp - pm)))
sqrt(mean((dl.dat$resp[dl.dat$time==52] - pm[dl.dat$time==52])^2))
mean(abs(dl.dat$resp[dl.dat$time==52] - pm[dl.dat$time==52]))

# plots for estimation:
dl.dat <- dl.dat %>% mutate(
  "mean"=pm, "lb"=lcl, "ub"=ucl
)


#longitudinal plot over dose:
ggplot(dl.dat,aes(x=time,group=factor(dose),color=factor(dose))) + 
  geom_line(aes(y = mean)) + geom_line(aes(y = lb), lty = 2) +
  geom_line(aes(y = ub), lty = 2) + geom_errorbar(aes(ymin=lb, ymax=ub)) +
  geom_point(data = msr, aes(x=time, y=mean, color=as.factor(dose))) +
  labs(x="Time",y="Response",color="Dose",
       title="Posterior mean (solid) and 2.5%, 97.5% quantiles (dashed)")

#plot for doses (for the final time)
dl.dat %>% dplyr::filter(time == 52) %>% 
  ggplot(data = .,aes(x=dose)) + geom_line(aes(y = mean), color="cornflowerblue") +
  geom_line(aes(y = ub), lty = 2, color="cornflowerblue") + 
  geom_line(aes(y = lb), lty = 2, color="cornflowerblue") + 
  geom_errorbar(aes(ymin=lb, ymax=ub), width=0.1, size=0.25) +
  geom_point(data = msr %>% dplyr::filter(time == 52), aes(x=dose, y=mean)) +
  scale_x_continuous(breaks = dose) +
  labs(x="Dose",y="Response",
       title="Posterior mean (solid) and 2.5%, 97.5% quantiles (dashed), time = 52")


##### diagnostic checks #####
library(mcmcse); library(coda); library(mcmcplots)
efsam.xi = efsam.mu = efsam.mu1 = efsam.mu2 = integer()
for(i in 1:ncol(sam.out$xi)){
  efsam.xi[i] = mcmcse::ess(sam.out$xi[,i])
}
for(i in 1:ncol(sam.out$mu)){
  efsam.mu[i] = mcmcse::ess(sam.out$mu[,i])
  efsam.mu1[i] = mcmcse::ess(sam.out$mu1[,i])
  efsam.mu2[i] = mcmcse::ess(sam.out$mu2[,i])
}
min(efsam.mu/nrow(sam.out$mu))
min(efsam.mu1/nrow(sam.out$mu1))
min(efsam.mu2/nrow(sam.out$mu2))
min(efsam.xi/nrow(sam.out$xi))
mcmcse::ess(sam.out$b0)/length(sam.out$b0)
mcmcse::ess(sam.out$b1)/length(sam.out$b1)
mcmcse::ess(sam.out$sig)/length(sam.out$sig)
mcmcse::ess(sam.out$tau)/length(sam.out$tau)
#checking traceplots:
xi_mat=matrix(sam.out$xi[,1:9],ncol = 9,
              dimnames = list(NULL,c(paste("xi[",1:9,"]",sep = ""))))
traplot(as.mcmc(xi_mat))
mu_mat=matrix(sam.out$mu[,1:9],ncol = 9,
              dimnames = list(NULL,c(paste("mu[",1:9,"]",sep = ""))))
traplot(as.mcmc(mu_mat))
mu1_mat=matrix(sam.out$mu1[,1:9],ncol = 9,
               dimnames = list(NULL,c(paste("mu1[",1:9,"]",sep = ""))))
traplot(as.mcmc(mu1_mat))
mu2_mat=matrix(sam.out$mu2[,c(1:3,5:7,9:11)],ncol = 9,
               dimnames=list(NULL,c(paste("mu2[",c(1:3,5:7,9:11),"]",sep=""))))
traplot(as.mcmc(mu2_mat))
traplot(as.mcmc(as.vector(sam.out$b0)),main="Traceplot for beta0")
traplot(as.mcmc(as.vector(sam.out$b1)),main="Traceplot for beta1")
traplot(as.mcmc(as.vector(sam.out$sig)),main="Traceplot for sig")
traplot(as.mcmc(as.vector(sam.out$tau)),main="Traceplot for tau")
