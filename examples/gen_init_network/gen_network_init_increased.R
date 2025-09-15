rm(list=ls())
my.seed=1
set.seed(my.seed)

library(ergm.tapered)
library(doRNG)
library(mfergm)
library(optimx)      # mfergm likelihood
library(R.utils)     # set time on ergm
library(doParallel)  # parallel loops using 'foreach'

#####################################################################
#                                                                   #
#     Create High Transitivity ERGM Params (using ERGM fitting)     #
#                                                                   #
#####################################################################

nsims       =  1000                              # number of networks simulated
n           =  50                               # number of nodes
# n           =  100                               # number of nodes
theta       =  c(-2,1,1,1) * c(2,2,1/n,1/n)      # starting theta value
mv_1        =  c(95.7479, 82.7075, 363.4233, 1.5*13.9101) # target mv (n=50)
# mv_1        =  c(393.0512, 341.0188, 3092.0576, 1.5*117.4754) # target mv (n=100)
niter       =  500

### optimized theta
# n = 50:   c(-4.008734051,  1.956079392, -0.008150082,  0.525968752)
# rerun: -3.73216764  1.79534689 -0.03263824  0.62993806 
# n = 100:  c(-3.0989464,   1.5654049,  -0.0625170,   0.5898891)
# n = 100 (from Mark): c(-3.022, 1.53, -0.0663, 0.60)



##################
#                #
#     Set-up     #
#                #
##################

# starting sim
sim <- initialize.network(theta, n, directed = FALSE)
x <- rbinom(n, 1, 0.5) # attributes
set.vertex.attribute(sim,"x",x) 

formula <- sim ~ edges + nodematch("x") + kstar(2) + triangles
names(mv_1) <- names(summary(formula))
names(theta) <- names(mv_1)

# function to optimize theta
optim_theta <- function(theta,sim,mv_1,nsims){
  a = foreach(i = 1:10, .combine = rbind) %dorng% {
    simulate_ergm.tapered(sim ~ edges+nodematch("x")+ kstar(2) + triangles,
                          tapering.centers=mv_1, tau=0.25/mv_1,
                          nsim = nsims,
                          # control=control.simulate.formula(MCMC.burnin=1000000, MCMC.interval=100000),
                          coef = theta,         
                          output = "stats"
    )
  }
  o <- colMeans(a)-mv_1
  # o2 <- c(3,3,0.5,1)*o*o
  # o2 <- c(1,1,1,3)*o*o
  o2 <- o*o
  message(sprintf("val = %f %f %f %f: %f", o[1],o[2],o[3],o[4], sqrt(sum(o2))))
  sqrt(sum(o2))
}


# iterations to optimize sim and theta
registerDoParallel(10)
for (iter in 1:niter)
{
  sim_fit_time = system.time({
    sim_fit = foreach(i = 1:20) %dorng% {
      skip_to_next <- FALSE
      ergm_sim = 
        tryCatch(withTimeout(
          ergm.tapered(formula, eval.loglik=FALSE, target.stats = mv_1,
                       control=control.ergm.tapered(
                         # parallel=4,
                         init=theta 
                         # ,MCMLE.MCMC.precision=0.001
                         # ,MCMC.burnin=1000000 
                         # ,MCMC.interval=10000
                       ) ),
          timeout = 60, onTimeout = "error"),
          error = function(e) {skip_to_next <<- NULL})
      ergm_sim
    }
  })
  print(sim_fit_time)
  
  degen = which(sapply(sim_fit, is.null))
  if (length(degen) > 0) {sim_fit = sim_fit[-degen]}
  
  sim <- sim_fit[[1]]$newnetwork
  
  theta_fit <- optim(par=theta, fn=optim_theta, sim=sim, mv_1=mv_1, nsims=nsims, control=list(maxit=50,abstol=2,trace=6))
  theta <- theta_fit$par
  
  names(theta) <- names(summary(formula))
  
  print(iter)
}

theta
summary(sim ~ edges + nodematch("x") + kstar(2) + triangles)
# save(sim, theta, file="sim_model1_n50_increased.RData")


dist = 5
i = 0
while (dist > 0.5) {
  i = i + 1
### test out sim and theta
set.seed(1)
test_sim <- simulate_ergm.tapered(sim ~ edges+nodematch("x")+ kstar(2) + triangles,
                             nsim = 1000,
                             tapering.centers=mv_1, tau=0.25/mv_1,
                             control=control.simulate.formula(MCMC.burnin=1000000, MCMC.interval=1000000),
                             coef = theta,
                             output = "stats")
                             
cbind(mv_1,colMeans(test_sim))
dist = sqrt(sum((mv_1 - colMeans(test_sim))^2))
print(i)
print(dist)
}
