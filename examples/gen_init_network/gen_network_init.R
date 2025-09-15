rm(list=ls())
my.seed=1
set.seed(my.seed)
setwd("~/ERGM")

library(ergm)
library(doRNG)
library(mfergm)
library(optimx)     # mfergm likelihood
library(R.utils)    # set time on ergm
library(doParallel) # parallel loops using 'foreach'

#####################################################################
#                                                                   #
#     Create High Transitivity ERGM Params (using ERGM fitting)     #
#                                                                   #
#####################################################################

niter       =  100
nsims       =  200                                        # number of networks simulated
n           =  200                                         # number of nodes
theta       =  c(-2,1,1,1) * c(2,2,1/n,1/n)               # true parameters for model 1
cores       =  10


##################
#                #
#     Set-up     #
#                #
##################

sim <- initialize.network(theta, n, directed = FALSE)
x <- rbinom(n, 1, 0.5) # attributes
set.vertex.attribute(sim, "x", x)
formula <- sim ~ edges + nodematch("x") + kstar(2) + triangles

summary_stats = summary_formula(sim ~ edges + nodematch("x") + kstar(2) + triangles)


registerDoParallel(cores)
system.time({
for (iter in 1:niter)
{
  # load("sim_model1_n50.RData")
  g_sim= foreach(i = 1:nsims) %dorng% {
    simulate(sim ~ edges + nodematch("x") + kstar(2) + triangles, 
             coef = theta
             ,control=control.simulate.formula(MCMC.burnin=1000000, MCMC.interval=1000000)
    )
  }
  sim = g_sim[[nsims]]
  sim_summary = summary_formula(sim ~ edges + nodematch("x") + kstar(2) + triangles)
  summary_stats = rbind(summary_stats, sim_summary)
}
})

save(sim, file="sim_model1_n100.RData")
save.image("gen_network_init_n100.RData")

set.seed(1)
g_sim_stats= foreach(i = 1:100000, .combine = rbind) %dorng% {
  simulate(sim ~ edges + nodematch("x") + kstar(2) + triangles, 
           coef = theta,
           output = "stats"
           ,control=control.simulate.formula(MCMC.burnin=1000000, MCMC.interval=1000000)
  )
}
colMeans(g_sim_stats)


