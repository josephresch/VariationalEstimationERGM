rm(list=ls())
my.seed=3141
set.seed(my.seed)

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

nsims       =  1000                               # number of networks simulated
n           =  100                               # number of nodes
mv_1 <- c(400.31, 349.88, 3254.42,     126.28)   # mean-value parameters standard model
theta       =  c(-2,1,1,1) * c(2,2,1/n,1/n)      # true parameters for model 2
##################
#                #
#     Set-up     #
#                #
##################

load("sim10.RData")
theta
set.seed(my.seed)
registerDoParallel(10)
g_sim_stats = foreach(i = 1:10, .combine = rbind) %dorng% {
 simulate(sim ~ edges + nodematch("x") + kstar(2) + triangles, 
                  nsim = 100*nsims,
                  coef = theta,
                  output = "stats",
		  control=control.simulate.formula(MCMC.burnin=1000000, MCMC.interval=10000)
                              )
}
apply(g_sim_stats,2,mean)
cbind(mv_1, apply(g_sim_stats,2,mean))
