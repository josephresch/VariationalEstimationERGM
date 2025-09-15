rm(list=ls())
my.seed=3141
set.seed(my.seed)

library(ergm.tapered)
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
mv_1 <- c(0.97*393.0512, 341.0188, 3092.0576, 1.5*117.4754)
theta <- c(-4.04102369, 1.92622055, -0.01309495, 0.40731208)
theta <- c(-4.04071794, 1.92361302, -0.01291308, 0.41322262)
##################
#                #
#     Set-up     #
#                #
##################

load("sim10.RData")
names(theta) <- names(summary(sim ~ edges + nodematch("x") + kstar(2) + triangles))
names(mv_1) <- names(summary(sim ~ edges + nodematch("x") + kstar(2) + triangles))
theta
set.seed(my.seed)
registerDoParallel(10)
g_sim_stats = foreach(i = 1:10, .combine = rbind) %dorng% {
 simulate_ergm.tapered(sim ~ edges + nodematch("x") + kstar(2) + triangles,
                  nsim = 100*nsims, tapering.centers=mv_1, tau=0.25/mv_1,
                  coef = theta,
                  output = "stats",
		  control=control.simulate.formula(MCMC.burnin=1000000, MCMC.interval=5000)
                              )
}
apply(g_sim_stats,2,mean)
cbind(mv_1, apply(g_sim_stats,2,mean))
save(g_sim_stats, file="g_sim_stats.RData")
