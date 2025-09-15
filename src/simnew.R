rm(list=ls())
my.seed=1
set.seed(my.seed)

library(ergm.tapered)
library(doRNG)
library(mfergm)
library(optimx)      # mfergm likelihood
library(R.utils)     # set time on ergm
library(doParallel)  # parallel loops using 'foreach'

load("sim_model1_n100_increased.RData")
ls()

mv_1        =  c(393.0512, 341.0188, 3092.0576, 1.5*117.4754)
theta = c(-3.0989464,   1.5654049,  -0.0625170,   0.5898891) 
theta = c(-3.01074016, 1.53137208, -0.06704853, 0.60209288)
theta = c(-3.022, 1.53, -0.0663, 0.60)

names(theta) <- names(summary(sim ~ edges+nodematch("x")+ kstar(2) + triangles))
names(mv_1) <- names(summary(sim ~ edges+nodematch("x")+ kstar(2) + triangles))

### test out sim and theta
test_sim <- simulate_ergm.tapered(sim ~ edges+nodematch("x")+ kstar(2) + triangles,
                             nsim = 100000,
                             tapering.centers=mv_1, tau=0.25/mv_1,
                             control=control.simulate.formula(parallel=10,MCMC.burnin=1000000, MCMC.interval=10000),
                             coef = theta,
                             output = "stats")

#                                         #1
#                theta           mv_1          
# edges       -3.0989464  379  393.0512  392.9802
# nodematch.x  1.5654049  331  341.0188  342.0167
# kstar2      -0.0625170 2837 3092.0576 3093.2289
# triangle     0.5898891  143  176.2131  176.9745

#                                     #2                             
# edges       -3.01074016  393.0512  393.6422
# nodematch.x  1.53137208  341.0188  341.1893
# kstar2      -0.06704853 3092.0576 3092.6556
# triangle     0.60209288  176.2131  176.5332

cbind(theta,summary(sim ~ edges+nodematch("x")+ kstar(2) + triangles),mv_1,colMeans(test_sim))
