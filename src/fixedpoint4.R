rm(list=ls())
my.seed=1
set.seed(my.seed)

pdf("model7_sim1.tapered.pdf")

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

nsims       =  1000                              # number of networks simulated
n           =  100                               # number of nodes
theta       =  c(-2,1,1,1) * c(2,2,1/n,1/n)      # true parameters for model 1
theta <- c(-1.842, 1.347, -0.154, 0.853) #perfect
theta <- c(-2.162, 1.457, -0.1314, 0.736)
theta <- c(-2.0967979, 1.4274966, -0.1387846, 0.7541693)
theta <- c(-2.1088249, 1.4353562, -0.1380418, 0.7504813)
theta <- c(-2.0766314, 1.4381472, -0.1398167, 0.7496554)
theta <- c(-1.9746393, 1.4805697, -0.1417278, 0.7249104)
theta <- c(-2.0028533, 1.5853063, -0.1655770, 0.7367452)
theta <- c(-1.8239798, 1.4320427, -0.1626898, 0.7308445)
theta <- c(-1.8793576, 1.4214192, -0.1594352, 0.7495866)
mv_1 <- c(394, 342, 3000, 180)
mv_1 <- c(358.16, 351.24, 1.0*2643.56, 1.0*123.69)
mv_1 <- c(400.31, 349.88, 3254.42, 1.3*126.28)
##################
#                #
#     Set-up     #
#                #
##################

sim <- initialize.network(theta, n, directed = FALSE)
x <- rbinom(n, 1, 0.5) # attributes
set.vertex.attribute(sim, # the name of the network object
                     "x", # the name we want to reference the variable by in that object
                     x # the value we are giving that variable
) 

load(file="sim2.RData")
formula <- sim ~ edges + nodematch("x") + kstar(2) + triangles
names(mv_1) <- names(summary(formula))
names(theta) <- names(mv_1)

if(F){
fit <- ergm.tapered(formula, eval.loglik=FALSE, target.stats=mv_1,
                    control=control.ergm.tapered(parallel=4,init=theta, MCMLE.MCMC.precision=0.001,MCMC.burnin=1000000, MCMC.interval=10000) )
sim <- fit$newnetwork
save(sim, file="sim2.RData")
}else{
load(file="sim2.RData")
#theta <- c(-4.000, 2.000, 0.010, 0.010)
sim <- simulate_ergm.tapered(sim ~ edges+nodematch("x")+ kstar(2) + triangles,
               tapering.centers=mv_1, tau=0.25/mv_1,
               control=control.simulate.formula(MCMC.burnin=1000000, MCMC.interval=10000),
               coef = theta)
}
cbind(theta,mv_1,summary(sim ~ edges+nodematch("x")+ kstar(2) + triangles))
pdf("sim.pdf")
plot(sim,col="x")
dev.off()

registerDoParallel(10)
fn <- function(theta,sim,mv_1,nsims){
  a = foreach(i = 1:10, .combine = rbind) %dorng% {
  simulate_ergm.tapered(sim ~ edges+nodematch("x")+ kstar(2) + triangles,
               tapering.centers=mv_1, tau=0.25/mv_1,
               nsim = nsims,
               control=control.simulate.formula(MCMC.burnin=1000000, MCMC.interval=10000),
               coef = theta,         
               output = "stats"
  )
  }
o <- colMeans(a)-mv_1
o2 <- c(3,3,0.5,1)*o*o
o2 <- o*o
o2 <- c(1,1,1,3)*o*o
message(sprintf("val = %f %f %f %f: %f", o[1],o[2],o[3],o[4], sqrt(sum(o2))))
sqrt(sum(o2))
}
theta
fn(theta,sim,mv_1,nsims)
fit <- optim(par=theta, fn=fn, sim=sim, mv_1=mv_1, nsims=nsims, control=list(maxit=50,abstol=2,trace=6))
fit
theta <- fit$par
names(theta) <- names(mv_1)
sim <- simulate_ergm.tapered(sim ~ edges+nodematch("x")+ kstar(2) + triangles,
               tapering.centers=mv_1, tau=0.25/mv_1,
               control=control.simulate.formula(MCMC.burnin=1000000, MCMC.interval=10000),
               coef = theta)
cbind(theta,mv_1,summary(sim ~ edges+nodematch("x")+ kstar(2) + triangles))
save(sim, theta, file="sim2.RData")
