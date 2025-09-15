rm(list=ls())
my.seed=2000
set.seed(my.seed)


library(ergm.tapered)
#library(mfergm)
library(optimx)     # mfergm likelihood
library(R.utils)    # set time on ergm
library(doParallel) # parallel loops using 'foreach'
library(doRNG) # parallel loops using 'foreach'

##################
#                #
#     Set-up     #
#                #
##################

#if(T){
nsims       =  10                               # number of networks simulated
n           =  2000                             # number of nodes
theta       =  c(-2,1,1,1) * c(2,2,1/n,1/n)     # true parameters for model

if(F){
registerDoParallel(5)
# Extract summary stats from 'g_sim'
g_sim_stats = foreach(i = 1:nsims, .combine = rbind) %dorng% {
  g <- readRDS(paste0("~/Drive/Students/JosephResch/OAVE/simdata/n2000/g_n2000_",i,".rds"))
  summary_formula(g ~ edges + nodematch("x") + kstar(2) + triangles)
}
attr(g_sim_stats,"rng") <- NULL
g_sim_stats

##################
#                #
#     MFERGM     #
#                #
##################

# Use "theta" on mfergm instead
tobs <- data.frame(matrix(NA, ncol = 4, nrow = nsims))
names(tobs) <- c("edges", "x", "kstar2", "triangles")
for (i in 1:nsims) 
{
  tobs[i,] = g_sim_stats[i,] / (c((n^2)/2, (n^2)/2, (n^3)/2 , (n^3)/2 ))  
}

# Initialize the parameters at MPLE
mcmc_params = foreach(i = 1:nsims, .combine = rbind) %dorng% {
  g <- readRDS(paste0("~/Drive/Students/JosephResch/OAVE/simdata/n2000/g_n2000_",i,".rds"))
  ergm(g ~ edges + nodematch("x") + kstar(2) + triangles
       )$coefficients
}
attr(mcmc_params,"rng") <- NULL

save.image("model1_n2000_mcmc.RData")
stopImplicitCluster()
}else{
load("model1_n2000_MPLE.RData")
init_params <- init_params[1:10,]

load("model1_n2000_mcmc.good.RData")
#attr(init_params,"rng") <- NULL
#attr(mcmc_params,"rng") <- NULL
}

str(init_params)
cbind(theta,t(init_params)) / c(2,2,1/n,1/n)
cbind(theta,t(mcmc_params)) / c(2,2,1/n,1/n)

cbind(theta,apply(init_params,2,mean)) / c(2,2,1/n,1/n)
cbind(theta,apply(mcmc_params,2,mean)) / c(2,2,1/n,1/n)

# MPLE
for(i in 1:ncol(init_params)){
  print( sprintf("%12s p-value: %f", colnames(init_params)[i],t.test(init_params[,i],mu=theta[i])$p.value))
}

# MFERGM
for(i in 1:ncol(init_params)){
  print( sprintf("%12s p-value: %f", colnames(init_params)[i],t.test(mcmc_params[,i],mu=theta[i])$p.value))
}

cbind(
sqrt(apply(init_params,2,var)),
sqrt(apply(mcmc_params,2,var))
)
cor(init_params, mcmc_params)
#save.image("model1_n2000_mcmc.RData")
