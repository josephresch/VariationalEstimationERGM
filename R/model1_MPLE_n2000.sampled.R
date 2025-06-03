rm(list=ls())
my.seed=2000
set.seed(my.seed)


library(ergm)
#library(ergm.tapered)
#library(mfergm)
library(optimx)     # mfergm likelihood
library(R.utils)    # set time on ergm
library(doParallel) # parallel loops using 'foreach'
library(doRNG) # parallel loops using 'foreach'

#####################################################################
#                                                                   #
#     Create High Transitivity ERGM Params (using ERGM fitting)     #
#                                                                   #
#####################################################################

nsims       =  200                              # number of networks simulated
n           =  2000                             # number of nodes
theta       =  c(-2,1,1,1) * c(2,2,1/n,1/n)     # true parameters for model


##################
#                #
#     Set-up     #
#                #
##################

if(F){
registerDoParallel(10)
## Initialize the parameters at MPLE
#samp_params_tapered = foreach(i = 1:nsims, .combine = rbind) %dorng% {
#  g <- readRDS(paste0("~/Downloads/n2000/g_n2000_",i,".rds"))
#  ergm.tapered(g ~ edges + nodematch("x") + kstar(2) + triangles, 
#               estimate = "MPLE", control=control.ergm.tapered(estimate.tapered.bias=FALSE))$coefficients
#}
#samp_params_tapered

## Extract summary stats from 'g_sim'
#g_sim_stats = foreach(i = 1:nsims, .combine = rbind) %dorng% {
#  g <- readRDS(paste0("~/Downloads/n2000/g_n2000_",i,".rds"))
#  summary_formula(g ~ edges + nodematch("x") + kstar(2) + triangles)
#}
#g_sim_stats

# Initialize the parameters at MPLE
samp_params = foreach(i = 1:nsims, .combine = rbind) %dorng% {
  g <- readRDS(paste0("~/Drive/Students/JosephResch/OAVE/simdata/n2000/g_n2000_",i,".rds"))
  ergm(g ~ edges + nodematch("x") + kstar(2) + triangles, 
       estimate = "MPLE", control=control.ergm(MPLE.samplesize=500000))$coefficients
}
save.image("model1_n2000_MPLE_sampled.RData")
stopImplicitCluster()
}else{
load("model1_n2000_MPLE_sampled.RData")
}

cbind(theta,t(samp_params)) / c(2,2,1/n,1/n)
#cbind(theta,t(samp_params_tapered)) / c(2,2,1/n,1/n)

cbind(theta,apply(samp_params,2,mean)) / c(2,2,1/n,1/n)
#cbind(theta,apply(samp_params_tapered,2,mean)) / c(2,2,1/n,1/n)

#MPLE
for(i in 1:ncol(samp_params)){
  print( sprintf("%12s p-value: %f", colnames(samp_params)[i],t.test(samp_params[,i],mu=theta[i])$p.value))
}

##tapered MPLE
#for(i in 1:ncol(samp_params)){
#  print( sprintf("%12s p-value: %f", colnames(samp_params)[i],t.test(samp_params_tapered[,i],mu=theta[i])$p.value))
#}

#cbind(
sqrt(apply(samp_params,2,var))
#sqrt(apply(samp_params_tapered,2,var))
#)

#save.image("model1_n2000_MPLE.RData")
