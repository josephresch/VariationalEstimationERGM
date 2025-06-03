rm(list=ls())
my.seed=19
set.seed(my.seed)
setwd("~/ERGM")

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
n           =  50                                 # number of nodes
mv_1        =  c(95.7479, 82.7075, 363.4233, 1.5*13.9101) # mean-value parameters for increased transitivity and tapered model
theta       =  c(-3.7322, 1.7953, -0.0326, 0.6299)  # increased trans theta
term_names  =  c("edges", "x", "kstar2", "triangles")
cores       =  10

##################
#                #
#     Set-up     #
#                #
##################

g <- initialize.network(theta, n, directed = FALSE)
x <- rbinom(n, 1, 0.5) # attributes
set.vertex.attribute(g, "x", x) 

formula <- g ~ edges + nodematch("x") + kstar(2) + triangles
names(mv_1) <- names(summary(formula))

load("sim_model1_n50_increased.RData")
set.seed(my.seed)
registerDoParallel(cores)


g_sim = foreach(i = 1:nsims) %dorng% {
  simulate_ergm.tapered(sim ~ edges + nodematch("x") + kstar(2) + triangles, 
                        tapering.centers=mv_1, tau=0.25/mv_1,
                        coef = theta
                        ,control=control.simulate.formula(MCMC.burnin=1000000, MCMC.interval=1000000)
  )
}

# Extract summary stats from 'g_sim'
g_sim_stats = foreach(i = 1:(nsims), .combine = rbind) %dorng% {
  as.vector(summary_formula(g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles))
}
cbind(mv_1, colMeans(g_sim_stats))

# pairs(g_sim_stats)

mv_s <- apply(g_sim_stats,2,mean)
mv_s
names(mv_s) <- c("edges", "x", "kstar2", "triangles")
cbind(mv_1, mv_s)
# mv_s <- mv_1
t.test(g_sim_stats[,1], mu=mv_1[1])
t.test(g_sim_stats[,2], mu=mv_1[2])
t.test(g_sim_stats[,3], mu=mv_1[3])
t.test(g_sim_stats[,4], mu=mv_1[4])

# Initialize the parameters at MPLE
init_params = foreach(i = 1:nsims, .combine = rbind) %dorng% {
  ergm.tapered(g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles, 
               estimate = "MPLE")$coefficients
}

##########################
#                        #
#     ERGM (Default)     #
#                        #
##########################

set.seed(my.seed)
### Compute ERGM Results
ergm_sim_list_tapered = foreach(i = 1:nsims) %dorng% {
  skip_to_next <- FALSE
  formula <- g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles
  ergm_sim = 
    tryCatch(withTimeout(
      ergm.tapered(formula, eval.loglik=FALSE,
                   control=control.ergm.tapered(init = init_params[i,],
                                                MCMC.burnin=100000,
                                                MCMC.interval=1000
                                                )),
      timeout = 5*60, onTimeout = "error"),
      error = function(e) {skip_to_next <<- NULL})
  ergm_sim
}

ergm_sim_estim_tapered = matrix(0, nrow = nsims, ncol = 4)
for(i in 1:nsims)
{
  if(!is.null(ergm_sim_list_tapered[[i]]))
  {
    est.params <- ergm_sim_list_tapered[[i]]$coefficients
    ergm_sim_estim_tapered[i,] <- est.params
  }
  else {ergm_sim_estim_tapered[i,] = c(NA, NA, NA, NA)}
}

ergm_result_tapered = apply(ergm_sim_estim_tapered, 2, mean, na.rm = T)
ergm_result_tapered

#### DEGENERACIES
degen_tp_ergm = which(!complete.cases(ergm_sim_estim_tapered))
degen_tp_ergm

ergm_sim_estim_tapered_filled = ergm_sim_estim_tapered

if(length(degen_tp_ergm) > 0 ){
  #degen_list_tapered = ergm_sim_list_tapered[[]]
  #ergm_degen_list_tapered = foreach(i = 1:length(degen_tp_ergm)) %dorng% {
  #  skip_to_next <- FALSE
  #  formula <- g_sim[[degen_tp_ergm[i]]] ~ edges + nodematch("x") + kstar(2) + triangles
  #  ergm_sim_tapered = tryCatch(withTimeout(ergm(formula, eval.loglik=FALSE,
  #                                        control=control.ergm(init=theta, MCMLE.confidence=0.95
  #                                       #   # MCMC.burnin=100000,
  #                                       #   # MCMC.interval=1000,
  #                                       #   # MCMC.samplesize = 5000,
  #                                          )
  #  ), timeout = 5*60, onTimeout = "error"),
  #  error = function(e) {skip_to_next <<- NULL})
  #  ergm_sim_tapered
  #}
  #ergm_degen_list_tapered
  
  # ergm_degen_estim_tapered = matrix(0, nrow = length(degen_tp_ergm), ncol = 4)
  for(i in 1:length(degen_tp_ergm))
  {
    #if(!is.null(ergm_degen_list_tapered[[i]]))
    #{
    # est.params <- ergm_degen_list_tapered[[i]]$coefficients
    # ergm_degen_estim_tapered[i,] <- est.params
    # ergm_sim_estim_tapered[degen_tp_ergm[i],] <- est.params
    ergm_sim_estim_tapered_filled[degen_tp_ergm[i],] <- init_params[degen_tp_ergm[i],]
    #}
    #else {ergm_degen_estim_tapered[i,] = c(NA, NA, NA, NA)}
  }
}



# compare mean-values
Sys.sleep(5)
ergm_mv_est_tapered = foreach(i = 1:nsims, .combine = rbind) %dorng% {
  simulate_ergm.tapered(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles,
                        nsim = 10, tapering.centers=mv_1, tau=0.25/mv_1,
                        control=control.simulate.formula(MCMC.burnin=100000, MCMC.interval=100000),
                        coef = ergm_sim_estim_tapered_filled[i,],         
                        output = "stats"
  )
}
colMeans(g_sim_stats)
ergm_mv_tapered = colMeans(ergm_mv_est_tapered)
ergm_mv_tapered


##################
#                #
#     MFERGM     #
#                #
##################

set.seed(my.seed)
# Use "theta" on mfergm instead
tobs <- data.frame(matrix(NA, ncol = 4, nrow = nsims))
names(tobs) <- c("edges", "x", "kstar2", "triangles")
for (i in 1:nsims) 
{
  tobs[i,] = g_sim_stats[i,] / (c((n^2)/2, (n^2)/2, (n^3)/2 , (n^3)/2 ))  
}

### Compute MFERGM Results
mfergm_estim <- NULL
# Double loop to reduce memory use
Sys.sleep(1)
for(j in 1:(nsims/100)){
  mfergm_estim_j = foreach(i = (100*(j-1)+(1:100)), .combine = rbind) %dorng% {
    pars <- init_params[i,] * c(.5,.5,n,n)
    addpars <- list(n=n, tobs = tobs[i,], x=x, ninit=1)
    cd.est <- tryCatch(optimx(pars, fn = loglikmf.model5, 
                              method = "BFGS", 
                              control = list(fnscale = -1), addpars = addpars), 
                       error = function(e) {return(c(NA, NA, NA, NA))})
    as.numeric(cd.est[1:4])
  }
  mfergm_estim <- rbind(mfergm_estim, mfergm_estim_j)
}
mfergm_estim = matrix(c(mfergm_estim[,1]*2, mfergm_estim[,2]*2, mfergm_estim[,3]/n, mfergm_estim[,4]/n),
                      nrow = nsims)
mfergm_result_med = apply(mfergm_estim, 2, median, na.rm = T) 
mfergm_result = apply(mfergm_estim, 2, mean, na.rm = T) 
cbind(theta,mfergm_result_med,mfergm_result)

#### DEGENERACIES
degen_tp_mfergm = which(!complete.cases(mfergm_estim))
degen_tp_mfergm

mfergm_estim_filled = mfergm_estim

if(length(degen_tp_mfergm) > 0){
  # Replace with MPLE
  for(i in 1:length(degen_tp_mfergm)){
    mfergm_estim_filled[degen_tp_mfergm[i],] <- init_params[degen_tp_mfergm[i],]
  }
}

# compare mean values
Sys.sleep(5)
mf_mv_est = foreach(i = 1:nsims, .combine = rbind) %dorng% {
  simulate_ergm.tapered(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles,
                        nsim = 10, tapering.centers=mv_1, tau=0.25/mv_1,
                        control=control.simulate.formula(MCMC.burnin=100000, MCMC.interval=100000),
                        coef =  mfergm_estim_filled[i,],    
                        output = "stats"
  )
}
mfergm_mv = colMeans(mf_mv_est)
mfergm_mv

# compare mean values
Sys.sleep(5)
MPLE_mv_est = foreach(i = 1:nsims, .combine = rbind) %dorng% {
  simulate_ergm.tapered(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles,
                        nsim = 10, tapering.centers=mv_1, tau=0.25/mv_1,
                        coef = init_params[i,],    
                        control=control.simulate.formula(MCMC.burnin=100000, MCMC.interval=100000),
                        output = "stats"
  )
}
MPLE_mv = colMeans(MPLE_mv_est)
MPLE_mv
MPLE_result = apply(init_params, 2, mean, na.rm = T) 

stopImplicitCluster()


mean_est <- cbind(theta, ergm_result_tapered, MPLE_result, mfergm_result)
rownames(mean_est) <- names(ergm_mv_tapered)
colnames(mean_est) <- c("true", "MCMC-MLE","MPLE","MFVLE")
mean_est

mean_mv_est <- cbind(round(colMeans(g_sim_stats), 2), round(ergm_mv_tapered,2), round(MPLE_mv,2), round(mfergm_mv, 2))
rownames(mean_mv_est) <- names(ergm_mv_tapered)
colnames(mean_mv_est) <- c("true", "MCMC-MLE","MPLE", "MFVLE")
mean_mv_est

# rm(ergm_sim_list_tapered)
# save.image("model1_n50_increased.RData")


outliers <- function(x, field = 1.5, na.rm = TRUE) 
{
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
  H <- field * IQR(x, na.rm = na.rm)
  which(x < (qnt[1] - H) | x > (qnt[2] + H))
}



# First natural params
complete_mf_results = mfergm_estim #[complete.cases(mfergm_estim),]
complete_tapered_results = ergm_sim_estim_tapered #[complete.cases(ergm_sim_estim_tapered),]
complete_MPLE_results = init_params #[complete.cases(init_params),]

bdrmse <- function(comp, theta, degen){
  outliers1 = union(outliers((comp[,1] - theta[1])^2), degen)
  outliers2 = union(outliers((comp[,2] - theta[2])^2), degen)
  outliers3 = union(outliers((comp[,3] - theta[3])^2), degen)
  outliers4 = union(outliers((comp[,4] - theta[4])^2), degen)
  rmse1 <- ((comp[,1] - theta[1])^2)
  rmse1[outliers1] <- max(rmse1[-outliers1])
  rmse2 <- ((comp[,2] - theta[2])^2)
  rmse2[outliers2] <- max(rmse2[-outliers2])
  rmse3 <- ((comp[,3] - theta[3])^2)
  rmse3[outliers3] <- max(rmse3[-outliers3])
  rmse4 <- ((comp[,4] - theta[4])^2)
  rmse4[outliers4] <- max(rmse4[-outliers4])
  c(sqrt(mean(rmse1)), sqrt(mean(rmse2)), sqrt(mean(rmse3)), sqrt(mean(rmse4)) ) 
}
bdmad <- function(comp, theta, degen){
  outliers1 = union(outliers(abs(comp[,1] - theta[1])),degen)
  outliers2 = union(outliers(abs(comp[,2] - theta[2])),degen)
  outliers3 = union(outliers(abs(comp[,3] - theta[3])),degen)
  outliers4 = union(outliers(abs(comp[,4] - theta[4])),degen)
  mad1 <- (abs(comp[,1] - theta[1]))
  mad1[outliers1] <- max(mad1[-outliers1])
  mad2 <- (abs(comp[,2] - theta[2]))
  mad2[outliers2] <- max(mad2[-outliers2])
  mad3 <- (abs(comp[,3] - theta[3]))
  mad3[outliers3] <- max(mad3[-outliers3])
  mad4 <- (abs(comp[,4] - theta[4]))
  mad4[outliers4] <- max(mad4[-outliers4])
  c((mean(mad1)), (mean(mad2)), (mean(mad3)), (mean(mad4)) ) 
}

# Degeneracy
length(intersect(which(is.na(mfergm_estim[,1])), which(is.na(ergm_sim_estim_tapered[,1]))))

# RMSE

MPLE_rmse <- bdrmse(complete_MPLE_results, theta,degen=NULL)
mf_rmse <- bdrmse(complete_mf_results, theta, degen_tp_mfergm)
tapered_rmse <- bdrmse(complete_tapered_results, theta, degen_tp_ergm)

# MAD
MPLE_mad <- bdmad(complete_MPLE_results, theta,degen=NULL)
mf_mad <- bdmad(complete_mf_results, theta, degen_tp_mfergm)
tapered_mad <- bdmad(complete_tapered_results, theta, degen_tp_ergm)

RMSE_natural_parameter <- cbind(round(tapered_rmse, 3), round(MPLE_rmse, 3), round(mf_rmse, 3))
rownames(RMSE_natural_parameter) <- names(mv_1)
colnames(RMSE_natural_parameter) <- c("MCMC-MLE","MPLE", "MFVLE")
# These are the RMSE of natural parameter estimates
RMSE_natural_parameter

MAD_natural_parameter <- cbind(round(tapered_mad, 3), round(MPLE_mad, 3), round(mf_mad, 3))
rownames(MAD_natural_parameter) <- names(mv_1)
colnames(MAD_natural_parameter) <- c("MCMC-MLE","MPLE", "MFVLE")
# These are the MAD of natural parameter estimates
MAD_natural_parameter


outliers1 = outliers((complete_tapered_results[,1] - theta[1])^2)
outliers2 = outliers((complete_tapered_results[,2] - theta[2])^2)
outliers3 = outliers((complete_tapered_results[,3] - theta[3])^2)
outliers4 = outliers((complete_tapered_results[,4] - theta[4])^2)
# plot((complete_mf_results[,1] - theta[1])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
# plot((complete_mf_results[,2] - theta[2])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
# plot((complete_mf_results[,3] - theta[3])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
# plot((complete_mf_results[,4] - theta[4])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
# plot(((complete_mf_results[,1] - theta[1])^2)[-outliers1], pch = 19, ylab = "", xlab = "Simulated Networks")
# plot(((complete_mf_results[,2] - theta[2])^2)[-outliers2], pch = 19, ylab = "", xlab = "Simulated Networks")
# plot(((complete_mf_results[,3] - theta[3])^2)[-outliers3], pch = 19, ylab = "", xlab = "Simulated Networks")
# plot(((complete_mf_results[,4] - theta[4])^2)[-outliers4], pch = 19, ylab = "", xlab = "Simulated Networks")


outliers1 = head(order((complete_mf_results[,1] - theta[1])^2, decreasing = T), n = 100)
outliers2 = head(order((complete_mf_results[,2] - theta[2])^2, decreasing = T), n = 100)
outliers3 = head(order((complete_mf_results[,3] - theta[3])^2, decreasing = T), n = 100)
outliers4 = head(order((complete_mf_results[,4] - theta[4])^2, decreasing = T), n = 100)
outliers_all = intersect(intersect(intersect(outliers1, outliers2), outliers3), outliers4)


# plot(density(complete_mf_results[,1]), main = "")
# plot(density(complete_mf_results[,2]), main = "")
# plot(density(complete_mf_results[,3]), main = "")
# plot(density(complete_mf_results[,4]), main = "")


# Now the mean values
complete_mf_results = mf_mv_est #[complete.cases(mf_mv_est),]
complete_tapered_results = ergm_mv_est_tapered #[complete.cases(ergm_mv_est_tapered),]
complete_MPLE_results = MPLE_mv_est #[complete.cases(MPLE_mv_est),]


# Degeneracy
length(intersect(which(is.na(mf_mv_est[,1])), which(is.na(ergm_mv_est_tapered[,1]))))

# MFERGM RMSE

outliers1 = outliers((complete_tapered_results[,1] - mv_s[1])^2)
outliers2 = outliers((complete_tapered_results[,2] - mv_s[2])^2)
outliers3 = outliers((complete_tapered_results[,3] - mv_s[3])^2)
outliers4 = outliers((complete_tapered_results[,4] - mv_s[4])^2)
# plot((complete_mf_results[,1] - mv_s[1])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
# plot((complete_mf_results[,2] - mv_s[2])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
# plot((complete_mf_results[,3] - mv_s[3])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
# plot((complete_mf_results[,4] - mv_s[4])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
# plot(((complete_mf_results[,1] - mv_s[1])^2)[-outliers1], pch = 19, ylab = "", xlab = "Simulated Networks")
# plot(((complete_mf_results[,2] - mv_s[2])^2)[-outliers2], pch = 19, ylab = "", xlab = "Simulated Networks")
# plot(((complete_mf_results[,3] - mv_s[3])^2)[-outliers3], pch = 19, ylab = "", xlab = "Simulated Networks")
# plot(((complete_mf_results[,4] - mv_s[4])^2)[-outliers4], pch = 19, ylab = "", xlab = "Simulated Networks")

# plot(density(complete_mf_results[,1]), main = "")
# plot(density(complete_mf_results[,2]), main = "")
# plot(density(complete_mf_results[,3]), main = "")
# plot(density(complete_mf_results[,4]), main = "")

# RMSE

MPLE_rmse <- bdrmse(complete_MPLE_results, mv_s, degen=NULL)
mf_rmse <- bdrmse(complete_mf_results, mv_s, degen_tp_mfergm)
tapered_rmse <- bdrmse(complete_tapered_results, mv_s, degen_tp_ergm)

# MAD
MPLE_mad <- bdmad(complete_MPLE_results, mv_s, degen=NULL)
mf_mad <- bdmad(complete_mf_results, mv_s, degen_tp_mfergm)
tapered_mad <- bdmad(complete_tapered_results, mv_s, degen_tp_ergm)

RMSE_mean_value_parameter <- cbind(round(tapered_rmse, 3), round(MPLE_rmse, 3), round(mf_rmse, 3))
rownames(RMSE_mean_value_parameter) <- names(mv_1)
colnames(RMSE_mean_value_parameter) <- c("MCMC-MLE","MPLE", "MFVLE")
# These are the RMSE of mean-value parameter estimates
RMSE_mean_value_parameter

MAD_mean_value_parameter <- cbind(round(tapered_mad, 3), round(MPLE_mad, 3), round(mf_mad, 3))
rownames(MAD_mean_value_parameter) <- names(mv_1)
colnames(MAD_mean_value_parameter) <- c("MCMC-MLE","MPLE", "MFVLE")
# These are the MAD of mean-value parameter estimates
MAD_mean_value_parameter

# plot(mfergm_estim[,4],ergm_sim_estim_tapered[,4])

# rm(ergm_sim_list_tapered)
# rm(g_sim)
# save.image("model1_n50_increased.RData")

# Summary
# These are the RMSE of natural parameter estimates
RMSE_natural_parameter
# These are the MAD of natural parameter estimates
MAD_natural_parameter
# These are the RMSE of mean-value parameter estimates
RMSE_mean_value_parameter
# These are the RMSE of mean-value parameter estimates
MAD_mean_value_parameter
