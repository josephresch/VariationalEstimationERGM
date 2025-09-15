rm(list=ls())
my.seed=108         # paper seed 108
set.seed(my.seed)
setwd("~/Desktop/Paper Submission/paper_results/")

library(ergm)
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

nsims       =  1000                                      # number of networks simulated
n           =  50                                        # number of nodes
mv_1        =  c(95.7479, 82.7075, 363.4233, 13.9101)    # long-run estimate
theta       =  c(-2,1,1,1) * c(2,2,1/n,1/n)              # true parameters for model 1
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
names(mv_1) <- names(summary(formula))    # just to get the name

# sim = g
load("sim_model1_n50.RData")             # long term simulated network

set.seed(my.seed)
registerDoParallel(cores)

# Simulate networks
g_sim = foreach(i = 1:nsims) %dorng% {
  simulate(sim ~ edges + nodematch("x") + kstar(2) + triangles, 
           coef = theta
           ,control=control.simulate.formula(MCMC.burnin=1000000, MCMC.interval=1000000)
  )
}

# Extract summary stats from 'g_sim'
g_sim_stats = foreach(i = 1:(nsims), .combine = rbind) %dorng% {
  as.vector(summary_formula(g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles))
}
apply(g_sim_stats,2,mean)
# mv_1 <- c(95.7479, 82.7075, 363.4233, 13.9101) # Long-run estimate


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
  ergm(g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles, 
       estimate = "MPLE")$coefficients
}

##########################
#                        #
#     ERGM (Default)     #
#                        #
##########################

set.seed(my.seed)
### Compute ERGM Results
ergm_sim_list = foreach(i = 1:nsims) %dorng% {
  skip_to_next <- FALSE
  formula <- g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles
  ergm_sim =
    tryCatch(withTimeout(
      ergm(formula, eval.loglik=FALSE,
           control=control.ergm(init = init_params[i,],
                                MCMC.burnin=100000,
                                MCMC.interval=1000
           )),
      timeout = 5*60, onTimeout = "error"),
      error = function(e) {skip_to_next <<- NULL})
  ergm_sim
}

ergm_sim_estim = matrix(0, nrow = nsims, ncol = 4)
for(i in 1:nsims)
{
  if(!is.null(ergm_sim_list[[i]]))
  {
    est.params <- ergm_sim_list[[i]]$coefficients
    ergm_sim_estim[i,] <- est.params
  }
  else {ergm_sim_estim[i,] = c(NA, NA, NA, NA)}
}

ergm_result = apply(ergm_sim_estim[complete.cases(ergm_sim_estim),], 2, mean)
ergm_result

#### Treat degenerately fitted ERGM's
degen_ergm = which(!complete.cases(ergm_sim_estim))
degen_ergm

ergm_sim_estim_filled = ergm_sim_estim

if(length(degen_ergm) > 0 ){
  degen_ergm_tapered = foreach(i = 1:length(degen_ergm)) %dorng% {
    skip_to_next <- FALSE
    formula <- g_sim[[degen_ergm[i]]] ~ edges + nodematch("x") + kstar(2) + triangles
    ergm_sim = 
      tryCatch(withTimeout(
        ergm.tapered(formula, eval.loglik=FALSE,
                     control=control.ergm.tapered(init = init_params[degen_ergm[i],],
                                                  MCMC.burnin=100000,
                                                  MCMC.interval=1000
                     )),
        timeout = 5*60, onTimeout = "error"),
        error = function(e) {skip_to_next <<- NULL})
    ergm_sim
  }
  
  # degen_ergm_tapered_estim = matrix(0, nrow = length(degen), ncol = 4)
  for(i in 1:length(degen_ergm))
  {
    ergm_sim_estim_filled[degen_ergm[i],] <- degen_ergm_tapered[[i]]$coefficients
    #if(!is.null(ergm_degen_list[[i]]))
    #{
    # est.params <- ergm_degen_list[[i]]$coefficients
    # degen_ergm_tapered_estim[i,] <- est.params
    # ergm_sim_estim[degen[i],] <- est.params
    #}
    #else {degen_ergm_tapered_estim[i,] = c(NA, NA, NA, NA)}
  }
}


Sys.sleep(5)
# compare mean-values
ergm_mv_est = foreach(i = 1:nsims, .combine = rbind) %dorng% {
  if (i %in% degen_ergm)
  {
    simulate_ergm.tapered(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles,
                          nsim = 10, 
                          tapering.centers=mv_1, tau=0.25/mv_1,
                          control=control.simulate.formula(MCMC.burnin=100000, MCMC.interval=100000),
                          coef = ergm_sim_estim_filled[i,],         
                          output = "stats")
  } else {
    simulate(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles,
             nsim = 10,
             control=control.simulate.formula(MCMC.burnin=100000, MCMC.interval=100000),
             coef = ergm_sim_estim_filled[i,],         
             output = "stats")
  }
}
colMeans(g_sim_stats)
ergm_mv = colMeans(ergm_mv_est)
ergm_mv



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

# Jitter intial params using mvn
init_mf_params = mvtnorm::rmvnorm(nsims, mean = colMeans(init_params), sigma = 10*cov(init_params))

### Compute MFERGM Results
mfergm_estim <- NULL
# Double loop to reduce memory use
Sys.sleep(5)
for(j in 1:(nsims/100)){
  mfergm_estim_j = foreach(i = (100*(j-1)+(1:100)), .combine = rbind) %dorng% {
    pars <- init_mf_params[i,] * c(.5,.5,n,n)
    addpars <- list(n=n, tobs = tobs[i,], x=x, ninit=1)
    cd.est <- tryCatch(optimx(pars, fn = loglikmf.model4, 
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
degen_mfergm = which(!complete.cases(mfergm_estim))
degen_mfergm

mfergm_estim_filled = mfergm_estim

if(length(degen_mfergm) > 0){
  # Replace with MPLE
  for(i in 1:length(degen_mfergm)){
    mfergm_estim_filled[degen_mfergm[i],] <- init_params[degen_mfergm[i],]
  }
}

# compare mean values
Sys.sleep(5)
set.seed(1)
mf_mv_est = foreach(i = 1:nsims, .combine = rbind) %dorng% {
  simulate(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles,
           nsim = 10,
           control=control.simulate.formula(MCMC.burnin=100000, MCMC.interval=100000),
           coef =  mfergm_estim_filled[i,],    
           output = "stats"
  )
}
mfergm_mv = colMeans(mf_mv_est)
mfergm_mv

# compare mean values
Sys.sleep(5)
set.seed(1)
MPLE_mv_est = foreach(i = 1:nsims, .combine = rbind) %dorng% {
  simulate(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles,
           nsim = 10,
           coef = init_params[i,],    
           control=control.simulate.formula(MCMC.burnin=100000, MCMC.interval=100000),
           output = "stats"
  )
}
MPLE_mv = colMeans(MPLE_mv_est)
MPLE_mv
MPLE_result = apply(init_params, 2, mean, na.rm = T) 

stopImplicitCluster()


mean_est <- cbind(theta, ergm_result, MPLE_result, mfergm_result)
rownames(mean_est) <- names(ergm_mv)
colnames(mean_est) <- c("true", "MCMC-MLE","MPLE","MFVLE")
mean_est

mean_mv_est <- cbind(round(colMeans(g_sim_stats), 2), round(ergm_mv,2), round(MPLE_mv,2), round(mfergm_mv, 2))
rownames(mean_mv_est) <- names(ergm_mv)
colnames(mean_mv_est) <- c("true", "MCMC-MLE","MPLE", "MFVLE")
mean_mv_est

# rm(ergm_sim_list_tapered)
# save.image("model1_n50.RData")


########################## Figure 1 and Figure 2 ########################## 
library(ggplot2)
### MCMLE Estimate Distributions
ggplot(data.frame(Edge = ergm_sim_estim[,1]),  
       aes(x=Edge))+
  geom_histogram(color="darkblue", fill="lightblue", bins = 20) + 
  ylim(0,175) +
  geom_vline(xintercept=-4, color ="orange", size=3, linetype="longdash") +
  ylab("Count") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))
  
ggplot(data.frame(Homophily = ergm_sim_estim[,2]), 
       aes(x=Homophily))+
  geom_histogram(color="darkblue", fill="lightblue", bins = 20) + 
  ylim(0,175) +
  geom_vline(xintercept=2, color ="orange", size=3, linetype="longdash") +
  ylab("Count") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data.frame(Two_Star = ergm_sim_estim[,3]), 
       aes(x=Two_Star))+
  geom_histogram(color="darkblue", fill="lightblue", bins = 20) + 
  ylim(0,175) +
  geom_vline(xintercept=.01, color ="orange", size=3, linetype="longdash") +
  ylab("Count") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data.frame(Triangle = ergm_sim_estim[,4]), 
       aes(x=Triangle))+
  geom_histogram(color="darkblue", fill="lightblue", bins = 20) + 
  ylim(0,175) +
  geom_vline(xintercept=.01, color ="orange", size=3, linetype="longdash") +
  ylab("Count") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))

# hist(ergm_sim_estim[,1], main = NULL, xlab = "Edge")
# hist(mfergm_estim[,1], main = NULL, xlab = "Triangle")

### MFERGM Estimate Distributions
ggplot(data.frame(Edge = mfergm_estim[,1][mfergm_estim[,1] < -1 & mfergm_estim[,1] > -7 ]), 
       aes(x=Edge))+
  geom_histogram(color="darkblue", fill="lemonchiffon", bins = 20) + 
  ylim(0,210) +
  geom_vline(xintercept=-4, color ="orange", size=3, linetype="longdash") +
  ylab("Count") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data.frame(Homophily = mfergm_estim[,2][mfergm_estim[,2] < 4 & mfergm_estim[,2] > 0 ]), 
       aes(x=Homophily))+
  geom_histogram(color="darkblue", fill="lemonchiffon", bins = 20) + 
  ylim(0,210) +
  geom_vline(xintercept=2, color ="orange", size=3, linetype="longdash") +
  ylab("Count") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data.frame(Two_Star = mfergm_estim[,3][mfergm_estim[,3] < 0.2 & mfergm_estim[,3] > -0.2 ]), 
       aes(x=Two_Star))+
  geom_histogram(color="darkblue", fill="lemonchiffon", bins = 20) + 
  ylim(0,210) +
  geom_vline(xintercept=.01, color ="orange", size=3, linetype="longdash") +
  ylab("Count") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data.frame(Triangle = mfergm_estim[,4][mfergm_estim[,4] < 0.5 & mfergm_estim[,4] > -0.5 ]), 
       aes(x=Triangle))+
  geom_histogram(color="darkblue", fill="lemonchiffon", bins = 20) + 
  ylim(0,210) +
  geom_vline(xintercept=.01, color ="orange", size=3, linetype="longdash") +
  ylab("Count") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))

### MPLE Estimate Distributions
ggplot(data.frame(Edge = init_params[,1]),  
       aes(x=Edge))+
  geom_histogram(color="darkblue", fill="lightgreen", bins = 20) + 
  ylim(0,175) +
  geom_vline(xintercept=-4, color ="orange", size=3, linetype="longdash") +
  ylab("Count") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data.frame(Homophily = init_params[,2]), 
       aes(x=Homophily))+
  geom_histogram(color="darkblue", fill="lightgreen", bins = 20) + 
  ylim(0,175) +
  geom_vline(xintercept=2, color ="orange", size=3, linetype="longdash") +
  ylab("Count") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data.frame(Two_Star = init_params[,3]), 
       aes(x=Two_Star))+
  geom_histogram(color="darkblue", fill="lightgreen", bins = 20) + 
  ylim(0,175) +
  geom_vline(xintercept=.01, color ="orange", size=3, linetype="longdash") +
  ylab("Count") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data.frame(Triangle = init_params[,4]), 
       aes(x=Triangle))+
  geom_histogram(color="darkblue", fill="lightgreen", bins = 20) + 
  ylim(0,175) +
  geom_vline(xintercept=.01, color ="orange", size=3, linetype="longdash") +
  ylab("Count") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))




# function for outliers
outliers <- function(x, field = 1.5, na.rm = TRUE) 
{
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
  H <- field * IQR(x, na.rm = na.rm)
  which(x < (qnt[1] - H) | x > (qnt[2] + H))
}

# replace degen and outlier fit values
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
  outliers1 = union(outliers((comp[,1] - theta[1])^2), degen)
  outliers2 = union(outliers((comp[,2] - theta[2])^2), degen)
  outliers3 = union(outliers((comp[,3] - theta[3])^2), degen)
  outliers4 = union(outliers((comp[,4] - theta[4])^2), degen)
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


# First natural params
complete_mf_results = mfergm_estim #[complete.cases(mfergm_estim),]
complete_ergm_results = ergm_sim_estim #[complete.cases(ergm_sim_estim),]
complete_MPLE_results = init_params #[complete.cases(init_params),]

# Degeneracy
length(intersect(which(is.na(mfergm_estim[,1])), which(is.na(ergm_sim_estim[,1]))))

# RMSE
MPLE_rmse <- bdrmse(complete_MPLE_results, theta,degen=NULL)
mf_rmse <- bdrmse(complete_mf_results, theta, degen_mfergm)
ergm_rmse <- bdrmse(complete_ergm_results, theta, degen_ergm)

# MAD
MPLE_mad <- bdmad(complete_MPLE_results, theta,degen=NULL)
mf_mad <- bdmad(complete_mf_results, theta, degen_mfergm)
ergm_mad <- bdmad(complete_ergm_results, theta, degen_ergm)

RMSE_natural_parameter <- cbind(round(ergm_rmse, 3), round(MPLE_rmse, 3), round(mf_rmse, 3))
rownames(RMSE_natural_parameter) <- names(mv_1)
colnames(RMSE_natural_parameter) <- c("MCMC-MLE","MPLE", "MFVLE")
# These are the RMSE of natural parameter estimates
RMSE_natural_parameter

MAD_natural_parameter <- cbind(round(ergm_mad, 3), round(MPLE_mad, 3), round(mf_mad, 3))
rownames(MAD_natural_parameter) <- names(mv_1)
colnames(MAD_natural_parameter) <- c("MCMC-MLE","MPLE", "MFVLE")
# These are the MAD of natural parameter estimates
MAD_natural_parameter


outliers1 = outliers((complete_ergm_results[,1] - theta[1])^2)
outliers2 = outliers((complete_ergm_results[,2] - theta[2])^2)
outliers3 = outliers((complete_ergm_results[,3] - theta[3])^2)
outliers4 = outliers((complete_ergm_results[,4] - theta[4])^2)
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
complete_ergm_results = ergm_mv_est #[complete.cases(ergm_mv_est),]
complete_MPLE_results = MPLE_mv_est #[complete.cases(MPLE_mv_est),]


# Degeneracy
length(intersect(which(is.na(mf_mv_est[,1])), which(is.na(ergm_mv_est[,1]))))

# MFERGM RMSE

outliers1 = outliers((complete_ergm_results[,1] - mv_s[1])^2)
outliers2 = outliers((complete_ergm_results[,2] - mv_s[2])^2)
outliers3 = outliers((complete_ergm_results[,3] - mv_s[3])^2)
outliers4 = outliers((complete_ergm_results[,4] - mv_s[4])^2)
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
mf_rmse <- bdrmse(complete_mf_results, mv_s, degen_mfergm)
ergm_rmse <- bdrmse(complete_ergm_results, mv_s, degen_ergm)

# MAD
MPLE_mad <- bdmad(complete_MPLE_results, mv_s, degen=NULL)
mf_mad <- bdmad(complete_mf_results, mv_s, degen_mfergm)
ergm_mad <- bdmad(complete_ergm_results, mv_s, degen_ergm)

RMSE_mean_value_parameter <- cbind(round(ergm_rmse, 3), round(MPLE_rmse, 3), round(mf_rmse, 3))
rownames(RMSE_mean_value_parameter) <- names(mv_1)
colnames(RMSE_mean_value_parameter) <- c("MCMC-MLE","MPLE", "MFVLE")
# These are the RMSE of mean-value parameter estimates
RMSE_mean_value_parameter

MAD_mean_value_parameter <- cbind(round(ergm_mad, 3), round(MPLE_mad, 3), round(mf_mad, 3))
rownames(MAD_mean_value_parameter) <- names(mv_1)
colnames(MAD_mean_value_parameter) <- c("MCMC-MLE","MPLE", "MFVLE")
# These are the MAD of mean-value parameter estimates
MAD_mean_value_parameter

# plot(mfergm_estim[,4],ergm_sim_estim[,4])

# rm(ergm_sim_list)
# rm(g_sim)

# Summary
length(degen_ergm) / nsims
length(degen_mfergm) / nsims
# These are the RMSE of natural parameter estimates
RMSE_natural_parameter
# These are the MAD of natural parameter estimates
MAD_natural_parameter
# These are the RMSE of mean-value parameter estimates
RMSE_mean_value_parameter
# These are the RMSE of mean-value parameter estimates
MAD_mean_value_parameter

