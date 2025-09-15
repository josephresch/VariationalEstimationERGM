rm(list=ls())
my.seed=1
set.seed(my.seed)

pdf("model4_sim1.tapered.pdf")

library(ergm.tapered)
library(mfergm)
library(optimx)     # mfergm likelihood
library(R.utils)    # set time on ergm
library(doParallel) # parallel loops using 'foreach'

#####################################################################
#                                                                   #
#     Create High Transitivity ERGM Params (using ERGM fitting)     #
#                                                                   #
#####################################################################

nsims       =  100                              # number of networks simulated
n           =  100                                # number of nodes
#theta       =  c(-3,2,1,3) * c(2,2,1/n,1/n)      # true parameters for model 2
theta       =  c(-2,1,1,1) * c(2,2,1/n,1/n)      # true parameters for model 1

##################
#                #
#     Set-up     #
#                #
##################

g <- initialize.network(theta, n, directed = FALSE)
x <- rbinom(n, 1, 0.5) # attributes
set.vertex.attribute(g, # the name of the network object
                     "x", # the name we want to reference the variable by in that object
                     x # the value we are giving that variable
) 

# Simulated networks 'g_sim' used for both ERGM and MFERGM
g_sim <- simulate_ergm.tapered(g ~ edges + nodematch("x") + kstar(2) + triangles, 
                  nsim = 10*nsims, tapering.centers=mv_i, tau=fit$tapering.coefficients,
                  coef = coef(fit),
		  control=control.simulate.formula(MCMC.burnin=100000, MCMC.interval=10000)
)
mv_1=apply(g_sim,2,mean))
formula <- g ~ edges + nodematch("x") + kstar(2) + triangles
theta

theta[4] <- theta[4]*10

# Simulated networks 'g_sim' used for both ERGM and MFERGM
g_sim <- simulate_ergm.tapered(g_sim[[nsim]] ~ edges + nodematch("x") + kstar(2) + triangles, 
                  nsim = 10*nsims, tapering.centers=mv_i, tau=fit$tapering.coefficients,
                  coef = coef(fit),
		  control=control.simulate.formula(MCMC.burnin=100000, MCMC.interval=10000)
                              )

# Extract summary stats from 'g_sim'
g_sim_stats = matrix(nrow = 1*nsims, ncol = 4)
for (i in 1:(1*nsims))
{
  g_sim_stats[i,] = summary_formula(g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles)
}

cbind(mv_1, apply(g_sim_stats,2,mean))
save.image("mfergm_params2_n100.RData")

# Initialize the parameters at MPLE
init_params = matrix(nrow = nsims, ncol = length(theta))
for (i in 1:nsims)
{
  init_params[i,] = ergm(g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles, 
                         estimate = "MPLE")$coefficients
}


##########################
#                        #
#     ERGM (Default)     #
#                        #
##########################

### Compute ERGM Results
#registerDoParallel(detectCores())
registerDoParallel(10)
ergm_sim_list = foreach(i = 1:nsims) %dopar% {
  # cat("***********************************\n")
  # cat(paste("estimating sample" ,i, "\n"))
  # cat("***********************************\n")
  skip_to_next <- FALSE
  formula <- g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles
  ergm_sim = 
    tryCatch(withTimeout(
                          ergm(formula, eval.loglik=FALSE, estimate="MPLE",
                               control=control.ergm(init = init_params[i,], MCMLE.confidence=0.95)),
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

ergm_result = apply(ergm_sim_estim, 2, median, na.rm = T)
ergm_result
theta
ergm_q = apply(ergm_sim_estim, 2, quantile, c(0.01,0.99), na.rm = T)

#### DEGENERACIES
degen = which(!complete.cases(ergm_sim_estim))
# 58 of em
degen

if(length(degen) > 0 ){
degen_list = ergm_sim_list[[]]
ergm_degen_list = foreach(i = 1:length(degen)) %dopar% {
  skip_to_next <- FALSE
  formula <- g_sim[[degen[i]]] ~ edges + nodematch("x") + kstar(2) + triangles
  ergm_sim = tryCatch(withTimeout(ergm(formula, eval.loglik=FALSE,
                                        control=control.ergm(init=theta, MCMLE.confidence=0.95
                                       #   # MCMC.burnin=100000,
                                       #   # MCMC.interval=1000,
                                       #   # MCMC.samplesize = 5000,
                                          )
  ), timeout = 5*60, onTimeout = "error"),
  error = function(e) {skip_to_next <<- NULL})
  ergm_sim
}

ergm_degen_estim = matrix(0, nrow = length(degen), ncol = 4)
for(i in 1:length(degen))
{
  if(!is.null(ergm_degen_list[[i]]))
  {
    est.params <- ergm_degen_list[[i]]$coefficients
    ergm_degen_estim[i,] <- est.params
    ergm_sim_estim[degen[i],] <- est.params
  }
  else {ergm_degen_estim[i,] = c(NA, NA, NA, NA)}
}
}

ergm_sim_estim
# compare mean-values
g3 = foreach(i = 1:nsims, .combine = rbind) %dopar% {
simulate_ergm.tapered(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles,
               nsim = 10,
               control=control.simulate.formula(MCMC.burnin=100000, MCMC.interval=10000),
               coef = ergm_sim_estim[i,],         
               output = "stats"
)
}
colMeans(g_sim_stats)
ergm_mv = colMeans(g3)

### Compute ERGM Results
ergm_sim_list_tapered = foreach(i = 1:nsims) %dopar% {
  # cat("***********************************\n")
  # cat(paste("estimating sample" ,i, "\n"))
  # cat("***********************************\n")
  skip_to_next <- FALSE
  formula <- g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles
  ergm_sim = 
    tryCatch(withTimeout(
                          ergm.tapered(formula, eval.loglik=FALSE,
                               control=control.ergm.tapered(init = init_params[i,])),
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

#### DEGENERACIES
degen_tapered = which(!complete.cases(ergm_sim_estim_tapered))
# 58 of em
degen_tapered

if(length(degen_tapered) > 0 ){
degen_list_tapered = ergm_sim_list_tapered[[]]
ergm_degen_list_tapered = foreach(i = 1:length(degen_tapered)) %dopar% {
  skip_to_next <- FALSE
  formula <- g_sim[[degen[i]]] ~ edges + nodematch("x") + kstar(2) + triangles
  ergm_sim_tapered = tryCatch(withTimeout(ergm(formula, eval.loglik=FALSE,
                                        control=control.ergm(init=theta, MCMLE.confidence=0.95
                                       #   # MCMC.burnin=100000,
                                       #   # MCMC.interval=1000,
                                       #   # MCMC.samplesize = 5000,
                                          )
  ), timeout = 5*60, onTimeout = "error"),
  error = function(e) {skip_to_next <<- NULL})
  ergm_sim_tapered
}
ergm_degen_list_tapered

ergm_degen_estim_tapered = matrix(0, nrow = length(degen), ncol = 4)
for(i in 1:length(degen_tapered))
{
  if(!is.null(ergm_degen_list_tapered[[i]]))
  {
    est.params <- ergm_degen_list_tapered[[i]]$coefficients
    ergm_degen_estim_tapered[i,] <- est.params
    ergm_sim_estim_tapered[degen_tapered[i],] <- est.params
  }
  else {ergm_degen_estim_tapered[i,] = c(NA, NA, NA, NA)}
}
}

ergm_result_tapered = apply(ergm_sim_estim_tapered, 2, median, na.rm = T)
ergm_result_tapered

ergm_sim_estim_tapered

# compare mean-values
g3_tapered = foreach(i = 1:nsims, .combine = rbind) %dopar% {
#g3_tapered <- simulate(g_sim[[1]] ~ edges+nodematch("x")+ kstar(2) + triangles,
simulate_ergm.tapered(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles,
               nsim = 10,
               control=control.simulate.formula(MCMC.burnin=100000, MCMC.interval=10000),
               coef = ergm_sim_estim_tapered[i,],         
               output = "stats"
)
}
colMeans(g_sim_stats)
ergm_mv_tapered = colMeans(g3_tapered)
ergm_mv_tapered

# Model statistics ‘kstar2’ and ‘triangle’ are not varying. This may indicate 
# that the observed data occupies an extreme point in the sample space or that 
# the estimation has reached a dead-end configuration.
# Warning: Model statistics ‘nodematch.x’ are linear combinations of some set of 
# preceding statistics at the current stage of the estimation. This may indicate 
# that the model is nonidentifiable.
# Error in ergm.MCMLE(init, nw, model, initialfit = (initialfit <- NULL),  : 
# Unconstrained MCMC sampling did not mix at all. Optimization cannot continue.
# In addition: Warning message:
# In ergm_MCMC_sample(s, control, theta = mcmc.init, verbose = max(verbose -  :
# Unable to reach target effective size in iterations alotted.


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

### Compute MFERGM Results
mfergm_estim = foreach(i = 1:nsims, .combine = rbind) %dopar% {
  pars <- init_params[i,] * c(.5,.5,n,n)
  addpars <- list(n=n, tobs = tobs[i,], x=x, ninit=1)
  cd.est <- tryCatch(optimx(pars, fn = loglikmf.model5, 
                            method = "BFGS", 
                            control = list(fnscale = -1), addpars = addpars), 
                     error = function(e) {return(c(NA, NA, NA, NA))})
  as.numeric(cd.est[1:4])
}
mfergm_estim = matrix(c(mfergm_estim[,1]*2, mfergm_estim[,2]*2, mfergm_estim[,3]/n, mfergm_estim[,4]/n),
               nrow = nsims)
mfergm_result_med = apply(mfergm_estim, 2, median, na.rm = T) 
mfergm_result = apply(mfergm_estim, 2, median, na.rm = T) 
cbind(theta,mfergm_result_med,mfergm_result)

#### DEGENERACIES
degen_mfergm = which(!complete.cases(mfergm_estim))
# 58 of em
degen_mfergm

if(length(degen_mfergm) > 0){
# Replace with MPLE
for(i in 1:length(degen_mfergm)){
  mfergm_estim[degen_mfergm[i],] <- init_params[degen_mfergm[i],]
}
}

mfergm_estim_bd <- mfergm_estim
for(i in 1:nsims){
 a = mfergm_estim_bd[i,] < ergm_q[1,]
 mfergm_estim_bd[i,a] <- ergm_q[1,a]
 a = mfergm_estim_bd[i,] > ergm_q[2,]
 mfergm_estim_bd[i,a] <- ergm_q[2,a]
}

# compare mean values
g5 = foreach(i = 1:nsims, .combine = rbind) %dopar% {
 simulate_ergm.tapered(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles,
               nsim = 10,
               control=control.simulate.formula(MCMC.burnin=100000, MCMC.interval=10000),
               coef = init_params[i,],    
               output = "stats"
)
}
MPLE_mv = colMeans(g5)
MPLE_result = apply(init_params, 2, median, na.rm = T) 

# compare mean values
g4 = foreach(i = 1:nsims, .combine = rbind) %dopar% {
 simulate(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles,
               nsim = 10,
               coef = mfergm_estim_bd[i,],    
               control=control.simulate.formula(MCMC.burnin=100000, MCMC.interval=10000),
               output = "stats"
)
}
colMeans(g_sim_stats)
mfergm_mv = colMeans(g4)
mfergm_mv

a <- cbind(theta, ergm_result, ergm_result_tapered, MPLE_result, mfergm_result)
rownames(a) <- names(ergm_mv)
colnames(a) <- c("true", "MCMC-MLE","Tapered","MPLE","MFVLE")
a

# compare mfergm result with "theta"
theta
ergm_result
ergm_result_tapered
mfergm_result


# results for report
round(theta, 3) 
round(colMeans(g_sim_stats), 2)
round(colMeans(init_params), 3)
round(ergm_result, 3) 
round(ergm_mv, 2) 
round(ergm_result_tapered, 3) 
round(ergm_mv_tapered, 2) 
round(mfergm_result, 3) 
round(mfergm_mv, 2) 
# / c(2,2,1/n,1/n)


a <- cbind(round(colMeans(g_sim_stats), 2), round(ergm_mv, 2), round(ergm_mv_tapered,2), round(MPLE_mv,2), round(mfergm_mv, 2))
rownames(a) <- names(ergm_mv)
colnames(a) <- c("true", "MCMC-MLE","Tapered","MPLE", "MFVLE")
a

# save.image("mfergm_params2_n100.RData")

hist(ergm_sim_estim[,4], main = NULL, xlab = "Triangle")
hist(ergm_sim_estim_tapered[,4], main = NULL, xlab = "Triangle")
hist(mfergm_estim[,4], main = NULL, xlab = "Triangle")


outliers <- function(x, field = 1.5, na.rm = TRUE) 
{
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
  H <- field * IQR(x, na.rm = na.rm)
  which(x < (qnt[1] - H) | x > (qnt[2] + H))
}



a  = is.na(ergm_sim_estim[,4])
b  = is.na(mfergm_estim[,4])
c = mfergm_estim[(!a)&(!b),]
c1 = c[,4]
d = ergm_sim_estim[(!a)&(!b),]
d1 = d[,4]
c1[c1 < -2] = d1[c1 < -2]
plot(c1~d1, col = (c1 == d1) + 1, pch = 16)
sum(c1 == d1)

f = c[-Reduce(union, list(outliers(c[,1], 0.6), outliers(c[,2], 0.3), outliers(c[,3], 1.5), outliers(c[,4], 1.5))),]



gtest_mf = matrix(nrow = nrow(c), ncol = 4)
for (i in 1:nrow(c))
{
  gtest_mf[i,] = simulate(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles, 
                          nsim = 1, 
                          coef = c[i,],        
                          output = "stats")
}
gtest_mc = matrix(nrow = nrow(d), ncol = 4)
for (i in 1:nrow(d))
{
  gtest_mc[i,] = simulate(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles, 
                          nsim = 1, 
                          coef = d[i,],        
                          output = "stats")
}
gtest_mc_tapered = matrix(nrow = nrow(d), ncol = 4)
for (i in 1:nrow(d))
{
  gtest_mc_tapered[i,] = as.vector(simulate_ergm.tapered(g_sim[[i]] ~ edges+nodematch("x")+ kstar(2) + triangles, 
                          nsim = 1, 
                          coef = d[i,],        
                          output = "stats"))
}
colMeans(gtest_mf)
colMeans(gtest_mc)
colMeans(gtest_mc_tapered)

summary(ergm_sim_estim)
summary(gtest_mf)
summary(gtest_mc)
summary(g_sim_stats)


gtest_mf_tri = gtest_mf[gtest_mf[,4]<60, 4]

hist(ergm_sim_estim[,1])
hist(ergm_sim_estim[,2])
hist(ergm_sim_estim[,3])
hist(ergm_sim_estim[,4])
hist(f[,1])
hist(f[,2])
hist(f[,3])
hist(f[,4])


gtest_mf_rm = gtest_mf[-Reduce(union, list(outliers(gtest_mf[,1], 0.7), outliers(gtest_mf[,2], 1), 
                                           outliers(gtest_mf[,3], 0.001), outliers(gtest_mf[,4], 0.001))),]
hist(gtest_mf[,1])
hist(gtest_mf[,2])
hist(gtest_mf[,3])
hist(gtest_mf[,4])

hist(ergm_sim_estim_tapered[,1])
hist(ergm_sim_estim_tapered[,2])
hist(ergm_sim_estim_tapered[,3])
hist(ergm_sim_estim_tapered[,4])







complete_mf_results = mfergm_estim[complete.cases(mfergm_estim),]
complete_mcmc_results = ergm_sim_estim[complete.cases(ergm_sim_estim),]
complete_tapered_results = ergm_sim_estim_tapered[complete.cases(ergm_sim_estim_tapered),]


# Degeneracy
length(intersect(which(is.na(mfergm_estim[,1])), which(is.na(ergm_sim_estim[,1]))))

# MFERGM RMSE

# mf
outliers1 = outliers((complete_mf_results[,1] - theta[1])^2)
outliers2 = outliers((complete_mf_results[,2] - theta[2])^2)
outliers3 = outliers((complete_mf_results[,3] - theta[3])^2)
outliers4 = outliers((complete_mf_results[,4] - theta[4])^2)
outliers_all = intersect(intersect(intersect(outliers1, outliers2), outliers3), outliers4)
mf_rmse1 = sqrt(sum(((complete_mf_results[,1] - theta[1])^2)[-outliers1]))
mf_rmse2 = sqrt(sum(((complete_mf_results[,2] - theta[2])^2)[-outliers2]))
mf_rmse3 = sqrt(sum(((complete_mf_results[,3] - theta[3])^2)[-outliers3]))
mf_rmse4 = sqrt(sum(((complete_mf_results[,4] - theta[4])^2)[-outliers4]))
c(mf_rmse1, mf_rmse2, mf_rmse3, mf_rmse4)

# mcmc
outliers1 = outliers((complete_mcmc_results[,1] - theta[1])^2)
outliers2 = outliers((complete_mcmc_results[,2] - theta[2])^2)
outliers3 = outliers((complete_mcmc_results[,3] - theta[3])^2)
outliers4 = outliers((complete_mcmc_results[,4] - theta[4])^2)
outliers_all = intersect(intersect(intersect(outliers1, outliers2), outliers3), outliers4)
mcmc_rmse1 = sqrt(sum(((complete_mcmc_results[,1] - theta[1])^2)[-outliers1]))
mcmc_rmse2 = sqrt(sum(((complete_mcmc_results[,2] - theta[2])^2)[-outliers2]))
mcmc_rmse3 = sqrt(sum(((complete_mcmc_results[,3] - theta[3])^2)[-outliers3]))
mcmc_rmse4 = sqrt(sum(((complete_mcmc_results[,4] - theta[4])^2)[-outliers4]))
c(mcmc_rmse1, mcmc_rmse2, mcmc_rmse3, mcmc_rmse4)

# tapered
outliers1 = outliers((complete_tapered_results[,1] - theta[1])^2)
outliers2 = outliers((complete_tapered_results[,2] - theta[2])^2)
outliers3 = outliers((complete_tapered_results[,3] - theta[3])^2)
outliers4 = outliers((complete_tapered_results[,4] - theta[4])^2)
outliers_all = intersect(intersect(intersect(outliers1, outliers2), outliers3), outliers4)
tapered_rmse1 = sqrt(sum(((complete_tapered_results[,1] - theta[1])^2)[-outliers1]))
tapered_rmse2 = sqrt(sum(((complete_tapered_results[,2] - theta[2])^2)[-outliers2]))
tapered_rmse3 = sqrt(sum(((complete_tapered_results[,3] - theta[3])^2)[-outliers3]))
tapered_rmse4 = sqrt(sum(((complete_tapered_results[,4] - theta[4])^2)[-outliers4]))
c(tapered_rmse1, tapered_rmse2, tapered_rmse3, tapered_rmse4)

plot((complete_mf_results[,1] - theta[1])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot((complete_mf_results[,2] - theta[2])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot((complete_mf_results[,3] - theta[3])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot((complete_mf_results[,4] - theta[4])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mf_results[,1] - theta[1])^2)[-outliers1], pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mf_results[,2] - theta[2])^2)[-outliers2], pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mf_results[,3] - theta[3])^2)[-outliers3], pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mf_results[,4] - theta[4])^2)[-outliers4], pch = 19, ylab = "", xlab = "Simulated Networks")


outliers1 = head(order((complete_mf_results[,1] - theta[1])^2, decreasing = T), n = 100)
outliers2 = head(order((complete_mf_results[,2] - theta[2])^2, decreasing = T), n = 100)
outliers3 = head(order((complete_mf_results[,3] - theta[3])^2, decreasing = T), n = 100)
outliers4 = head(order((complete_mf_results[,4] - theta[4])^2, decreasing = T), n = 100)
outliers_all = intersect(intersect(intersect(outliers1, outliers2), outliers3), outliers4)


plot(density(complete_mf_results[,1]), main = "")
plot(density(complete_mf_results[,2]), main = "")
plot(density(complete_mf_results[,3]), main = "")
plot(density(complete_mf_results[,4]), main = "")



# MAD

mcmc_mad = matrix(nrow = nrow(complete_mcmc_results), ncol = 4)
for (i in 1:nrow(complete_mcmc_results))
{
  mcmc_mad[i,] = abs(complete_mcmc_results[i,] - theta)
}
round(apply(mcmc_mad, 2, median), 3)

tapered_mad = matrix(nrow = nrow(complete_tapered_results), ncol = 4)
for (i in 1:nrow(complete_tapered_results))
{
  tapered_mad[i,] = abs(complete_tapered_results[i,] - theta)
}
round(apply(tapered_mad, 2, median), 3)

mf_mad = matrix(nrow = nrow(complete_mf_results), ncol = 4)
for (i in 1:nrow(complete_mf_results))
{
  mf_mad[i,] = abs(complete_mf_results[i,] - theta)
}
round(apply(mf_mad, 2, median), 3)

stopImplicitCluster()

save.image("mfergm_params2_n1000.RData")
