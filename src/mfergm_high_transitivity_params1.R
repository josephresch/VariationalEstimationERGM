rm(list=ls())
set.seed(1000)

library(ergm)
library(mfergm)
library(optimx)     # mfergm likelihood
library(R.utils)    # set time on ergm
library(doParallel) # parallel loops using 'foreach'







#####################################################################
#                                                                   #
#     Create High Transitivity ERGM Params (using ERGM fitting)     #
#                                                                   #
#####################################################################

ht_factor   =  1.5
nsims       =  1000                             # number simulated
ht_theta    =  c(NA, NA, NA, NA)
estim.table =  matrix(NA, nrow = nsims, ncol = 4)
n           =  50
theta       =  c(-2,1,1,1) * c(2,2,1/n,1/n)      # true parameters for model
degeneracies = c()

g <- initialize.network(theta, n, directed = FALSE)
x <- rbinom(n,1,.5) # attributes
set.vertex.attribute(g, # the name of the network object
                     "x", # the name we want to reference the variable by in that object
                     x # the value we are giving that variable
) 


g1 <- simulate(g ~ edges + nodematch("x")+ kstar(2) + triangles, 
               nsim = nsims,
               coef = theta
)


mv_params1 = t(sapply(1:nsims, function(i) 
  summary_formula(g1[[i]] ~ edges + nodematch("x") + kstar(2) + triangles)))
mv_params1_avg = colMeans(mv_params1)


while (class(ht_theta) == "logical" | sum(is.na(estim.table[,1])) > nsims*0.05)
{
  mv_params2_avg = mv_params1_avg * c(1, 1, 1, ht_factor)
  
  ##### ERGM
  registerDoParallel(14)
  ergm_list_ht = foreach(i = 1:nsims) %dopar% {
    skip_to_next <- FALSE
    formula <- g1[[i]] ~ edges + nodematch("x") + kstar(2) + triangles
    ergm_it = tryCatch(withTimeout(ergm(formula,
                                        target.stats = mv_params2_avg), timeout = 1000, onTimeout = "error"), 
                       error = function(e) {skip_to_next <<- NULL})
    ergm_it
  }
  stopImplicitCluster()
  
  estim.table = matrix(0, nrow = nsims, ncol = 4)
  for(i in 1:nsims)
  {
    if(!is.null(ergm_list_ht[[i]]))
    {
      est.params <- ergm_list_ht[[i]]$coefficients
      estim.table[i,] <- est.params
    }
    else {estim.table[i,] = c(NA, NA, NA, NA)}
  }
  
  degeneracies = c(degeneracies, sum(is.na(estim.table[,1])))
  ht_factor   =  ht_factor - 0.05 # factor of increased transitivity for triangle term
  
  # here, we change the mean-value parameterization back to ergm coefs
  ht_theta = apply(estim.table, 2, median, na.rm = TRUE)
}




######################################################################################################################################

##################
#                #
#     Set-up     #
#                #
##################

set.seed(1000)

g <- initialize.network(ht_theta, n, directed = FALSE)
x <- rbinom(n, 1, 0.5) # attributes
set.vertex.attribute(g, # the name of the network object
                     "x", # the name we want to reference the variable by in that object
                     x # the value we are giving that variable
) 

# Simulated networks 'g_sim' using 'ht_theta'
g_sim <- simulate(g ~ edges + nodematch("x")+ kstar(2) + triangles, 
                  nsim = nsims,
                  coef = ht_theta
)


# Extract summary stats from 'g_sim'
g_sim_stats = t(sapply(1:nsims, function(i) 
  summary(g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles)))

colMeans(g_sim_stats)
mv_params2_avg

# Initialize the parameters at MPLE
init_params = t(sapply(1:nsims, function(i) 
  ergm(g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles, 
       estimate = "MPLE")$coefficients))



##########################
#                        #
#     ERGM (Default)     #
#                        #
##########################

### Compute ERGM Results
registerDoParallel(14)
ergm_sim_list = foreach(i = 1:nsims) %dopar% {
  skip_to_next <- FALSE
  formula <- g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles
  ergm_sim = tryCatch(withTimeout(ergm(formula,
                                       # control=control.ergm(
                                       #   MCMC.burnin=100000,
                                       #   MCMC.interval=1000,
                                       #   init = init_params[i,])
                                       target.stats = colMeans(g_sim_stats)
  ), 
  timeout = 1000, onTimeout = "error"), 
  error = function(e) {skip_to_next <<- NULL})
  ergm_sim
}
stopImplicitCluster()


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

# compare mean-values
g3 <- simulate(g ~ edges + nodematch("x")+ kstar(2) + triangles, 
               nsim = nsims,
               coef = ergm_result,
               output = "stats"
)
colMeans(g_sim_stats)
ergm_mv = colMeans(g3)
ergm_mv



##################
#                #
#     MFERGM     #
#                #
##################

# Use "ht_theta" on mfergm instead
tobs <- data.frame(matrix(NA, ncol = 4, nrow = nsims))
names(tobs) <- c("edges", "x", "kstar2", "triangles")
for (i in 1:nsims) 
{
  tobs[i,] = g_sim_stats[i,]/(c((n^2)/2, (
    n^2)/2, (n^3)/2 , (n^3)/2 ))  
}

### Compute MFERGM Results
registerDoParallel(14)
mfergm_estim = foreach(i = 1:nsims, .combine = rbind) %dopar% {
  pars <- init_params[i,] * c(.5,.5,n,n)
  addpars <- list(n=n, tobs = tobs[i,], x=x, ninit=1)
  cd.est <- tryCatch(optimx(pars, fn = loglikmf.model4, 
                            method = "BFGS", 
                            control = list(fnscale = -1), addpars = addpars), 
                     error = function(e) {return(c(NA, NA, NA, NA))})
  as.numeric(cd.est[1:4])
}
stopImplicitCluster()
mfergm_result = apply(mfergm_estim, 2, median, na.rm = T) * c(2,2,1/n,1/n)

# compare mean values
g4 <- simulate(g ~ edges+nodematch("x")+ kstar(2) + triangles,
               nsim = nsims,
               coef = mfergm_result,    
               output = "stats"
)

colMeans(g_sim_stats)
mfergm_mv = colMeans(g4)
mfergm_mv

# compare mfergm result with "ht_theta"
ht_theta
ergm_result
mfergm_result


# results for report
round(mv_params1_avg, 2)
round(mv_params2_avg, 2)
round(ht_theta, 3)
round(apply(init_params, 2, median, na.rm = T), 3)
round(ergm_result, 3)
round(ergm_mv, 2) 
round(mfergm_result, 3)
round(mfergm_mv, 2) 




# save.image("mfergm_ht_params1_n50.RData")







degen_mf_results = mfergm_estim[!complete.cases(mfergm_estim),]
degen_mcmc_results = ergm_sim_estim[!complete.cases(ergm_sim_estim),]


# Degeneracy
length(intersect(which(is.na(mfergm_estim[,1])), which(is.na(ergm_sim_estim[,1]))))

# MFERGM RMSE
outliers <- function(x, na.rm = TRUE) 
{
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  which(x < (qnt[1] - H) | x > (qnt[2] + H))
}

# mf
outliers1 = outliers((complete_mf_results[,1] - ht_theta[1])^2)
outliers2 = outliers((complete_mf_results[,2] - ht_theta[2])^2)
outliers3 = outliers((complete_mf_results[,3] - ht_theta[3])^2)
outliers4 = outliers((complete_mf_results[,4] - ht_theta[4])^2)
outliers = intersect(intersect(intersect(outliers1, outliers2), outliers3), outliers4)
mf_rmse1 = sqrt(sum(((complete_mf_results[,1] - ht_theta[1])^2)[-outliers1]))
mf_rmse2 = sqrt(sum(((complete_mf_results[,2] - ht_theta[2])^2)[-outliers2]))
mf_rmse3 = sqrt(sum(((complete_mf_results[,3] - ht_theta[3])^2)[-outliers3]))
mf_rmse4 = sqrt(sum(((complete_mf_results[,4] - ht_theta[4])^2)[-outliers4]))
c(mf_rmse1, mf_rmse2, mf_rmse3, mf_rmse4)

# mcmc
outliers1 = outliers((complete_mcmc_results[,1] - ht_theta[1])^2)
outliers2 = outliers((complete_mcmc_results[,2] - ht_theta[2])^2)
outliers3 = outliers((complete_mcmc_results[,3] - ht_theta[3])^2)
outliers4 = outliers((complete_mcmc_results[,4] - ht_theta[4])^2)
outliers = intersect(intersect(intersect(outliers1, outliers2), outliers3), outliers4)
mcmc_rmse1 = sqrt(sum(((complete_mcmc_results[,1] - ht_theta[1])^2)[-outliers1]))
mcmc_rmse2 = sqrt(sum(((complete_mcmc_results[,2] - ht_theta[2])^2)[-outliers2]))
mcmc_rmse3 = sqrt(sum(((complete_mcmc_results[,3] - ht_theta[3])^2)[-outliers3]))
mcmc_rmse4 = sqrt(sum(((complete_mcmc_results[,4] - ht_theta[4])^2)[-outliers4]))
c(mcmc_rmse1, mcmc_rmse2, mcmc_rmse3, mcmc_rmse4)

plot((complete_mfergm_results[,1] - ht_theta[1])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot((complete_mfergm_results[,2] - ht_theta[2])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot((complete_mfergm_results[,3] - ht_theta[3])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot((complete_mfergm_results[,4] - ht_theta[4])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mfergm_results[,1] - ht_theta[1])^2)[-outliers1], pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mfergm_results[,2] - ht_theta[2])^2)[-outliers2], pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mfergm_results[,3] - ht_theta[3])^2)[-outliers3], pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mfergm_results[,4] - ht_theta[4])^2)[-outliers4], pch = 19, ylab = "", xlab = "Simulated Networks")


outliers1 = head(order((complete_mfergm_results[,1] - ht_theta[1])^2, decreasing = T), n = 100)
outliers2 = head(order((complete_mfergm_results[,2] - ht_theta[2])^2, decreasing = T), n = 100)
outliers3 = head(order((complete_mfergm_results[,3] - ht_theta[3])^2, decreasing = T), n = 100)
outliers4 = head(order((complete_mfergm_results[,4] - ht_theta[4])^2, decreasing = T), n = 100)
outliers = intersect(intersect(intersect(outliers1, outliers2), outliers3), outliers4)


plot(density(complete_mfergm_results[,1]), main = "")
plot(density(complete_mfergm_results[,2]), main = "")
plot(density(complete_mfergm_results[,3]), main = "")
plot(density(complete_mfergm_results[,4]), main = "")



# MAD

mcmc_mad = matrix(nrow = nrow(complete_mcmc_results), ncol = 4)
for (i in 1:nrow(complete_mcmc_results))
{
  mcmc_mad[i,] = abs(complete_mcmc_results[i,] - ht_theta)
}
round(apply(mcmc_mad, 2, median), 3)

mf_mad = matrix(nrow = nrow(complete_mf_results), ncol = 4)
for (i in 1:nrow(complete_mf_results))
{
  mf_mad[i,] = abs(complete_mf_results[i,] - ht_theta)
}
round(apply(mf_mad, 2, median), 3)

