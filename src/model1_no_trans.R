rm(list=ls())
set.seed(1000)

# install.packages("rootSolve")
# install.packages("~/handcock/mfergm-master", repos = NULL, type="source")


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

nsims       =  1000                              # number of networks simulated
n           =  200                                # number of nodes
theta       =  c(-2,1,1,0) * c(2,2,1/n,1/n)      # true parameters for model


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
g_sim <- simulate(g ~ edges + nodematch("x") + kstar(2) + triangles, 
                  nsim = nsims,
                  coef = theta
)

# Extract summary stats from 'g_sim'
g_sim_stats = matrix(nrow = nsims, ncol = 4)
for (i in 1:nsims)
{
  g_sim_stats[i,] = summary_formula(g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles)
}
colMeans(g_sim_stats)
# both homophily and triangle set to zero
# n = 50
# 84.036  73.246 278.005   8.928
# 97.061  84.308 378.036  14.472
# n = 100
# 344.668  300.963 2381.074   77.241
# 400.312  349.880 3254.415  126.281
# n = 200
# 1361.655  1181.495 18405.026   585.892
# 1574.424  1362.559 24736.702   920.997
# only triangle set to zero
# n = 50
# 96.617  83.958 373.970  14.048
# n = 100
# 397.715  347.622 3210.465  123.067
# n = 200
# 1561.451  1350.792 24323.540   891.638


# asdf_index = g_sim_stats[,1] > mean(g_sim_stats[,1]) - sd(g_sim_stats[,1]) & g_sim_stats[,1] < mean(g_sim_stats[,1]) + sd(g_sim_stats[,1]) &
#              g_sim_stats[,2] > mean(g_sim_stats[,2]) - sd(g_sim_stats[,2]) & g_sim_stats[,2] < mean(g_sim_stats[,2]) + sd(g_sim_stats[,2]) &
#              g_sim_stats[,3] > mean(g_sim_stats[,3]) - sd(g_sim_stats[,3]) & g_sim_stats[,3] < mean(g_sim_stats[,3]) + sd(g_sim_stats[,3]) &
#              g_sim_stats[,4] > mean(g_sim_stats[,4]) - sd(g_sim_stats[,4]) & g_sim_stats[,4] < mean(g_sim_stats[,4]) + sd(g_sim_stats[,4])
# 
# asdf = g_sim_stats[asdf_index,]
# g_sim_stats = asdf
# nsims = nrow(asdf)
# g_sim = g_sim[which(asdf_index)]



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
registerDoParallel(detectCores())
ergm_sim_list = foreach(i = 1:nsims) %dopar% {
  # cat("***********************************\n")
  # cat(paste("estimating sample" ,i, "\n"))
  # cat("***********************************\n")
  skip_to_next <- FALSE
  formula <- g_sim[[i]] ~ edges + nodematch("x") + kstar(2) + triangles
  ergm_sim = 
    tryCatch(withTimeout(
                          ergm(formula,
                               control=control.ergm(init = init_params[i,])),
               timeout = 10000, onTimeout = "error"),
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
g3 <- simulate(g ~ edges+nodematch("x")+ kstar(2) + triangles,
               nsim = nsims,
               coef = ergm_result,         
               output = "stats"
)
colMeans(g_sim_stats)
ergm_mv = colMeans(g3)


#### DEGENERACIES
# degen = which(!complete.cases(ergm_sim_estim))
# # 58 of em
# 
# degen_list = ergm_sim_list[[]]
# registerDoParallel(14)
# ergm_degen_list = foreach(i = 1:length(degen)) %dopar% {
#   skip_to_next <- FALSE
#   formula <- g_sim[[degen[i]]] ~ edges + nodematch("x") + kstar(2) + triangles
#   ergm_sim = tryCatch(withTimeout(ergm(formula,
#                                         control=control.ergm(
#                                        #   # MCMC.burnin=100000,
#                                        #   # MCMC.interval=1000,
#                                        #   # MCMC.samplesize = 5000,
#                                           )
#   ), timeout = 1000, onTimeout = "error"), 
#   error = function(e) {skip_to_next <<- NULL})
#   ergm_sim
# }
# stopImplicitCluster()
# 
# ergm_degen_estim = matrix(0, nrow = length(degen), ncol = 4)
# for(i in 1:length(degen))
# {
#   if(!is.null(ergm_degen_list[[i]]))
#   {
#     est.params <- ergm_degen_list[[i]]$coefficients
#     ergm_degen_estim[i,] <- est.params
#   }
#   else {ergm_degen_estim[i,] = c(NA, NA, NA, NA)}
# }


# Model statistics ‘kstar2’ and ‘triangle’ are not varying. This may indicate 
# that the observed data occupies an extreme point in the sample space or that 
# the estimation has reached a dead-end configuration.
# Warning: Model statistics ‘nodematch.x’ are linear combinations of some set of 
# preceding statistics at the current stage of the estimation. This may indicate 
# that the model is nonidentifiable.


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
registerDoParallel(detectCores())
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
mfergm_estim = matrix(c(mfergm_estim[,1]*2, mfergm_estim[,2]*2, mfergm_estim[,3]/n, mfergm_estim[,4]/n),
               nrow = nsims)
mfergm_result = apply(mfergm_estim, 2, median, na.rm = T) 

# compare mean values
g4 <- simulate(g ~ edges+nodematch("x")+ kstar(2) + triangles,
               nsim = nsims,
               coef = mfergm_result,    
               output = "stats"
)
colMeans(g_sim_stats)
mfergm_mv = colMeans(g4)
mfergm_mv

# compare mfergm result with "theta"
theta
ergm_result
mfergm_result


# results for report
round(theta, 3) 
round(colMeans(g_sim_stats), 2)
round(colMeans(init_params), 3)
round(ergm_result, 3) 
round(ergm_mv, 2) 
round(mfergm_result, 3) 
round(mfergm_mv, 2) 
# / c(2,2,1/n,1/n)



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
  gtest_mf[i,] = simulate(g ~ edges+nodematch("x")+ kstar(2) + triangles, 
                          nsim = 1, 
                          coef = c[i,],        
                          output = "stats")
}
gtest_mc = matrix(nrow = nrow(d), ncol = 4)
for (i in 1:nrow(d))
{
  gtest_mc[i,] = simulate(g ~ edges+nodematch("x")+ kstar(2) + triangles, 
                          nsim = 1, 
                          coef = d[i,],        
                          output = "stats")
}
colMeans(gtest_mf)
colMeans(gtest_mc)

summary(ergm_sim_estim)
summary(gtest)
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
me(gtest_mf[,1])
hist(gtest_mf[,2])
hist(gtest_mf[,3])
hist(gtest_mf[,4])










complete_mf_results = mfergm_estim[complete.cases(mfergm_estim),]
complete_mcmc_results = ergm_sim_estim[complete.cases(ergm_sim_estim),]


# Degeneracy
length(intersect(which(is.na(mfergm_estim[,1])), which(is.na(ergm_sim_estim[,1]))))

# MFERGM RMSE

# mf
outliers1 = outliers((complete_mf_results[,1] - theta[1])^2)
outliers2 = outliers((complete_mf_results[,2] - theta[2])^2)
outliers3 = outliers((complete_mf_results[,3] - theta[3])^2)
outliers4 = outliers((complete_mf_results[,4] - theta[4])^2)
outliers = intersect(intersect(intersect(outliers1, outliers2), outliers3), outliers4)
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
outliers = intersect(intersect(intersect(outliers1, outliers2), outliers3), outliers4)
mcmc_rmse1 = sqrt(sum(((complete_mcmc_results[,1] - theta[1])^2)[-outliers1]))
mcmc_rmse2 = sqrt(sum(((complete_mcmc_results[,2] - theta[2])^2)[-outliers2]))
mcmc_rmse3 = sqrt(sum(((complete_mcmc_results[,3] - theta[3])^2)[-outliers3]))
mcmc_rmse4 = sqrt(sum(((complete_mcmc_results[,4] - theta[4])^2)[-outliers4]))
c(mcmc_rmse1, mcmc_rmse2, mcmc_rmse3, mcmc_rmse4)

plot((complete_mfergm_results[,1] - theta[1])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot((complete_mfergm_results[,2] - theta[2])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot((complete_mfergm_results[,3] - theta[3])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot((complete_mfergm_results[,4] - theta[4])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mfergm_results[,1] - theta[1])^2)[-outliers1], pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mfergm_results[,2] - theta[2])^2)[-outliers2], pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mfergm_results[,3] - theta[3])^2)[-outliers3], pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mfergm_results[,4] - theta[4])^2)[-outliers4], pch = 19, ylab = "", xlab = "Simulated Networks")


outliers1 = head(order((complete_mfergm_results[,1] - theta[1])^2, decreasing = T), n = 100)
outliers2 = head(order((complete_mfergm_results[,2] - theta[2])^2, decreasing = T), n = 100)
outliers3 = head(order((complete_mfergm_results[,3] - theta[3])^2, decreasing = T), n = 100)
outliers4 = head(order((complete_mfergm_results[,4] - theta[4])^2, decreasing = T), n = 100)
outliers = intersect(intersect(intersect(outliers1, outliers2), outliers3), outliers4)


plot(density(complete_mfergm_results[,1]), main = "")
plot(density(complete_mfergm_results[,2]), main = "")
plot(density(complete_mfergm_results[,3]), main = "")
plot(density(complete_mfergm_results[,4]), main = "")



# MAD

mcmc_mad = matrix(nrow = nrow(complete_mcmc_results), ncol = 4)
for (i in 1:nrow(complete_mcmc_results))
{
  mcmc_mad[i,] = abs(complete_mcmc_results[i,] - theta)
}
round(apply(mcmc_mad, 2, median), 3)

mf_mad = matrix(nrow = nrow(complete_mf_results), ncol = 4)
for (i in 1:nrow(complete_mf_results))
{
  mf_mad[i,] = abs(complete_mf_results[i,] - theta)
}
round(apply(mf_mad, 2, median), 3)

