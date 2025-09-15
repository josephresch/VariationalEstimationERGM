rm(list=ls())
my.seed=1
set.seed(my.seed)

pdf("model9_sim1.tapered.pdf")

library(ergm.tapered)
library(doRNG)
library(mfergm)
library(optimx)     # mfergm likelihood
library(R.utils)    # set time on ergm
library(doParallel) # parallel loops using 'foreach'


if(T){
load("mfergm_params10_n1000_1.RData")
v1 <- mfergm_estim
v2 <- ergm_sim_estim_tapered
v3 <- init_params
v4 <- degen_mfergm
v5 <- degen_tapered
v6 <- mf_mv_est
v7 <- ergm_mv_est_tapered
v8 <- MPLE_mv_est
v9 <- mv_s/5

load("mfergm_params10_n1000_2.RData")
v1 <- rbind(v1,mfergm_estim)
v2 <- rbind(v2,ergm_sim_estim_tapered)
v3 <- rbind(v3,init_params)
v4 <- c(v4,degen_mfergm+200)
v5 <- c(v5,degen_tapered+200)
v6 <- rbind(v6,mf_mv_est)
v7 <- rbind(v7,ergm_mv_est_tapered)
v8 <- rbind(v8,MPLE_mv_est)
v9 <- v9 + mv_s/5

load("mfergm_params10_n1000_3.RData")
v1 <- rbind(v1,mfergm_estim)
v2 <- rbind(v2,ergm_sim_estim_tapered)
v3 <- rbind(v3,init_params)
v4 <- c(v4,degen_mfergm+2*200)
v5 <- c(v5,degen_tapered+2*200)
v6 <- rbind(v6,mf_mv_est)
v7 <- rbind(v7,ergm_mv_est_tapered)
v8 <- rbind(v8,MPLE_mv_est)
v9 <- v9 + mv_s/5

load("mfergm_params10_n1000_4.RData")
v1 <- rbind(v1,mfergm_estim)
v2 <- rbind(v2,ergm_sim_estim_tapered)
v3 <- rbind(v3,init_params)
v4 <- c(v4,degen_mfergm+3*200)
v5 <- c(v5,degen_tapered+3*200)
v6 <- rbind(v6,mf_mv_est)
v7 <- rbind(v7,ergm_mv_est_tapered)
v8 <- rbind(v8,MPLE_mv_est)
v9 <- v9 + mv_s/5

load("mfergm_params10_n1000_5.RData")
v1 <- rbind(v1,mfergm_estim)
v2 <- rbind(v2,ergm_sim_estim_tapered)
v3 <- rbind(v3,init_params)
v4 <- c(v4,degen_mfergm+4*200)
v5 <- c(v5,degen_tapered+4*200)
v6 <- rbind(v6,mf_mv_est)
v7 <- rbind(v7,ergm_mv_est_tapered)
v8 <- rbind(v8,MPLE_mv_est)
v9 <- v9 + mv_s/5

mfergm_estim <- v1
ergm_sim_estim_tapered <- v2
init_params <- v3
degen_mfergm <- v4
degen_tapered <- v5
mf_mv_est <- v6
ergm_mv_est_tapered <- v7
MPLE_mv_est <- v8
}else{
load("mfergm_params10_n1000.RData")
}

# First natural params
complete_mf_results = mfergm_estim[complete.cases(mfergm_estim),]
complete_tapered_results = ergm_sim_estim_tapered[complete.cases(ergm_sim_estim_tapered),]
complete_MPLE_results = init_params[complete.cases(init_params),]

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
mf_rmse <- bdrmse(complete_mf_results, theta, degen_mfergm)
tapered_rmse <- bdrmse(complete_tapered_results, theta, degen_tapered)

# MAD
MPLE_mad <- bdmad(complete_MPLE_results, theta,degen=NULL)
mf_mad <- bdmad(complete_mf_results, theta, degen_mfergm)
tapered_mad <- bdmad(complete_tapered_results, theta, degen_tapered)

RMSE_natural_parameter <- cbind(round(tapered_rmse, 3), round(MPLE_rmse, 3), round(mf_rmse, 3))
rownames(RMSE_natural_parameter) <- names(mv_s)
colnames(RMSE_natural_parameter) <- c("MCMC-MLE","MPLE", "MFVLE")
# These are the RMSE of natural parameter estimates
RMSE_natural_parameter

MAD_natural_parameter <- cbind(round(tapered_mad, 3), round(MPLE_mad, 3), round(mf_mad, 3))
rownames(MAD_natural_parameter) <- names(mv_s)
colnames(MAD_natural_parameter) <- c("MCMC-MLE","MPLE", "MFVLE")
# These are the MAD of natural parameter estimates
MAD_natural_parameter


outliers1 = outliers((complete_tapered_results[,1] - theta[1])^2)
outliers2 = outliers((complete_tapered_results[,2] - theta[2])^2)
outliers3 = outliers((complete_tapered_results[,3] - theta[3])^2)
outliers4 = outliers((complete_tapered_results[,4] - theta[4])^2)
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


# Now the mean values
complete_mf_results = mf_mv_est[complete.cases(mf_mv_est),]
complete_tapered_results = ergm_mv_est_tapered[complete.cases(ergm_mv_est_tapered),]
complete_MPLE_results = MPLE_mv_est[complete.cases(MPLE_mv_est),]


# Degeneracy
length(intersect(which(is.na(mf_mv_est[,1])), which(is.na(ergm_mv_est_tapered[,1]))))

# MFERGM RMSE

outliers1 = outliers((complete_tapered_results[,1] - mv_s[1])^2)
outliers2 = outliers((complete_tapered_results[,2] - mv_s[2])^2)
outliers3 = outliers((complete_tapered_results[,3] - mv_s[3])^2)
outliers4 = outliers((complete_tapered_results[,4] - mv_s[4])^2)
plot((complete_mf_results[,1] - mv_s[1])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot((complete_mf_results[,2] - mv_s[2])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot((complete_mf_results[,3] - mv_s[3])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot((complete_mf_results[,4] - mv_s[4])^2, pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mf_results[,1] - mv_s[1])^2)[-outliers1], pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mf_results[,2] - mv_s[2])^2)[-outliers2], pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mf_results[,3] - mv_s[3])^2)[-outliers3], pch = 19, ylab = "", xlab = "Simulated Networks")
plot(((complete_mf_results[,4] - mv_s[4])^2)[-outliers4], pch = 19, ylab = "", xlab = "Simulated Networks")

plot(density(complete_mf_results[,1]), main = "")
plot(density(complete_mf_results[,2]), main = "")
plot(density(complete_mf_results[,3]), main = "")
plot(density(complete_mf_results[,4]), main = "")

# RMSE

MPLE_rmse <- bdrmse(complete_MPLE_results, mv_s, degen=NULL)
mf_rmse <- bdrmse(complete_mf_results, mv_s, degen_mfergm)
tapered_rmse <- bdrmse(complete_tapered_results, mv_s, degen_tapered)

# MAD
MPLE_mad <- bdmad(complete_MPLE_results, mv_s, degen=NULL)
mf_mad <- bdmad(complete_mf_results, mv_s, degen_mfergm)
tapered_mad <- bdmad(complete_tapered_results, mv_s, degen_tapered)

RMSE_mean_value_parameter <- cbind(round(tapered_rmse, 3), round(MPLE_rmse, 3), round(mf_rmse, 3))
rownames(RMSE_mean_value_parameter) <- names(mv_s)
colnames(RMSE_mean_value_parameter) <- c("MCMC-MLE","MPLE", "MFVLE")
# These are the RMSE of mean-value parameter estimates
RMSE_mean_value_parameter

MAD_mean_value_parameter <- cbind(round(tapered_mad, 3), round(MPLE_mad, 3), round(mf_mad, 3))
rownames(MAD_mean_value_parameter) <- names(mv_s)
colnames(MAD_mean_value_parameter) <- c("MCMC-MLE","MPLE", "MFVLE")
# These are the MAD of mean-value parameter estimates
MAD_mean_value_parameter

plot(mfergm_estim[,4],ergm_sim_estim_tapered[,4])

mean_est <- cbind(theta, ergm_result_tapered, MPLE_result, mfergm_result)
rownames(mean_est) <- names(mv_s)
colnames(mean_est) <- c("true", "MCMC-MLE","MPLE","MFVLE")
mean_est

mean_mv_est <- cbind(round(colMeans(g_sim_stats), 2), round(ergm_mv_tapered,2), round(MPLE_mv,2), round(mfergm_mv, 2))
rownames(mean_mv_est) <- names(mv_s)
colnames(mean_mv_est) <- c("true", "MCMC-MLE","MPLE", "MFVLE")
mean_mv_est

# Summary
length(degen_tapered) / nsims
length(degen_mfergm) / nsims
mean_est
mean_mv_est
# These are the RMSE of natural parameter estimates
RMSE_natural_parameter
# These are the MAD of natural parameter estimates
MAD_natural_parameter
# These are the RMSE of mean-value parameter estimates
RMSE_mean_value_parameter
# These are the RMSE of mean-value parameter estimates
MAD_mean_value_parameter

rm(g.sim)
rm(ergm_sim_list_tapered)
save.image("mfergm_params10_n1000.RData")
