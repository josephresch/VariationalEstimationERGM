#### DEGENERACIES
length(which(!complete.cases(ergm_sim_estim)))
length(which(!complete.cases(mfergm_estim)))

degen = which(!complete.cases(ergm_sim_estim))
# 58 of em

degen_list = ergm_sim_list[[]]
registerDoParallel(14)
ergm_degen_list = foreach(i = 1:length(degen)) %dopar% {
  skip_to_next <- FALSE
  formula <- g_sim[[degen[i]]] ~ edges + nodematch("x") + kstar(2) + triangles
  ergm_sim = tryCatch(withTimeout(ergm(formula,
                                       control=control.ergm(
                                         #   # MCMC.burnin=100000,
                                         #   # MCMC.interval=1000,
                                         #   # MCMC.samplesize = 5000,
                                       )
  ), timeout = 1000, onTimeout = "error"),
  error = function(e) {skip_to_next <<- NULL})
  ergm_sim
}
stopImplicitCluster()

ergm_degen_estim = matrix(0, nrow = length(degen), ncol = 4)
for(i in 1:length(degen))
{
  if(!is.null(ergm_degen_list[[i]]))
  {
    est.params <- ergm_degen_list[[i]]$coefficients
    ergm_degen_estim[i,] <- est.params
  }
  else {ergm_degen_estim[i,] = c(NA, NA, NA, NA)}
}