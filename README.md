# Comparison of ERGM Likelihood Inference Methods (With Discussion on ERGM Variational Inference)

This repository includes R implementation of a comparison framework to compare likelihood inference methods of Exponential-family Random Graph Models (ERGMs). In particular, we compared the two leading methods (MCMC, Pseudolikelihood) to a newly proposed mean-field variational inference (Mele and Zhu, 2023) using a framework that:
- Compares likelihood estimation using mean value parameterization for ERGMs, thereby comparing the network statistics that are produced (Handcock, 2003)
- Stress tests inference methods by manually increasing network edge-dependency (transivity) by manipulating the mean value parameters

These findings lead to a discussion on the effectiveness of variational inference on ERGM, and how network transitivity plays a part. Our methods and findings are compiled in the following manuscript, which is currently in review.

**Resch, J., Handcock, M. S.** 
*Assessing Variational Likelihood Estimation in Network Formation Models*. *In Review*.  

## R Code

Implementations are provided fit ERGMs on microeconomics networks with varying sizes of 50, 100, 200, and 2,000. For each network size, we compare the standard model with the manually increased transitivity model using methods:
- Markov Chain Monte Carlo Maxmimum Likelihood Estimation (MCMC-MLE)
- Maximum Pseudolikelihood Estimation (MPLE)
- Mean-Field Varitional Likelihood Estimation (MFVLE)
