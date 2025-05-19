################################################################################
#                                                                              #  
#        Simulations post-clustering inference with Data Fission               #      
#       Control of type I error rate when spurious clusters are built          #
#                     Benjamin Hivert - 04/04/2025                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
#
# Theses simulations have for objectives to study Type I error rate in a setting 
# where more clusters than the true numbers of component in the mixture are built.
# To ensure a fair comparision between marginal and conditional data fission, 
# clustering will be applied only on one component of the mixture, leading to two
# spurious clusters. Statistical hypothesis testing will be applied on this spurious
# clusters.
# Moreover, data fissions will be applied with empirical estimate of corresponding
# covariances. 
#
#------------------------------------------------------------------------------#

# Library 

library(pbapply)
library(dplyr)
library(mvtnorm)

#------------------------------------------------------------------------------#
#                   Functions used to simulate data                            #
#------------------------------------------------------------------------------#

# Marginal variance of a Gaussian Mixture
mixture_variance <- function(pi, mu_list, sigma_list) {
  # Function that computes the marginal variance of a Gaussian mixtures models of
  # parameters :
  #pi: a vector of length G that contains mixing proportion 
  # mu_list: a list of length G that contains vector of length p that are the means of each components
  #sigma_list: a list of G that contains matrix of dimension p x p that are the component-specific covariance
  #
  #return a matrix of dimension p x p that contains the marginal covariance of the mixture 
  
  G <- length(pi)  #Number of component in the mixture
  d <- length(mu_list[[1]])  # Dimension
  mu_Z <- Reduce("+", Map("*", pi, mu_list))  # Marginal mean of the mixture
  
  # Initialisation of the covariance matrix
  Sigma <- matrix(0, nrow = d, ncol = d)
  
  # Marginal variance computation
  for (g in 1:G) {
    Sigma <- Sigma + pi[g] * (sigma_list[[g]] + tcrossprod(mu_list[[g]]))
  }
  Sigma <- Sigma - tcrossprod(mu_Z)
  
  return(Sigma)
}

#------------------------------------------------------------------------------#

# Global function to perform one simulation 
sim_fun <- function(n, 
                    pi,
                    list_mu, 
                    list_Sigma, 
                    tau){
  
  # Function that simulates a bivariate G-component mixture models with:
  #pi: a vector of length G of mixing proportion
  #list_mu: a list of length G containing component specific means (vectors of length p)
  #list_Sigma: a list of length G containing component specific variances (matrix of length p x p)
  #n: the sample size
  #tau: tuning parameter for data fission
  #
  #return a list of two data.frames:
  #pval: A data.frame containing results of post-clustering inference (pvalues and ARI)
  #Covariance: A data.frame containing results about the data fission (marginal and conditional covariance between X1 and X2)
  
  G <- length(pi)
  Z <- sample(1:G, n, replace = TRUE, prob = pi)
  X <- matrix(NA, nrow = n, ncol = 2)
  for (g in 1:G){
    X[Z==g,] <- rmvnorm(n = sum(Z==g), 
                        mean = list_mu[[g]], 
                        sigma = list_Sigma[[g]])
  }
  
  ##############################################################################
  #                        Marginal Data Fission                              #
  
  margCov_emp <- cov(X)
  
  W <- rmvnorm(n, sigma = margCov_emp)
  
  X1 <- X + tau*W
  X2 <- X - (1/tau) * W
  
  #Clustering 
  cl <- kmeans(X1[Z==1,], centers = 2, nstart = 100)$cluster

  #Test
  pval <- apply(X2[Z==1,], 2, function(x){t.test(x~cl)$p.value})
  
  
  #Results
  res_margFission_pval <- data.frame(pvalues = pval,
                                     Fission = "Marginal",
                                     Variable = c("X1", 'X2'), 
                                     n = n,
                                     tau = tau)
 
  
  ##############################################################################
  #                            Conditional Fission                             #
  
  W_cond <- matrix(NA, nrow = n, ncol = 2)
  for(g in 1:G){
    W_cond[Z==g] <- rmvnorm(sum(Z == g), sigma = cov(X[Z==g, ]))
  }
  
  X1_cond <- X + tau * W_cond
  X2_cond <- X - (1 / tau) * W_cond
  
  #Clustering
  cl_cond <- kmeans(X1_cond[Z==1,], centers = 2, nstart = 100)$cluster

  #Test
  pval_cond <- apply(X2_cond[Z==1,], 2, function(x){t.test(x~cl_cond)$p.value})
  
  
  #Results
  res_condFission_pval <- data.frame(pvalues = pval_cond,
                                     Fission = "Conditional",
                                     Variable = c("X1", 'X2'),
                                     n =n,
                                     tau = tau)
  

  
  return(list(pval = rbind.data.frame(res_margFission_pval, res_condFission_pval)))
}

#------------------------------------------------------------------------------#
# Simulation function for varying sampple size

sim_n <- function(n_vec, 
                  pi,
                  list_mu, 
                  list_Sigma,
                  tau){
  #n_vec: a vector of length N of different sample size
  
  temp <- lapply(n_vec, function(n){
    res_temp <- sim_fun(n = n,
                        pi = pi,
                        list_mu = list_mu, 
                        list_Sigma = list_Sigma,
                        tau = tau)
    return(res_temp)})
  sim_cov_list <- lapply(temp, function(l){l$Covariance})
  sim_pvalues_list <- lapply(temp, function(l){l$pval})
  return(list(Covariance = do.call("rbind.data.frame", sim_cov_list),
              pval = do.call("rbind.data.frame", sim_pvalues_list)))
}

#------------------------------------------------------------------------------#
# Simulation function for varying sampple size

sim_n_tau <- function(n_vec,
                      pi,
                      list_mu,
                      list_Sigma, 
                      tau_vec){
  #n_vec: a vector of length N of different sample size
  #tau_vec: a vector of length T of different values for the tunning parameter tau
  
  temp <- lapply(tau_vec, function(t){
    res_temp <- sim_n(n_vec = n_vec,
                      pi = pi,
                      list_mu = list_mu,
                      list_Sigma = list_Sigma,
                      tau = t)
    return(res_temp)})
  sim_cov_list <- lapply(temp, function(l){l$Covariance})
  sim_pvalues_list <- lapply(temp, function(l){l$pval})
  return(list(Covariance = do.call("rbind.data.frame", sim_cov_list),
              pval = do.call("rbind.data.frame", sim_pvalues_list)))
}


#------------------------------------------------------------------------------#
#                          Parameters definition                               #
#------------------------------------------------------------------------------#

# Parameters of interest
sample_size <- c(50, 100, 250, 500, 1000, 5000, 10000)
tau_grid <- c(seq(0.1, 1, by = .05), seq(1.5, 5, by = 0.5))

# Parameters of the mixtures (2 components)
pi <- c(0.5, 0.5)
mu <- list(c(0,5),
           c(0,0))
Sigma <- list(diag(1,2),
              cbind(c(1, 0.5),
                    c(0.5, 1)))

# Run simulation
sim <- sim_n_tau(n_vec = sample_size,
                 pi = pi,
                 list_mu = mu,
                 list_Sigma = Sigma,
                 tau_vec = tau_grid)

# Save ouput 
slar_taskid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
write.csv(sim$pval, file = paste0("results/MarginalConditionalFission_typeI_", slar_taskid, ".csv"))

