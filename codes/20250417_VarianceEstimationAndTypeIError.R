################################################################################
#                                                                              #  
#           Simulations post-clustering inference with Data fission            #      
#        Link between biase variance and Type I error rate inflation           #
#                       Benjamin Hivert - 17/04/2025                           #  
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
#
# Theses simulations have for objectives to verifiy the analytical expression of 
# the type I error rate as a function of the biase in estimation of variance
# for Gaussian data fission.
#
#------------------------------------------------------------------------------#

# Library
library(dplyr)

#------------------------------------------------------------------------------#
#                   Functions used to simulate data                            #
#------------------------------------------------------------------------------#

sim_fun <- function(n,sigma, sigma_hat, tau){
  # Function that simulate an univariate realisation of a Gaussian distribution 
  # with mean 0 and variance sigma^2
  # Parameters :
  # n: a numerical value defining the sample size
  # sigma: a numerical value defining the true standard deviation of the Gaussian distribution
  # sigma_hat: a numerical value used in data fission a an estimate of the Gaussian standard deviation 
  # tau: the tunning parameter used in data fission 

    X <- rnorm(n, mean = 0, sd = sigma)
    Z <- rnorm(n, mean = 0, sd = sigma_hat)
    X1 <- X + tau*Z
    X2 <- X - (1/tau)*Z
    
    cl <- kmeans(X1, centers = 2)$cluster
    test <- t.test(X2~cl)
    return(data.frame(pvalues = test$p.value,
                      sigma = sigma, 
                      sigma_hat = sigma_hat, 
                      n = n))
}


sim_fun_sigma_grid <- function(n,sigma, sigma_hat, tau){
  # Function that apply the univariate simulation over a grid of sigma_hat
  # Parameters :
  # n: a numerical value defining the sample size
  # sigma: a numerical value defining the true standard deviation of the Gaussian distribution
  # sigma_hat: a vector of values used in data fission a estimates of the Gaussian standard deviation 
  # tau: the tunning parameter used in data fission 
  res_temp <- pblapply(sigma_hat, function(s){sim_fun(n, sigma, sigma_hat = s, tau)})
  res <- do.call("rbind.data.frame", res_temp)
  return(res)
}

sim_fun_sigma <- function(n,sigma, sigma_hat, tau){
  # Function that apply the univariate simulation over a grid of sigma_hat and a grid of sigma
  # Parameters :
  # n: a numerical value defining the sample size
  # sigma: a vector of values defining the true standard deviation of the Gaussian distribution
  # sigma_hat: a vector of values used in data fission a estimates of the Gaussian standard deviation 
  # tau: the tunnining parameter used in data fission 
  res_temp <- pblapply(sigma, function(s){sim_fun_sigma_grid(n, sigma = s, sigma_hat, tau)})
  res <- do.call("rbind.data.frame", res_temp)
  return(res)
}

sim_fun_n <- function(n_grid,sigma, sigma_hat, tau){
  # Function that apply the univariate simulation over a grid of sigma_hat and a grid of sigma for a grid of sample size
  # Parameters :
  # n: a vector of values defining the sample size
  # sigma: a vector of values defining the true standard deviation of the Gaussian distribution
  # sigma_hat: a vector of values used in data fission a estimates of the Gaussian standard deviation 
  # tau: the tunnining parameter used in data fission 
  res_temp <- pblapply(n_grid, function(n){sim_fun_sigma_grid(n=n, sigma = sigma, sigma_hat, tau)})
  res <- do.call("rbind.data.frame", res_temp)
  return(res)
}

#------------------------------------------------------------------------------#
#                            Parameter definitions                             #
#------------------------------------------------------------------------------#

n <- 100
n_grid <- c(50, 100, 200, 500, 1000)

sigma <- c(0.1, 0.5, 1, 2)
sigma_grid <- sort(c(seq(0, 4, length.out = 50), 2, 0.1))
sigma_grid2 <- seq(0.8, 1.8, length = 50)
nsimu <- 1000
tau <- .4

#------------------------------------------------------------------------------#
#                                Run simulations                               #
#------------------------------------------------------------------------------#

pblapply(1:nsimu, function(ns){
  path <- paste0("raw_results/sim_VarianceEstimationAndTypeIError/results_", ns, ".csv")
  res_temp <- sim_fun_sigma(n = n,
                sigma = sigma, 
                sigma_hat = sigma_grid,
                tau = tau)
  write.csv(res_temp, file = path, row.names = FALSE)
}, cl = 7) 


pblapply(1:nsimu, function(ns){
  path <- paste0("raw_results/sim_VarianceEstimationAndTypeIErrorSampleSize/results_", ns, ".csv")
  res_temp <- sim_fun_n(n_grid = n_grid,
                            sigma = 1, 
                            sigma_hat = sigma_grid2,
                            tau = tau)
  write.csv(res_temp, file = path, row.names = FALSE)
}, cl = 7) 
