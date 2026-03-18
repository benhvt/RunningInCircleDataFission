###############################################################################
# 
#   Simulations binomial negatives 
#
###############################################################################

library(dplyr)
library(mvtnorm)
library(countsplit)
library(pbapply)

# Simulation functions 

sim_fun <- function(n, 
                    rho, 
                    theta, 
                    error, 
                    tau){
  
  theta_hat <- theta*error
  # Simulate NB correlated data
  p <- length(theta)
  sigma2 <- log(1 + 1 / theta)
  
  # corrélation latente
  R <- matrix(rho, p, p)
  diag(R) <- 1
  
  # matrice de covariance latente
  Sigma <- diag(sqrt(sigma2)) %*% R %*% diag(sqrt(sigma2))
  
  Z <- rmvnorm(n, mean = rep(0, p), sigma = Sigma)
  Lambda <- sweep(exp(Z), 2, mu, "*")
  
  Y <- matrix(rpois(n * p, Lambda), n, p)
  colnames(Y) <- paste0("NB_", 1:p)
  
  Y_Gauss <- rmvnorm(n, mean = mu, sigma = Sigma)
  
  # Perform thinning 
  
  # Wrong theta
  split <- countsplit(Y, 
                      overdisps = theta_hat, 
                      epsilon = rep(tau, 2))
  
  km_thin <- kmeans(split[[1]], centers = 2, nstart = 100)$cluster
  pval_thin <- pbsapply(1:ncol(Y), function(s){
    wilcox.test(split[[2]][,s]~km_thin)$p.value
  })
  
  # Oracle
  split_oracle <- countsplit(Y, 
                             overdisps = theta, 
                             epsilon = rep(tau, 2))
  
  km_thin_oracle <- kmeans(split_oracle[[1]], centers = 2, nstart = 100)$cluster
  pval_thin_oracle <- pbsapply(1:ncol(Y), function(s){
    wilcox.test(split_oracle[[2]][,s]~km_thin_oracle)$p.value
  })
  
  # Perform data_fission
  # Wrong variance
  sigma2_hat <- log(1 + 1 / theta_hat)
  Sigma_hat <- diag(sqrt(sigma2_hat)) %*% R %*% diag(sqrt(sigma2_hat))
  
  Z_fiss <- rmvnorm(n, 
                    sigma = Sigma_hat)
  Y_Gauss_Fiss1 <- Y_Gauss + tau*Z_fiss
  Y_Gauss_Fiss2 <- Y_Gauss - (1/tau)*Z_fiss
  
  km_fiss <- kmeans(Y_Gauss_Fiss1, centers = 2, nstart = 100)$cluster
  pval_fiss <- pbsapply(1:ncol(Y_Gauss_Fiss2), function(s){
    wilcox.test(Y_Gauss_Fiss2[,s]~km_fiss)$p.value
  })
  
  Z_fiss_oracle <- rmvnorm(n, 
                           sigma = Sigma)
  Y_Gauss_Fiss1_oracle <- Y_Gauss + tau*Z_fiss_oracle
  Y_Gauss_Fiss2_oracle <- Y_Gauss - (1/tau)*Z_fiss_oracle
  
  km_fiss_oracle <- kmeans(Y_Gauss_Fiss1_oracle,
                           centers = 2, 
                           nstart = 100)$cluster
  pval_fiss_oracle <- pbsapply(1:ncol(Y_Gauss_Fiss2_oracle), function(s){
    wilcox.test(Y_Gauss_Fiss2_oracle[,s]~km_fiss_oracle)$p.value
  })
  
  res <- data.frame(pval = c(pval_thin, 
                             pval_thin_oracle,
                             pval_fiss,
                             pval_fiss_oracle),
                    Variable = rep(paste0("Var", 1:p), 4),
                    Method = rep(c("NB", "Gauss"), each = 2*p),
                    Estimation = rep(rep(c("Wrong", "Oracle"), each = p), 2),
                    Theta = theta,
                    Error = error,
                    Rho = rho, 
                    n = n)
  return(res)
}

sim_fun_cor <- function(n, 
                        rho_grid, 
                        theta, 
                        error, 
                        tau){
  res_cor <- lapply(1:length(rho_grid), function(r){
    return(sim_fun(n=n,
                   rho = rho_grid[r],
                   theta = theta,
                   error = error,
                   tau = tau))
  })
  return(do.call("rbind.data.frame", res_cor))
}

sim_fun_cor_theta <- function(n, 
                              rho_grid, 
                              theta, 
                              error_grid, 
                              tau){
  # browser()
  res_theta <- lapply(1:length(error_grid), function(t){
    return(sim_fun_cor(n=n,
                       rho = rho_grid,
                       theta = theta,
                       error = error_grid[t],
                       tau = tau))
  })
  return(do.call("rbind.data.frame", res_theta))
}



# Param def 
p <- 50
mu <- runif(p, 2, 20)
# theta <- runif(p, 1, 10)
theta <- rep(10, p)
# theta_hat <- c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100)
theta_hat <- c(0.001, 0.01, 0.05, 0.1, 0.5, 0.75, 0.9, 0.95, 
               1, 1.01, 1.1, 1.25, 1.5, 2, 2.5, 5, 10)


rho <- c(0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9)
nsim <- 500

res <- sim_fun_cor_theta(n = 100, 
                         rho_grid = rho, 
                         theta = theta,
                         error_grid = theta_hat,
                         tau = .5)


# Save ouput 
slar_taskid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
write.csv(res, file = paste0("results/CorrelationAndBiais_", slar_taskid, ".csv"))