################################################################################
#                                                                              #  
#        Simulations post-clustering inference with Data thinning              #      
#                     Negative Binomial Mixture setting                        #
#                       Benjamin Hivert - 17/04/2025                           #  
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
#
# Theses simulations have for objectives to study the behaviour of both 
# conditional and marginal data thinning when applied for mixture of negative
# binomial distributions. 
#
#------------------------------------------------------------------------------#

# Library 

library(datathin)
library(ggplot2)
library(pbapply)
library(dplyr)
library(latex2exp)
library(patchwork)
theme_set(theme_bw())


#------------------------------------------------------------------------------#
#                   Functions used to simulate data                            #
#------------------------------------------------------------------------------#


order_cluster <- function(x, cl){
  # x : the variable to test (where the cluster should be ordered) 
  # cl : A three clusters partitions of the data
  # return: a vector containing the two clusters originating from the same components
  
  df <- data.frame(x = x,
                   Cluster = as.factor(cl))
  
  ord_mean <- df %>% group_by(Cluster) %>% summarise(M = mean(x))
  clDiff <- abs(c(ord_mean$M[1]-ord_mean$M[2],
                  ord_mean$M[1]-ord_mean$M[3],
                  ord_mean$M[2] - ord_mean$M[3]))
  
  cl.ord <- which.min(clDiff)
  if(cl.ord == 1){
    cl1 <- 1
    cl2 <- 2 
  }else if(cl.ord == 2){
    cl1 <- 1
    cl2 <- 3
  }else{
    cl1 <- 2
    cl2 <- 3
  }
  return(c(cl1, cl2))
}

sim_fun <- function(n, 
                    probs, 
                    size, 
                    seed = 310123){
  
  # Function that simulates a two components negative binomial mixture with parameters :
  #n: sample size
  # probs: a list of two values giving each componenents probabilities
  # size: a list of two values giving the dispersion parameters
  #
  #return a data frame with Type I error rate of both marginal and conditional data thinninhg
  
  Z <-  rep(1:2, each = n/2)
  
  x1 <- rnbinom(n=n/2, prob = probs[1], size = size[1])
  x2 <- rnbinom(n=n/2, prob = probs[2], size = size[2])
  y <- rnbinom(n=n, prob = .5, size = 5)
  X <- cbind(c(x1,x2), 
             y)
  
################################################################################
#                           Marginal  thinning                                 #
  
  marginal_theta <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[,p], mu = mean(X[,p]))
  })
  
  res_marginal <- countsplit::countsplit(X, overdisps = marginal_theta)
  clust_marginal <-kmeans(as.matrix(res_marginal[[1]]), centers = 3)$cluster
  clToTest_marginal <- order_cluster(as.matrix(res_marginal[[1]])[,1], 
                                     clust_marginal)
  
  test_marginal <- wilcox.test(as.matrix(res_marginal[[2]])[clust_marginal == clToTest_marginal[1],1], 
                               as.matrix(res_marginal[[2]])[clust_marginal == clToTest_marginal[2],1])$p.value
  
################################################################################
#                           Conditionnal thinning                              #
  
  conditional_theta_Z1 <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[Z==1,p], mu = mean(X[Z==1,p]))
  })
  
  conditional_theta_Z2 <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[Z==2,p], mu = mean(X[Z==2,p]))
  })

  conditional_theta <- list(conditional_theta_Z1, 
                            conditional_theta_Z2)
  res_conditional_temp <- lapply(1:2, function(cl){
    res_datathin <- countsplit::countsplit(X[Z==cl, ], 
                                           overdisps = conditional_theta[[cl]])
  })
  
  Xtrain_cond <- rbind(as.matrix(res_conditional_temp[[1]][[1]]),
                       as.matrix(res_conditional_temp[[2]][[1]]))
  Xtest_cond <- rbind(as.matrix(res_conditional_temp[[1]][[2]]),
                      as.matrix(res_conditional_temp[[2]][[2]]))
  
  
  clust_cond <- kmeans(Xtrain_cond, centers=3)$cluster
  clToTest_cond <- order_cluster(Xtrain_cond[,1], clust_cond)
  test_cond <- wilcox.test(Xtest_cond[clust_cond==clToTest_cond[1],1],
                           Xtest_cond[clust_cond==clToTest_cond[2],1])$p.value 
  
  return(data.frame(pvalues = c(test_marginal, test_cond),
                    Fission = c("Marginal thinning",
                                "Conditional thinning")))
}

#------------------------------------------------------------------------------#
#                          Parameters definition                               #
#------------------------------------------------------------------------------#

n <- 100
probs <- c(0.5, 0.4)
size <- c(5, 40)
nsimu <- 1000

# ---------------------------------------------------------------------------- #
#                                  Run simulations                             #
#------------------------------------------------------------------------------#

res_nb <- pblapply(1:nsimu, function(ns){
  set.seed(310123*ns)
  path <- paste("raw_results/sim_TypeIThinningNegBin/results_", ns, ".csv")
  temp <- sim_fun(n =n,
          probs = probs,
          size = size)
  write.csv(temp, file = path, row.names = FALSE)
}, cl = 5)

typeI_nb <- do.call("rbind.data.frame", res_nb)

