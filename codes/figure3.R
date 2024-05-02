# --------------------------------- Figure 3 --------------------------------- #
library(dplyr)

#-- Functions 
source("utils.R")

# kernel <- function(u, h) {
#   a <- h *sqrt(3)
#   return(ifelse(abs(u) < a, 0.5/a, 0))
# }

# ID of task

slar_taskid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#Gaussian
kernel <- function(u, h) {
  return(dnorm(u, sd=h))
}

weightvar <- function(x,obs, h = NULL){
  if (is.null(h)) {
    h <- sd(x) * (4/(3*n))^(1/5)
  }
  weigth <- kernel(x-obs, h)
  return(modi::weighted.var(obs, weigth))
}

local_var <- function(x){
  n <- length(x)
  h0 <- ks::hpi(x, deriv.order = 1)
  pilot <- ks::kde(x, h = h0, eval.points=x)$estimate
  lambda <- exp(mean(log(pilot)))
  hi <- h0*sqrt(lambda/pilot)
  variance_wk <- sapply(1:n, FUN = function(i){weightvar(x[i], obs = x, h = hi[i])})
  return(variance_wk)
}

# compute_hi <- function(x, obs){
#   # browser()
#   e2 <- (x-obs)^2
#   spread_knns <- (cumsum(sort(e2))/(1:length(e2)-1))[-length(e2)][-1]
#   y <- spread_knns
#   n <- 1:length(spread_knns)
#   bkp_ind <- floor(segmented::segmented(lm(y~n))$psi[,"Est."])
#   spread_sqrt <- sqrt(spread_knns[bkp_ind])
#   
#   hi <- abs(x-obs[order(e2)[bkp_ind+1]])
#   return(hi)
# }

local_var_CP <- function(x){
  variance_wk <- sapply(1:length(x), FUN = function(i){
    e2 <- (x[i]-x)^2 #Calcul des distances entre une observations et le reste
    spread_knns <- (cumsum(sort(e2))/(1:length(e2)-1))[-length(e2)][-1]
    y <- spread_knns
    xx <- 1:length(spread_knns)
    bkp_ind <- floor(segmented::segmented(lm(y~xx))$psi[,"Est."])
    spread_sqrt <- sqrt(spread_knns[bkp_ind])
    hi <- abs(x[i]-x[order(e2)[bkp_ind+1]])
    return(weightvar(x[i], obs = x, h = hi))})
  return(variance_wk)
}

sim_fun <- function(seed, delta, n, sd, tau){
  # First Version of the local variance estimation
  # browser()
  set.seed(250124*seed)
  X <- c(rnorm(n/2, mean = 0, sd = sd),
         rnorm(n/2, mean = delta, sd = sd))
  variance <- local_var(X)
  Z <- sapply(variance, function(s){rnorm(1, mean = 0, sd = sqrt(s))})
  fX  <- X + tau*Z 
  gX <- X - (1/tau)*Z 
  
  cl <- kmeans(fX, centers = 3, nstart = 100)$cluster
  df <- data.frame(gX = gX, 
                   Cluster = as.factor(cl))
  
  clToTest <- order_cluster(gX, cl)
  res.ttest <- try(t.test(gX[cl ==clToTest[1]], gX[cl==clToTest[2]]))
  if(class(res.ttest) == "try-error"){
    pval <- NA}
  else{
    pval <- res.ttest$p.value}
  
  res.ttestPower <- try(t.test(gX[cl ==clToTest[1]], gX[cl==which(levels(as.factor(cl)) %in% clToTest == F)]))
  
  if(class(res.ttestPower) == "try-error"){
    pvalPower <- NA}
  else{
    pvalPower <- res.ttestPower$p.value}
  
  #Second version of the variance estimaiton
  variance_CP <- local_var_CP(X)
  Z_CP <- sapply(variance_CP, function(s){rnorm(1, mean = 0, sd = sqrt(s))})
  fX_CP  <- X + tau*Z_CP
  gX_CP <- X - (1/tau)*Z_CP

  cl_CP <- kmeans(fX_CP, centers = 3, nstart = 100)$cluster
  df_CP <- data.frame(gX = gX_CP,
                      Cluster = as.factor(cl_CP))

  clToTest_CP <- order_cluster(gX_CP, cl_CP)
  res.ttest_CP <- t.test(gX_CP[cl_CP ==clToTest_CP[1]], gX_CP[cl_CP==clToTest_CP[2]])
  if(class(res.ttest_CP) == "try-error"){
    pvalCP <- NA}
  else{
    pvalCP <- res.ttest_CP$p.value}

  res.ttest_CPPower <- try(t.test(gX_CP[cl_CP ==clToTest_CP[1]], gX_CP[cl_CP==which(levels(as.factor(cl_CP)) %in% clToTest_CP == F)]))
  
  if(class(res.ttest_CPPower) == "try-error"){
    pvalPower_CP <- NA}
  else{
    pvalPower_CP <- res.ttest_CPPower$p.value}
  
  return(data.frame(pval = c(pval, pvalCP),
                    pvalPower = c(pvalPower, pvalPower_CP),
                    Variance_hat = c(mean(variance), mean(variance_CP)),
                    Variance_hat_comp1 = c(mean(variance[1:(n/2)]),
                                           mean(variance_CP[1:(n/2)])),
                    Variance_hat_comp2 = c(mean(variance[((n/2)+1):n]),
                                           mean(variance_CP[((n/2)+1):n])),
                    sdVariance_hat = c(sd(variance), sd(variance_CP)),
                    Median_Variance = c(median(variance), median(variance_CP)),
                    Median_Variance_comp1 = c(median(variance[1:(n/2)]), median(variance_CP[1:(n/2)])),
                    Median_Variance_comp2 = c(median(variance[((n/2)+1):n]), median(variance_CP[((n/2)+1):n])),
                    Q1_Variance = c(quantile(variance, 0.25), quantile(variance_CP, 0.25)),
                    Q3_Variance = c(quantile(variance, 0.75), quantile(variance_CP, 0.75)),
                    Method = c("Individual hi", "Change Point"),
                    delta = delta,
                    sigma = sd))
}

apply_sim_delta <- function(seed, delta_grid, n, sd, tau){
  res_delta <- lapply(delta_grid, function(d){
    res_temp <- sim_fun(seed = seed, delta = d, n=n, sd = sd, tau = tau)
    return(res_temp)})
  temp <- do.call("rbind.data.frame", res_delta)
  # temp$delta <- delta_grid
  return(temp)
}

apply_sim_delta_sigma <- function(seed = seed, delta_grid, n, sd_grid, tau){
  res_sigma <- lapply(sd_grid, function(s){
    res_temp <- apply_sim_delta(seed = seed, delta_grid = delta_grid, n = n, sd = s, tau = tau)
    # res_temp$sigma <- s
    return(res_temp)
  })
  temp <- do.call("rbind.data.frame", res_sigma)
  return(temp)
}

#--Param 
n <- 100
delta_grid <- c(seq(0, 3, length.out = 25), seq(3.5, 100, length.out = 25))

tau <- .4
sigma <- c(0.1, 0.5, 1, 2) 

res <- apply_sim_delta_sigma(seed = slar_taskid, delta_grid = delta_grid, n = n, sd_grid = sigma, tau = tau)

#-- Save Results 

write.csv(res, paste0("results/figure3/results_", slar_taskid, ".csv"), row.names = FALSE)
