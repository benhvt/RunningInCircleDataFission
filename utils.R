library(dplyr)

km_fun <- function(x, K){
  km <- kmeans(x, centers = K, nstart = 1)
  return(as.factor(km$cluster))
}

intra_cluster_fission <- function(X, tau=1, cl_ref, sigma_c = NULL){
  #cl_ref: A reference clustering used to compute intra-cluster variance
  # sigma_c: a list containing intra-cluster cov if knwon
  if (!is.factor(cl_ref)){
    cl_ref <- as.factor(cl_ref)
  }
  if (is.null(sigma_c)){
    sigma_c <- lapply(levels(cl_ref), function(c){
      cov(X[cl_ref == c,])
    })
  }
  intra_cluster_res <- lapply(levels(cl_ref), function(c){
    
    fission_c <- data_fission(as.matrix(X)[cl_ref==c,], Sigma = sigma_c[[c]], tau = tau)
  })
  
  fX <- lapply(intra_cluster_res, function(l){
    return(l$fX)
  })
  
  gX <- lapply(intra_cluster_res, function(l){
    return(l$gX)
  })
  
  return(list(fX =  do.call('rbind',fX),
              gX =  do.call('rbind',gX)))
}

order_cluster <- function(x, cl){
  # x : the variable to test (where the cluster should be ordered) 
  # cl : A three clusters partitions of the data
  
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

compute_typeI <- function(n, tau, sigma, sigma_hat, alpha){
  qu_t <- qt(alpha/2, df = n-2, lower.tail = FALSE)
  qu_n <- qnorm(alpha/2, lower.tail = F)
  cor_fg <- (sigma^2-sigma_hat^2)/sqrt((sigma^2 + (tau^2)*sigma_hat^2)*(sigma^2 + (1/tau^2)*sigma_hat^2))  
  mean_t <- sqrt(n)*sqrt((2/pi)*cor_fg^2)/sqrt(1-(2/pi)*cor_fg^2)
  
  return(pnorm(qu_n, mean = mean_t, sd= 1, lower.tail = FALSE) + pnorm(-qu_n, mean = mean_t, sd= 1, lower.tail = TRUE))
}


# Non parametric variance estimation 
kernel <- function(u, h) {
  a <- h *sqrt(3)
  return(ifelse(abs(u) < a, 0.5/a, 0))
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

