################################################################################
#                                                                              #
#                       Prepared results to make figure                        #
#                         Benjamin Hivert - 17/04/2025                         #
#                                                                              #
################################################################################

# This R script contains code to load all results from simulations studies and 
# application to prepare the mains results to make figures

#------------------------------------------------------------------------------#
#                                 Figures                                      #
#               Behaviour of marginal and conditional data fission             # 
#------------------------------------------------------------------------------#

# Results under the ideal scenario

pvalues_files <- list.files(path = "raw_results/sim_MarginalConditionalFission/",
                            pattern = "_pvalues",
                            full.names = TRUE)

pvaluesResults_list <- lapply(pvalues_files, read.csv)

pvaluesResults <- do.call("rbind.data.frame", pvaluesResults_list)

write.csv(pvaluesResults, file = "results/IdealScenario_GaussianFission.csv",
          row.names = FALSE)

# Results under the adverse scenario 

typeI_files <- list.files(path = "raw_results/sim_TypeI_MarginalConditionalFission/",
                          pattern = "_typeI",
                          full.names = TRUE)

typeIResults_list <- lapply(typeI_files, read.csv)

typesIResults <- do.call("rbind.data.frame", typeIResults_list)
write.csv(typesIResults, file = "results/AdverseScenario_GaussianFission.csv",
          row.names = FALSE)


#------------------------------------------------------------------------------#
#                                 Figures                                      #
#       Bias in the estimation of the covariance & Type I error rate           # 
#------------------------------------------------------------------------------#
 
# Function of the original variance 
bias_results_files <- list.files(path = "raw_results/sim_VarianceEstimationAndTypeIError/",
                                 full.names = TRUE)
bias_list <- lapply(bias_results_files, read.csv)

bias_results <- do.call("rbind.data.frame", bias_list)
write.csv(bias_results, file = "results/VarianceEstimationAndTypeIError.csv",
          row.names = FALSE)


# Function of the sample size 

bias_results_n_files <- list.files(path = "raw_results/sim_VarianceEstimationAndTypeIErrorSampleSize/",
                                 full.names = TRUE)
bias_n_list <- lapply(bias_results_n_files, read.csv)

bias_results_n <- do.call("rbind.data.frame", bias_n_list)
write.csv(bias_results_n, file = "results/VarianceEstimationAndTypeIErrorSampleSize.csv",
          row.names = FALSE)

#------------------------------------------------------------------------------#
#                                 Figures                                      #
#               Data thinning for negative binomial mixture                    # 
#------------------------------------------------------------------------------#

negbin_files <- list.files(path = "raw_results/sim_TypeIThinningNegBin/",
                           full.names = TRUE)
negbin_list <- lapply(negbin_files, read.csv)

negbin_results <- do.call("rbind.data.frame", negbin_list)
write.csv(negbin_results, file = "results/TypeIThinningNegBin.csv",
          row.names = FALSE)
