# We've already seen that bias correction works, when correcting with the true radial bias. Of course, this won't really by known.
# However, second-order MRV might help. This script attempts to simulate true biases, as well as biases assessed assuming
# second order MRV in an oracle setting, where we know the index of regular variation, the SOMRV parameter, and the bias of intermediate
# projection.
#
# Jimmy Butler
# March 2026

set.seed(45678)

library(dplyr)
library(mvtnorm)
library(mnormt)
library(argparse)
library(matrixStats)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')

# parsing command line specified parameters
parser <- ArgumentParser(description = "Args for Oracle SOMRV Bias Correction.")
parser$add_argument("--save_df_path", type = "character", default = NULL, 
                    help = "File path to save the csv comparing true vs bootstrap biases.")

args <- parser$parse_args()

# simulation parameters
ns <- c(1000, 5000, 10000, 50000)
C <- 5
xi <- 1/4
gamma1 <- 1/2
gamma2 <- 2/3
rho <- -1/2

# number of theta to use
n_coords <- 50
n_monte_carlo <- 5000

pn <- function(n){ C/n } 
dists <- c('bivt')
lbs <- c(0,0)
ubs <- c(200, 200)

# initialize a list to accumulate results across distributions and n values
all_results_list <- list()

# loop over the distributions
for (i in 1:length(dists)) {
    dist <- dists[i]
    sampling_func <- loadSamplingFunction(dist)
    
    print(paste0('Starting new distribution: ', dist))
    
    for (k in 1:length(ns)) {
        n <- ns[k]
        p <- pn(n)
        q1 <- n^(-gamma1)
        q2 <- n^(-gamma2)

        # compute the true radial biases via monte carlo simulation, to compare with the bootstrap estimates
        p_isoline <- drawIsoline(dist=dist, numCoords=n_coords, gridUbs=ubs, gridLbs=lbs, prob=p)
        p_thetas_true <- atan2(p_isoline[,2], p_isoline[,1])
        rp_thetas_true <- sqrt(p_isoline[,1]**2 + p_isoline[,2]**2)
        rp_hat_theta_mat <- matrix(NA, nrow=n_monte_carlo, ncol=length(p_thetas_true))

        print(paste0('Starting n: ', n, ' | Computing True MC Radial Bias...'))
        for (m in 1:n_monte_carlo) {
            dat_mc <- sampling_func(n)
            isoline_ext_estimate <- drawExtremeIsoline(dat=dat_mc, p=p, n_coords=n_coords, grid_lbs=lbs, gamma=gamma1, xi=xi)
            rp_hat_theta_mat[m,] <- sqrt(isoline_ext_estimate[,1]**2 + isoline_ext_estimate[,2]**2)
        }
        
        true_biases <- colMeans(rp_hat_theta_mat) - rp_thetas_true
        
        # the second order MRV bias ratio
        M <- ((q2/p)^(xi))*(((q1/p)^(rho) - 1)/((q1/q2)^(rho) - 1))

        q2_isoline <- drawIsoline(dist=dist, numCoords=n_coords, gridUbs=ubs, gridLbs=lbs, prob=q2)
        q2_thetas_true <- atan2(q2_isoline[,2], q2_isoline[,1])
        rq2_thetas_true <- sqrt(q2_isoline[,1]**2 + q2_isoline[,2]**2)
        rq2_hat_theta_mat <- matrix(NA, nrow=n_monte_carlo, ncol=length(q2_thetas_true))
        
        # compute monte carlo bias of intermediate projection
        for (m in 1:n_monte_carlo) {
            dat_mc <- sampling_func(n)
            isoline_q2_estimate <- drawExtremeIsoline(dat=dat_mc, p=q2, n_coords=n_coords, grid_lbs=lbs, gamma=gamma1, xi=xi)
            rq2_hat_theta_mat[m,] <- sqrt(isoline_q2_estimate[,1]**2 + isoline_q2_estimate[,2]**2)
        }

        somrv_oracle_biases <- M * (colMeans(rq2_hat_theta_mat) - rq2_thetas_true)
        
        # store the current iteration's results in a long dataframe
        current_results <- data.frame(
            dist = dist,
            n = n,
            theta = p_thetas_true,
            true_bias = true_biases,
            somrv_bias_correction = somrv_oracle_biases
        )
        
        # append to the list
        all_results_list[[length(all_results_list) + 1]] <- current_results
    }
}

# bind all iterations into one master dataframe
final_df <- bind_rows(all_results_list)

# save to CSV if a path was provided
if (!is.null(args$save_df_path)) {
    write.csv(final_df, file = paste0(args$save_df_path, 'somrv_oracle_bias.csv'), row.names = FALSE)
    print(paste0("Results successfully saved to: ", args$save_df_path))
} else {
    print("No save path provided. Results calculated but not saved.")
}