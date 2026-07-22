# Script to compare bootstrap-estimated radial biases against true Monte Carlo expected biases
# for extreme isoline level sets. This is in response to the simulations of knownindex_oracleline_asymmetric_strategy1.R
# which show that, even with no sample splitting at all, bias correcting the radius estimates at each theta
# with the true bias leads to the correct coverage (if anything, we are overcovered, which is preferable..)
#
# Jimmy Butler
# March 2026

set.seed(45678)

library(dplyr)
library(mvtnorm)
library(mnormt)
library(foreach)
library(doSNOW)
library(parallel)
library(argparse)
library(matrixStats)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')

# parsing command line specified parameters
parser <- ArgumentParser(description = "Args for Bootstrap Radial Bias Exploration.")
parser$add_argument("--n_cores", type = "integer", default = 64, 
                    help = "Number of cores for parallel processing.")
parser$add_argument("--save_df_path", type = "character", default = NULL, 
                    help = "File path to save the csv comparing true vs bootstrap biases.")

args <- parser$parse_args()

# simulation parameters
ns <- c(1000, 5000, 10000, 50000)
C <- 5
xi <- 1/4
gamma <- 1/2

# number of theta to use
n_coords <- 50
n_monte_carlo <- 2000
B_boot <- 2000
n_sims <- 500

pn <- function(n){ C/n } 
dists <- c('bivt')
lbs <- c(0,0)
ubs <- c(200, 200)

n_cores <- args$n_cores

# Function to compute the bootstrap estimates of the radial bias, at each angle theta
estimateBootstrapRadialBias <- function(dat, p, B, n_coords, gamma, xi, lbs) {
    # estimate the "true" isoline in the bootstrap world
    base_isoline <- drawExtremeIsoline(dat=dat, p=p, n_coords=n_coords, grid_lbs=lbs, gamma=gamma, xi=xi)
    # true radii in the bootstrap world
    base_r <- sqrt(base_isoline[,1]^2 + base_isoline[,2]^2)
    thetas <- atan2(base_isoline[,2], base_isoline[,1])
    
    n_dat <- nrow(dat)
    boot_r_mat <- matrix(NA, nrow=B, ncol=n_coords)
    
    # bootstrap loop
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        boot_dat <- dat[idx, , drop=FALSE]
        
        # estimate extreme isoline using the bootstrap data, and compute bootstrap estimates of radii
        boot_isoline <- drawExtremeIsoline(dat=boot_dat, p=p, n_coords=n_coords, grid_lbs=lbs, gamma=gamma, xi=xi)
        boot_r_mat[k, ] <- sqrt(boot_isoline[,1]^2 + boot_isoline[,2]^2)
    }
    
    # compute average radii in bootstrap world, and compare to truth in bootstrap world
    boot_expected_r <- colMeans(boot_r_mat)
    boot_bias <- boot_expected_r - base_r
    
    return(list(thetas = thetas, 
                base_r = base_r, 
                boot_expected_r = boot_expected_r, 
                boot_bias = boot_bias))
}

# loop over the distributions
for (i in 1:length(dists)) {
    dist <- dists[i]
    sampling_func <- loadSamplingFunction(dist)
    
    print(paste0('Starting new distribution: ', dist))
    
    for (k in 1:length(ns)) {
        n <- ns[k]
        p <- pn(n)

        # compute the true radial biases via monte carlo simulation, to compare with the bootstrap estimates
        isoline <- drawIsoline(dist=dist, numCoords=n_coords, gridUbs=ubs, gridLbs=lbs, prob=p)
        thetas_true <- atan2(isoline[,2], isoline[,1])
        r_thetas_true <- sqrt(isoline[,1]**2 + isoline[,2]**2)
        r_hat_theta_mat <- matrix(NA, nrow=n_monte_carlo, ncol=length(thetas_true))

        print(paste0('Starting n: ', n, ' | Computing True MC Radial Bias...'))
        for (m in 1:n_monte_carlo) {
            dat_mc <- sampling_func(n)
            isoline_ext_estimate <- drawExtremeIsoline(dat=dat_mc, p=p, n_coords=n_coords, grid_lbs=lbs, gamma=gamma, xi=xi)
            r_hat_theta_mat[m,] <- sqrt(isoline_ext_estimate[,1]**2 + isoline_ext_estimate[,2]**2)
        }
        
        true_b_thetas <- colMeans(r_hat_theta_mat) - r_thetas_true

        # set up parallel computations
        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        clusterEvalQ(clust, {
                source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')
        })
            
        pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
        progress <- function(iter) setTxtProgressBar(pb, iter)
        opts <- list(progress = progress)

        print(paste0('Evaluating bootstrap radial bias over ', n_sims, ' simulations...'))
        
        # define function to be executed in parallel
        parallelizedCode <- function(ind, true_b_thetas) {
            dat <- sampling_func(n)
            
            # compute bootstrap bias for this specific sample
            boot_res <- estimateBootstrapRadialBias(dat=dat, p=p, B=B_boot, n_coords=n_coords, 
                                                    gamma=gamma, xi=xi, lbs=lbs)
            
            # create a dataframe comparing true vs bootstrap bias point-by-point
            res_df <- data.frame(
                dist = dist,
                p = p,
                n = n,
                sim_id = ind,
                theta = boot_res$thetas,
                expected_r_bootstrapworld = boot_res$boot_expected_r,
                true_r_bootstrapworld = boot_res$base_r,
                true_bias = true_b_thetas,
                boot_bias = boot_res$boot_bias
            )
            
            return(res_df)
        }

        # execute parallel loop
        samp_df_list <- foreach(l = 1:n_sims, 
                .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l, true_b_thetas)

        close(pb)
        stopCluster(clust)

        # Append concise data frames to the CSV
        coverage_df <- do.call(rbind, samp_df_list)
        save_fname <- paste0(args$save_df_path, 'radial_bias_bootstrap_comparison_gamma0.5_knownindex.csv')
        file_exists <- file.exists(save_fname)
            
        write.table(coverage_df, 
                    file = save_fname, 
                    sep = ",", 
                    row.names = FALSE, 
                    col.names = !file_exists, 
                    append = file_exists)     
    }
}