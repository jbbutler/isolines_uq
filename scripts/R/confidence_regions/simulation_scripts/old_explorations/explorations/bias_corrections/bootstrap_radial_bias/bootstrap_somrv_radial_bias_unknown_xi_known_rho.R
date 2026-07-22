# This script builds upon the success of oracle_somrv_radial_bias.R, where instead of assuming we know the bias from 
# q1 to q2 (the intermediate projections), we now estimate it via a non-parametric bootstrap for each simulated set of
# data. Averaging these estimates across the draws of data, we should hopefully still recover the bias output of
# oracle_somrv_radial_bias.R, and since that output looked VERY close to the true bias, our bias corrected estimators would
# be ~nearly~ unbiased.
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
parser <- ArgumentParser(description = "Args for Empirical SOMRV Radial Bias Exploration.")
parser$add_argument("--n_cores", type = "integer", default = 64, 
                    help = "Number of cores for parallel processing.")
parser$add_argument("--save_df_path", type = "character", default = NULL, 
                    help = "File path to save the csv comparing true vs SOMRV biases.")

args <- parser$parse_args()

# simulation parameters
ns <- c(1000, 5000, 10000, 50000)
C <- 5
rho <- -1/2 # assume we know the second-order MRV parameter (Oracle)
gamma1 <- 1/2
gamma2 <- 2/3

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

# Function to compute the Empirical SOMRV estimates of the radial bias
estimateSomrvBootstrapBias <- function(dat, p, q1, q2, B, n_coords, gamma1, gamma2, rho, lbs) {

    base_radii <- sqrt(dat[,1]^2 + dat[,2]^2)
    base_xi_hat <- estimate_xi_hill(base_radii, gamma1)

    # generate true q2 isoline in the bootstrap world
    isoline_q2_emp <- drawEmpiricalIsoline(dat=dat, n_coords=n_coords, grid_lbs=lbs, p=q2)
    rq2_emp <- sqrt(isoline_q2_emp[,1]^2 + isoline_q2_emp[,2]^2)
    thetas <- atan2(isoline_q2_emp[,2], isoline_q2_emp[,1])
    
    n_dat <- nrow(dat)
    boot_rq2_proj_mat <- matrix(NA, nrow=B, ncol=n_coords)
    
    # compute expectation of projection estimator from q1 to q2
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        boot_dat <- dat[idx, , drop=FALSE]
        
        # re-estimate xi for the bootstrap sample to propagate uncertainty
        boot_radii <- sqrt(boot_dat[,1]^2 + boot_dat[,2]^2)
        
        # project from q1 to q2 using the bootstrap sample and the bootstrap xi
        boot_isoline_q2 <- drawExtremeIsoline(dat=boot_dat, p=q2, n_coords=n_coords, grid_lbs=lbs, gamma=gamma1, xi=base_xi_hat)
        boot_rq2_proj_mat[k, ] <- sqrt(boot_isoline_q2[,1]^2 + boot_isoline_q2[,2]^2)
    }
    
    # calculate estimator bias in bootstrap world
    expected_rq2_boot <- colMeans(boot_rq2_proj_mat)
    in_sample_bias_est <- expected_rq2_boot - rq2_emp
    
    # compute the multiplier M using the base sample's xi_hat
    M <- ((q2/p)^base_xi_hat) * (((q1/p)^rho - 1) / ((q1/q2)^rho - 1))
    
    # extrapolate to extreme target p
    somrv_bias_est <- M * in_sample_bias_est
    
    return(list(thetas = thetas, 
                rq2_emp = rq2_emp, 
                expected_rq2_boot = expected_rq2_boot,
                in_sample_bias_est = in_sample_bias_est,
                somrv_bias_est = somrv_bias_est,
                base_xi_hat = base_xi_hat))
}

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

        # Compute the true radial biases of estimator from q1 to p via monte carlo simulation
        isoline <- drawIsoline(dist=dist, numCoords=n_coords, gridUbs=ubs, gridLbs=lbs, prob=p)
        thetas_true <- atan2(isoline[,2], isoline[,1])
        r_thetas_true <- sqrt(isoline[,1]^2 + isoline[,2]^2)
        r_hat_theta_mat <- matrix(NA, nrow=n_monte_carlo, ncol=length(thetas_true))

        print(paste0('Starting n: ', n, ' | Computing True MC Naive Radial Bias...'))
        for (m in 1:n_monte_carlo) {
            dat_mc <- sampling_func(n)
            
            # 3. Estimate xi for each MC draw to establish the true baseline behavior
            mc_radii <- sqrt(dat_mc[,1]^2 + dat_mc[,2]^2)
            mc_xi_hat <- estimate_xi_hill(mc_radii, gamma1)
            
            # naive projection from q1 out to p
            isoline_ext_estimate <- drawExtremeIsoline(dat=dat_mc, p=p, n_coords=n_coords, grid_lbs=lbs, gamma=gamma1, xi=mc_xi_hat)
            r_hat_theta_mat[m,] <- sqrt(isoline_ext_estimate[,1]^2 + isoline_ext_estimate[,2]^2)
        }
        
        # true bias calculation
        true_bias <- colMeans(r_hat_theta_mat) - r_thetas_true

        # set up parallel computations
        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        clusterEvalQ(clust, {
                source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')
        })
            
        pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
        progress <- function(iter) setTxtProgressBar(pb, iter)
        opts <- list(progress = progress)

        print(paste0('Evaluating Empirical SOMRV radial bias over ', n_sims, ' simulations...'))
        
        # define function to be executed in parallel
        parallelizedCode <- function(ind, true_bias, q1, q2) {
            dat <- sampling_func(n)
            
            # compute SOMRV bias estimate for this specific sample
            somrv_res <- estimateSomrvBootstrapBias(dat=dat, p=p, q1=q1, q2=q2, B=B_boot, n_coords=n_coords, 
                                                    gamma1=gamma1, gamma2=gamma2, rho=rho, lbs=lbs)
            
            # create a dataframe comparing true naive bias vs empirical SOMRV estimate
            res_df <- data.frame(
                dist = dist,
                p = p,
                n = n,
                sim_id = ind,
                theta = somrv_res$thetas,
                true_bias = true_bias,
                somrv_bias_est = somrv_res$somrv_bias_est,
                in_sample_bias_est = somrv_res$in_sample_bias_est,
                xi_hat = somrv_res$base_xi_hat
            )
            
            return(res_df)
        }

        # execute parallel loop
        samp_df_list <- foreach(l = 1:n_sims, 
                .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l, true_bias, q1, q2)

        close(pb)
        stopCluster(clust)

        # Append concise data frames to the CSV
        coverage_df <- do.call(rbind, samp_df_list)
        save_fname <- paste0(args$save_df_path, 'somrv_simulation3_nobootxi.csv')
        file_exists <- file.exists(save_fname)
            
        write.table(coverage_df, 
                    file = save_fname, 
                    sep = ",", 
                    row.names = FALSE, 
                    col.names = !file_exists, 
                    append = file_exists)     
    }
}