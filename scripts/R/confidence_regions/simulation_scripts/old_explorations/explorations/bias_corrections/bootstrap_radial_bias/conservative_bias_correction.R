# This script computes the Conservative Worst-Case Spatial Envelope using the 
# Bias Ratio trick in the limit as rho -> 0. It completely bypasses the bootstrap
# and the estimation of A(q1) or rho.
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
parser <- ArgumentParser(description = "Args for Conservative Envelope Bound Exploration.")
parser$add_argument("--n_cores", type = "integer", default = 64, 
                    help = "Number of cores for parallel processing.")
parser$add_argument("--save_df_path", type = "character", default = NULL, 
                    help = "File path to save the csv comparing true radii to envelope bounds.")

args <- parser$parse_args()

# simulation parameters
ns <- c(1000, 5000, 10000, 50000, 100000, 500000)
C <- 5
gamma1 <- 1/2
gamma2 <- 2/3

# number of theta to use
n_coords <- 50
n_sims <- 500

pn <- function(n){ C/n } 
dists <- c('bivt')
lbs <- c(0,0)
ubs <- c(200, 200)

n_cores <- args$n_cores

# Function to compute the Conservative Envelope Bounds
estimateConservativeEnvelope <- function(dat, p, q1, q2, n_coords, gamma1, lbs) {

    base_radii <- sqrt(dat[,1]^2 + dat[,2]^2)
    base_xi_hat <- estimate_xi_hill(base_radii, gamma1)

    # 1. Calculate the First-Order Naive Bound at p (Inner or Outer depending on bias)
    isoline_p_naive <- drawExtremeIsoline(dat=dat, p=p, n_coords=n_coords, grid_lbs=lbs, gamma=gamma1, xi=base_xi_hat)
    Rp_first_order <- sqrt(isoline_p_naive[,1]^2 + isoline_p_naive[,2]^2)
    thetas <- atan2(isoline_p_naive[,2], isoline_p_naive[,1])

    # 2. Empirical Anchor at q2
    isoline_q2_emp <- drawEmpiricalIsoline(dat=dat, n_coords=n_coords, grid_lbs=lbs, p=q2)
    Rq2_emp <- sqrt(isoline_q2_emp[,1]^2 + isoline_q2_emp[,2]^2)
    
    # 3. Naive Projection at q2
    isoline_q2_proj <- drawExtremeIsoline(dat=dat, p=q2, n_coords=n_coords, grid_lbs=lbs, gamma=gamma1, xi=base_xi_hat)
    Rq2_proj <- sqrt(isoline_q2_proj[,1]^2 + isoline_q2_proj[,2]^2)
    
    # 4. In-Sample Gap at q2
    in_sample_gap <- Rq2_proj - Rq2_emp
    
    # 5. Worst-Case Multiplier (rho -> 0 Limit)
    log_ratio <- log(q1 / p) / log(q1 / q2)
    M_limit <- ((q2 / p)^base_xi_hat) * log_ratio
    
    # 6. Maximal Second-Order Bound at p
    Rp_worst_case <- Rp_first_order - (M_limit * in_sample_gap)
    
    # 7. Construct Envelope Boundaries
    # The true boundary is mathematically bounded between the zero-correction and maximal-correction projections
    r_lower <- pmin(Rp_first_order, Rp_worst_case)
    r_upper <- pmax(Rp_first_order, Rp_worst_case)
    
    return(list(thetas = thetas, 
                r_naive = Rp_first_order,
                r_lower = r_lower,
                r_upper = r_upper,
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

        # Compute the TRUE mathematical isoline radii at probability p
        isoline_true <- drawIsoline(dist=dist, numCoords=n_coords, gridUbs=ubs, gridLbs=lbs, prob=p)
        thetas_true <- atan2(isoline_true[,2], isoline_true[,1])
        r_thetas_true <- sqrt(isoline_true[,1]^2 + isoline_true[,2]^2)

        # set up parallel computations
        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        clusterEvalQ(clust, {
                source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')
        })
            
        pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
        progress <- function(iter) setTxtProgressBar(pb, iter)
        opts <- list(progress = progress)

        print(paste0('Evaluating Conservative Envelope Bounds over ', n_sims, ' simulations...'))
        
        # define function to be executed in parallel
        parallelizedCode <- function(ind, r_thetas_true, q1, q2) {
            dat <- sampling_func(n)
            
            # compute Conservative Envelope Bounds for this specific sample
            env_res <- estimateConservativeEnvelope(dat=dat, p=p, q1=q1, q2=q2, n_coords=n_coords, 
                                                    gamma1=gamma1, lbs=lbs)
            
            # create a dataframe saving true radius vs estimated bounds
            res_df <- data.frame(
                dist = dist,
                p = p,
                n = n,
                sim_id = ind,
                theta = env_res$thetas,
                r_true = r_thetas_true,
                r_naive = env_res$r_naive,
                r_lower = env_res$r_lower,
                r_upper = env_res$r_upper,
                xi_hat = env_res$base_xi_hat
            )
            
            return(res_df)
        }

        # execute parallel loop
        samp_df_list <- foreach(l = 1:n_sims, 
                .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l, r_thetas_true, q1, q2)

        close(pb)
        stopCluster(clust)

        # Append concise data frames to the CSV
        coverage_df <- do.call(rbind, samp_df_list)
        save_fname <- paste0(args$save_df_path, 'conservative_envelope_simulation.csv')
        file_exists <- file.exists(save_fname)
            
        write.table(coverage_df, 
                    file = save_fname, 
                    sep = ",", 
                    row.names = FALSE, 
                    col.names = !file_exists, 
                    append = file_exists)     
    }
}