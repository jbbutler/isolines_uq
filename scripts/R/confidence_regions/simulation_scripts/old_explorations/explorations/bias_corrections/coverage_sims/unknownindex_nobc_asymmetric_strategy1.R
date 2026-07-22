# Script to compute coverage rates of confidence tube procedure, drawing data from the bivariate t but
# doing no transforms, no sample splitting, and no bias correction. 
# NOW ESTIMATING XI DYNAMICALLY VIA THE HILL ESTIMATOR.
# Mostly serving as a basis of comparison with bias correction strategies to see how much better things are.
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
parser <- ArgumentParser(description = "Args for NO BIAS CORRECTION (Unknown Xi) survival function simulations.")
parser$add_argument("--save_full_path", type = "character", default = NULL, required = FALSE,
                    help = "Optional file path to save full details about simulation (data sample, etc.).")
parser$add_argument("--n_cores", type = "integer", default = 64, 
                    help = "Number of cores for parallel processing.")
parser$add_argument("--save_df_path", type = "character", default = NULL, 
                    help = "Optional file path to save just the csv of coverage results for each simulation.")

args <- parser$parse_args()

if (is.null(args$save_full_path) && is.null(args$save_df_path)) {
    stop("Error: Both --save_full_path and --save_df_path are NULL. 
          \n You must provide at least one output destination so results are saved.")
}

# simulation parameters
ns <- c(1000, 5000, 10000, 50000)
alphas <- c(0.05, 0.1, 0.01) 
C <- 5

gamma <- 1/2 
B <- 2000
n_sims <- 500
n_coords <- 100

pn <- function(n){ C/n } 
dists <- c('bivt')
lbs <- c(0,0)
ubs <- c(200, 200)

n_cores <- args$n_cores
save_full_path <- args$save_full_path

# a function to compute an extreme region given true radial biases (with unknown xi)
computeExtremeRegion <- function(dat, alphas, p, B, gamma, n_coords, lbs=c(0,0), progress_bar=FALSE) {

    construction_dat <- dat
    
    # dynamically estimate xi via Hill estimator
    base_radii <- sqrt(construction_dat[,1]^2 + construction_dat[,2]^2)
    xi_hat <- estimate_xi_hill(base_radii, gamma)

    # construct the naive isoline to take the supremum over
    ext_region_pts <- drawExtremeIsoline(dat=construction_dat, p=p, n_coords=n_coords, grid_lbs=lbs, gamma=gamma, xi=xi_hat)
    theta <- atan2(ext_region_pts[,2], ext_region_pts[,1])  
    
    rp <- sqrt(rowSums(ext_region_pts^2))
    inv_cos <- 1/cos(theta)
    inv_sin <- 1/sin(theta)

    proj_x <- construction_dat[,1] %o% inv_cos
    proj_y <- construction_dat[,2] %o% inv_sin
    Z_static <- pmin(proj_x, proj_y)

    match_x <- outer(construction_dat[,1], ext_region_pts[,1], '>')
    match_y <- outer(construction_dat[,2], ext_region_pts[,2], '>')
    hit_static <- (match_x & match_y)*1

    n_dat <- nrow(construction_dat)
    qn <- n_dat^(-gamma)
    
    boot_draws_Wplus <- rep(NA, B)
    boot_draws_Wminus <- rep(NA, B)

    if (progress_bar) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        Z_boot <- Z_static[idx, , drop=FALSE]
        hit_boot <- hit_static[idx, , drop=FALSE]
        emp_probs <- colMeans(hit_boot)
        
        rq <- colQuantiles(Z_boot, probs=1-qn, type=1, drop=TRUE)
        tail_probs <- ((rq/rp)^(1/xi_hat))*qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)
        
        # Taking the max/min over the entire defined region minus the base p
        boot_draws_Wplus[k] <- max(final_probs - p)
        boot_draws_Wminus[k] <- min(final_probs - p)
        
        if (progress_bar) utils::setTxtProgressBar(pb, k)
    }
    if (progress_bar) close(pb)

    c_plus_estimates <- list()
    c_minus_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]       
        c_plus_estimate <- as.numeric(quantile(boot_draws_Wplus, probs = 1-(alpha/2)))
        c_minus_estimate <- as.numeric(quantile(boot_draws_Wminus, probs = (alpha/2)))
        c_plus_estimates[as.character(alpha)] <- c_plus_estimate
        c_minus_estimates[as.character(alpha)] <- c_minus_estimate
    }

    survFunc <- function(pts) {
        survProb <- vectorizedBlendedSurv(pts, construction_dat, gamma, xi_hat)
        return(survProb)
    }

    res_lst <- list()
    res_lst$dat <- dat
    res_lst$transformed <- FALSE
    res_lst$c_plus_estimates <- c_plus_estimates
    res_lst$c_minus_estimates <- c_minus_estimates
    res_lst$p <- p
    res_lst$gamma <- gamma
    res_lst$xi_hat <- xi_hat
    res_lst$survFunc <- survFunc

    return(res_lst)
}

# loop over the distributions
for (i in 1:length(dists)) {
    dist <- dists[i]
    sampling_func <- loadSamplingFunction(dist)
    
    if (!is.null(save_full_path)) {
        distpath <- paste0(save_full_path, dist, '/')
        if (!dir.exists(distpath)) dir.create(distpath, recursive = TRUE)
    }
    print(paste0('Starting new distribution: ', dist))
    
    for (k in 1:length(ns)) {
        n <- ns[k]
        p <- pn(n)

        # compute true p-isoline
        isoline <- drawIsoline(dist=dist, numCoords=n_coords, gridUbs=ubs, gridLbs=lbs, prob=p)

        thetas <- atan2(isoline[,2], isoline[,1])
        print(paste0('Starting new n: ', n))
        
        # create the n-specific directory BEFORE the parallel loop
        n_dirpath <- NULL
        if (!is.null(save_full_path)) {
            n_dirpath <- paste0(distpath, 'n', n, '/')
            if (!dir.exists(n_dirpath)) dir.create(n_dirpath, recursive = TRUE)
        }

        # set up parallel computations
        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        clusterEvalQ(clust, {
                source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')
        })
            
        pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
        progress <- function(iter) setTxtProgressBar(pb, iter)
        opts <- list(progress = progress)

        # define function to be executed in parallel
        parallelizedCode <- function(ind) {
            dat <- sampling_func(n)
            
            # compute extreme region estimating xi dynamically inside the function
            regions <- computeExtremeRegion(dat=dat, alphas=alphas, p=p, B=B, 
                                            gamma=gamma, n_coords=n_coords)

            is_covereds <- evaluateCoverage(regions, isoline)
            regions$is_covereds <- is_covereds
            
            # save directly to disk from the worker node
            if (!is.null(n_dirpath)) {
                saveRDS(regions, file=paste0(n_dirpath, 'simulation', ind, '.RData'))
            }
            
            # return only the concise row needed for the CSV
            res_df <- data.frame(
                dist = dist,
                p = p,
                n = n,
                alpha = names(regions$is_covereds),
                is_covered = unlist(regions$is_covereds, use.names=FALSE),
                c_plus_estimate = unlist(regions$c_plus_estimates[names(regions$is_covereds)]),
                c_minus_estimate = unlist(regions$c_minus_estimates[names(regions$is_covereds)]),
                xi_hat = regions$xi_hat  # Saved to track shape estimator variance
            )
            
            return(res_df)
        }

        # execute parallel loop, returning a list of simple data frames.
        samp_df_list <- foreach(l = 1:n_sims, 
                .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l)

        close(pb)
        stopCluster(clust)

        # save CSV if save_df_path argument is provided
        if (!is.null(args$save_df_path)) {
            # combine concise data frame rows returned from parallel workers
            coverage_df <- do.call(rbind, samp_df_list)

            # Updated output file name to reflect unknown xi
            save_fname <- paste0(args$save_df_path, 'extremal_coverage_nobc_unknownxi.csv')
            file_exists <- file.exists(save_fname)
                
            # write to CSV using write.table for append support
            write.table(coverage_df, 
                        file = save_fname, 
                        sep = ",", 
                        row.names = FALSE, 
                        col.names = !file_exists, 
                        append = file_exists)     
        }
    }
}