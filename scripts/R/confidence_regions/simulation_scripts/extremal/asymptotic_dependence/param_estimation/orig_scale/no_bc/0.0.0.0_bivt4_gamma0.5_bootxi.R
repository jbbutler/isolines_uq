# Script to compute coverage for confidence tubes in the style of Mammen and Polonik (2013), 
# where tube bounds are determined ONLY by variance in the extreme value estimator.
# No structural bias correction is applied. 
#
# **MODIFIED: This version RE-ESTIMATES xi via the Hill estimator inside every bootstrap loop, in contrast to exp0.1
#
# Jimmy Butler
# June 2026

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
parser <- ArgumentParser(description = "Args for Variance-Only (No Bias Correction) extremal simulations with xi re-estimation.")
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

ns <- c(1000, 5000, 10000, 50000, 100000)
alphas <- c(0.05, 0.1, 0.01)
C <- 5

gamma1 <- 0.5
B <- 5000
n_sims <- 500
n_coords <- 50

pn <- function(n){ C/n } 
dists <- c('bivt')
lbs <- c(0,0)
ubs <- c(200, 200)

n_cores <- args$n_cores
save_full_path <- args$save_full_path

# function that constructs confidence tube(s) using ONLY variance, without bias correction.
computeExtremeRegion <- function(dat, alphas, p, B, gamma1=1/2, ncoords=50, lbs=c(0,0), progress_bar=FALSE, verbose=FALSE) {
    
    n_dat <- nrow(dat)
    q1 <- n_dat^(-gamma1)

    if (verbose) print('Estimating original xi via Hill estimator.')
    base_radii <- sqrt(dat[,1]^2 + dat[,2]^2)
    xi_hat <- estimate_xi_hill(base_radii, gamma1)

    # 1. Isoline Estimates (Using original xi_hat)
    isoline_p <- drawExtremeIsoline(dat=dat, p=p, n_coords=ncoords, grid_lbs=lbs, gamma=gamma1, xi=xi_hat)
    
    if (verbose) print('Estimating pure variance parameters via Centered Bootstrap (with xi re-estimation).')

    # Set up hit matrix for the bootstrap loop directly over the p-isoline
    thetas <- atan2(isoline_p[,2], isoline_p[,1])
    rp <- sqrt(rowSums(isoline_p^2))
    
    # SAFELY construct Z_static avoiding NaN generation on the axes
    Z_static <- matrix(NA, nrow = n_dat, ncol = length(thetas))
    pos_orthant_mask <- (dat[,1] > 0) & (dat[,2] > 0)
    
    ind_0 <- which(thetas == 0)
    if (length(ind_0) > 0) Z_static[, ind_0] <- dat[,1] * pos_orthant_mask
    
    ind_90 <- which(thetas == pi/2)
    if (length(ind_90) > 0) Z_static[, ind_90] <- dat[,2] * pos_orthant_mask
    
    inds_no_ax <- which(!(thetas == 0 | thetas == pi/2))
    if (length(inds_no_ax) > 0) {
        Z_static[, inds_no_ax] <- pmin(
            dat[,1] %o% (1/cos(thetas[inds_no_ax])), 
            dat[,2] %o% (1/sin(thetas[inds_no_ax]))
        )
    }

    hit_static <- (outer(dat[,1], isoline_p[,1], '>') & 
                   outer(dat[,2], isoline_p[,2], '>')) * 1.0

    qn <- n_dat^(-gamma1)
    
    boot_probs_mat <- matrix(NA, nrow=B, ncol=ncoords)

    if (progress_bar) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        
        boot_base_radii <- base_radii[idx]
        xi_boot <- estimate_xi_hill(boot_base_radii, gamma1)
        
        # Tabulate trick for massive memory savings
        w <- tabulate(idx, nbins = n_dat)
        emp_probs <- as.numeric(w %*% hit_static) / n_dat
        
        Z_boot <- Z_static[idx, , drop=FALSE]
        rq <- colQuantiles(Z_boot, probs=1-qn, type=1, drop=TRUE)

        rm(Z_boot) # explicit garbage collection
        
        # NOTE: Using the newly re-estimated xi_boot for projection
        tail_probs <- ((rq/rp)^(1/xi_boot))*qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)
        
        boot_probs_mat[k, ] <- final_probs
        
        if (progress_bar) utils::setTxtProgressBar(pb, k)
    }
    if (progress_bar) close(pb)
    
    rm(Z_static, hit_static) # Final cleanup

    # 5. Bootstrap Centering 
    mean_boot_probs <- colMeans(boot_probs_mat)
    boot_devs <- sweep(boot_probs_mat, MARGIN = 2, STATS = mean_boot_probs, FUN = "-")
    
    boot_draws_Wplus <- apply(boot_devs, 1, max)
    boot_draws_Wminus <- apply(boot_devs, 1, min)

    c_plus_estimates <- list()
    c_minus_estimates <- list()
    
    # 6. Bind Tube Bounds
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]        
        c_plus_raw <- as.numeric(quantile(boot_draws_Wplus, probs = 1-(alpha/2)))
        c_minus_raw <- as.numeric(quantile(boot_draws_Wminus, probs = (alpha/2)))
        
        # Since there is no bias correction, the adjusted estimates are exactly the raw estimates
        c_plus_estimates[as.character(alpha)] <- c_plus_raw
        c_minus_estimates[as.character(alpha)] <- c_minus_raw
    }

    survFunc <- function(pts) {
        survProb <- vectorizedBlendedSurv(pts, dat, gamma1, xi_hat)
        return(survProb)
    }

    # Compile the final result list
    res_lst <- list(
        dat = dat,
        p = p,
        gamma1 = gamma1,
        xi_hat = xi_hat,
        thetas = thetas,
        r_naive = rp,
        c_plus_estimates = c_plus_estimates,
        c_minus_estimates = c_minus_estimates,
        survFunc = survFunc,
        transformed = FALSE
    )

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
            
            # use the Variance-Only bounding function
            regions <- computeExtremeRegion(
                dat=dat, alphas=alphas, p=p, B=B, 
                gamma1=gamma1, ncoords=n_coords, lbs=lbs
            )

            is_covereds <- evaluateCoverage(regions, isoline)
            regions$is_covereds <- is_covereds
            
            # save directly to disk from the worker node
            if (!is.null(n_dirpath)) {
                saveRDS(regions, file=paste0(n_dirpath, 'simulation', ind, '.RData'))
            }
            
            # return comprehensive row needed for the CSV tracking
            res_df <- data.frame(
                dist = dist,
                p = p,
                n = n,
                alpha = names(regions$is_covereds),
                is_covered = unlist(regions$is_covereds, use.names=FALSE),
                c_plus_adj = unlist(regions$c_plus_estimates[names(regions$is_covereds)]),
                c_minus_adj = unlist(regions$c_minus_estimates[names(regions$is_covereds)]),
                xi_hat = regions$xi_hat  
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

            # Updated output filename so it doesn't overwrite your control script output
            save_fname <- paste0(args$save_df_path, 'exp0.2_ctrl.csv')
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