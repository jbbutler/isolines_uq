# Script to evaluate the coverage properties of the confidence tube method using the 
# TRUE EXPECTED BIAS, computed via Monte Carlo simulation before the bootstrap phase.
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
parser <- ArgumentParser(description = "Args for TRUE EXPECTED BIAS survival function simulations.")
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

gamma <- 2/3 
B <- 2000
n_sims <- 500
M_mc <- 1000 # number of Monte Carlo iterations to compute the true expectation

pn <- function(n){ C/n } 
dists <- c('bivt') 
lbs_lst <- list('bivt'=c(-2,-2), 'bivgauss'=c(-2,-2), 'karachi'=c(50, 0))
ubs_lst <- list('bivt'=c(200, 200), 'bivgauss'=c(200,200), 'karachi'=c(140,100))

n_cores <- args$n_cores
save_full_path <- args$save_full_path

computeExtremeRegionBC_GivenBounds <- function(dat, alphas, p, B, gamma, b1, b2, ncoords=100, n_lines=10, xi=1, lbs=c(0,0), 
                                               sup_region=TRUE, transformed=TRUE, progress_bar=FALSE, verbose=FALSE) {
    # function that constructs confidence tube(s) given a particular dataset and bounds b1 and b2.
    if (transformed) {
        # define a function for the transform
        transform <- function(pts) {
            transformed_pts <- 1/(1-pt(pts, df=4)) - 1
            return(transformed_pts)
        }

        # define a function for the inverse transform
        inv_transform <- function(pts) {
            orig_pts <- qt(pts/(1+pts), df=4)
            return(orig_pts)
        }

        # transform the data. We will now only work with the transformed data.
        transformed_dat_X1 <- transform(dat[,1])
        transformed_dat_X2 <- transform(dat[,2])
        construction_dat <- data.frame(X1=transformed_dat_X1, X2=transformed_dat_X2)

    } else {
        construction_dat <- dat
    }

    # Construct the supremum region B using the provided b1 and b2 bounds
    if (sup_region) {
        if (verbose) print('Constructing region using provided true expected bounds.')
        
        p_seq <- seq(b1, b2, length.out = n_lines)
        # generate points for all isolines in the sequence and stack them
        region_pts_list <- lapply(p_seq, function(prob) {
            drawExtremeIsoline(construction_dat, prob, n_coords=ncoords, grid_lbs=lbs, gamma, xi)
        })
        ext_region_pts <- do.call(rbind, region_pts_list)
        
    } else {
        ext_region_pts <- drawExtremeIsoline(construction_dat, p, n_coords=ncoords, grid_lbs=lbs, gamma, xi)
    }

    if (verbose) print('Estimating tube size parameters via bootstrap.')

    # Set up hit matrix for the bootstrap loop over the constructed region
    theta <- atan2(ext_region_pts[,2], ext_region_pts[,1])
    rp <- sqrt(rowSums(ext_region_pts^2))
    
    n_dat <- nrow(construction_dat)
    
    # SAFELY construct Z_static avoiding NaN generation on the axes
    Z_static <- matrix(NA, nrow = n_dat, ncol = length(theta))
    pos_orthant_mask <- (construction_dat[,1] > 0) & (construction_dat[,2] > 0)
    
    ind_0 <- which(theta == 0)
    if (length(ind_0) > 0) Z_static[, ind_0] <- construction_dat[,1] * pos_orthant_mask
    
    ind_90 <- which(theta == pi/2)
    if (length(ind_90) > 0) Z_static[, ind_90] <- construction_dat[,2] * pos_orthant_mask
    
    inds_no_ax <- which(!(theta == 0 | theta == pi/2))
    if (length(inds_no_ax) > 0) {
        Z_static[, inds_no_ax] <- pmin(
            construction_dat[,1] %o% (1/cos(theta[inds_no_ax])), 
            construction_dat[,2] %o% (1/sin(theta[inds_no_ax]))
        )
    }

    hit_static <- (outer(construction_dat[,1], ext_region_pts[,1], '>') & 
                   outer(construction_dat[,2], ext_region_pts[,2], '>')) * 1.0

    qn <- n_dat^(-gamma)
    
    boot_draws_Wplus <- rep(NA, B)
    boot_draws_Wminus <- rep(NA, B)

    if (progress_bar) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        
        # Tabulate trick for massive memory savings
        w <- tabulate(idx, nbins = n_dat)
        emp_probs <- as.numeric(w %*% hit_static) / n_dat
        
        Z_boot <- Z_static[idx, , drop=FALSE]
        rq <- colQuantiles(Z_boot, probs=1-qn, type=1, drop=TRUE)

        rm(Z_boot) # explicit garbage collection
        
        tail_probs <- ((rq/rp)^(1/xi))*qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)
        
        # Taking the max/min over the entire defined region minus the base p
        boot_draws_Wplus[k] <- max(final_probs - p)
        boot_draws_Wminus[k] <- min(final_probs - p)
        
        if (progress_bar) utils::setTxtProgressBar(pb, k)
    }
    if (progress_bar) close(pb)
    
    rm(Z_static, hit_static) # Final cleanup

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
        survProb <- vectorizedBlendedSurv(pts, construction_dat, gamma, xi)
        return(survProb)
    }

    res_lst <- list()
    res_lst$dat <- dat
    res_lst$transformed <- transformed
    res_lst$c_plus_estimates <- c_plus_estimates
    res_lst$c_minus_estimates <- c_minus_estimates
    res_lst$p <- p
    res_lst$gamma <- gamma
    res_lst$xi <- xi
    res_lst$survFunc <- survFunc

    if (sup_region) {
        res_lst$b1 <- b1
        res_lst$b2 <- b2
    }

    if (transformed) {
        res_lst$transform_func1 <- function(pts) transform(pts)
        res_lst$transform_func2 <- function(pts) transform(pts)
        res_lst$inv_transform_func1 <- function(pts) inv_transform(pts)
        res_lst$inv_transform_func2 <- function(pts) inv_transform(pts)  
    }

    return(res_lst)
}

# function to compute the true relative bias bounds (Axis Safe & Memory Optimized)
computeTrueExpectedBiasBounds <- function(sampling_func, n, p, gamma, true_isoline, xi=1, M=1000) {
    n_angles <- nrow(true_isoline)
    prob_estimates_matrix <- matrix(NA, nrow = M, ncol = n_angles)
    
    for (m in 1:M) {
        # draw independent sample of size n
        dat <- sampling_func(n)
        
        # transform marginals, assuming you know the marginals
        apply_transform <- function(pts) {
            1/(1-pt(pts, df=4)) - 1
        }
        
        # transform the TRUE isoline into this sample's empirical Pareto space
        trans_isoline <- cbind(
            apply_transform(true_isoline[,1]),
            apply_transform(true_isoline[,2])
        )
        
        construction_dat <- data.frame(
            X1 = apply_transform(dat[,1]),
            X2 = apply_transform(dat[,2])
        )
        
        # evaluate the extremal survival function
        theta_true <- atan2(trans_isoline[,2], trans_isoline[,1])
        rp_true <- sqrt(rowSums(trans_isoline^2))
        
        n_dat <- nrow(construction_dat)
        
        # SAFELY construct Z_static_true avoiding NaN generation on the axes
        Z_static_true <- matrix(NA, nrow = n_dat, ncol = length(theta_true))
        pos_orthant_mask <- (construction_dat[,1] > 0) & (construction_dat[,2] > 0)
        
        ind_0 <- which(theta_true == 0)
        if (length(ind_0) > 0) Z_static_true[, ind_0] <- construction_dat[,1] * pos_orthant_mask
        
        ind_90 <- which(theta_true == pi/2)
        if (length(ind_90) > 0) Z_static_true[, ind_90] <- construction_dat[,2] * pos_orthant_mask
        
        inds_no_ax <- which(!(theta_true == 0 | theta_true == pi/2))
        if (length(inds_no_ax) > 0) {
            Z_static_true[, inds_no_ax] <- pmin(
                construction_dat[,1] %o% (1/cos(theta_true[inds_no_ax])), 
                construction_dat[,2] %o% (1/sin(theta_true[inds_no_ax]))
            )
        }

        hit_static_true <- (outer(construction_dat[,1], trans_isoline[,1], '>') & 
                            outer(construction_dat[,2], trans_isoline[,2], '>')) * 1.0
        
        qn <- n^(-gamma)
        emp_probs_true <- colMeans(hit_static_true)
        rq_true <- colQuantiles(Z_static_true, probs=1-qn, type=1, drop=TRUE)
        tail_probs_true <- ((rq_true/rp_true)^(1/xi)) * qn
        
        # store the estimated probabilities along the true isoline for this iteration
        prob_estimates_matrix[m, ] <- ifelse(emp_probs_true >= qn, emp_probs_true, tail_probs_true)
        
        # garbage collection: free the matrices immediately before the next M loop
        rm(Z_static_true, hit_static_true, trans_isoline, construction_dat, dat)
    }
    
    # compute the expectation over all M independent samples
    expected_probs <- colMeans(prob_estimates_matrix)
    
    # calculate true relative bias b(r_theta, theta)
    biases <- (expected_probs / p) - 1
    
    # extract the bounding scalars anchoring to 0
    b1 <- p / (1 + max(max(biases), 0))
    b2 <- p / (1 + min(min(biases), 0))
    
    return(list(b1 = b1, b2 = b2, expected_biases = biases))
}

# loop over the distributions
for (i in 1:length(dists)) {
    dist <- dists[i]
    lbs <- lbs_lst[[dist]]
    ubs <- ubs_lst[[dist]]
    
    sampling_func <- loadSamplingFunction(dist)
    
    if (!is.null(save_full_path)) {
        distpath <- paste0(save_full_path, dist, '/')
        if (!dir.exists(distpath)) dir.create(distpath, recursive = TRUE)
    }
    print(paste0('Starting new distribution: ', dist))
    
    for (k in 1:length(ns)) {
        n <- ns[k]
        p <- pn(n)
        
        # compute the TRUE p-isoline for this specific distribution
        isoline <- drawIsoline(dist=dist, numCoords=200, gridUbs=ubs, gridLbs=lbs, prob=p)
        print(paste0('Starting new n: ', n))
        
        # create the n-specific directory BEFORE the parallel loop
        n_dirpath <- NULL
        if (!is.null(save_full_path)) {
            n_dirpath <- paste0(distpath, 'n', n, '/')
            if (!dir.exists(n_dirpath)) dir.create(n_dirpath, recursive = TRUE)
        }
        
        # pre-compute the bias using Monte Carlo
        print(paste0('Calculating True Expected Bias using M=', M_mc, ' iterations...'))
        true_bias_bounds <- computeTrueExpectedBiasBounds(sampling_func=sampling_func, n=n, p=p, 
                                                          gamma=gamma, true_isoline=isoline, M=M_mc)
        print(sprintf("True Expected Bounds for n=%d: b1=%g, b2=%g", n, true_bias_bounds$b1, true_bias_bounds$b2))

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
        parallelizedCode <- function(ind, fixed_b1, fixed_b2) {
            dat <- sampling_func(n)
            
            # use the new bounded function
            regions <- computeExtremeRegionBC_GivenBounds(dat=dat, alphas=alphas, p=p, B=B, 
                                                          gamma=gamma, b1=fixed_b1, b2=fixed_b2, ncoords=50, 
                                                          sup_region=TRUE, transformed=TRUE, progress_bar=FALSE)

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
                c_minus_estimate = unlist(regions$c_minus_estimates[names(regions$is_covereds)])
            )
            
            return(res_df)
        }

        # execute parallel loop, returning a list of simple data frames.
        samp_df_list <- foreach(l = 1:n_sims, 
                .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l, true_bias_bounds$b1, true_bias_bounds$b2)

        close(pb)
        stopCluster(clust)

        # save CSV if save_df_path argument is provided
        if (!is.null(args$save_df_path)) {
            # combine concise data frame rows returned from parallel workers
            coverage_df <- do.call(rbind, samp_df_list)

            save_fname <- paste0(args$save_df_path, 'extremal_coverage_true_bias_knownmarginals_strategy1.csv')
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