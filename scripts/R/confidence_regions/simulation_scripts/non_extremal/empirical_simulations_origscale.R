# Script to evaluate the coverage properties of the confidence tube method using the empirical survival function,
# by simulating data from known distributions with different distributional characteristics.
#
# MODIFIED: also computes tube SIZE (absolute + relative radial width) for each simulated tube
# by drawing the inner/outer wall isolines and reusing the level-p isoline as the reference.
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

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/confidenceRegions.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')

# parsing command line specified parameters
parser <- ArgumentParser(description = "Args for empirical survival function simulations (Asymmetric).")
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
ns <- c(1000, 5000, 10000, 50000, 100000)
alphas <- c(0.05, 0.1, 0.01)
C <- 5
p_funcs <- list()

# add in central functions
p_funcs[['central']] <- list()
p_funcs[['central']][['0.05']] <- function(n) 0.05
# add in intermediate functions
p_funcs[['intermediate']] <- list()
p_funcs[['intermediate']][['1_3']] <- function(n) n^(-1/3)
p_funcs[['intermediate']][['1_2']] <- function(n) n^(-1/2)
p_funcs[['intermediate']][['2_3']] <- function(n) n^(-2/3)
p_funcs[['extremal']] <- list()
# add in extremal functions
p_funcs[['extremal']][[as.character(C)]] <- function(n) C/n

dists <- c('bivt', 'bivgauss', 'karachi') # distributions to simulate from
# list of lower left hand corners of space, depending on distribution
lbs_lst <- list('bivt'=c(0,0), 'bivgauss'=c(0,0), 'karachi'=c(50, 0))
# list of upper right hand corners of space, depending on distribution
ubs_lst <- list('bivt'=c(200, 200), 'bivgauss'=c(200,200), 'karachi'=c(140,100))

n_sims <- 500 # number of simulations for coverage statistics
B <- 5000 # number of bootstrap iterations for each confidence tube
n_coords <- 50

# command line parameters
n_cores <- args$n_cores
save_full_path <- args$save_full_path


computeEmpiricalRegion <- function(dat, alphas, p, B, lbs=c(0,0), ncoords=n_coords, verbose=FALSE) {
    
    n_dat <- nrow(dat)
    
    # Draw the base isoline
    emp_isoline <- drawEmpiricalIsoline(dat=dat, n_coords=ncoords, grid_lbs=lbs, p)
    
    # Static hit matrix for rapid resampling
    # Determines if data point i exceeds isoline coordinate j
    hit_static <- (outer(dat[,1], emp_isoline[,1], '>') & 
                   outer(dat[,2], emp_isoline[,2], '>')) * 1.0
                   
    boot_probs_mat <- matrix(NA, nrow=B, ncol=nrow(emp_isoline))

    if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    
    # Fast Bootstrap Loop using a tabulate trick
    for (k in 1:B) {
        # Sample just the indices
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        
        # Tabulate occurrences of each index
        w <- tabulate(idx, nbins = n_dat)
        
        # Matrix multiply the weights by the static hit matrix to get empirical probs
        emp_probs <- as.numeric(w %*% hit_static) / n_dat
        
        boot_probs_mat[k, ] <- emp_probs
        
        if (verbose) utils::setTxtProgressBar(pb, k)
    }
    if (verbose) close(pb)
    
    # Bootstrap Centering (using mean instead of p)
    mean_boot_probs <- colMeans(boot_probs_mat)
    boot_devs <- sweep(boot_probs_mat, MARGIN = 2, STATS = mean_boot_probs, FUN = "-")
    
    # Extract global max and min deviations across the curve for each draw
    boot_draws_Wplus <- apply(boot_devs, 1, max)
    boot_draws_Wminus <- apply(boot_devs, 1, min)
    
    c_plus_estimates <- list()
    c_minus_estimates <- list()
    
    # 5. Build Final Asymmetric Bounds
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]        
        c_plus_estimates[[as.character(alpha)]] <- as.numeric(quantile(boot_draws_Wplus, probs=1-(alpha/2), type=1))
        c_minus_estimates[[as.character(alpha)]] <- as.numeric(quantile(boot_draws_Wminus, probs=alpha/2, type=1))
    }

    survFunc <- function(pts) {
        survProb <- vectorizedEmpSurv(pts, dat)
        return(survProb)
    }

    r_ref <- sqrt((emp_isoline[,1] - lbs[1])^2 + (emp_isoline[,2] - lbs[2])^2)

    size_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha  <- alphas[i]
        akey   <- as.character(alpha)
        cplus  <- as.numeric(c_plus_estimates[[akey]])
        cminus <- as.numeric(c_minus_estimates[[akey]])

        # Unboundedness check: c_minus <= -p  <=>  outer survival level p + c_minus <= 0,
        # an unattainable level, so the outer wall never closes.
        if (cminus <= -p) {
            size_estimates[[akey]] <- data.frame(
                is_unbounded     = TRUE,
                width_abs_median = NA_real_,
                width_abs_max    = NA_real_,
                width_rel_median = NA_real_,
                width_rel_max    = NA_real_
            )
        } else {
            s_hi <- p + cplus     # inner wall: higher survival -> smaller radius
            s_lo <- p + cminus    # outer wall: lower  survival -> larger  radius

            inner_iso <- drawEmpiricalIsoline(dat=dat, n_coords=ncoords, grid_lbs=lbs, s_hi)
            outer_iso <- drawEmpiricalIsoline(dat=dat, n_coords=ncoords, grid_lbs=lbs, s_lo)

            r_inner <- sqrt((inner_iso[,1] - lbs[1])^2 + (inner_iso[,2] - lbs[2])^2)
            r_outer <- sqrt((outer_iso[,1] - lbs[1])^2 + (outer_iso[,2] - lbs[2])^2)

            sz <- computeTubeSize(r_inner, r_outer, r_ref)
            size_estimates[[akey]] <- cbind(data.frame(is_unbounded = FALSE), sz)
        }
    }

    # Compile result object
    res_lst <- list(
        dat = dat,
        c_plus_estimates = c_plus_estimates,
        c_minus_estimates = c_minus_estimates,
        p = p,
        survFunc = survFunc,
        emp_isoline = emp_isoline,
        size_estimates = size_estimates,
        transformed=FALSE
    )

    return(res_lst)
}

# part of loop that will be parallelized: the computation of c for each simulation
parallelizedCode <- function(ind) {

    dat <- sampling_func(n)
    
    # Swap out to use the newly defined optimized, asymmetric function
    regions <- computeEmpiricalRegion(
        dat=dat, alphas=alphas, p=p, B=B, 
        lbs=lbs, ncoords=n_coords, verbose=FALSE
    )

    # evaluate coverage for each region
    is_covereds <- evaluateCoverage(regions, isoline)
    regions$is_covereds <- is_covereds
    
    return(regions)
}

# loop over the distributions
for (i in 1:length(dists)) {
    dist <- dists[i]
    lbs <- lbs_lst[[dist]]
    ubs <- ubs_lst[[dist]]
    
    sampling_func <- loadSamplingFunction(dist)
    print(paste0('Starting new distribution: ', dist))

    # loop over regimes
    for (regime in names(p_funcs)) {
        
        for (label in names(p_funcs[[regime]])) {
            
            current_p_func <- p_funcs[[regime]][[label]]
            # define path specific to Dist -> Regime -> Label
            if (!is.null(args$save_full_path)) {
                save_path_full <- paste0(args$save_full_path, dist, '/', regime, '/', label, '/')
                if (!dir.exists(save_path_full)) dir.create(save_path_full, recursive = TRUE)
            }

            print(paste0('  Running Regime: ', regime, ' | Label: ', label))
            
            for (k in 1:length(ns)) {
                n <- ns[k]
                p <- current_p_func(n)
                isoline <- drawIsoline(dist=dist, numCoords=n_coords, gridUbs=ubs, gridLbs=lbs, prob=p)

                print(paste0('Starting new n: ', n))
    
                clust <- makeSOCKcluster(n_cores)
                registerDoSNOW(clust)
                clusterEvalQ(clust, {
                    source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')
                })
                
                pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
                progress <- function(n) setTxtProgressBar(pb ,n)
                opts <- list(progress = progress)
    
                samp_regions <- foreach(l = 1:n_sims, 
                        .options.snow = opts, 
                        .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l)
    
                close(pb)
                stopCluster(clust)
    
                # only save .RData if save_full_path is provided
                if (!is.null(args$save_full_path)) {
                    saveRDS(samp_regions, file=paste0(save_path_full, 'n', n, '_tubes.RData'))
                }
                # save CSV if save_df_path argument is provided
                if (!is.null(args$save_df_path)) {
                    
                    coverage_df <- do.call(rbind, lapply(samp_regions, function(res) {

                    alpha_names <- names(res$is_covereds)
                    # stack the per-alpha size rows in the same order as alpha_names
                    size_df <- do.call(rbind, res$size_estimates[alpha_names])

                    data.frame(
                      dist = dist,
                      regime = regime,
                      label = label,
                      n = n,
                      p = p,
                      alpha = alpha_names,
                      is_covered = unlist(res$is_covereds, use.names=FALSE),
                      c_plus_estimate = unlist(res$c_plus_estimates[alpha_names]),
                      c_minus_estimate = unlist(res$c_minus_estimates[alpha_names]),
                      is_unbounded  = size_df$is_unbounded,
                      width_abs_median = size_df$width_abs_median,
                      width_abs_max    = size_df$width_abs_max,
                      width_rel_median = size_df$width_rel_median,
                      width_rel_max    = size_df$width_rel_max
                    )
                  }))
          
                  save_fname <- paste0(args$save_df_path, 'empirical_simulations_origscale.csv')
                  file_exists <- file.exists(save_fname)
          
                  write.table(coverage_df, 
                          file = save_fname, 
                          sep = ",", 
                          row.names = FALSE, 
                          col.names = !file_exists, 
                          append = file_exists)
                }
            }
        }
    }
}