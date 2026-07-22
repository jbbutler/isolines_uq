# Script to evaluate the coverage properties of the confidence tube method using the empirical survival function,
# by simulating data from known distributions with different distributional characteristics.
# Specifically using the methodology which splits the sample in 3, one group for isoline estimation
# one group for bootstrapping, and a separate group for evaluating the survival function..
#
# Jimmy Butler
# February 2025

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
source('~/isolines_uq/scripts/R/confidence_regions/modules/confidenceRegions.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')

# parsing command line specified parameters
parser <- ArgumentParser(description = "Args for empirical survival function simulations.")
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
ns <- 3*c(1000, 5000, 10000, 50000) # sample sizes to test, 3 times because sample will be split in thirds..
alphas <- c(0.05, 0.1, 0.01) # alphas for 1-alpha tubes
C <- 15
gamma <- 2/3
pn <- function(n){ C/n } # exceedance probabilities
dists <- c('bivt') # distributions to simulate from
# list of lower left hand corners of space, depending on distribution
lbs_lst <- list('bivt'=c(-2,-2), 'bivgauss'=c(-2,-2), 'karachi'=c(50, 0))
# list of upper right hand corners of space, depending on distribution
# note: this is for determining which points are on the isoline, not actually used to draw tubes
# need to be big enough for root-finding algorithm
ubs_lst <- list('bivt'=c(200, 200), 'bivgauss'=c(200,200), 'karachi'=c(140,100))

n_sims <- 500 # number of simulations for coverage statistics
B <- 2000 # number of bootstrap iterations for each confidence tube

# command line parameters
n_cores <- args$n_cores
save_full_path <- args$save_full_path

# part of loop that will be parallelized: the computation of c for each simulation
parallelizedCode <- function(ind) {

    # sample dataset
    dat <- sampling_func(n)
    
    # compute confidence regions for different alpha, for a particular p
    regions <- computeExtremeRegionThreewayIsoBootSurvTwoCs(dat=dat,
                                 alphas=alphas,
                                 p=p,
                                 B=B,
                                 gamma=gamma,
                                 xi=1,
                                 verbose=FALSE)

    transformed_isoline <- data.frame(X1=regions$transform_func1(isoline[,1]),
                                  X2=regions$transform_func2(isoline[,2]))

    # evaluate coverage for each region
    is_covereds <- evaluateCoverage(regions, transformed_isoline)
    # add into the list
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

# loop over the distributions
for (i in 1:length(dists)) {
    dist <- dists[i]
    lbs <- lbs_lst[[dist]]
    ubs <- ubs_lst[[dist]]
    
    sampling_func <- loadSamplingFunction(dist)
    
    # only create dist directory if save_full_path is provided
    distpath <- NULL
    if (!is.null(save_full_path)) {
        distpath <- paste0(save_full_path, dist, '/')
        if (!dir.exists(distpath)) dir.create(distpath, recursive = TRUE)
    }
    
    print(paste0('Starting new distribution: ', dist))
    
    for (k in 1:length(ns)) {
        n <- ns[k]
        p <- pn(n)
        
        # create the true p-isoline for this distribution
        isoline <- drawIsoline(dist=dist, numCoords=200, gridUbs=ubs, gridLbs=lbs, prob=p)
        print(paste0('Starting new n: ', n))

        # create the n-specific directory BEFORE the parallel loop
        n_dirpath <- NULL
        if (!is.null(distpath)) {
            n_dirpath <- paste0(distpath, 'n', n, '/')
            if (!dir.exists(n_dirpath)) dir.create(n_dirpath, recursive = TRUE)
        }

        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        
        # we need to make sure the workers load the required scripts
        clusterEvalQ(clust, {
                source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')
        })
            
        pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
        progress <- function(n_prog) setTxtProgressBar(pb ,n_prog)
        opts <- list(progress = progress)

        # Execute parallel loop, relying on foreach's auto-export for global variables
        samp_df_list <- foreach(l = 1:n_sims, 
                .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l)

        close(pb)
        stopCluster(clust)

        # save CSV if save_df_path argument is provided
        if (!is.null(args$save_df_path)) {
                
            # combine concise data frame rows returned from parallel workers
            coverage_df <- do.call(rbind, samp_df_list)

            # determine if we need to write headers (only if file doesn't exist)
            save_fname <- paste0(args$save_df_path, 'extremal_coverage_threeway_iso_2000boot_surv_C15_twoCs.csv')
            file_exists <- file.exists(save_fname)
                
            # write to CSV using write.table for append support
            write.table(coverage_df, 
                        file = save_fname, 
                        sep = ",", 
                        row.names = FALSE, 
                        col.names = !file_exists, # Write headers only if new file
                        append = file_exists)     # Append if file exists
        }
    }
}